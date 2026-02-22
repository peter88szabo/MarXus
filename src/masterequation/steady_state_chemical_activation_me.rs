use super::collisional_relaxation::{
    band_limits, compute_alpha_cm1, compute_collision_frequency_s_inv,
};
use super::state_index::GlobalLayout;
use super::reaction_network::{
    ChemicalActivationDefinition, CollisionKernelImplementation, MasterEquationSettings,
    MicrocanonicalProvider, SolutionResults, WellDefinition,
};
use crate::numeric::iterative_solvers::{solve_bicgstab_left_jacobi_dense, BiCgStabDiagnostics};
use crate::numeric::ldlt_solvers::{solve_symmetric_indefinite_ldlt_bunch_kaufman, LdltDiagnostics};
use crate::numeric::linear_algebra::{
    cholesky_solve_spd_with_diagnostics, CholeskyDiagnostics, DenseMatrix, DiagonalScale,
};

pub enum LinearSolveMethod {
    CholeskySpd,
    LdltSymmetricIndefinite,
    BiCgStab,
}

pub struct MasterEquationSolveDiagnostics {
    pub collision_conservation_max_abs: f64,
    pub transformed_symmetry_relative_frobenius: f64,
    pub solve_method: LinearSolveMethod,
    pub cholesky: Option<CholeskyDiagnostics>,
    pub ldlt: Option<LdltDiagnostics>,
    pub bicgstab: Option<BiCgStabDiagnostics>,
}

#[derive(Default)]
struct CollisionDiagnostics {
    max_column_sum_abs: f64,
}

fn absolute_energy_cm1(well: &WellDefinition, local_grain: usize) -> f64 {
    let abs_grain = (local_grain as isize) + well.alignment_offset_in_grains;
    (abs_grain as f64) * well.energy_grain_width_cm1
}

fn map_grain_by_aligned_energy(
    from_well: &WellDefinition,
    from_grain: usize,
    to_well: &WellDefinition,
) -> Result<Option<usize>, String> {
    let e_abs = absolute_energy_cm1(from_well, from_grain);
    let delta_to = to_well.energy_grain_width_cm1;
    if delta_to <= 0.0 {
        return Err("Non-positive energy_grain_width_cm1 in target well.".into());
    }

    let exact = (e_abs / delta_to) - (to_well.alignment_offset_in_grains as f64);
    let rounded = exact.round();
    let to_grain_isize = rounded as isize;
    if to_grain_isize < 0 {
        return Ok(None);
    }
    let to_grain = to_grain_isize as usize;

    let e_back = absolute_energy_cm1(to_well, to_grain);
    let mismatch = (e_abs - e_back).abs();
    let tol = 0.5 * delta_to + 1e-10 * delta_to;
    if mismatch > tol {
        let e_min = absolute_energy_cm1(to_well, to_well.lowest_included_grain_index);
        let e_max = absolute_energy_cm1(
            to_well,
            to_well.one_past_highest_included_grain_index.saturating_sub(1),
        );
        let lo = e_min.min(e_max) - delta_to;
        let hi = e_min.max(e_max) + delta_to;
        if e_abs >= lo && e_abs <= hi {
            return Err(format!(
                "Energy alignment mismatch between wells (ΔE differs or offsets inconsistent): E_abs={:.6} cm^-1 cannot be mapped within tolerance {:.6} cm^-1.",
                e_abs, tol
            ));
        }
        return Ok(None);
    }

    if to_grain < to_well.lowest_included_grain_index
        || to_grain >= to_well.one_past_highest_included_grain_index
    {
        return Ok(None);
    }

    Ok(Some(to_grain))
}

/// Main engine: assembles and solves the steady-state chemical activation master equation.
pub struct MasterEquationEngine {
    wells: Vec<WellDefinition>,
    settings: MasterEquationSettings,
}

impl MasterEquationEngine {
    pub fn new(wells: Vec<WellDefinition>, settings: MasterEquationSettings) -> Self {
        Self { wells, settings }
    }

    pub fn solve_steady_state_chemical_activation_with_diagnostics(
        &self,
        micro: &dyn MicrocanonicalProvider,
        activation: ChemicalActivationDefinition,
    ) -> Result<(SolutionResults, MasterEquationSolveDiagnostics), String> {
        let mut collision_diag = CollisionDiagnostics::default();
        let global_layout = GlobalLayout::from_wells(&self.wells)?;

        // 1) Assemble the raw operator L (collisions + sinks + inter-well couplings)
        let mut raw_operator = DenseMatrix::zeros(global_layout.total_state_count);

        // 1a) Collisions within each well
        for (well_index, well) in self.wells.iter().enumerate() {
            self.add_collision_terms_for_well(
                micro,
                &global_layout,
                well_index,
                well,
                &mut raw_operator,
                &mut collision_diag,
            )?;
        }

        // 1b) Reaction sinks (diagonal loss from all channels)
        self.add_reaction_sinks(micro, &global_layout, &mut raw_operator)?;

        // 1c) Inter-well couplings (off-diagonal population transfer for internal channels)
        self.add_interwell_couplings(micro, &global_layout, &mut raw_operator)?;

        // 2) Build chemical activation source vector s (raw, then normalized)
        let source_vector = self.build_activation_source(micro, &global_layout, &activation)?;

        // 3) Similarity transform weights W (diagonal)
        let similarity_scale = self.build_similarity_scale(micro, &global_layout)?;

        // Transform matrix A = W L W^{-1}, rhs b = W^{-1} s
        let transformed_operator = raw_operator.similarity_transform(&similarity_scale);
        let transformed_rhs = similarity_scale.inverse_apply_to_vector(&source_vector);

        // We solve (-A) q = b.
        let negative_transformed = transformed_operator.scaled(-1.0);
        let symmetry_rel = negative_transformed.symmetry_relative_frobenius();

        // 4) Solve using SPD Cholesky if possible, else fall back.
        let mut diag = MasterEquationSolveDiagnostics {
            collision_conservation_max_abs: collision_diag.max_column_sum_abs,
            transformed_symmetry_relative_frobenius: symmetry_rel,
            solve_method: LinearSolveMethod::CholeskySpd,
            cholesky: None,
            ldlt: None,
            bicgstab: None,
        };

        let transformed_solution = match cholesky_solve_spd_with_diagnostics(
            &negative_transformed,
            &transformed_rhs,
        ) {
            Ok((x, chol)) => {
                diag.solve_method = LinearSolveMethod::CholeskySpd;
                diag.cholesky = Some(chol);
                x
            }
            Err(chol_err) => {
                // If matrix is close to symmetric, try LDLT; otherwise fall back to BiCGSTAB.
                let sym_tol = 1e-12;
                if symmetry_rel <= sym_tol {
                    match solve_symmetric_indefinite_ldlt_bunch_kaufman(
                        &negative_transformed,
                        &transformed_rhs,
                    ) {
                        Ok((x, ldlt)) => {
                            diag.solve_method = LinearSolveMethod::LdltSymmetricIndefinite;
                            diag.ldlt = Some(ldlt);
                            x
                        }
                        Err(ldlt_err) => {
                            // Final fallback: BiCGSTAB with left Jacobi preconditioner.
                            let (x, bicg) = solve_bicgstab_left_jacobi_dense(
                                &negative_transformed,
                                &transformed_rhs,
                                1e-10,
                                8000,
                            )
                            .map_err(|e| {
                                format!(
                                    "All solvers failed.\nCholesky: {chol_err}\nLDLT: {ldlt_err}\nBiCGSTAB: {e}"
                                )
                            })?;
                            diag.solve_method = LinearSolveMethod::BiCgStab;
                            diag.bicgstab = Some(bicg);
                            x
                        }
                    }
                } else {
                    let (x, bicg) = solve_bicgstab_left_jacobi_dense(
                        &negative_transformed,
                        &transformed_rhs,
                        1e-10,
                        8000,
                    )
                    .map_err(|e| format!("Cholesky failed: {chol_err}\nBiCGSTAB failed: {e}"))?;
                    diag.solve_method = LinearSolveMethod::BiCgStab;
                    diag.bicgstab = Some(bicg);
                    x
                }
            }
        };

        // 5) Recover p = W q
        let steady_state_population = similarity_scale.apply_to_vector(&transformed_solution);

        // 6) Compute macroscopic (averaged) rate constants
        let macro_rates =
            self.compute_macroscopic_rates(micro, &global_layout, &steady_state_population)?;

        Ok((macro_rates, diag))
    }

    /// Solve the steady-state chemical activation problem.
    ///
    /// Returns:
    /// - steady_state_population: p (stacked over all wells/grains)
    /// - macroscopic_rate_constants: per (well, channel) averaged rates k(T,p)
    /// - total_outgoing_rate_constant: sum of outgoing channels (+ stabilization leakage if truncation is used)
    pub fn solve_steady_state_chemical_activation(
        &self,
        micro: &dyn MicrocanonicalProvider,
        activation: ChemicalActivationDefinition,
    ) -> Result<SolutionResults, String> {
        let (res, _diag) =
            self.solve_steady_state_chemical_activation_with_diagnostics(micro, activation)?;
        Ok(res)
    }

    /// Adds collision terms for one well using a band-limited kernel with detailed balance.
    ///
    /// Kernel construction (implemented here):
    /// - For a given well and grain j, we define a *proposal* weight for i within a band |i-j|<=bw:
    ///
    ///   w(j->i) = exp( -|i-j| * ΔE / alpha(T) ) * equilibrium_correction(j,i)
    ///
    /// - Detailed balance is enforced using an equilibrium-like factor:
    ///
    ///   equilibrium_correction(j,i) = [ρ(i)/ρ(j)] * exp( -(i-j)*ΔE / (k_B T) )   for i>j
    ///                               = 1                                       for i<=j
    ///
    /// - Then probabilities are normalized:
    ///
    ///   P(j->i) = w(j->i) / Σ_i w(j->i)
    ///
    /// Collision operator contribution (continuous-time Markov chain):
    ///
    ///   L(i,j) += z * P(j->i)     for i != j
    ///   L(j,j) += z * (P(j->j) - 1) = -z * Σ_{i!=j} P(j->i)
    ///
    fn add_collision_terms_for_well(
        &self,
        micro: &dyn MicrocanonicalProvider,
        layout: &GlobalLayout,
        well_index: usize,
        well: &WellDefinition,
        operator: &mut DenseMatrix,
        diag: &mut CollisionDiagnostics,
    ) -> Result<(), String> {
        let temperature = self.settings.temperature_kelvin;
        let pressure = self.settings.pressure_torr;

        let alpha_cm1 = compute_alpha_cm1(
            well.collision_params.alpha_at_1000K_cm1,
            well.collision_params.alpha_temperature_exponent,
            temperature,
        );

        let collision_frequency_s_inv = compute_collision_frequency_s_inv(
            &well.collision_params,
            &self.settings,
            temperature,
            pressure,
        )?;

        let local_start = well.lowest_included_grain_index;
        let local_end_exclusive = well.one_past_highest_included_grain_index;
        let band = self.settings.collision_band_half_width;

        let grain_width = well.energy_grain_width_cm1;
        let boltzmann = self.settings.boltzmann_constant_wavenumber_per_kelvin;

        match self.settings.collision_kernel_implementation {
            CollisionKernelImplementation::Mess => {
                // MESS-like kernel: per-source normalization *including* a ΔE=0 self-weight.
                //
                // Implementation detail: include `target_grain == source_grain` with weight 1.0
                // so that p_same = 1 / nfac and Σ_{i≠j} p_move = 1 - p_same.
                for source_grain in local_start..local_end_exclusive {
                    let (target_min, target_max_exclusive) =
                        band_limits(source_grain, local_start, local_end_exclusive, band);

                    let rho_source = micro.density_of_states(well_index, source_grain);
                    if rho_source <= 0.0 {
                        return Err(format!(
                            "Non-positive density of states at well={}, grain={}",
                            well.well_name, source_grain
                        ));
                    }

                    // Build unnormalized weights w(source->target)
                    let mut unnormalized_weights: Vec<(usize, f64)> = Vec::new();
                    for target_grain in target_min..target_max_exclusive {
                        let delta_grains = if target_grain >= source_grain {
                            (target_grain - source_grain) as i64
                        } else {
                            -((source_grain - target_grain) as i64)
                        };

                        // exp(-|ΔE|/alpha)
                        let step_energy_cm1 = (delta_grains.abs() as f64) * grain_width;
                        let exponential_down = (-step_energy_cm1 / alpha_cm1).exp();

                        // enforce detailed balance for upward moves via equilibrium-like factor
                        let detailed_balance_factor = if target_grain > source_grain {
                            let rho_target = micro.density_of_states(well_index, target_grain);
                            if rho_target <= 0.0 {
                                return Err(format!(
                                    "Non-positive density of states at well={}, grain={}",
                                    well.well_name, target_grain
                                ));
                            }
                            let upward_energy_cm1 =
                                ((target_grain - source_grain) as f64) * grain_width;
                            // exp( -ΔE / (k_B T) ) * ρ(target)/ρ(source)
                            (rho_target / rho_source)
                                * (-upward_energy_cm1 / (boltzmann * temperature)).exp()
                        } else {
                            1.0
                        };

                        let weight = exponential_down * detailed_balance_factor;
                        unnormalized_weights.push((target_grain, weight.max(0.0)));
                    }

                    let sum_weights: f64 = unnormalized_weights.iter().map(|(_, w)| *w).sum();
                    if sum_weights <= 0.0 || !sum_weights.is_finite() {
                        return Err(format!(
                            "Collision kernel normalization failed at well={}, grain={}",
                            well.well_name, source_grain
                        ));
                    }

                    let global_source = layout.global_index_of(well_index, source_grain)?;
                    let mut offdiag_prob_sum = 0.0;

                    for (target_grain, weight) in unnormalized_weights {
                        let probability = weight / sum_weights;
                        let global_target = layout.global_index_of(well_index, target_grain)?;

                        if global_target != global_source {
                            operator.add(
                                global_target,
                                global_source,
                                collision_frequency_s_inv * probability,
                            );
                            offdiag_prob_sum += probability;
                        }
                    }

                    let p_same = (1.0 - offdiag_prob_sum).max(0.0);
                    operator.add(
                        global_source,
                        global_source,
                        collision_frequency_s_inv * (p_same - 1.0),
                    );

                    // Column-sum diagnostic: Σ_i L(i,j) should be 0 for collision-only operator.
                    let residual = (p_same + offdiag_prob_sum) - 1.0;
                    diag.max_column_sum_abs =
                        diag.max_column_sum_abs.max((collision_frequency_s_inv * residual).abs());
                }
            }
            CollisionKernelImplementation::Spd => {
                // MarXus SPD kernel: do not include a ΔE=0 weight; instead treat the raw band
                // weights as per-collision move probabilities, and apply a global cap in the
                // reactive region by scaling the effective collision frequency.
                //
                // This matches the single-well `CollisionKernelModel::ExponentialBanded` logic.
                let mut weights_by_source: Vec<Vec<(usize, f64)>> =
                    vec![Vec::new(); local_end_exclusive.saturating_sub(local_start)];
                let mut out_raw: Vec<f64> = vec![0.0; weights_by_source.len()];

                let mut loss_rates: Vec<f64> = vec![0.0; weights_by_source.len()];
                for (idx, source_grain) in (local_start..local_end_exclusive).enumerate() {
                    let mut sum_k = 0.0;
                    for ch in 0..well.channels.len() {
                        sum_k += micro
                            .microcanonical_rate(well_index, ch, source_grain)
                            .max(0.0);
                    }
                    loss_rates[idx] = sum_k;
                }

                for (idx, source_grain) in (local_start..local_end_exclusive).enumerate() {
                    let (target_min, target_max_exclusive) =
                        band_limits(source_grain, local_start, local_end_exclusive, band);

                    let rho_source = micro.density_of_states(well_index, source_grain);
                    if rho_source <= 0.0 {
                        return Err(format!(
                            "Non-positive density of states at well={}, grain={}",
                            well.well_name, source_grain
                        ));
                    }

                    for target_grain in target_min..target_max_exclusive {
                        if target_grain == source_grain {
                            continue;
                        }
                        let d = if target_grain > source_grain {
                            target_grain - source_grain
                        } else {
                            source_grain - target_grain
                        };
                        let step_energy_cm1 = (d as f64) * grain_width;
                        let base = (-step_energy_cm1 / alpha_cm1).exp();

                        let weight = if target_grain > source_grain {
                            let rho_target = micro.density_of_states(well_index, target_grain);
                            if rho_target <= 0.0 {
                                return Err(format!(
                                    "Non-positive density of states at well={}, grain={}",
                                    well.well_name, target_grain
                                ));
                            }
                            base * (rho_target / rho_source)
                                * (-step_energy_cm1 / (boltzmann * temperature)).exp()
                        } else {
                            base
                        };

                        let w = weight.max(0.0);
                        weights_by_source[idx].push((target_grain, w));
                        out_raw[idx] += w;
                    }
                }

                let max_out_reactive = out_raw
                    .iter()
                    .copied()
                    .zip(loss_rates.iter().copied())
                    .filter(|(_, k)| *k > 0.0)
                    .map(|(out, _)| out)
                    .fold(0.0_f64, |a, b| a.max(b));
                let max_out_all = out_raw.iter().copied().fold(0.0_f64, |a, b| a.max(b));
                let out_cap = if max_out_reactive > 0.0 {
                    max_out_reactive
                } else {
                    max_out_all
                };
                let omega_eff = if out_cap > 1.0 {
                    collision_frequency_s_inv / out_cap
                } else {
                    collision_frequency_s_inv
                };

                for (idx, source_grain) in (local_start..local_end_exclusive).enumerate() {
                    let global_source = layout.global_index_of(well_index, source_grain)?;
                    let mut col_sum = -omega_eff * out_raw[idx];
                    for (target_grain, w) in &weights_by_source[idx] {
                        if *w <= 0.0 {
                            continue;
                        }
                        let global_target = layout.global_index_of(well_index, *target_grain)?;
                        operator.add(global_target, global_source, omega_eff * (*w));
                        col_sum += omega_eff * (*w);
                    }

                    // diagonal: -ω_eff * Σ_{i≠j} p_move(j->i)
                    operator.add(global_source, global_source, -omega_eff * out_raw[idx]);
                    diag.max_column_sum_abs = diag.max_column_sum_abs.max(col_sum.abs());
                }
            }
        }

        Ok(())
    }

    /// Adds diagonal reaction sinks:
    ///
    /// For each state (well, grain):
    ///   L(state,state) -= Σ_c k_c(E)
    ///
    /// Outgoing-channel thresholding is applied here: if a channel is outgoing and k < threshold,
    /// it is set to zero in the sink sum.
    fn add_reaction_sinks(
        &self,
        micro: &dyn MicrocanonicalProvider,
        layout: &GlobalLayout,
        operator: &mut DenseMatrix,
    ) -> Result<(), String> {
        let out_thresh = self.settings.outgoing_rate_threshold;
        let skip_internal = self.settings.enforce_interwell_detailed_balance;

        for (well_index, well) in self.wells.iter().enumerate() {
            for local_grain in
                well.lowest_included_grain_index..well.one_past_highest_included_grain_index
            {
                let global_state = layout.global_index_of(well_index, local_grain)?;

                let mut total_loss_rate = 0.0;

                // k(E)=0 for nonreactive grains
                    if local_grain >= well.nonreactive_grain_count {
                        for (channel_index, channel) in well.channels.iter().enumerate() {
                            if skip_internal && channel.connected_well_index.is_some() {
                                continue;
                            }
                            let mut k =
                                micro.microcanonical_rate(well_index, channel_index, local_grain);
                            if channel.connected_well_index.is_none() && k < out_thresh {
                                k = 0.0;
                            }
                            total_loss_rate += k;
                        }
                }

                operator.add(global_state, global_state, -total_loss_rate);
            }
        }

        Ok(())
    }

    /// Adds inter-well coupling terms for internal channels (population transfer).
    ///
    /// For each well w, grain i, and each internal channel to well u:
    ///   let i' = i + offset(w) - offset(u)
    /// if i' is inside [lowest(u), upbound(u)):
    ///   L( (u,i'), (w,i) ) += k_{w->u}(i)
    ///
    /// Internal detailed-balance correction is *not* implemented here (by design).
    /// If you want it, we can add a preprocessing pass that rescales paired rates
    /// using your exact correction logic.
    fn add_interwell_couplings(
        &self,
        micro: &dyn MicrocanonicalProvider,
        layout: &GlobalLayout,
        operator: &mut DenseMatrix,
    ) -> Result<(), String> {
        let internal_thresh = self.settings.internal_rate_threshold;
        let enforce_db = self.settings.enforce_interwell_detailed_balance;
        let temperature = self.settings.temperature_kelvin;
        let boltzmann = self.settings.boltzmann_constant_wavenumber_per_kelvin;

        if !enforce_db {
            for (from_well_index, from_well) in self.wells.iter().enumerate() {
                for from_grain in from_well.lowest_included_grain_index
                    ..from_well.one_past_highest_included_grain_index
                {
                    if from_grain < from_well.nonreactive_grain_count {
                        continue;
                    }

                    let global_from = layout.global_index_of(from_well_index, from_grain)?;

                    for (channel_index, channel) in from_well.channels.iter().enumerate() {
                        let Some(to_well_index) = channel.connected_well_index else {
                            continue;
                        };

                        let k_forward =
                            micro.microcanonical_rate(from_well_index, channel_index, from_grain);
                        if k_forward < internal_thresh {
                            continue;
                        }

                        let to_well = &self.wells[to_well_index];
                        let Some(to_grain) =
                            map_grain_by_aligned_energy(from_well, from_grain, to_well)?
                        else {
                            continue;
                        };

                        let global_to = layout.global_index_of(to_well_index, to_grain)?;
                        operator.add(global_to, global_from, k_forward);
                    }
                }
            }

            return Ok(());
        }

        // Enforced detailed balance:
        // - Require unique internal channels between wells (one per direction).
        // - Use energy-based grain mapping (supports differing ΔE with commensurate offsets).
        // - Replace internal sinks by explicit bidirectional transfers (including diagonal losses).

        let mut unique_to: Vec<std::collections::HashMap<usize, usize>> =
            vec![std::collections::HashMap::new(); self.wells.len()];
        for (w, well) in self.wells.iter().enumerate() {
            for (ch, channel) in well.channels.iter().enumerate() {
                let Some(to) = channel.connected_well_index else {
                    continue;
                };
                if unique_to[w].insert(to, ch).is_some() {
                    return Err(format!(
                        "Multiple internal channels from well {} to well {} are not supported when enforce_interwell_detailed_balance=true.",
                        w, to
                    ));
                }
            }
        }

        let mut links: Vec<(usize, usize, usize, usize)> = Vec::new();
        for from in 0..self.wells.len() {
            for (&to, &ch_from_to) in unique_to[from].iter() {
                if to <= from {
                    continue;
                }
                let Some(&ch_to_from) = unique_to[to].get(&from) else {
                    return Err(format!(
                        "Missing reverse internal channel for wells {} <-> {} with enforce_interwell_detailed_balance=true.",
                        from, to
                    ));
                };
                links.push((from, to, ch_from_to, ch_to_from));
            }
        }

        for (w_i, w_j, ch_i_to_j, ch_j_to_i) in links {
            let well_i = &self.wells[w_i];
            let well_j = &self.wells[w_j];

            for grain_i in well_i.lowest_included_grain_index..well_i.one_past_highest_included_grain_index
            {
                if grain_i < well_i.nonreactive_grain_count {
                    continue;
                }

                let Some(grain_j) = map_grain_by_aligned_energy(well_i, grain_i, well_j)? else {
                    continue;
                };
                if grain_j < well_j.nonreactive_grain_count {
                    continue;
                }

                let k_ij = micro
                    .microcanonical_rate(w_i, ch_i_to_j, grain_i)
                    .max(0.0);
                let k_ji = micro
                    .microcanonical_rate(w_j, ch_j_to_i, grain_j)
                    .max(0.0);

                if k_ij < internal_thresh || k_ji < internal_thresh {
                    continue;
                }

                let rho_i = micro.density_of_states(w_i, grain_i);
                let rho_j = micro.density_of_states(w_j, grain_j);
                if rho_i <= 0.0 || rho_j <= 0.0 {
                    return Err("Non-positive density of states in inter-well detailed balance.".into());
                }

                let e_i = absolute_energy_cm1(well_i, grain_i);
                let e_j = absolute_energy_cm1(well_j, grain_j);
                let w_i_eq = rho_i * (-e_i / (boltzmann * temperature)).exp();
                let w_j_eq = rho_j * (-e_j / (boltzmann * temperature)).exp();
                if w_i_eq <= 0.0 || w_j_eq <= 0.0 || !w_i_eq.is_finite() || !w_j_eq.is_finite() {
                    return Err("Invalid equilibrium weights in inter-well detailed balance.".into());
                }

                let g = (k_ij * k_ji).sqrt();
                let ratio = (w_j_eq / w_i_eq).sqrt();
                let k_ij_corr = g * ratio;
                let k_ji_corr = g / ratio;

                let gi = layout.global_index_of(w_i, grain_i)?;
                let gj = layout.global_index_of(w_j, grain_j)?;

                operator.add(gj, gi, k_ij_corr);
                operator.add(gi, gi, -k_ij_corr);

                operator.add(gi, gj, k_ji_corr);
                operator.add(gj, gj, -k_ji_corr);
            }
        }

        Ok(())
    }

    /// Build the chemical activation source vector (normalized).
    ///
    /// Implemented injection distribution (per your doc’s CA concept):
    ///
    ///   s_raw(w*, i) = ρ(i) * exp( - i ΔE / (k_B T) ) * k_recomb(i)
    ///   s(w*, i)     = s_raw / Σ_i s_raw
    ///
    /// All other wells: s = 0.
    fn build_activation_source(
        &self,
        micro: &dyn MicrocanonicalProvider,
        layout: &GlobalLayout,
        activation: &ChemicalActivationDefinition,
    ) -> Result<Vec<f64>, String> {
        let mut source = vec![0.0; layout.total_state_count];

        let w_star = activation.activated_well_index;
        if w_star >= self.wells.len() {
            return Err("Activated well index out of range".into());
        }
        let well = &self.wells[w_star];
        let c_star = activation.recombination_channel_index;
        if c_star >= well.channels.len() {
            return Err("Recombination channel index out of range".into());
        }

        let temperature = self.settings.temperature_kelvin;
        let boltzmann = self.settings.boltzmann_constant_wavenumber_per_kelvin;

        let mut normalization_sum = 0.0;

        for local_grain in
            well.lowest_included_grain_index..well.one_past_highest_included_grain_index
        {
            let rho = micro.density_of_states(w_star, local_grain);
            if rho <= 0.0 {
                continue;
            }

            let k_recomb = if local_grain >= well.nonreactive_grain_count {
                micro
                    .microcanonical_rate(w_star, c_star, local_grain)
                    .max(0.0)
            } else {
                0.0
            };

            let energy_cm1 = absolute_energy_cm1(well, local_grain);
            let boltzmann_weight = (-energy_cm1 / (boltzmann * temperature)).exp();

            let raw = rho * boltzmann_weight * k_recomb;

            let global_idx = layout.global_index_of(w_star, local_grain)?;
            source[global_idx] = raw;
            normalization_sum += raw;
        }

        if normalization_sum <= 0.0 {
            return Err(
                "Chemical activation source normalization is zero; check k(E) and rho(E) inputs."
                    .into(),
            );
        }

        for x in &mut source {
            *x /= normalization_sum;
        }

        Ok(source)
    }

    /// Build similarity transform scale W (diagonal).
    ///
    /// We use an equilibrium-like weight:
    ///
    ///   W(state) = sqrt( ρ(i) * exp( -E_abs / (k_B T) ) )
    ///
    /// where E_abs uses the alignment offset:
    ///
    ///   E_abs = (i + offset_w) * ΔE
    fn build_similarity_scale(
        &self,
        micro: &dyn MicrocanonicalProvider,
        layout: &GlobalLayout,
    ) -> Result<DiagonalScale, String> {
        let mut diag = vec![0.0; layout.total_state_count];
        let temperature = self.settings.temperature_kelvin;
        let boltzmann = self.settings.boltzmann_constant_wavenumber_per_kelvin;

        for (well_index, well) in self.wells.iter().enumerate() {
            let grain_width = well.energy_grain_width_cm1;
            for local_grain in
                well.lowest_included_grain_index..well.one_past_highest_included_grain_index
            {
                let rho = micro.density_of_states(well_index, local_grain);
                if rho <= 0.0 {
                    return Err(format!(
                        "Non-positive density of states in similarity scale at well={}, grain={}",
                        well.well_name, local_grain
                    ));
                }

                let abs_grain = (local_grain as isize) + well.alignment_offset_in_grains;
                let abs_energy_cm1 = (abs_grain as f64) * grain_width;
                let abs_boltzmann = (-abs_energy_cm1 / (boltzmann * temperature)).exp();

                let global_idx = layout.global_index_of(well_index, local_grain)?;
                diag[global_idx] = (rho * abs_boltzmann).sqrt();
            }
        }

        Ok(DiagonalScale { diagonal: diag })
    }

    /// Compute macroscopic rates by population-weighted microcanonical averaging:
    ///
    /// Let population vector be p(state). Normalize:
    ///   P_total = Σ_state p(state)
    ///
    /// For each well w and channel c:
    ///   k_w,c(T,p) = (1/P_total) Σ_{i in included grains} p(w,i) * k_w,c(i)
    ///
    /// Also compute total outgoing:
    ///   k_out_total = Σ_{outgoing channels} k_w,c + (optional truncation leakage, not included in this simplified kernel)
    fn compute_macroscopic_rates(
        &self,
        micro: &dyn MicrocanonicalProvider,
        layout: &GlobalLayout,
        population: &[f64],
    ) -> Result<SolutionResults, String> {
        if population.len() != layout.total_state_count {
            return Err("Population vector length mismatch".into());
        }

        let total_population: f64 = population.iter().sum();
        if total_population <= 0.0 {
            return Err("Total population is non-positive after solve".into());
        }

        let mut per_well_per_channel_rates: Vec<Vec<f64>> = Vec::new();
        let mut total_outgoing = 0.0;

        for (well_index, well) in self.wells.iter().enumerate() {
            let mut channel_rates = vec![0.0; well.channels.len()];

            for local_grain in
                well.lowest_included_grain_index..well.one_past_highest_included_grain_index
            {
                let global_idx = layout.global_index_of(well_index, local_grain)?;
                let weight = population[global_idx] / total_population;

                if local_grain < well.nonreactive_grain_count {
                    continue;
                }

                for (channel_index, channel) in well.channels.iter().enumerate() {
                    let k = micro
                        .microcanonical_rate(well_index, channel_index, local_grain)
                        .max(0.0);
                    channel_rates[channel_index] += weight * k;

                    if channel.connected_well_index.is_none() {
                        total_outgoing += weight * k;
                    }
                }
            }

            per_well_per_channel_rates.push(channel_rates);
        }

        Ok(SolutionResults {
            steady_state_population: population.to_vec(),
            per_well_per_channel_rates,
            total_outgoing_rate_constant: total_outgoing,
        })
    }
}
