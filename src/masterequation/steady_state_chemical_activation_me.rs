use super::collisional_relaxation::{
    band_limits, compute_alpha_cm1, compute_collision_frequency_s_inv,
};
use super::state_index::GlobalLayout;
use super::reaction_network::{
    ChemicalActivationDefinition, CollisionKernelImplementation, MasterEquationSettings,
    MicrocanonicalProvider, SolutionResults, WellDefinition,
};
use crate::numeric::linear_algebra::{cholesky_solve_spd, DenseMatrix, DiagonalScale};

/// Main engine: assembles and solves the steady-state chemical activation master equation.
pub struct MasterEquationEngine {
    wells: Vec<WellDefinition>,
    settings: MasterEquationSettings,
}

impl MasterEquationEngine {
    pub fn new(wells: Vec<WellDefinition>, settings: MasterEquationSettings) -> Self {
        Self { wells, settings }
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

        // 4) Solve (-A) q = b using Cholesky (requires SPD)
        let negative_transformed = transformed_operator.scaled(-1.0);

        let transformed_solution = cholesky_solve_spd(&negative_transformed, &transformed_rhs)
            .map_err(|e| format!("Cholesky solve failed: {e}"))?;

        // 5) Recover p = W q
        let steady_state_population = similarity_scale.apply_to_vector(&transformed_solution);

        // 6) Compute macroscopic (averaged) rate constants
        let macro_rates =
            self.compute_macroscopic_rates(micro, &global_layout, &steady_state_population)?;

        Ok(macro_rates)
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
                    for (target_grain, w) in &weights_by_source[idx] {
                        if *w <= 0.0 {
                            continue;
                        }
                        let global_target = layout.global_index_of(well_index, *target_grain)?;
                        operator.add(global_target, global_source, omega_eff * (*w));
                    }

                    // diagonal: -ω_eff * Σ_{i≠j} p_move(j->i)
                    operator.add(global_source, global_source, -omega_eff * out_raw[idx]);
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

        for (well_index, well) in self.wells.iter().enumerate() {
            for local_grain in
                well.lowest_included_grain_index..well.one_past_highest_included_grain_index
            {
                let global_state = layout.global_index_of(well_index, local_grain)?;

                let mut total_loss_rate = 0.0;

                // k(E)=0 for nonreactive grains
                    if local_grain >= well.nonreactive_grain_count {
                        for (channel_index, channel) in well.channels.iter().enumerate() {
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

                    let aligned_to_grain_isize = (from_grain as isize)
                        + from_well.alignment_offset_in_grains
                        - to_well.alignment_offset_in_grains;

                    if aligned_to_grain_isize < 0 {
                        continue;
                    }
                    let to_grain = aligned_to_grain_isize as usize;

                    if to_grain < to_well.lowest_included_grain_index
                        || to_grain >= to_well.one_past_highest_included_grain_index
                    {
                        continue;
                    }

                    let global_to = layout.global_index_of(to_well_index, to_grain)?;
                    operator.add(global_to, global_from, k_forward);
                }
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
        let grain_width = well.energy_grain_width_cm1;

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

            let energy_cm1 = (local_grain as f64) * grain_width;
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
