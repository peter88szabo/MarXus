use super::chemical_activation_source::build_normalized_source_distribution;
use super::tridiagonal_solvers::{solve_spd_symmetric_tridiagonal, solve_tridiagonal_general};
use super::energy_grained_me::{
    CollisionKernelModel, EnergyGrid, MicrocanonicalGridData, OlzmannMasterEquationSettings,
    SourceConstructionChoice, SteadyStateChemicalActivationResult, ThermalDissociationResult,
};
use crate::numeric::linear_algebra::{cholesky_solve_spd, DenseMatrix};

/// Main solver: single intermediate, stepladder collisions, steady-state CA.
pub struct OlzmannStepladderMasterEquationSolver {
    pub energy_grid: EnergyGrid,
    pub settings: OlzmannMasterEquationSettings,
}

impl OlzmannStepladderMasterEquationSolver {
    /// Solve J N_ss = R_form * F using similarity-transformed symmetric tridiagonal system.
    pub fn solve_steady_state_with_source(
        &self,
        micro: &dyn MicrocanonicalGridData,
        formation_flux: f64,
        normalized_source_distribution: &[f64],
    ) -> Result<SteadyStateChemicalActivationResult, String> {
        self.validate_source_distribution(normalized_source_distribution)?;

        // MESS-normalization kernel is generally non-symmetric; solve the untransformed
        // banded system directly with an iterative solver.
        if let CollisionKernelModel::ExponentialBandedMess {
            mean_downstep_wavenumber,
            exponent_cutoff,
        } = &self.settings.collision_kernel_model
        {
            let (diag, upper, lower) = self.build_untransformed_banded_exponential_operator_mess(
                micro,
                *mean_downstep_wavenumber,
                *exponent_cutoff,
            )?;
            let rhs = normalized_source_distribution
                .iter()
                .map(|f_i| formation_flux * f_i)
                .collect::<Vec<f64>>();
            let populations = solve_banded_bicgstab_left_jacobi(&diag, &upper, &lower, &rhs, 1e-10, 4000)
                .or_else(|_| solve_banded_gauss_seidel(&diag, &upper, &lower, &rhs, 1e-12, 40000))?;
            return self.postprocess_steady_state_populations(micro, populations);
        }

        let similarity_weights = self.build_similarity_weights(micro)?;
        let rhs_transformed = normalized_source_distribution
            .iter()
            .zip(similarity_weights.iter())
            .map(|(f_i, w_i)| (formation_flux * f_i) / w_i)
            .collect::<Vec<f64>>();

        let x_transformed = match &self.settings.collision_kernel_model {
            CollisionKernelModel::Stepladder => {
                let (diag_a, offdiag_a) =
                    self.build_transformed_symmetric_tridiagonal_operator(micro, &similarity_weights)?;
                match solve_spd_symmetric_tridiagonal(&diag_a, &offdiag_a, &rhs_transformed) {
                    Ok(x) => x,
                    Err(_) => solve_tridiagonal_general(
                        &offdiag_a,
                        &diag_a,
                        &offdiag_a,
                        &rhs_transformed,
                    )?,
                }
            }
            CollisionKernelModel::ExponentialBanded {
                mean_downstep_wavenumber,
                exponent_cutoff,
            } => {
                let (diag_a, bands) = self
                    .build_transformed_symmetric_banded_exponential_operator(
                        micro,
                        &similarity_weights,
                        *mean_downstep_wavenumber,
                        *exponent_cutoff,
                    )?;
                solve_spd_symmetric_banded_cholesky(&diag_a, &bands, &rhs_transformed)
                    .or_else(|e| {
                        // Fallback: expand to dense and try dense Cholesky.
                        // This is slower but can rescue borderline cases without stalling.
                        solve_spd_symmetric_banded_dense_fallback(&diag_a, &bands, &rhs_transformed)
                            .map_err(|e2| format!("{e}; dense fallback: {e2}"))
                    })?
            }
            CollisionKernelModel::ExponentialBandedMess {
                ..
            } => unreachable!("ExponentialBandedMess is solved in untransformed space above."),
        };

        let steady_state_populations = x_transformed
            .iter()
            .zip(similarity_weights.iter())
            .map(|(x_i, w_i)| x_i * w_i)
            .collect::<Vec<f64>>();

        self.postprocess_steady_state_populations(micro, steady_state_populations)
    }

    fn postprocess_steady_state_populations(
        &self,
        micro: &dyn MicrocanonicalGridData,
        steady_state_populations: Vec<f64>,
    ) -> Result<SteadyStateChemicalActivationResult, String> {
        let total_population: f64 = steady_state_populations.iter().sum();
        if total_population <= 0.0 {
            return Err(
                "Steady-state population sum is non-positive; check inputs / SPD assumptions."
                    .into(),
            );
        }

        let normalized_distribution = steady_state_populations
            .iter()
            .map(|n_i| n_i / total_population)
            .collect::<Vec<f64>>();

        let n_channels = micro.number_of_unimolecular_channels();
        let mut channel_k_ca = vec![0.0; n_channels];
        let mut total_unimolecular_loss = 0.0;

        for (i, weight) in normalized_distribution.iter().copied().enumerate() {
            let mut sum_loss_here = 0.0;

            for (ch, k_ca) in channel_k_ca.iter_mut().enumerate().take(n_channels) {
                let k = micro.microcanonical_rate_for_channel(ch, i)
                    .max(0.0);
                *k_ca += weight * k;
                sum_loss_here += k;
            }

            total_unimolecular_loss += weight * sum_loss_here;
        }

        Ok(SteadyStateChemicalActivationResult {
            steady_state_bin_populations: steady_state_populations,
            normalized_steady_state_distribution: normalized_distribution,
            chemically_activated_rate_constants: channel_k_ca,
            total_unimolecular_loss_rate_constant: total_unimolecular_loss,
        })
    }

    /// Build normalized F from a source-choice variant and solve.
    pub fn solve_steady_state_with_source_choice(
        &self,
        micro: &dyn MicrocanonicalGridData,
        formation_flux: f64,
        source_choice: &SourceConstructionChoice,
    ) -> Result<SteadyStateChemicalActivationResult, String> {
        let f =
            build_normalized_source_distribution(&self.energy_grid, &self.settings, source_choice)?;
        self.solve_steady_state_with_source(micro, formation_flux, &f)
    }

    /// Solve steady state with an equilibrium-like source:
    ///   F_i ∝ ρ(E_i) exp(-E_i / k_B T)
    ///
    /// This provides a no-reactants steady-state mode that still uses
    /// the full collision/reaction operator J.
    pub fn solve_steady_state_with_equilibrium_source(
        &self,
        micro: &dyn MicrocanonicalGridData,
        formation_flux: f64,
    ) -> Result<SteadyStateChemicalActivationResult, String> {
        let n = self.energy_grid.number_of_bins;
        let t = self.settings.temperature_kelvin;
        let k_b = self.settings.boltzmann_constant_wavenumber_per_kelvin;

        let mut source = vec![0.0; n];
        for (i, value) in source.iter_mut().enumerate() {
            let rho = micro.density_of_states(i);
            if rho <= 0.0 {
                return Err(format!("Non-positive density of states ρ at bin {}.", i));
            }
            let e = self.energy_grid.energy_at_bin(i);
            *value = rho * (-e / (k_b * t)).exp();
        }
        let z: f64 = source.iter().sum();
        if z <= 0.0 || !z.is_finite() {
            return Err("Equilibrium source normalization failed.".into());
        }
        for x in &mut source {
            *x /= z;
        }

        self.solve_steady_state_with_source(micro, formation_flux, &source)
    }

    /// Solve thermal unimolecular dissociation without chemical-activation source feeding.
    ///
    /// Uses canonical population weights:
    ///   n_i ∝ ρ(E_i) exp(-E_i / k_B T)
    /// then averages microcanonical rates over that distribution.
    pub fn solve_thermal_dissociation(
        &self,
        micro: &dyn MicrocanonicalGridData,
    ) -> Result<ThermalDissociationResult, String> {
        let n = self.energy_grid.number_of_bins;
        let t = self.settings.temperature_kelvin;
        let k_b = self.settings.boltzmann_constant_wavenumber_per_kelvin;

        let mut thermal = vec![0.0; n];
        for (i, value) in thermal.iter_mut().enumerate() {
            let rho = micro.density_of_states(i);
            if rho <= 0.0 {
                return Err(format!("Non-positive density of states ρ at bin {}.", i));
            }
            let e = self.energy_grid.energy_at_bin(i);
            *value = rho * (-e / (k_b * t)).exp();
        }

        let z: f64 = thermal.iter().sum();
        if z <= 0.0 || !z.is_finite() {
            return Err("Thermal normalization failed (non-positive partition sum).".into());
        }
        for x in &mut thermal {
            *x /= z;
        }

        let n_channels = micro.number_of_unimolecular_channels();
        let mut channel_rates = vec![0.0; n_channels];
        let mut total = 0.0;
        for (i, w) in thermal.iter().copied().enumerate() {
            let mut sum_loss_here = 0.0;
            for (ch, k_th) in channel_rates.iter_mut().enumerate().take(n_channels) {
                let k = micro.microcanonical_rate_for_channel(ch, i)
                    .max(0.0);
                *k_th += w * k;
                sum_loss_here += k;
            }
            total += w * sum_loss_here;
        }

        Ok(ThermalDissociationResult {
            normalized_thermal_distribution: thermal,
            thermal_rate_constants: channel_rates,
            total_unimolecular_loss_rate_constant: total,
        })
    }

    fn validate_source_distribution(&self, source: &[f64]) -> Result<(), String> {
        if source.len() != self.energy_grid.number_of_bins {
            return Err("Source distribution length does not match energy grid.".into());
        }
        if source.iter().any(|x| *x < 0.0) {
            return Err("Source distribution has negative entries.".into());
        }
        let sum: f64 = source.iter().sum();
        if (sum - 1.0).abs() > 1e-10 {
            return Err(format!(
                "Source distribution must be normalized (sum=1). Got sum={}",
                sum
            ));
        }

        Ok(())
    }

    fn build_similarity_weights(
        &self,
        micro: &dyn MicrocanonicalGridData,
    ) -> Result<Vec<f64>, String> {
        let n = self.energy_grid.number_of_bins;
        let t = self.settings.temperature_kelvin;
        let k_b = self.settings.boltzmann_constant_wavenumber_per_kelvin;

        let mut weights = vec![0.0; n];
        for (i, value) in weights.iter_mut().enumerate() {
            let rho = micro.density_of_states(i);
            if rho <= 0.0 {
                return Err(format!("Non-positive density of states ρ at bin {}.", i));
            }

            let energy = self.energy_grid.energy_at_bin(i);
            let boltz = (-energy / (k_b * t)).exp();
            *value = (rho * boltz).sqrt();
        }

        Ok(weights)
    }


    fn build_transformed_symmetric_tridiagonal_operator(
        &self,
        micro: &dyn MicrocanonicalGridData,
        similarity_weights: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), String> {
        let n = self.energy_grid.number_of_bins;
        let omega = self.settings.collision_frequency_per_second;
        let k_capture = self.settings.pseudo_first_order_capture_loss_per_second;

        let base_down = self.settings.stepladder_base_downward_probability;
        if !(0.0..=1.0).contains(&base_down) {
            return Err("stepladder_base_downward_probability must be in [0,1].".into());
        }

        let t = self.settings.temperature_kelvin;
        let k_b = self.settings.boltzmann_constant_wavenumber_per_kelvin;

        let loss_rates = self.build_loss_rates(micro);

        let mut diag_a = vec![0.0; n];
        let mut offdiag_a = vec![0.0; n.saturating_sub(1)];

        let delta_e = self.energy_grid.bin_width_wavenumber;
        let boltz_step = (-delta_e / (k_b * t)).exp();

        for source in 0..n {
            let mut p_down = if source > 0 { base_down } else { 0.0 };

            let mut p_up = if source + 1 < n {
                let rho_here = micro.density_of_states(source);
                let rho_up = micro.density_of_states(source + 1);
                if rho_here <= 0.0 || rho_up <= 0.0 {
                    return Err("Non-positive density of states in stepladder ratio.".into());
                }

                let ratio_r = (rho_up / rho_here) * boltz_step;
                base_down * ratio_r
            } else {
                0.0
            };

            let sum = p_down + p_up;
            if sum > 1.0 {
                p_down /= sum;
                p_up /= sum;
            }

            let j_diag = omega * (p_down + p_up) + loss_rates[source] + k_capture;
            diag_a[source] = j_diag;

            // Build the (similarity-transformed) symmetric off-diagonal A_{i+1,i}.
            //
            // We assemble the raw operator J in "target,row; source,col" form:
            //   J_{i+1,i} = -ω p_up(i)
            //
            // Similarity transform A = W J W^{-1} gives:
            //   A_{i+1,i} = J_{i+1,i} * (W_{i+1} / W_i)
            //
            // With detailed balance, this is symmetric, so we only need one direction.
            if source + 1 < n && p_up > 0.0 {
                let upper = source + 1;
                let j_off = -omega * p_up;
                offdiag_a[source] = j_off * (similarity_weights[upper] / similarity_weights[source]);
            }
        }

        if diag_a.iter().any(|d| *d <= 0.0) {
            return Err(
                "Transformed operator has non-positive diagonal entries; SPD likely violated."
                    .into(),
            );
        }

        Ok((diag_a, offdiag_a))
    }

    fn build_transformed_symmetric_banded_exponential_operator(
        &self,
        micro: &dyn MicrocanonicalGridData,
        similarity_weights: &[f64],
        mean_downstep_wavenumber: f64,
        exponent_cutoff: f64,
    ) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
        let n = self.energy_grid.number_of_bins;
        let omega = self.settings.collision_frequency_per_second;
        let k_capture = self.settings.pseudo_first_order_capture_loss_per_second;
        let t = self.settings.temperature_kelvin;
        let k_b = self.settings.boltzmann_constant_wavenumber_per_kelvin;
        let delta_e = self.energy_grid.bin_width_wavenumber;
        if mean_downstep_wavenumber <= 0.0 {
            return Err(
                "mean_downstep_wavenumber must be > 0 for exponential banded kernel.".into(),
            );
        }
        if exponent_cutoff <= 0.0 {
            return Err("exponent_cutoff must be > 0 for exponential banded kernel.".into());
        }
        let half_bandwidth = ((exponent_cutoff * mean_downstep_wavenumber / delta_e).ceil()
            as usize)
            .max(1)
            .min(64)
            .min(n.saturating_sub(1));
        let mut q_up = (0..half_bandwidth)
            .map(|d| vec![0.0; n.saturating_sub(d + 1)])
            .collect::<Vec<Vec<f64>>>();
        let mut q_down = (0..half_bandwidth)
            .map(|d| vec![0.0; n.saturating_sub(d + 1)])
            .collect::<Vec<Vec<f64>>>();
        let mut out_raw = vec![0.0; n];

        for i in 0..n {
            let rho_i = micro.density_of_states(i);
            if rho_i <= 0.0 {
                return Err(format!("Non-positive density of states ρ at bin {}.", i));
            }
            let max_target = (i + half_bandwidth).min(n - 1);
            for j in (i + 1)..=max_target {
                let offset = j - i;
                let d_e = (offset as f64) * delta_e;
                let qd = (-d_e / mean_downstep_wavenumber).exp();
                let rho_j = micro.density_of_states(j);
                if rho_j <= 0.0 {
                    return Err(format!("Non-positive density of states ρ at bin {}.", j));
                }
                let qu = qd * (rho_j / rho_i) * (-d_e / (k_b * t)).exp();
                q_up[offset - 1][i] = qu;
                q_down[offset - 1][i] = qd;
                out_raw[i] += qu;
                out_raw[j] += qd;
            }
        }

        let loss_rates = self.build_loss_rates(micro);

        // Scale the raw exponential weights with a *single global factor* so that, in the
        // chemically relevant region (where k(E) > 0), the total leave-probability per
        // collision is <= 1.
        //
        // Using the global maximum across *all* grains can over-dampen collisions due to
        // extreme density-of-states ratios far below the reactive threshold. This heuristic
        // chooses the cap from the reactive region to better match typical master-equation
        // practice for falloff calculations.
        let max_out_reactive = out_raw
            .iter()
            .copied()
            .zip(loss_rates.iter().copied())
            .filter(|(_, k)| *k > 0.0)
            .map(|(out, _)| out)
            .fold(0.0_f64, |a, b| a.max(b));
        let max_out_all = out_raw
            .iter()
            .copied()
            .fold(0.0_f64, |a, b| a.max(b));
        let out_cap = if max_out_reactive > 0.0 {
            max_out_reactive
        } else {
            max_out_all
        };
        let omega_eff = if out_cap > 1.0 { omega / out_cap } else { omega };

        let mut diag = vec![0.0; n];
        let mut bands = (0..half_bandwidth)
            .map(|d| vec![0.0; n.saturating_sub(d + 1)])
            .collect::<Vec<Vec<f64>>>();

        for i in 0..n {
            diag[i] += omega_eff * out_raw[i] + loss_rates[i] + k_capture;
        }

        for d in 1..=half_bandwidth {
            for i in 0..(n - d) {
                let j = i + d;
                let p_up = q_up[d - 1][i];
                let p_down = q_down[d - 1][i];
                if p_up <= 0.0 && p_down <= 0.0 {
                    continue;
                }
                let a_ij_from_down =
                    -omega_eff * p_down * (similarity_weights[j] / similarity_weights[i]);
                let a_ij_from_up =
                    -omega_eff * p_up * (similarity_weights[i] / similarity_weights[j]);
                bands[d - 1][i] = 0.5 * (a_ij_from_down + a_ij_from_up);
            }
        }

        if diag.iter().any(|x| *x <= 0.0 || !x.is_finite()) {
            return Err("Banded exponential operator has non-positive/invalid diagonal.".into());
        }

        Ok((diag, bands))
    }

    fn build_untransformed_banded_exponential_operator_mess(
        &self,
        micro: &dyn MicrocanonicalGridData,
        mean_downstep_wavenumber: f64,
        exponent_cutoff: f64,
    ) -> Result<(Vec<f64>, Vec<Vec<f64>>, Vec<Vec<f64>>), String> {
        // This follows the normalization logic in:
        //   utils/MESS/src/libmess/new_mess.cc (Well::set, ENERGY TRANSFER KERNEL)
        //
        // Note: MarXus uses a single effective bath; MESS supports mixtures via kernel_fraction.
        // Here we implement the single-bath case (kernel_fraction = 1).
        let n = self.energy_grid.number_of_bins;
        let omega = self.settings.collision_frequency_per_second;
        let k_capture = self.settings.pseudo_first_order_capture_loss_per_second;
        let t = self.settings.temperature_kelvin;
        let k_b = self.settings.boltzmann_constant_wavenumber_per_kelvin;
        let delta_e = self.energy_grid.bin_width_wavenumber;

        if mean_downstep_wavenumber <= 0.0 {
            return Err(
                "mean_downstep_wavenumber must be > 0 for exponential banded kernel.".into(),
            );
        }
        if exponent_cutoff <= 0.0 {
            return Err("exponent_cutoff must be > 0 for exponential banded kernel.".into());
        }
        if delta_e <= 0.0 {
            return Err("energy grid bin width must be > 0.".into());
        }

        let half_bandwidth = ((exponent_cutoff * mean_downstep_wavenumber / delta_e).ceil()
            as usize)
            .max(1)
            .min(64)
            .min(n.saturating_sub(1));

        let loss_rates = self.build_loss_rates(micro);

        // energy_transfer[d] = exp(-d ΔE / <ΔE_down>), with energy_transfer[0]=1.
        let mut energy_transfer = vec![0.0; half_bandwidth + 1];
        energy_transfer[0] = 1.0;
        for d in 1..=half_bandwidth {
            let d_e = (d as f64) * delta_e;
            energy_transfer[d] = (-d_e / mean_downstep_wavenumber).exp();
        }

        let mut upper = (0..half_bandwidth)
            .map(|d| vec![0.0; n.saturating_sub(d + 1)])
            .collect::<Vec<Vec<f64>>>();
        let mut lower = (0..half_bandwidth)
            .map(|d| vec![0.0; n.saturating_sub(d + 1)])
            .collect::<Vec<Vec<f64>>>();
        // MESS-style banded exponential kernel (default block) assembled from the MESS source:
        //   utils/MESS/src/libmess/mess_dd.cc and utils/MESS/src/libmess/new_mess.cc
        //
        // MESS defines an internal kernel matrix K with:
        // - off-diagonals negative (generator form)
        // - row sums equal to zero (diagonal built by subtracting the off-diagonals)
        //
        // Our steady-state balance system uses the transpose convention "target,row; source,col"
        // (column sums are zero for a pure collision generator). Therefore we build:
        //   J_coll = ω * K^T
        //
        // Additionally, MESS' energy bins are indexed from high-to-low energy in parts of the
        // code path (see `ener -= energy_step()` loops), whereas MarXus indexes low-to-high.
        // To match MESS' kernel algebra, we build K in a reversed index `m = n-1-i` and map
        // back to MarXus indices.
        let mut diag = loss_rates
            .iter()
            .map(|k_loss| k_loss + k_capture)
            .collect::<Vec<f64>>();

        let map_mess_to_marxus = |m: usize| -> usize { n - 1 - m };

        let mut diag_k = vec![0.0_f64; n]; // diagonal of K in MESS indexing (row-sum form)

        for e_m in 0..n {
            let e_i = map_mess_to_marxus(e_m);
            let rho_e = micro.density_of_states(e_i);
            if rho_e <= 0.0 {
                return Err(format!("Non-positive density of states ρ at bin {}.", e_i));
            }

            // nfac = energy_transfer[0] + Σ (K_raw(e_m, j_m) + K_raw(j_m, e_m))
            // with K_raw(e_m, j_m)=energy_transfer[d] and
            //      K_raw(j_m, e_m)=energy_transfer[d] * ρ(e)/ρ(j) / exp(d ΔE / kT)
            let mut nfac = energy_transfer[0];
            for d in 1..=half_bandwidth {
                let j_m = e_m + d;
                if j_m >= n {
                    break;
                }
                let j_i = map_mess_to_marxus(j_m);
                let rho_j = micro.density_of_states(j_i);
                if rho_j <= 0.0 {
                    return Err(format!("Non-positive density of states ρ at bin {}.", j_i));
                }

                let thermal = ((d as f64) * delta_e / (k_b * t)).exp();
                let w = energy_transfer[d];
                let w_back = w * rho_e / rho_j / thermal.max(1e-300);
                nfac += w + w_back;
            }

            if !(nfac > 0.0) || !nfac.is_finite() {
                return Err(format!(
                    "Invalid MESS-like kernel normalization at bin {}.",
                    e_i
                ));
            }

            // Fill banded entries for J_coll = ω * K^T.
            for d in 1..=half_bandwidth {
                let j_m = e_m + d;
                if j_m >= n {
                    break;
                }
                let j_i = map_mess_to_marxus(j_m);
                let rho_j = micro.density_of_states(j_i);
                if rho_j <= 0.0 {
                    return Err(format!("Non-positive density of states ρ at bin {}.", j_i));
                }

                let thermal = ((d as f64) * delta_e / (k_b * t)).exp();
                let w = energy_transfer[d];

                // K(e_m, j_m) and K(j_m, e_m) after division by -nfac (negative values).
                let k_e_j = -w / nfac;
                let k_j_e = -(w * rho_e / rho_j / thermal.max(1e-300)) / nfac;

                // Diagonal of K (row-sum form): K(r,r) = -Σ_{c≠r} K(r,c).
                diag_k[e_m] -= k_e_j;
                diag_k[j_m] -= k_j_e;

                // Map into MarXus indices with transpose: J_coll = ω * K^T.
                //
                // J_coll(row=j, col=e) = ω * K(e, j)  and  J_coll(row=e, col=j) = ω * K(j, e).
                let row1 = j_i;
                let col1 = e_i;
                let val1 = omega * k_e_j;

                let row2 = e_i;
                let col2 = j_i;
                let val2 = omega * k_j_e;

                // Store into banded representation.
                if row1 < col1 {
                    let off = col1 - row1;
                    if off <= half_bandwidth {
                        upper[off - 1][row1] = val1;
                    }
                } else {
                    let off = row1 - col1;
                    if off <= half_bandwidth {
                        lower[off - 1][col1] = val1;
                    }
                }

                if row2 < col2 {
                    let off = col2 - row2;
                    if off <= half_bandwidth {
                        upper[off - 1][row2] = val2;
                    }
                } else {
                    let off = row2 - col2;
                    if off <= half_bandwidth {
                        lower[off - 1][col2] = val2;
                    }
                }
            }
        }

        // Add collision diagonal (same under transpose), mapped back to MarXus indices.
        for e_m in 0..n {
            let e_i = map_mess_to_marxus(e_m);
            diag[e_i] += omega * diag_k[e_m];
        }

        if diag.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err("Invalid diagonal in MESS-like untransformed operator.".into());
        }

        Ok((diag, upper, lower))
    }

    fn build_loss_rates(&self, micro: &dyn MicrocanonicalGridData) -> Vec<f64> {
        let n = self.energy_grid.number_of_bins;
        let n_channels = micro.number_of_unimolecular_channels();
        let mut loss_rates = vec![0.0; n];
        for (i, value) in loss_rates.iter_mut().enumerate() {
            let mut sum_k = 0.0;
            for ch in 0..n_channels {
                sum_k += micro.microcanonical_rate_for_channel(ch, i)
                    .max(0.0);
            }
            *value = sum_k;
        }
        loss_rates
    }
}

fn banded_matvec_general(
    diag: &[f64],
    upper: &[Vec<f64>],
    lower: &[Vec<f64>],
    x: &[f64],
) -> Vec<f64> {
    let n = diag.len();
    let mut y = vec![0.0; n];
    for i in 0..n {
        y[i] = diag[i] * x[i];
    }
    let bw = upper.len().min(lower.len());
    for d in 1..=bw {
        for i in 0..n.saturating_sub(d) {
            let j = i + d;
            y[i] += upper[d - 1][i] * x[j];
            y[j] += lower[d - 1][i] * x[i];
        }
    }
    y
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn l2_norm(v: &[f64]) -> f64 {
    v.iter().copied().map(|x| x * x).sum::<f64>().sqrt()
}

fn solve_banded_bicgstab_left_jacobi(
    diag: &[f64],
    upper: &[Vec<f64>],
    lower: &[Vec<f64>],
    rhs: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<Vec<f64>, String> {
    let n = diag.len();
    if rhs.len() != n {
        return Err("RHS length mismatch in BiCGSTAB solve.".into());
    }
    if n == 0 {
        return Err("Empty system in BiCGSTAB solve.".into());
    }
    if tol <= 0.0 || !tol.is_finite() {
        return Err("Invalid tolerance in BiCGSTAB solve.".into());
    }
    if diag.iter().any(|d| !d.is_finite() || *d == 0.0) {
        return Err("Invalid diagonal in BiCGSTAB solve.".into());
    }

    // Left Jacobi preconditioner: solve A' x = b', A' = M^{-1}A, b' = M^{-1}b with M=diag(A).
    let mut b_prime = vec![0.0; n];
    for i in 0..n {
        b_prime[i] = rhs[i] / diag[i];
    }
    let b_norm = l2_norm(&b_prime).max(1e-300);

    let mut x = vec![0.0; n];
    let mut r = b_prime.clone(); // x=0 => r=b'
    let r0 = r.clone();

    let mut p = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut s = vec![0.0; n];
    let mut t = vec![0.0; n];

    let mut rho_old = 1.0_f64;
    let mut alpha = 1.0_f64;
    let mut omega = 1.0_f64;

    if l2_norm(&r) / b_norm <= tol {
        return Ok(x);
    }

    for _iter in 0..max_iter {
        let rho_new = dot(&r0, &r);
        if rho_new.abs() < 1e-300 {
            return Err("BiCGSTAB breakdown (rho ~ 0).".into());
        }

        let beta = (rho_new / rho_old) * (alpha / omega);
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        // v = A' p = M^{-1}(A p)
        let ap = banded_matvec_general(diag, upper, lower, &p);
        for i in 0..n {
            v[i] = ap[i] / diag[i];
        }

        let denom = dot(&r0, &v);
        if denom.abs() < 1e-300 {
            return Err("BiCGSTAB breakdown (r0·v ~ 0).".into());
        }
        alpha = rho_new / denom;

        for i in 0..n {
            s[i] = r[i] - alpha * v[i];
        }
        if l2_norm(&s) / b_norm <= tol {
            for i in 0..n {
                x[i] += alpha * p[i];
            }
            return Ok(x);
        }

        // t = A' s
        let as_vec = banded_matvec_general(diag, upper, lower, &s);
        for i in 0..n {
            t[i] = as_vec[i] / diag[i];
        }

        let t_dot_t = dot(&t, &t);
        if t_dot_t.abs() < 1e-300 {
            return Err("BiCGSTAB breakdown (t·t ~ 0).".into());
        }
        omega = dot(&t, &s) / t_dot_t;
        if omega.abs() < 1e-300 {
            return Err("BiCGSTAB breakdown (omega ~ 0).".into());
        }

        for i in 0..n {
            x[i] += alpha * p[i] + omega * s[i];
        }
        for i in 0..n {
            r[i] = s[i] - omega * t[i];
        }

        if l2_norm(&r) / b_norm <= tol {
            return Ok(x);
        }

        rho_old = rho_new;
    }

    Err("BiCGSTAB did not converge within max_iter.".into())
}

fn solve_banded_gauss_seidel(
    diag: &[f64],
    upper: &[Vec<f64>],
    lower: &[Vec<f64>],
    rhs: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<Vec<f64>, String> {
    let n = diag.len();
    if rhs.len() != n {
        return Err("RHS length mismatch in Gauss-Seidel solve.".into());
    }
    if n == 0 {
        return Err("Empty system in Gauss-Seidel solve.".into());
    }
    if diag.iter().any(|d| !d.is_finite() || *d == 0.0) {
        return Err("Invalid diagonal in Gauss-Seidel solve.".into());
    }

    let bw = upper.len().min(lower.len());
    let rhs_norm = l2_norm(rhs).max(1e-300);
    let mut x = vec![0.0; n];

    for _iter in 0..max_iter {
        let x_old = x.clone();
        for i in 0..n {
            let mut sum = rhs[i];

            // columns < i: (row i, col i-d) stored in lower[d-1][i-d]
            let dmax = bw.min(i);
            for d in 1..=dmax {
                let col = i - d;
                sum -= lower[d - 1][col] * x[col];
            }

            // columns > i: (row i, col i+d) stored in upper[d-1][i]
            let dmax = bw.min(n - 1 - i);
            for d in 1..=dmax {
                let col = i + d;
                sum -= upper[d - 1][i] * x_old[col];
            }

            x[i] = sum / diag[i];
        }

        let ax = banded_matvec_general(diag, upper, lower, &x);
        let mut r = vec![0.0; n];
        for i in 0..n {
            r[i] = ax[i] - rhs[i];
        }
        if l2_norm(&r) <= tol * rhs_norm {
            return Ok(x);
        }
    }

    Err("Gauss-Seidel did not converge within max_iter.".into())
}

fn solve_spd_symmetric_banded_cholesky(
    diag: &[f64],
    bands: &[Vec<f64>],
    rhs: &[f64],
) -> Result<Vec<f64>, String> {
    let n = diag.len();
    if rhs.len() != n {
        return Err("RHS length mismatch in banded Cholesky solve.".into());
    }
    if n == 0 {
        return Err("Empty system in banded Cholesky solve.".into());
    }
    if diag.iter().any(|d| !d.is_finite() || *d <= 0.0) {
        return Err("Invalid diagonal in banded Cholesky solve.".into());
    }
    if bands.iter().any(|b| b.iter().any(|x| !x.is_finite())) {
        return Err("Invalid band entries in banded Cholesky solve.".into());
    }

    // A is represented by:
    // - diag[i] = A_{i,i}
    // - bands[d-1][i] = A_{i, i+d} = A_{i+d, i} for d=1..bw
    let bw = bands.len();
    for (d, b) in bands.iter().enumerate() {
        let expected = n.saturating_sub(d + 1);
        if b.len() != expected {
            return Err("Band length mismatch in banded Cholesky solve.".into());
        }
    }

    // Cholesky factor L (lower-triangular banded), stored similarly:
    // - l_diag[i] = L_{i,i}
    // - l_sub[d-1][i] = L_{i+d, i}
    let mut l_diag = vec![0.0; n];
    let mut l_sub = (0..bw)
        .map(|d| vec![0.0; n.saturating_sub(d + 1)])
        .collect::<Vec<Vec<f64>>>();

    for i in 0..n {
        let k_min = i.saturating_sub(bw);

        // Compute L_{i,i}
        let mut sum = diag[i];
        for k in k_min..i {
            let d = i - k;
            if d == 0 || d > bw {
                continue;
            }
            let lik = l_sub[d - 1][k];
            sum -= lik * lik;
        }
        if !(sum > 0.0) || !sum.is_finite() {
            return Err(format!(
                "Banded Cholesky failed at i={} (non-SPD pivot).",
                i
            ));
        }
        l_diag[i] = sum.sqrt();

        // Compute column i below diagonal within bandwidth
        let j_max = (i + bw).min(n.saturating_sub(1));
        for j in (i + 1)..=j_max {
            let dj = j - i;
            let mut a_ji = bands[dj - 1][i];

            // subtract Σ_k L_{j,k} L_{i,k}, with k in intersection of bands
            let k_start = k_min.max(j.saturating_sub(bw));
            for k in k_start..i {
                let dik = i - k;
                let djk = j - k;
                if dik == 0 || djk == 0 || dik > bw || djk > bw {
                    continue;
                }
                let lik = l_sub[dik - 1][k];
                let ljk = l_sub[djk - 1][k];
                a_ji -= ljk * lik;
            }

            l_sub[dj - 1][i] = a_ji / l_diag[i];
        }
    }

    // Forward solve: L y = rhs
    let mut y = vec![0.0; n];
    for i in 0..n {
        let k_min = i.saturating_sub(bw);
        let mut sum = rhs[i];
        for k in k_min..i {
            let d = i - k;
            if d == 0 || d > bw {
                continue;
            }
            sum -= l_sub[d - 1][k] * y[k];
        }
        y[i] = sum / l_diag[i];
    }

    // Back solve: L^T x = y
    let mut x = vec![0.0; n];
    for i_rev in 0..n {
        let i = n - 1 - i_rev;
        let j_max = (i + bw).min(n.saturating_sub(1));
        let mut sum = y[i];
        for j in (i + 1)..=j_max {
            let d = j - i;
            sum -= l_sub[d - 1][i] * x[j];
        }
        x[i] = sum / l_diag[i];
    }

    Ok(x)
}

fn solve_spd_symmetric_banded_dense_fallback(
    diag: &[f64],
    bands: &[Vec<f64>],
    rhs: &[f64],
) -> Result<Vec<f64>, String> {
    let n = diag.len();
    if rhs.len() != n {
        return Err("RHS length mismatch in dense fallback solve.".into());
    }
    let mut a = DenseMatrix::zeros(n);
    for i in 0..n {
        a.set(i, i, diag[i]);
    }
    for d in 1..=bands.len() {
        let b = &bands[d - 1];
        for i in 0..b.len() {
            let j = i + d;
            let v = b[i];
            a.set(i, j, v);
            a.set(j, i, v);
        }
    }
    cholesky_solve_spd(&a, rhs).map_err(|e| format!("Dense Cholesky fallback failed: {e}"))
}
