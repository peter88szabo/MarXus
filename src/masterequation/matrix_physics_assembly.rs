use super::collisional_relaxation::{
    band_limits, compute_alpha_cm1, compute_collision_frequency_s_inv,
};
use super::reaction_network::{
    CollisionKernelImplementation, MasterEquationSettings, MicrocanonicalProvider, WellDefinition,
};
use super::state_index::GlobalLayout;
use crate::numeric::linear_algebra::DenseMatrix;

pub struct OperatorAssemblyDiagnostics {
    pub collision_conservation_max_abs: f64,
}

#[derive(Default)]
struct CollisionDiagnostics {
    max_column_sum_abs: f64,
}

pub(crate) fn absolute_energy_cm1(well: &WellDefinition, local_grain: usize) -> f64 {
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
            to_well
                .one_past_highest_included_grain_index
                .saturating_sub(1),
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

pub fn assemble_raw_operator(
    wells: &[WellDefinition],
    settings: &MasterEquationSettings,
    micro: &dyn MicrocanonicalProvider,
    layout: &GlobalLayout,
) -> Result<(DenseMatrix, OperatorAssemblyDiagnostics), String> {
    let mut collision_diag = CollisionDiagnostics::default();
    let mut raw_operator = DenseMatrix::zeros(layout.total_state_count);

    for (well_index, well) in wells.iter().enumerate() {
        add_collision_terms_for_well(
            settings,
            micro,
            layout,
            well_index,
            well,
            &mut raw_operator,
            &mut collision_diag,
        )?;
    }

    add_reaction_sinks(wells, settings, micro, layout, &mut raw_operator)?;
    add_interwell_couplings(wells, settings, micro, layout, &mut raw_operator)?;

    Ok((
        raw_operator,
        OperatorAssemblyDiagnostics {
            collision_conservation_max_abs: collision_diag.max_column_sum_abs,
        },
    ))
}

fn add_collision_terms_for_well(
    settings: &MasterEquationSettings,
    micro: &dyn MicrocanonicalProvider,
    layout: &GlobalLayout,
    well_index: usize,
    well: &WellDefinition,
    operator: &mut DenseMatrix,
    diag: &mut CollisionDiagnostics,
) -> Result<(), String> {
    let temperature = settings.temperature_kelvin;
    let pressure = settings.pressure_torr;

    let alpha_cm1 = compute_alpha_cm1(
        well.collision_params.alpha_at_1000K_cm1,
        well.collision_params.alpha_temperature_exponent,
        temperature,
    );

    let collision_frequency_s_inv =
        compute_collision_frequency_s_inv(&well.collision_params, temperature, pressure)?;

    let local_start = well.lowest_included_grain_index;
    let local_end_exclusive = well.one_past_highest_included_grain_index;
    let band = settings.collision_band_half_width;

    let grain_width = well.energy_grain_width_cm1;
    let boltzmann = settings.boltzmann_constant_wavenumber_per_kelvin;

    match settings.collision_kernel_implementation {
        CollisionKernelImplementation::Mess => {
            // Exponential-down model with detailed balance and exact normalization
            // (Robertson, Comprehensive Chemical Kinetics 43 (2019), eqs. 4.4, 4.7, 4.11, 4.16):
            //   deactivating (j <= i): P(j|i) = A_i exp(-(E_i - E_j)/alpha)
            //   activating   (j >  i): P(j|i) = A_j (rho_j/rho_i) exp(-(E_j - E_i)(1/alpha + 1/kT))
            //   sum_j P(j|i) = 1.
            // The normalization equations are upper triangular in A_i and are solved by back
            // substitution from the top grain, where only deactivating collisions exist.
            let n_local = local_end_exclusive.saturating_sub(local_start);
            let beta = 1.0 / (boltzmann * temperature);

            let mut rho = vec![0.0; n_local];
            for (idx, grain) in (local_start..local_end_exclusive).enumerate() {
                let r = micro.density_of_states(well_index, grain);
                if r <= 0.0 {
                    return Err(format!(
                        "Non-positive density of states at well={}, grain={}",
                        well.well_name, grain
                    ));
                }
                rho[idx] = r;
            }

            let deactivating = |steps: usize| (-(steps as f64) * grain_width / alpha_cm1).exp();
            let activating = |steps: usize| {
                (-(steps as f64) * grain_width * (1.0 / alpha_cm1 + beta)).exp()
            };

            // For sparse low-energy states the back substitution can return non-positive
            // coefficients (Robertson, CCK 43, p. 294). Levels below half the lowest reaction
            // threshold E0 of the well stay at their equilibrium populations, so they all share the
            // normalization of the grain above them; detailed balance is unaffected because the
            // activating probabilities are built from the same coefficients. A breakdown above
            // E0/2 (or in a well without any reactive channel) is reported as an error.
            let threshold_idx = (0..n_local).find(|&idx| {
                let grain = local_start + idx;
                (0..well.channels.len())
                    .any(|ch| micro.microcanonical_rate(well_index, ch, grain) > 0.0)
            });
            let cut_idx = threshold_idx.map(|t| t / 2).unwrap_or(0);

            let mut a_norm = vec![0.0; n_local];
            for idx in (0..n_local).rev() {
                if idx < cut_idx && idx + 1 < n_local {
                    a_norm[idx] = a_norm[idx + 1];
                    continue;
                }
                let grain = local_start + idx;
                let (target_min, target_max_exclusive) =
                    band_limits(grain, local_start, local_end_exclusive, band);

                let down_sum: f64 = (target_min..=grain).map(|t| deactivating(grain - t)).sum();
                let up_sum: f64 = ((grain + 1)..target_max_exclusive)
                    .map(|t| {
                        let t_idx = t - local_start;
                        a_norm[t_idx] * (rho[t_idx] / rho[idx]) * activating(t - grain)
                    })
                    .sum();

                let remainder = 1.0 - up_sum;
                if !(remainder > 0.0) || !down_sum.is_finite() {
                    return Err(format!(
                        "Exponential-down normalization failed at well={}, grain={}: activating \
                         collisions already carry a probability of {:.6}.",
                        well.well_name, grain, up_sum
                    ));
                }
                a_norm[idx] = remainder / down_sum;
            }

            for idx in 0..n_local {
                let source_grain = local_start + idx;
                let (target_min, target_max_exclusive) =
                    band_limits(source_grain, local_start, local_end_exclusive, band);
                let global_source = layout.global_index_of(well_index, source_grain)?;

                let mut out_probability = 0.0;
                for target_grain in target_min..target_max_exclusive {
                    if target_grain == source_grain {
                        continue;
                    }
                    let t_idx = target_grain - local_start;
                    let probability = if target_grain < source_grain {
                        a_norm[idx] * deactivating(source_grain - target_grain)
                    } else {
                        a_norm[t_idx]
                            * (rho[t_idx] / rho[idx])
                            * activating(target_grain - source_grain)
                    };
                    let global_target = layout.global_index_of(well_index, target_grain)?;
                    operator.add(
                        global_target,
                        global_source,
                        collision_frequency_s_inv * probability,
                    );
                    out_probability += probability;
                }
                operator.add(
                    global_source,
                    global_source,
                    -collision_frequency_s_inv * out_probability,
                );

                // Elastic (j = i) probability is A_i; the column is normalized to 1 by construction.
                let residual = a_norm[idx] + out_probability - 1.0;
                diag.max_column_sum_abs = diag
                    .max_column_sum_abs
                    .max((collision_frequency_s_inv * residual).abs());
            }
        }
        CollisionKernelImplementation::Spd => {
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

                operator.add(global_source, global_source, -omega_eff * out_raw[idx]);
                diag.max_column_sum_abs = diag.max_column_sum_abs.max(col_sum.abs());
            }
        }
    }

    Ok(())
}

fn add_reaction_sinks(
    wells: &[WellDefinition],
    settings: &MasterEquationSettings,
    micro: &dyn MicrocanonicalProvider,
    layout: &GlobalLayout,
    operator: &mut DenseMatrix,
) -> Result<(), String> {
    let out_thresh = settings.outgoing_rate_threshold;
    let skip_internal = settings.enforce_interwell_detailed_balance;

    for (well_index, well) in wells.iter().enumerate() {
        for local_grain in
            well.lowest_included_grain_index..well.one_past_highest_included_grain_index
        {
            let global_state = layout.global_index_of(well_index, local_grain)?;

            let mut total_loss_rate = 0.0;
            if local_grain >= well.nonreactive_grain_count {
                for (channel_index, channel) in well.channels.iter().enumerate() {
                    if skip_internal && channel.connected_well_index.is_some() {
                        continue;
                    }
                    let mut k = micro.microcanonical_rate(well_index, channel_index, local_grain);
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

fn add_interwell_couplings(
    wells: &[WellDefinition],
    settings: &MasterEquationSettings,
    micro: &dyn MicrocanonicalProvider,
    layout: &GlobalLayout,
    operator: &mut DenseMatrix,
) -> Result<(), String> {
    let internal_thresh = settings.internal_rate_threshold;
    let enforce_db = settings.enforce_interwell_detailed_balance;
    let temperature = settings.temperature_kelvin;
    let boltzmann = settings.boltzmann_constant_wavenumber_per_kelvin;

    if !enforce_db {
        for (from_well_index, from_well) in wells.iter().enumerate() {
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

                    let to_well = &wells[to_well_index];
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

    let mut unique_to: Vec<std::collections::HashMap<usize, usize>> =
        vec![std::collections::HashMap::new(); wells.len()];
    for (w, well) in wells.iter().enumerate() {
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
    for from in 0..wells.len() {
        for (&to, &ch_from_to) in &unique_to[from] {
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
        let well_i = &wells[w_i];
        let well_j = &wells[w_j];

        for grain_i in
            well_i.lowest_included_grain_index..well_i.one_past_highest_included_grain_index
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

            let k_ij = micro.microcanonical_rate(w_i, ch_i_to_j, grain_i).max(0.0);
            let k_ji = micro.microcanonical_rate(w_j, ch_j_to_i, grain_j).max(0.0);

            if k_ij < internal_thresh || k_ji < internal_thresh {
                continue;
            }

            let rho_i = micro.density_of_states(w_i, grain_i);
            let rho_j = micro.density_of_states(w_j, grain_j);
            if rho_i <= 0.0 || rho_j <= 0.0 {
                return Err(
                    "Non-positive density of states in inter-well detailed balance.".into(),
                );
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::KB_CM;
    use crate::masterequation::microcanonical_builder::{
        ArrayMicrocanonicalProvider, MicrocanonicalNetworkData,
    };
    use crate::masterequation::reaction_network::{
        CollisionKernelImplementation, CollisionModelParams, MultiwellLinearSolver,
        ReactionChannel,
    };

    fn one_well(
        with_reaction: bool,
    ) -> (Vec<WellDefinition>, MasterEquationSettings, MicrocanonicalNetworkData, Vec<f64>) {
        let n_grains = 60;
        let channels = if with_reaction {
            vec![ReactionChannel { name: "P".to_string(), connected_well_index: None }]
        } else {
            vec![]
        };
        let well = WellDefinition {
            well_name: "W".to_string(),
            energy_grain_width_cm1: 20.0,
            lowest_included_grain_index: 0,
            one_past_highest_included_grain_index: n_grains,
            alignment_offset_in_grains: 0,
            nonreactive_grain_count: 0,
            collision_params: CollisionModelParams {
                lennard_jones_sigma_angstrom: 5.0,
                lennard_jones_epsilon_kelvin: 300.0,
                reduced_mass_amu: 25.0,
                alpha_at_1000K_cm1: 300.0,
                alpha_temperature_exponent: 0.85,
            },
            channels: channels.clone(),
        };
        let settings = MasterEquationSettings {
            temperature_kelvin: 500.0,
            pressure_torr: 10.0,
            boltzmann_constant_wavenumber_per_kelvin: KB_CM,
            collision_band_half_width: 20,
            collision_kernel_implementation: CollisionKernelImplementation::Mess,
            outgoing_rate_threshold: 0.0,
            internal_rate_threshold: 0.0,
            enforce_interwell_detailed_balance: false,
            linear_solver: MultiwellLinearSolver::Direct,
            krylov_tolerance: 1e-12,
            krylov_max_iter: 100,
            gmres_restart: 10,
        };
        let rho: Vec<f64> = (0..n_grains).map(|i| (1.0 + 0.05 * i as f64).powi(10)).collect();
        // Reaction threshold at grain 40 (k = 0 below), so E0/2 corresponds to grain 20.
        let k: Vec<f64> = (0..n_grains)
            .map(|i| if i < 40 { 0.0 } else { 1.0e7 * (i - 39) as f64 })
            .collect();
        let data = MicrocanonicalNetworkData {
            rho_by_well: vec![rho.clone()],
            k_by_well_by_channel: vec![if with_reaction { vec![k] } else { vec![] }],
            channels_by_well: vec![channels],
        };
        (vec![well], settings, data, rho)
    }

    #[test]
    fn exponential_down_kernel_is_normalized_and_obeys_detailed_balance() {
        // Robertson, Comprehensive Chemical Kinetics 43 (2019), Ch. 4, eqs. 4.4 and 4.6:
        // every collision column must conserve population and R_ij f_j = R_ji f_i
        // with the Boltzmann distribution f_i = rho_i exp(-E_i/kT).
        let (wells, settings, data, rho) = one_well(true);
        let micro = ArrayMicrocanonicalProvider::new(&data);
        let layout = GlobalLayout::from_wells(&wells).unwrap();
        let (r_full, _) = assemble_raw_operator(&wells, &settings, &micro, &layout).unwrap();
        let n = rho.len();
        let d_e = wells[0].energy_grain_width_cm1;
        let t = settings.temperature_kelvin;

        // Remove the reactive loss from the diagonal to isolate the collision operator.
        let k_total = |i: usize| micro.microcanonical_rate(0, 0, i);
        let r = |i: usize, j: usize| if i == j { r_full.get(i, j) + k_total(i) } else { r_full.get(i, j) };

        let f: Vec<f64> =
            (0..n).map(|i| rho[i] * (-(i as f64) * d_e / (KB_CM * t)).exp()).collect();
        let z = r(1, 0).abs().max(r(0, 0).abs());
        for j in 0..n {
            let col: f64 = (0..n).map(|i| r(i, j)).sum();
            assert!(col.abs() < 1e-9 * z, "column {j} sums to {col:e}");
            for i in 0..n {
                if i == j {
                    continue;
                }
                let lhs = r(i, j) * f[j];
                let rhs = r(j, i) * f[i];
                let scale = lhs.abs().max(rhs.abs());
                if scale > 0.0 {
                    assert!(
                        ((lhs - rhs) / scale).abs() < 1e-10,
                        "detailed balance violated for ({i},{j}): {lhs:e} vs {rhs:e}"
                    );
                }
            }
        }
    }

    #[test]
    fn normalization_breakdown_without_reaction_threshold_is_an_error() {
        // Without a reactive channel there is no E0 below which the equilibrium populations
        // may share one normalization, so a breakdown of the back substitution must be reported.
        let (wells, settings, data, _) = one_well(false);
        let micro = ArrayMicrocanonicalProvider::new(&data);
        let layout = GlobalLayout::from_wells(&wells).unwrap();
        let result = assemble_raw_operator(&wells, &settings, &micro, &layout);
        assert!(result.is_err(), "expected a normalization error");
    }
}
