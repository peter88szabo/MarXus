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
        compute_collision_frequency_s_inv(&well.collision_params, settings, temperature, pressure)?;

    let local_start = well.lowest_included_grain_index;
    let local_end_exclusive = well.one_past_highest_included_grain_index;
    let band = settings.collision_band_half_width;

    let grain_width = well.energy_grain_width_cm1;
    let boltzmann = settings.boltzmann_constant_wavenumber_per_kelvin;

    match settings.collision_kernel_implementation {
        CollisionKernelImplementation::Mess => {
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

                let mut unnormalized_weights: Vec<(usize, f64)> = Vec::new();
                for target_grain in target_min..target_max_exclusive {
                    let delta_grains = if target_grain >= source_grain {
                        (target_grain - source_grain) as i64
                    } else {
                        -((source_grain - target_grain) as i64)
                    };

                    let step_energy_cm1 = (delta_grains.abs() as f64) * grain_width;
                    let exponential_down = (-step_energy_cm1 / alpha_cm1).exp();

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

                let residual = (p_same + offdiag_prob_sum) - 1.0;
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
