use super::energy_grained_me::{EnergyGrid, OlzmannMasterEquationSettings, SourceConstructionChoice};

/// Build a normalized source distribution F over the energy grid.
pub fn build_normalized_source_distribution(
    energy_grid: &EnergyGrid,
    settings: &OlzmannMasterEquationSettings,
    source_choice: &SourceConstructionChoice,
) -> Result<Vec<f64>, String> {
    let mut raw = match source_choice {
        SourceConstructionChoice::ThermalTransitionStateSumOfStates {
            threshold_wavenumber,
            transition_state_sum_of_states,
        } => build_source_thermal_ts_sum_of_states(
            energy_grid,
            settings,
            *threshold_wavenumber,
            transition_state_sum_of_states,
        )?,

        SourceConstructionChoice::NonthermalShiftedConvolution {
            threshold_wavenumber,
            reactant_a_distribution,
            reactant_b_distribution,
        } => build_source_nonthermal_shifted_convolution(
            energy_grid,
            *threshold_wavenumber,
            reactant_a_distribution,
            reactant_b_distribution,
        )?,

        SourceConstructionChoice::ShiftApproximation {
            reactant_distribution,
            reaction_energy_wavenumber,
            partner_mean_energy_wavenumber,
        } => build_source_shift_approximation(
            energy_grid,
            reactant_distribution,
            *reaction_energy_wavenumber,
            *partner_mean_energy_wavenumber,
        )?,
    };

    normalize_nonnegative_vector_in_place(&mut raw)?;
    Ok(raw)
}

fn build_source_thermal_ts_sum_of_states(
    energy_grid: &EnergyGrid,
    settings: &OlzmannMasterEquationSettings,
    threshold_wavenumber: f64,
    transition_state_sum_of_states: &[f64],
) -> Result<Vec<f64>, String> {
    if transition_state_sum_of_states.is_empty() {
        return Err("Transition-state sum-of-states array must be non-empty.".into());
    }
    if transition_state_sum_of_states.iter().any(|x| *x < 0.0) {
        return Err("Transition-state sum-of-states must be nonnegative.".into());
    }

    let n = energy_grid.number_of_bins;
    let offset_bins = energy_grid.energy_shift_to_bin_shift(threshold_wavenumber)?;
    if offset_bins < 0 {
        return Err("Negative threshold energy shift is not supported for this source.".into());
    }
    let offset_bins = offset_bins as usize;

    let t = settings.temperature_kelvin;
    let k_b = settings.boltzmann_constant_wavenumber_per_kelvin;
    let delta_e = energy_grid.bin_width_wavenumber;

    let mut raw = vec![0.0; n];
    for (i, value) in raw.iter_mut().enumerate().skip(offset_bins) {
        let k = i - offset_bins;
        if k >= transition_state_sum_of_states.len() {
            break;
        }
        let w_ts = transition_state_sum_of_states[k];
        let epsilon = (k as f64) * delta_e;
        *value = w_ts * (-epsilon / (k_b * t)).exp();
    }

    Ok(raw)
}

fn build_source_nonthermal_shifted_convolution(
    energy_grid: &EnergyGrid,
    threshold_wavenumber: f64,
    reactant_a_distribution: &[f64],
    reactant_b_distribution: &[f64],
) -> Result<Vec<f64>, String> {
    let n = energy_grid.number_of_bins;

    validate_normalized_distribution("reactant_a_distribution", reactant_a_distribution, n)?;
    validate_normalized_distribution("reactant_b_distribution", reactant_b_distribution, n)?;

    let offset_bins = energy_grid.energy_shift_to_bin_shift(threshold_wavenumber)?;
    if offset_bins < 0 {
        return Err("Negative threshold energy shift is not supported for this source.".into());
    }
    let offset_bins = offset_bins as usize;

    let mut raw = vec![0.0; n];

    for (i, value) in raw.iter_mut().enumerate().skip(offset_bins) {
        let k = i - offset_bins;

        let mut sum = 0.0;
        let jmax = k.min(n - 1);
        for j in 0..=jmax {
            let b_index = k.checked_sub(j).unwrap_or(usize::MAX);
            if b_index >= n {
                continue;
            }
            sum += reactant_a_distribution[j] * reactant_b_distribution[b_index];
        }
        *value = sum.max(0.0);
    }

    Ok(raw)
}

fn build_source_shift_approximation(
    energy_grid: &EnergyGrid,
    reactant_distribution: &[f64],
    reaction_energy_wavenumber: f64,
    partner_mean_energy_wavenumber: f64,
) -> Result<Vec<f64>, String> {
    let n = energy_grid.number_of_bins;
    validate_normalized_distribution("reactant_distribution", reactant_distribution, n)?;

    let shift_wavenumber = reaction_energy_wavenumber - partner_mean_energy_wavenumber;
    let shift_bins = energy_grid.energy_shift_to_bin_shift(shift_wavenumber)?;

    let mut raw = vec![0.0; n];

    for (i, value) in raw.iter_mut().enumerate() {
        let src = (i as isize) + shift_bins;
        if src >= 0 && src < n as isize {
            *value = reactant_distribution[src as usize];
        }
    }

    Ok(raw)
}

fn validate_normalized_distribution(
    name: &str,
    dist: &[f64],
    expected_len: usize,
) -> Result<(), String> {
    if dist.len() != expected_len {
        return Err(format!(
            "{} must have length {} (same as energy grid). Got {}.",
            name,
            expected_len,
            dist.len()
        ));
    }
    if dist.iter().any(|x| *x < 0.0) {
        return Err(format!("{} contains negative entries.", name));
    }
    let sum: f64 = dist.iter().sum();
    if (sum - 1.0).abs() > 1e-10 {
        return Err(format!(
            "{} must be normalized (sum=1). Got sum={}",
            name, sum
        ));
    }

    Ok(())
}

fn normalize_nonnegative_vector_in_place(vec: &mut [f64]) -> Result<(), String> {
    if vec.iter().any(|x| *x < 0.0) {
        return Err(
            "Vector contains negative entries; cannot normalize as a probability distribution."
                .into(),
        );
    }

    let sum: f64 = vec.iter().sum();
    if sum <= 0.0 || !sum.is_finite() {
        return Err("Vector has zero (or invalid) total weight; cannot normalize.".into());
    }

    for x in vec.iter_mut() {
        *x /= sum;
    }

    Ok(())
}
