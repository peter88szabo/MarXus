use std::collections::HashMap;

use MarXus::masterequation::energy_grained_me::{
    CollisionKernelModel, EnergyGrid, OlzmannMasterEquationSettings, SourceConstructionChoice,
};
use MarXus::masterequation::energy_grained_steady_state::OlzmannStepladderMasterEquationSolver;
use MarXus::masterequation::placeholder_microcanonical::PlaceholderMicrocanonicalData;

fn main() -> Result<(), String> {
    let input = include_str!("steadystate_me_input.txt");
    let cfg = parse_input(input)?;

    let grid = EnergyGrid {
        number_of_bins: get_usize(&cfg, "parameters", "number_of_bins")?,
        bin_width_wavenumber: get_f64(&cfg, "parameters", "bin_width_wavenumber")?,
        energy_origin_wavenumber: get_f64(&cfg, "parameters", "energy_origin_wavenumber")?,
    };

    let settings = OlzmannMasterEquationSettings {
        temperature_kelvin: get_f64(&cfg, "parameters", "temperature_kelvin")?,
        boltzmann_constant_wavenumber_per_kelvin: get_f64(
            &cfg,
            "parameters",
            "boltzmann_constant_wavenumber_per_kelvin",
        )?,
        collision_frequency_per_second: get_f64(
            &cfg,
            "parameters",
            "collision_frequency_per_second",
        )?,
        pseudo_first_order_capture_loss_per_second: get_f64(
            &cfg,
            "parameters",
            "pseudo_first_order_capture_loss_per_second",
        )?,
        stepladder_base_downward_probability: get_f64(
            &cfg,
            "parameters",
            "stepladder_base_downward_probability",
        )?,
        collision_kernel_model: CollisionKernelModel::Stepladder,
    };

    let formation_flux = get_f64(&cfg, "parameters", "formation_flux")?;

    let provider = get_string(&cfg, "microcanonical", "provider")?;
    if provider != "placeholder" {
        return Err(format!(
            "Unsupported microcanonical provider '{}'. Only 'placeholder' is supported in this example.",
            provider
        ));
    }

    let channel_count = get_usize(&cfg, "microcanonical", "placeholder_channel_count")?;
    let micro = PlaceholderMicrocanonicalData { channel_count };

    let solver = OlzmannStepladderMasterEquationSolver {
        energy_grid: grid.clone(),
        settings,
    };

    let active = get_string(&cfg, "source", "active")?;
    let source = match active.as_str() {
        "thermal_ts" => SourceConstructionChoice::ThermalTransitionStateSumOfStates {
            threshold_wavenumber: get_f64(&cfg, "source_thermal_ts", "threshold_wavenumber")?,
            transition_state_sum_of_states: get_vec_f64(
                &cfg,
                "source_thermal_ts",
                "transition_state_sum_of_states",
            )?,
        },
        "convolution" => SourceConstructionChoice::NonthermalShiftedConvolution {
            threshold_wavenumber: get_f64(&cfg, "source_convolution", "threshold_wavenumber")?,
            reactant_a_distribution: fit_distribution_to_grid(
                get_vec_f64(&cfg, "source_convolution", "reactant_a_distribution")?,
                grid.number_of_bins,
            )?,
            reactant_b_distribution: fit_distribution_to_grid(
                get_vec_f64(&cfg, "source_convolution", "reactant_b_distribution")?,
                grid.number_of_bins,
            )?,
        },
        "shift" => SourceConstructionChoice::ShiftApproximation {
            reaction_energy_wavenumber: get_f64(
                &cfg,
                "source_shift",
                "reaction_energy_wavenumber",
            )?,
            partner_mean_energy_wavenumber: get_f64(
                &cfg,
                "source_shift",
                "partner_mean_energy_wavenumber",
            )?,
            reactant_distribution: fit_distribution_to_grid(
                get_vec_f64(&cfg, "source_shift", "reactant_distribution")?,
                grid.number_of_bins,
            )?,
        },
        _ => {
            return Err(format!(
                "Invalid source.active '{}'. Use thermal_ts | convolution | shift.",
                active
            ));
        }
    };

    let result = solver.solve_steady_state_with_source_choice(&micro, formation_flux, &source)?;

    println!("Active source: {}", active);
    println!(
        "Total loss = {:.6e}",
        result.total_unimolecular_loss_rate_constant
    );
    println!(
        "Channel k = {:?}",
        result.chemically_activated_rate_constants
    );

    Ok(())
}

type InputMap = HashMap<(String, String), String>;

fn parse_input(input: &str) -> Result<InputMap, String> {
    let mut map: InputMap = HashMap::new();
    let mut section = String::new();

    let mut current_key: Option<String> = None;
    let mut current_value = String::new();

    for (idx, raw_line) in input.lines().enumerate() {
        let line_no = idx + 1;

        let mut line = raw_line.trim().to_string();
        if let Some(hash) = line.find('#') {
            line = line[..hash].trim().to_string();
        }
        if line.is_empty() {
            continue;
        }

        if let Some(key) = &current_key {
            current_value.push(' ');
            current_value.push_str(&line);
            if line.contains(']') {
                map.insert(
                    (section.clone(), key.clone()),
                    current_value.trim().to_string(),
                );
                current_key = None;
                current_value.clear();
            }
            continue;
        }

        if line.starts_with('[') && line.ends_with(']') {
            section = line[1..line.len() - 1].trim().to_string();
            if section.is_empty() {
                return Err(format!("Empty section name at line {}", line_no));
            }
            continue;
        }

        let parts: Vec<&str> = line.splitn(2, '=').collect();
        if parts.len() != 2 {
            return Err(format!("Expected key=value at line {}", line_no));
        }

        if section.is_empty() {
            return Err(format!("Key outside any section at line {}", line_no));
        }

        let key = parts[0].trim().to_string();
        let value = parts[1].trim().to_string();

        if key.is_empty() {
            return Err(format!("Empty key at line {}", line_no));
        }

        if value.starts_with('[') && !value.contains(']') {
            current_key = Some(key);
            current_value = value;
            continue;
        }

        map.insert((section.clone(), key), value);
    }

    if current_key.is_some() {
        return Err("Unterminated list value in input.".into());
    }

    Ok(map)
}

fn get_raw<'a>(map: &'a InputMap, section: &str, key: &str) -> Result<&'a str, String> {
    map.get(&(section.to_string(), key.to_string()))
        .map(|s| s.as_str())
        .ok_or_else(|| format!("Missing {}.{}", section, key))
}

fn get_string(map: &InputMap, section: &str, key: &str) -> Result<String, String> {
    Ok(get_raw(map, section, key)?.trim().to_string())
}

fn get_usize(map: &InputMap, section: &str, key: &str) -> Result<usize, String> {
    get_raw(map, section, key)?
        .parse::<usize>()
        .map_err(|_| format!("Invalid usize for {}.{}", section, key))
}

fn get_f64(map: &InputMap, section: &str, key: &str) -> Result<f64, String> {
    get_raw(map, section, key)?
        .parse::<f64>()
        .map_err(|_| format!("Invalid f64 for {}.{}", section, key))
}

fn get_vec_f64(map: &InputMap, section: &str, key: &str) -> Result<Vec<f64>, String> {
    let raw = get_raw(map, section, key)?.trim();
    let mut s = raw.to_string();
    if s.starts_with('[') && s.ends_with(']') {
        s = s[1..s.len() - 1].to_string();
    }

    if s.trim().is_empty() {
        return Ok(Vec::new());
    }

    let mut out = Vec::new();
    for token in s.split(',') {
        let t = token.trim();
        if t.is_empty() {
            continue;
        }
        let v = t
            .parse::<f64>()
            .map_err(|_| format!("Invalid list number for {}.{}: '{}'", section, key, t))?;
        out.push(v);
    }

    Ok(out)
}

fn fit_distribution_to_grid(values: Vec<f64>, n: usize) -> Result<Vec<f64>, String> {
    if values.is_empty() {
        return Err("Distribution list is empty.".into());
    }
    if values.iter().any(|x| *x < 0.0 || !x.is_finite()) {
        return Err("Distribution has negative or invalid entries.".into());
    }
    if values.len() > n {
        return Err(format!(
            "Distribution length {} exceeds number_of_bins {}.",
            values.len(),
            n
        ));
    }

    let mut out = values;
    if out.len() < n {
        out.resize(n, 0.0);
    }

    let sum: f64 = out.iter().sum();
    if sum <= 0.0 {
        return Err("Distribution sum must be positive.".into());
    }

    for x in &mut out {
        *x /= sum;
    }

    Ok(out)
}
