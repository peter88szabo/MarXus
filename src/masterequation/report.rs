//! General reporting helpers for master-equation calculations.
//!
//! Keep these printers physics-facing: print energies, barriers, T/P conditions,
//! high-pressure-limit reference rates, and steady-state master-equation outputs.

use std::collections::{BTreeSet, HashMap};

use crate::masterequation::high_pressure_limit::HighPressureThermoResult;
use crate::masterequation::reaction_network::{SolutionResults, WellDefinition};

pub struct ReactionEnergetics {
    pub well_name: String,
    pub product_name: String,
    pub barrier_name: String,
    pub barrier_label: String,
    pub well_dh0_kcal_mol: f64,
    pub product_dh0_kcal_mol: f64,
    pub barrier_dh0_kcal_mol: f64,
}

impl ReactionEnergetics {
    pub fn delta_e_kcal_mol(&self) -> f64 {
        self.product_dh0_kcal_mol - self.well_dh0_kcal_mol
    }
    pub fn ea_forward_kcal_mol(&self) -> f64 {
        self.barrier_dh0_kcal_mol - self.well_dh0_kcal_mol
    }
    pub fn ea_reverse_kcal_mol(&self) -> f64 {
        self.barrier_dh0_kcal_mol - self.product_dh0_kcal_mol
    }
}

pub fn print_reaction_header(temperature_kelvin: f64, pressure_atm: f64, e: &ReactionEnergetics) {
    println!(
        "Steady-state dissociation ({} -> {}, barrier {})",
        e.well_name, e.product_name, e.barrier_name
    );
    println!(
        "T = {:.1} K, P = {:.3} atm",
        temperature_kelvin, pressure_atm
    );
    println!("Energetics (kcal/mol):");
    println!("E({}) = {:.3}", e.well_name, e.well_dh0_kcal_mol);
    println!("E({}) = {:.3}", e.product_name, e.product_dh0_kcal_mol);
    println!("E({}) = {:.3}", e.barrier_name, e.barrier_dh0_kcal_mol);
    println!(
        "DeltaE({} -> {}) = {:.3}",
        e.well_name,
        e.product_name,
        e.delta_e_kcal_mol()
    );
    println!(
        "Ea_forward({} -> {}) = {:.3}",
        e.well_name,
        e.product_name,
        e.ea_forward_kcal_mol()
    );
    println!(
        "Ea_reverse({} -> {}) = {:.3}",
        e.product_name,
        e.well_name,
        e.ea_reverse_kcal_mol()
    );
    println!("TS structure energetics:");
    println!("TS label: {}", e.barrier_label);
    println!(
        "E(TS) = E({}) = {:.3} kcal/mol",
        e.barrier_name, e.barrier_dh0_kcal_mol
    );
    println!(
        "E(TS)-E({}) = {:.3} kcal/mol",
        e.well_name,
        e.ea_forward_kcal_mol()
    );
    println!(
        "E(TS)-E({}) = {:.3} kcal/mol",
        e.product_name,
        e.ea_reverse_kcal_mol()
    );
}

pub fn print_high_pressure_limit(
    e: &ReactionEnergetics,
    forward: &HighPressureThermoResult,
    reverse: &HighPressureThermoResult,
    tunneling_model: &str,
    tunneling_kappa: f64,
) {
    println!("High-pressure TST:");
    println!(
        "k_inf({} -> {}) = {:.6e} s^-1  [DeltaG‡={:.3} kcal/mol, kappa={:.3}]",
        e.well_name, e.product_name, forward.rate_constant, forward.d_g_kcal_mol, tunneling_kappa
    );
    println!(
        "Tunneling model: {} (with stability guard)",
        tunneling_model
    );
    println!(
        "Thermo terms: dZPE={:.3} kcal/mol, dH0={:.3} kcal/mol, prefactor={:.6e}",
        forward.d_zpe_kcal_mol, forward.d_h0_kcal_mol, forward.prefactor
    );
    println!(
        "k_inf({} -> {}) = {:.6e} cm^3 molecule^-1 s^-1  [DeltaG‡={:.3} kcal/mol]",
        e.product_name, e.well_name, reverse.rate_constant, reverse.d_g_kcal_mol
    );
}

pub fn print_source_threshold(threshold_wavenumber: f64) {
    let e0_kcal = threshold_wavenumber * 2.859_144e-3;
    println!(
        "TS source threshold: E0 = {:.3} cm^-1 ({:.3} kcal/mol)",
        threshold_wavenumber, e0_kcal
    );
}

pub fn print_steady_state_summary(total_loss: f64, channel_rates: &[f64]) {
    println!("Mode: chemical activation");
    println!("Total loss = {:.6e}", total_loss);
    if let Some(first) = channel_rates.first().copied() {
        println!("Channel-0 k = {:.6e} s^-1", first);
    }
    println!("Channel k = {:?}", channel_rates);
}

pub fn print_species_species_table(
    temperature_kelvin: f64,
    pressure_atm: f64,
    e: &ReactionEnergetics,
    hp_forward: f64,
    hp_reverse: f64,
    one_atm_forward: f64,
    one_atm_reverse: f64,
) {
    println!();
    println!(
        "______________________________________________________________________________________"
    );
    println!("Species-Species Rate Tables:");
    println!();
    println!("Temperature = {:.0} K", temperature_kelvin);
    println!();
    println!("High Pressure Rate Coefficients:");
    println!("From\\To{:>14}{:>13}", e.well_name, e.product_name);
    println!("{:<6}{:>14}{:>13.6e}", e.well_name, "***", hp_forward);
    println!("{:<6}{:>14.6e}{:>13}", e.product_name, hp_reverse, "***");
    println!();
    println!(
        "Temperature = {:.0} K    Pressure = {:.3} atm",
        temperature_kelvin, pressure_atm
    );
    println!("From\\To{:>14}{:>13}", e.well_name, e.product_name);
    println!(
        "{:<6}{:>14.6e}{:>13.6e}",
        e.well_name, one_atm_forward, one_atm_forward
    );
    println!(
        "{:<6}{:>14.6e}{:>13.6e}",
        e.product_name, one_atm_reverse, one_atm_reverse
    );
}

pub fn print_species_species_rate_matrix(
    temperature_kelvin: f64,
    pressure_label: Option<&str>,
    species_names: &[String],
    matrix: &[Vec<Option<f64>>],
) -> Result<(), String> {
    if matrix.len() != species_names.len() {
        return Err("Rate matrix row count must match species_names length.".into());
    }
    if matrix.iter().any(|row| row.len() != species_names.len()) {
        return Err("Rate matrix must be square with size species_names.len().".into());
    }

    let name_width = species_names
        .iter()
        .map(|s| s.len())
        .max()
        .unwrap_or(4)
        .max("From\\To".len())
        .min(24);

    // Fixed column width to keep alignment stable.
    let col_width = 14usize;

    println!();
    println!(
        "______________________________________________________________________________________"
    );
    println!("Species-Species Rate Tables:");
    println!();
    println!("Temperature = {:.0} K", temperature_kelvin);
    println!();

    if pressure_label.is_none() {
        println!("High Pressure Rate Coefficients:");
    } else {
        println!(
            "Temperature = {:.0} K    Pressure = {}",
            temperature_kelvin,
            pressure_label.unwrap()
        );
    }

    // Header row
    print!("{:<name_width$}", "From\\To", name_width = name_width);
    for name in species_names {
        print!("{:>col_width$}", name, col_width = col_width);
    }
    println!();

    for (i, from) in species_names.iter().enumerate() {
        print!("{:<name_width$}", from, name_width = name_width);
        for j in 0..species_names.len() {
            if i == j {
                print!("{:>col_width$}", "***", col_width = col_width);
                continue;
            }
            match matrix[i][j] {
                Some(v) => print!("{:>col_width$.6e}", v, col_width = col_width),
                None => print!("{:>col_width$}", "***", col_width = col_width),
            }
        }
        println!();
    }

    Ok(())
}

/// Build a full "species-species" macroscopic rate matrix from a multiwell steady-state solution.
///
/// What this includes
/// ------------------
/// - Rows ("From"): each well in `wells`
/// - Columns ("To"): all wells + all sink endpoints that appear in the channel list
/// - Entries: `solution.per_well_per_channel_rates[w][c]`, mapped to the corresponding target
///
/// Notes
/// -----
/// - This reports *macroscopic* rates (population-weighted averages) from the ME steady-state
///   solution, not microcanonical k(E).
/// - Pure bimolecular association rates (e.g. "R -> G2") are not produced by the ME solve itself,
///   because the reactant pool is not represented as an energy-grained well. Those must be
///   computed separately (tight TST or PST capture) and inserted into a reporting table by the caller.
pub fn build_species_species_rate_matrix_from_solution(
    wells: &[WellDefinition],
    solution: &SolutionResults,
) -> Result<(Vec<String>, Vec<Vec<Option<f64>>>), String> {
    if wells.len() != solution.per_well_per_channel_rates.len() {
        return Err("Solution per_well_per_channel_rates length mismatch with wells.".into());
    }

    // 1) Collect all species names that can appear as destinations.
    let mut names: BTreeSet<String> = BTreeSet::new();
    for w in wells {
        names.insert(w.well_name.clone());
    }
    for w in wells {
        for ch in &w.channels {
            if let Some(j) = ch.connected_well_index {
                if let Some(wj) = wells.get(j) {
                    names.insert(wj.well_name.clone());
                }
            } else if let Some((_from, to)) = ch.name.split_once("_to_") {
                names.insert(to.to_string());
            }
        }
    }

    let species_names: Vec<String> = names.into_iter().collect();
    let index_of: HashMap<String, usize> = species_names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.clone(), i))
        .collect();

    // 2) Fill matrix.
    let n = species_names.len();
    let mut matrix: Vec<Vec<Option<f64>>> = vec![vec![None; n]; n];

    for (well_index, well) in wells.iter().enumerate() {
        let from_i = *index_of
            .get(&well.well_name)
            .ok_or_else(|| "Missing well name in species index map.".to_string())?;

        if well.channels.len() != solution.per_well_per_channel_rates[well_index].len() {
            return Err(format!(
                "Solution channel rate count mismatch for well '{}'.",
                well.well_name
            ));
        }

        for (channel_index, ch) in well.channels.iter().enumerate() {
            let k = solution.per_well_per_channel_rates[well_index][channel_index];
            if !(k.is_finite() && k > 0.0) {
                continue;
            }

            let to_name = if let Some(j) = ch.connected_well_index {
                wells
                    .get(j)
                    .ok_or_else(|| "Invalid connected_well_index in channel.".to_string())?
                    .well_name
                    .clone()
            } else if let Some((_from, to)) = ch.name.split_once("_to_") {
                to.to_string()
            } else {
                ch.name.clone()
            };

            let to_j = match index_of.get(&to_name) {
                Some(v) => *v,
                None => continue,
            };
            matrix[from_i][to_j] = Some(k);
        }
    }

    Ok((species_names, matrix))
}
