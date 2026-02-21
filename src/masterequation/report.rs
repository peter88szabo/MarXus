//! General reporting helpers for master-equation calculations.
//!
//! Keep these printers physics-facing: print energies, barriers, T/P conditions,
//! high-pressure-limit reference rates, and steady-state master-equation outputs.

use crate::masterequation::high_pressure_limit::HighPressureThermoResult;

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
    println!("T = {:.1} K, P = {:.3} atm", temperature_kelvin, pressure_atm);
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
    println!("Tunneling model: {} (with stability guard)", tunneling_model);
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
    println!("______________________________________________________________________________________");
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
