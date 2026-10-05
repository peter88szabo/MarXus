//! Steady-state chemical-activation calculation for an input deck in the MESS input format.
//!
//!   cargo run --release --example chemical_activation_from_deck [deck.inp]
//!
//! The deck's `Reactant` (a Bimolecular species) forms the wells through its barriers; the source is
//! thermal (Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014), eqs. 7 and 9). The results table is
//! written for the intermediate steady state (absorbing barrier 10 k_BT below the lowest threshold of
//! each well) and for the final steady state; see `chemical_activation_driver.rs` for the columns.

use MarXus::masterequation::chemical_activation_driver::{
    run_chemical_activation, write_results_table, ChemicalActivationRun, SourceSpecification,
};
use MarXus::masterequation::chemical_activation_from_mess_input::{chemical_activation_model_from_mess, MessNetworkSettings};
use MarXus::masterequation::chemical_activation_network::{AbsorbingBarrier, ChemicalActivationOptions, SteadyState};
use MarXus::masterequation::chemical_activation_steady_state::LinearSolver;
use MarXus::masterequation::mess_input::parse_mess_input_file;

fn main() -> Result<(), String> {
    let path = std::env::args().nth(1).unwrap_or_else(|| "examples/c2h3_chemical_activation.inp".to_string());
    let deck = parse_mess_input_file(&path)?;
    let model = chemical_activation_model_from_mess(&deck, &MessNetworkSettings::default())?;
    if model.entrance_channels.is_empty() {
        return Err(format!("{path}: no barrier connects the Reactant of the deck to a well."));
    }

    let network = &model.network;
    println!("# deck {path}: grain {:.3} cm-1", network.grain_width_cm1);
    for well in &network.wells {
        println!(
            "#   well {:<8} grains {:>6} (absolute {} .. {}), channels: {}",
            well.name,
            well.grain_count(),
            well.bottom_offset_grains,
            well.bottom_offset_grains + well.grain_count() as isize - 1,
            well.channels.iter().map(|c| c.name.as_str()).collect::<Vec<_>>().join(", ")
        );
    }

    let mut stdout = std::io::stdout();
    for (label, steady_state) in [
        ("intermediate steady state", SteadyState::Intermediate { barrier: AbsorbingBarrier::default() }),
        ("final steady state", SteadyState::Final),
    ] {
        let run = ChemicalActivationRun {
            temperatures_kelvin: model.temperatures_kelvin.clone(),
            pressures_torr: model.pressures_torr.clone(),
            options: ChemicalActivationOptions { collision_model: model.collision_model, steady_state },
            solver: LinearSolver::BandedCholesky,
            source: SourceSpecification::ThermalEntrance { channels: model.entrance_channels.clone() },
            tolerance: 1e-8,
        };
        println!("\n# {label}");
        let results = run_chemical_activation(network, &run)?;
        write_results_table(network, &results, &mut stdout).map_err(|e| e.to_string())?;
    }
    Ok(())
}
