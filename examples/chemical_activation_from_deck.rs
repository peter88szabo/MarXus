//! Steady-state chemical-activation calculation for an input deck in the MESS input format.
//!
//!   cargo run --release --example chemical_activation_from_deck -- [deck.inp] [reactant]
//!       [--barrier-kt X] [--steady-state intermediate|final|eigenvalue|both|all] [--eigen-solver inverse|full|lapack]
//!       [--sum-rule-tolerance X]
//!
//! --barrier-kt X   absorbing barrier X k_BT below the lowest threshold of each well (default 10; a
//!                  smaller value for wells that are shallow compared with 10 k_BT plus their thermal width;
//!                  the stabilization then depends on this choice)
//! --steady-state   which solution to compute: the two steady states (`both`, default), one of them, the
//!                  eigenvalue analysis (`eigenvalue`), or all three (`all`)
//! --eigen-solver   eigenvalue analysis by inverse iteration with the banded Cholesky factor (`inverse`,
//!                  default), by the full Householder/QL decomposition (`full`), or by the full decomposition
//!                  of LAPACK DSYEVD (`lapack`, needs the `openblas` build feature, on by default)
//! --sum-rule-tolerance X   relative deviation |lambda_1 - k_uni| / k_uni of the eigenvalue analysis above which
//!                  a warning is printed (default 1.5e-2); the result is kept
//!
//! Eigenvalue analysis (Olzmann's solution method, `chemical_activation_eigen.rs`) of J without absorbing
//! barrier: k_uni(T, p) = sum_j k_j^th + k_c[D], the specific rate coefficients averaged over the
//! normalized thermal eigenvector (Gonzalez-Garcia, Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010),
//! text after eq. 12), equal to the lowest eigenvalue lambda_1 (eq. 12), which is printed beside it; a
//! deviation between the two above the tolerance is reported as a warning (on stderr and in the table).
//! The output explains this with the reference. For a single well formed through one entrance channel, the
//! association rate coefficient follows by detailed balance, k(R -> W, T, p) = k_uni(T, p) K(T), with the
//! equilibrium constant K = k_inf,assoc/k_inf,diss of the high-pressure rate coefficients of the entrance
//! channel; no absorbing barrier is involved.
//!
//! The deck's `Reactant` (a Bimolecular species) forms the wells through its barriers; the source is
//! thermal (Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014), eqs. 7 and 9). The results table is
//! written for the intermediate steady state (absorbing barrier 10 k_BT below the lowest threshold of
//! each well) and for the final steady state; see `chemical_activation_driver.rs` for the columns.
//! For the intermediate steady state the bimolecular rate coefficients k(R -> X) = k_inf Phi_X follow
//! (Pilling, Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003), eq. 44): k_inf is the high-pressure rate
//! coefficient of the Reactant forming the wells, X a stabilized well, a product channel or the
//! bimolecular sink of a well.

use MarXus::masterequation::chemical_activation_driver::{
    run_chemical_activation, run_thermal_rate_coefficients, write_results_table, write_thermal_table,
    ChemicalActivationRun, SourceSpecification,
};
use MarXus::masterequation::chemical_activation_eigen::{EigenSolver, DEFAULT_SUM_RULE_TOLERANCE};
use MarXus::masterequation::chemical_activation_from_mess_input::{chemical_activation_model_from_mess, MessNetworkSettings};
use MarXus::masterequation::chemical_activation_network::{
    AbsorbingBarrier, ChannelDestination, ChemicalActivationOptions, SteadyState,
};
use MarXus::masterequation::chemical_activation_steady_state::LinearSolver;
use MarXus::masterequation::mess_input::parse_mess_input_file;

fn main() -> Result<(), String> {
    // Positional arguments: deck, reactant; options: see the module documentation.
    let mut positional = Vec::new();
    let mut barrier_kt = 10.0;
    let mut selection = "both".to_string();
    let mut eigen_solver = EigenSolver::default();
    let mut sum_rule_tolerance = DEFAULT_SUM_RULE_TOLERANCE;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--barrier-kt" => {
                let value = args.next().ok_or("--barrier-kt needs a value")?;
                barrier_kt = value.parse::<f64>().map_err(|_| format!("--barrier-kt: invalid number '{value}'"))?;
            }
            "--steady-state" => selection = args.next().ok_or("--steady-state needs a value")?,
            "--eigen-solver" => {
                eigen_solver = match args.next().ok_or("--eigen-solver needs a value")?.as_str() {
                    "inverse" => EigenSolver::InverseIteration,
                    "full" => EigenSolver::FullDecomposition,
                    "lapack" => EigenSolver::FullDecompositionLapack,
                    other => return Err(format!("--eigen-solver: unknown choice '{other}' (inverse, full, lapack)")),
                }
            }
            "--sum-rule-tolerance" => {
                let value = args.next().ok_or("--sum-rule-tolerance needs a value")?;
                sum_rule_tolerance =
                    value.parse::<f64>().map_err(|_| format!("--sum-rule-tolerance: invalid number '{value}'"))?;
            }
            _ => positional.push(arg),
        }
    }
    let path = positional.first().cloned().unwrap_or_else(|| "examples/c2h3_chemical_activation.inp".to_string());
    let mut deck = parse_mess_input_file(&path)?;
    // Optional second positional argument: the Bimolecular species that forms the wells (overrides `Reactant`).
    if let Some(reactant) = positional.get(1) {
        deck.global.reactant_name = Some(reactant.clone());
    }
    let intermediate = ("intermediate steady state", SteadyState::Intermediate {
        barrier: AbsorbingBarrier::BelowLowestThreshold { kt_multiple: barrier_kt },
    });
    let (solutions, eigenvalue_analysis) = match selection.as_str() {
        "both" => (vec![intermediate, ("final steady state", SteadyState::Final)], false),
        "intermediate" => (vec![intermediate], false),
        "final" => (vec![("final steady state", SteadyState::Final)], false),
        "eigenvalue" => (vec![], true),
        "all" => (vec![intermediate, ("final steady state", SteadyState::Final)], true),
        other => return Err(format!("--steady-state: unknown choice '{other}' (intermediate, final, eigenvalue, both, all)")),
    };
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
    for (label, steady_state) in solutions {
        let is_intermediate = matches!(steady_state, SteadyState::Intermediate { .. });
        let run = ChemicalActivationRun {
            temperatures_kelvin: model.temperatures_kelvin.clone(),
            pressures_torr: model.pressures_torr.clone(),
            options: ChemicalActivationOptions { collision_model: model.collision_model, steady_state: steady_state.clone() },
            solver: LinearSolver::BandedCholesky,
            source: SourceSpecification::ThermalEntrance { channels: model.entrance_channels.clone() },
            tolerance: 1e-8,
        };
        match &steady_state {
            SteadyState::Intermediate { .. } => println!("\n# {label} (absorbing barrier {barrier_kt} k_BT below the lowest threshold)"),
            _ => println!("\n# {label}"),
        }
        // Each (T, p) separately, so that a condition without a valid steady state (e.g. the final steady
        // state of a deep well without a sink at low temperature, or an absorbing barrier below the well
        // bottom at high temperature) is reported and the other conditions are still computed.
        let mut results = Vec::new();
        for &t in &run.temperatures_kelvin {
            for &p in &run.pressures_torr {
                let single = ChemicalActivationRun { temperatures_kelvin: vec![t], pressures_torr: vec![p], ..run.clone() };
                match run_chemical_activation(network, &single) {
                    Ok(mut r) => results.append(&mut r),
                    Err(e) => println!("# not available: {e}"),
                }
            }
        }
        if results.is_empty() {
            continue;
        }
        write_results_table(network, &results, &mut stdout).map_err(|e| e.to_string())?;

        // Bimolecular rate coefficients of the Reactant, k(R -> X) = k_inf Phi_X (cm3/s).
        if let (true, Some(k_inf), Some(reactant)) =
            (is_intermediate, &model.entrance_high_pressure_rate, deck.global.reactant_name.as_ref())
        {
            println!("\n# bimolecular rate coefficients of {reactant} [cm3/s], {label}");
            let mut header = vec!["T[K]".to_string(), "P[Torr]".to_string(), "k_inf".to_string()];
            for well in &network.wells {
                header.push(format!("k({reactant}->{})", well.name));
            }
            for well in &network.wells {
                for channel in &well.channels {
                    if let ChannelDestination::Products { name } = &channel.destination {
                        if name != reactant {
                            header.push(format!("k({reactant}->{name} via {})", channel.name));
                        }
                    }
                }
                if well.bimolecular_sink_s_inv > 0.0 {
                    header.push(format!("k({reactant}->sink of {})", well.name));
                }
            }
            println!("{}", header.join(","));
            for r in &results {
                let k = k_inf.rate_cm3_s(r.conditions.temperature_kelvin);
                let mut row = vec![
                    format!("{}", r.conditions.temperature_kelvin),
                    format!("{}", r.conditions.pressure_torr),
                    format!("{k:.6e}"),
                ];
                for w in &r.result.wells {
                    row.push(format!("{:.6e}", k * w.stabilization_yield));
                }
                for (w, well) in network.wells.iter().enumerate() {
                    for c in r.result.channels.iter().filter(|c| c.well == w) {
                        if let ChannelDestination::Products { name } = &c.destination {
                            if name != reactant {
                                row.push(format!("{:.6e}", k * c.flux));
                            }
                        }
                    }
                    if well.bimolecular_sink_s_inv > 0.0 {
                        row.push(format!("{:.6e}", k * r.result.wells[w].bimolecular_sink_yield));
                    }
                }
                println!("{}", row.join(","));
            }
        }
    }

    if eigenvalue_analysis {
        println!(
            "\n# eigenvalue analysis ({eigen_solver:?}) of J without absorbing barrier; sum-rule tolerance \
             {sum_rule_tolerance:e}"
        );
        let mut results = Vec::new();
        for &t in &model.temperatures_kelvin {
            for &p in &model.pressures_torr {
                match run_thermal_rate_coefficients(network, &[t], &[p], model.collision_model, eigen_solver, sum_rule_tolerance) {
                    Ok(mut r) => {
                        for warning in r.iter().filter_map(|c| c.thermal.warning.as_ref()) {
                            eprintln!("warning: T = {t} K, p = {p} Torr: {warning}");
                        }
                        results.append(&mut r)
                    }
                    Err(e) => println!("# not available: {e}"),
                }
            }
        }
        if !results.is_empty() {
            write_thermal_table(network, &results, &mut stdout).map_err(|e| e.to_string())?;
        }
        // Association by detailed balance for a single well with one entrance channel.
        if let (1, [(w, c)], Some(k_inf), Some(reactant)) = (
            network.wells.len(),
            model.entrance_channels.as_slice(),
            &model.entrance_high_pressure_rate,
            deck.global.reactant_name.as_ref(),
        ) {
            println!("\n# bimolecular rate coefficients of {reactant} [cm3/s] by detailed balance, eigenvalue analysis");
            println!("T[K],P[Torr],k_inf_assoc[cm3/s],k_inf_diss[1/s],k_uni[1/s],k({reactant}->{})", network.wells[*w].name);
            for r in &results {
                let k_assoc_inf = k_inf.rate_cm3_s(r.conditions.temperature_kelvin);
                let entrance = r.thermal.channels.iter().find(|ch| ch.well == *w && ch.channel == *c).unwrap();
                let k_diss_inf = entrance.high_pressure_rate_s_inv;
                println!(
                    "{},{},{:.6e},{:.6e},{:.6e},{:.6e}",
                    r.conditions.temperature_kelvin,
                    r.conditions.pressure_torr,
                    k_assoc_inf,
                    k_diss_inf,
                    r.thermal.k_uni_s_inv,
                    r.thermal.k_uni_s_inv * k_assoc_inf / k_diss_inf
                );
            }
        }
    }
    Ok(())
}
