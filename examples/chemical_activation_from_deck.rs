//! Master-equation calculation for an input deck in the MESS input format.
//!
//!   cargo run --release --example chemical_activation_from_deck -- [deck.inp] [reactant]
//!       [--method steady-state|cse] [--steady-state intermediate|final|both] [--barrier-kt X]
//!       [--eigen-solver inverse|full|lapack] [--sum-rule-tolerance X] [--tunneling exact-eckart|mess-eckart]
//!
//! Two solution methods (`masterequation::solution_method`), chosen in the `MarXus ... End` block of the deck
//! header (`mess_input.rs`) or on the command line, which overrides the deck:
//!
//! 1. Steady state, J N = F (Gonzalez-Garcia, Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010) [GO10],
//!    eqs. 7, 8), in two versions (GO10 Sec. 3.2):
//!    - intermediate steady state: absorbing barrier X k_BT below the lowest threshold of each well;
//!    - final steady state: no absorbing barrier. It includes the thermal rate coefficients of the same J
//!      (`chemical_activation_eigen.rs`): k_uni(T, p) = sum_j k_j^th + k_c[D], the specific rate
//!      coefficients averaged over the normalized eigenvector of the lowest eigenvalue, the thermal
//!      steady-state population (GO10, text after eq. 12), and lambda_1 (eq. 12) beside it. A deviation
//!      between the two above the tolerance is reported as a warning (on stderr and in the table); the
//!      output explains this with the reference. For a single well formed through one entrance channel,
//!      the association rate coefficient follows by detailed balance, k(R -> W, T, p) = k_uni(T, p) K(T),
//!      with K = k_inf,assoc/k_inf,diss of the high-pressure rate coefficients of the entrance channel.
//! 2. Phenomenological rate coefficients from the chemically significant eigenvalues (CSE: Miller,
//!    Klippenstein, J. Phys. Chem. A 110, 10528 (2006); Georgievskii et al., J. Phys. Chem. A 117, 12146
//!    (2013); `chemically_significant_eigenvalues.rs`). It needs all eigenpairs: LAPACK unless the full
//!    Householder/QL decomposition is asked for; inverse iteration is refused.
//!
//! Options (deck keyword in the MarXus block in brackets):
//! --method               `steady-state` (default) or `cse` [Method SteadyState | CSE]
//! --steady-state         versions of the steady-state method: `intermediate`, `final` or `both` (default)
//!                        [SteadyState Intermediate | Final | Both]
//! --barrier-kt X         absorbing barrier of the intermediate steady state (default 10; a smaller value for
//!                        wells that are shallow compared with 10 k_BT plus their thermal width; the
//!                        stabilization then depends on this choice) [AbsorbingBarrierBelowThreshold[kT] X]
//! --eigen-solver         thermal eigenpair of the final steady state by inverse iteration with the banded
//!                        Cholesky factor (`inverse`, default), by the full Householder/QL decomposition
//!                        (`full`), or by LAPACK DSYEVD (`lapack`, needs the `openblas` build feature, on by
//!                        default); for the CSE method `lapack` (default) or `full`
//!                        [EigenSolver InverseIteration | FullDecomposition | Lapack]
//! --sum-rule-tolerance X relative deviation |lambda_1 - k_uni| / k_uni of the thermal eigenpair above which a
//!                        warning is printed (default 1.5e-2); the result is kept [SumRuleTolerance X]
//! --tunneling            transmission model of `Tunneling Eckart` blocks: the exact Eckart probability
//!                        (Miller 1979, eq. 8; `exact-eckart`, default) or the MESS semiclassical model
//!                        (`mess-eckart`, `tunneling::mess_eckart_tunneling`), to reproduce MESS results
//! A setting that the selected solution does not use is reported as a note, not refused.
//!
//! The deck's `Reactant` (a Bimolecular species) forms the wells through its barriers; the source is
//! thermal (Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014), eqs. 7 and 9). The results table is
//! written for each steady state; see `chemical_activation_driver.rs` for the columns.
//! For the intermediate steady state the bimolecular rate coefficients k(R -> X) = k_inf Phi_X follow
//! (Pilling, Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003), eq. 44): k_inf is the high-pressure rate
//! coefficient of the Reactant forming the wells, X a stabilized well, a product channel or the
//! bimolecular sink of a well.

use MarXus::masterequation::chemical_activation_driver::{
    run_chemical_activation, run_phenomenological_rates, run_thermal_rate_coefficients, write_phenomenological_tables,
    write_results_table, write_thermal_table, ChemicalActivationRun, SourceSpecification,
};
use MarXus::masterequation::chemical_activation_from_mess_input::{
    chemical_activation_model_from_mess, EckartTunnelingModel, MessNetworkSettings,
};
use MarXus::masterequation::chemical_activation_network::{
    AbsorbingBarrier, ChannelDestination, ChemicalActivationOptions, SteadyState,
};
use MarXus::masterequation::chemical_activation_steady_state::LinearSolver;
use MarXus::masterequation::mess_input::parse_mess_input_file;
use MarXus::masterequation::solution_method::{
    eigen_solver_from_keyword, Solution, SolutionMethod, SolutionSettings, SteadyStateVersions,
};

fn main() -> Result<(), String> {
    // Positional arguments: deck, reactant; options: see the module documentation.
    let mut positional = Vec::new();
    let mut command_line = SolutionSettings::default();
    let mut settings = MessNetworkSettings::default();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let mut value = || args.next().ok_or(format!("{arg} needs a value"));
        let number = |v: String| v.parse::<f64>().map_err(|_| format!("{arg}: invalid number '{v}'"));
        match arg.as_str() {
            "--method" => command_line.method = Some(SolutionMethod::from_keyword(&value()?)?),
            "--steady-state" => command_line.steady_state = Some(SteadyStateVersions::from_keyword(&value()?)?),
            "--barrier-kt" => command_line.absorbing_barrier_kt = Some(number(value()?)?),
            "--eigen-solver" => command_line.eigen_solver = Some(eigen_solver_from_keyword(&value()?)?),
            "--sum-rule-tolerance" => command_line.sum_rule_tolerance = Some(number(value()?)?),
            "--tunneling" => {
                settings.eckart_tunneling = match value()?.as_str() {
                    "exact-eckart" => EckartTunnelingModel::Exact,
                    "mess-eckart" => EckartTunnelingModel::Mess,
                    other => return Err(format!("--tunneling: unknown model '{other}' (exact-eckart, mess-eckart)")),
                }
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
    let resolved = deck.global.solution.overridden_by(&command_line).resolve()?;
    // The steady states to solve, the settings of the thermal eigenpair of the final steady state, and the
    // eigen-solver of the CSE method.
    let (solutions, thermal, cse) = match resolved.solution {
        Solution::SteadyState { intermediate_absorbing_barrier_kt, final_steady_state } => {
            let mut solutions = Vec::new();
            if let Some(kt) = intermediate_absorbing_barrier_kt {
                let barrier = AbsorbingBarrier::BelowLowestThreshold { kt_multiple: kt };
                solutions.push(("intermediate steady state", SteadyState::Intermediate { barrier }));
            }
            if final_steady_state.is_some() {
                solutions.push(("final steady state", SteadyState::Final));
            }
            (solutions, final_steady_state, None)
        }
        Solution::ChemicallySignificantEigenvalues { eigen_solver } => (Vec::new(), None, Some(eigen_solver)),
    };
    let model = chemical_activation_model_from_mess(&deck, &settings)?;
    if model.entrance_channels.is_empty() {
        return Err(format!("{path}: no barrier connects the Reactant of the deck to a well."));
    }

    let network = &model.network;
    println!("# deck {path}: grain {:.3} cm-1", network.grain_width_cm1);
    let method = match resolved.solution {
        Solution::SteadyState { intermediate_absorbing_barrier_kt, final_steady_state } => {
            let mut versions = Vec::new();
            if let Some(kt) = intermediate_absorbing_barrier_kt {
                versions.push(format!("intermediate (absorbing barrier {kt} k_BT below the lowest threshold)"));
            }
            if let Some(t) = final_steady_state {
                versions.push(format!(
                    "final, with its thermal eigenpair by {:?} (sum-rule tolerance {:e})",
                    t.eigen_solver, t.sum_rule_tolerance
                ));
            }
            format!("steady state J N = F (GO10 eqs. 7, 8): {}", versions.join("; "))
        }
        Solution::ChemicallySignificantEigenvalues { eigen_solver } => {
            format!("chemically significant eigenvalues (MK06; G13), all eigenpairs by {eigen_solver:?}")
        }
    };
    println!("# solution method: {method}");
    for note in &resolved.unused_settings {
        println!("# note: {note}");
        eprintln!("note: {note}");
    }
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
            SteadyState::Intermediate { barrier: AbsorbingBarrier::BelowLowestThreshold { kt_multiple } } => {
                println!("\n# {label} (absorbing barrier {kt_multiple} k_BT below the lowest threshold)")
            }
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

    // Thermal rate coefficients of the final steady state: lowest eigenpair of its J (GO10 eq. 12).
    if let Some(thermal) = thermal {
        println!(
            "\n# thermal rate coefficients of the final steady state: lowest eigenpair of J (GO10 eq. 12) by {:?}; \
             sum-rule tolerance {:e}",
            thermal.eigen_solver, thermal.sum_rule_tolerance
        );
        let mut results = Vec::new();
        for &t in &model.temperatures_kelvin {
            for &p in &model.pressures_torr {
                match run_thermal_rate_coefficients(
                    network,
                    &[t],
                    &[p],
                    model.collision_model,
                    thermal.eigen_solver,
                    thermal.sum_rule_tolerance,
                ) {
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
            println!("\n# bimolecular rate coefficients of {reactant} [cm3/s] by detailed balance from the thermal k_uni of the final steady state");
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
    if let Some(solver) = cse {
        println!("\n# CSE analysis ({solver:?})");
        let capture = |t: f64| model.entrance_high_pressure_rate.as_ref().map_or(f64::NAN, |k| k.rate_cm3_s(t));
        let reactant = deck.global.reactant_name.as_deref().filter(|_| model.entrance_high_pressure_rate.is_some());
        for &t in &model.temperatures_kelvin {
            for &p in &model.pressures_torr {
                match run_phenomenological_rates(network, &[t], &[p], model.collision_model, solver, reactant.map(|r| (r, &capture as &dyn Fn(f64) -> f64))) {
                    Ok(results) => {
                        for warning in results.iter().flat_map(|r| r.rates.warnings.iter()) {
                            eprintln!("warning: T = {t} K, p = {p} Torr: {warning}");
                        }
                        write_phenomenological_tables(&results, &mut stdout).map_err(|e| e.to_string())?;
                    }
                    Err(e) => println!("# not available: {e}"),
                }
            }
        }
    }
    Ok(())
}
