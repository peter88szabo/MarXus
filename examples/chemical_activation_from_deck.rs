//! Master-equation calculation for an input deck in the MESS input format.
//!
//!   cargo run --release --example chemical_activation_from_deck -- [deck.inp] [reactant]
//!       [--method steady-state|cse] [--steady-state intermediate|final|both] [--barrier-kt X]
//!       [--eigen-solver inverse|full|lapack] [--sum-rule-tolerance X] [--tunneling exact-eckart|mess-eckart]
//!       [--csv FILE] [--ncore N] [--integrator rodas4|rodas3|ros4|ros3|ros2] [--initial pulse|continuous]
//!       [--time-range T1 T2] [--times-per-decade N] [--integration-tolerance X]
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
//! --method               `steady-state` (default), `cse` or `time-integration` [Method SteadyState | CSE | TimeIntegration]
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
//! --csv FILE             write the machine-readable tables (CSV with titled blocks) to FILE, and every quantity
//!                        group of the report (yields, rate coefficients) to FILE_tables.csv
//! --integrator           time integration: Rosenbrock method, default rodas4 [Integrator]
//! --initial              time integration: `pulse` (default; N(0) = F) or `continuous` (formation R F from N = 0)
//!                        [InitialState]
//! --time-range T1 T2     time integration: first and last output time in s, default 1e-12 1e2 [TimeRange[s]]
//! --times-per-decade N   time integration: output times per decade, default 4 [TimesPerDecade]
//! --integration-tolerance X   time integration: relative tolerance, default 1e-6 [IntegrationTolerance]
//! --ncore N              number of cores: the conditions (T, p) are computed in parallel (rayon), in batches of
//!                        up to N at a time; overrides NCores of the deck; without either, RAYON_NUM_THREADS,
//!                        otherwise all logical cores [NCores]
//! A setting that the selected solution does not use is reported as a note, not refused.
//!
//! Output: a human-readable report on stdout (`report_sections.rs`, `report_tables.rs`):
//! - RUN SETTINGS (method, solvers, conditions, grid, collisions, tunneling, source) and CHEMICAL NETWORK
//!   (wells, channels with their rate models, sinks), then the energetics in kcal/mol relative to the Reactant;
//! - for every solution, its quantities (yields in %, rate coefficients) as tables by temperature (rows:
//!   pressure), by pressure (rows: temperature) and temperature-pressure tables (one per quantity);
//! - steady state: the yields and k_ca of each steady state, the thermal rate coefficients and branching of the
//!   final steady state, the thermal fate of each well, and chemical activation and thermal reaction separately
//!   (prompt, through the stabilized wells) and together (final steady state);
//! - CSE: species-to-species tables per condition, the rate coefficients from every species, the prompt branching
//!   of the reactant, the thermal fate of each well and the long-time yields (direct + through the wells).
//!
//! The deck's `Reactant` (a Bimolecular species) forms the wells through its barriers; the source is
//! thermal (Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014), eqs. 7 and 9). The machine-readable table
//! of each steady state is described in `chemical_activation_driver.rs`.
//! For the intermediate steady state the bimolecular rate coefficients k(R -> X) = k_inf Phi_X follow
//! (Pilling, Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003), eq. 44): k_inf is the high-pressure rate
//! coefficient of the Reactant forming the wells, X a stabilized well, a product channel or the
//! bimolecular sink of a well.

use std::io::Write;

use MarXus::constants::KB_CM;
use MarXus::masterequation::chemical_activation_driver::{
    run_chemical_activation, run_phenomenological_rates, run_thermal_rate_coefficients,
    run_thermal_well_fates, write_phenomenological_tables, write_results_table,
    write_thermal_table, ChemicalActivationRun, SourceSpecification, ThermalConditionResult,
};
use MarXus::masterequation::chemical_activation_from_mess_input::{
    chemical_activation_model_from_mess, EckartTunnelingModel, MessNetworkSettings,
};
use MarXus::masterequation::chemical_activation_network::{
    AbsorbingBarrier, ChannelDestination, ChemicalActivationOptions, CollisionModel, Conditions,
    SteadyState,
};
use MarXus::masterequation::chemical_activation_operator::low_energy_reservoirs;
use MarXus::masterequation::chemical_activation_sources::thermal_entrance_source;
use MarXus::masterequation::chemical_activation_steady_state::LinearSolver;
use MarXus::masterequation::direct_time_integration::{
    integrate_master_equation, log_spaced_times, InitialState, TimeIntegrationSettings,
};
use MarXus::masterequation::mess_input::parse_mess_input_file;
use MarXus::masterequation::parallel_conditions::ConditionPool;
use MarXus::masterequation::report_sections::{
    cse_groups, steady_state_groups, thermal_groups, time_integration_groups, well_fate_groups,
    write_cse_species_tables, write_energetics, write_groups, write_network_summary,
    write_time_evolution_tables,
};
use MarXus::masterequation::report_tables::{write_groups_csv, Quantity, QuantityGroup};
use MarXus::numeric::lapack_interface::set_blas_threads;

/// The groups with `section: ` before their titles (the titles of the tables file).
fn prefixed(section: &str, groups: &[QuantityGroup]) -> Vec<QuantityGroup> {
    groups
        .iter()
        .map(|g| QuantityGroup {
            title: format!("{section}: {}", g.title),
            quantities: g.quantities.clone(),
        })
        .collect()
}
use MarXus::masterequation::solution_method::{
    eigen_solver_from_keyword, initial_state_from_keyword, integrator_from_keyword, Solution,
    SolutionMethod, SolutionSettings, STEADY_STATE_KEYWORD_REPLACED,
};

fn main() -> Result<(), String> {
    // Positional arguments: deck, reactant; options: see the module documentation.
    let mut positional = Vec::new();
    let mut command_line = SolutionSettings::default();
    let mut settings = MessNetworkSettings::default();
    let mut csv_path: Option<String> = None;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let mut value = || args.next().ok_or(format!("{arg} needs a value"));
        let number = |v: String| v.parse::<f64>().map_err(|_| format!("{arg}: invalid number '{v}'"));
        match arg.as_str() {
            "--method" => command_line.method = Some(SolutionMethod::from_keyword(&value()?)?),
            "--steady-state" => {
                return Err(format!(
                    "--steady-state no longer exists: {STEADY_STATE_KEYWORD_REPLACED}."
                ))
            }
            "--barrier-kt" => command_line.absorbing_barrier_kt = Some(number(value()?)?),
            "--eigen-solver" => command_line.eigen_solver = Some(eigen_solver_from_keyword(&value()?)?),
            "--sum-rule-tolerance" => command_line.sum_rule_tolerance = Some(number(value()?)?),
            "--integrator" => command_line.integrator = Some(integrator_from_keyword(&value()?)?),
            "--initial" => {
                command_line.initial_state = Some(initial_state_from_keyword(&value()?)?)
            }
            "--time-range" => {
                let first = number(value()?)?;
                let last = number(value()?)?;
                command_line.time_range_s = Some((first, last));
            }
            "--times-per-decade" => {
                let v = value()?;
                command_line.times_per_decade = Some(
                    v.parse::<usize>()
                        .map_err(|_| format!("--times-per-decade: invalid number '{v}'"))?,
                );
            }
            "--integration-tolerance" => {
                command_line.integration_tolerance = Some(number(value()?)?)
            }
            "--csv" => csv_path = Some(value()?),
            "--ncore" => {
                let v = value()?;
                command_line.cores = Some(
                    v.parse::<usize>()
                        .map_err(|_| format!("--ncore: invalid number '{v}'"))?,
                );
            }
            "--tunneling" => {
                settings.eckart_tunneling = match value()?.as_str() {
                    "exact-eckart" => EckartTunnelingModel::Exact,
                    "mess-eckart" => EckartTunnelingModel::Mess,
                    other => return Err(format!("--tunneling: unknown model '{other}' (exact-eckart, mess-eckart)")),
                }
            }
            other if other.starts_with("--") => {
                return Err(format!(
                    "unknown option '{other}' (--method, --barrier-kt, --eigen-solver, --sum-rule-tolerance, \
                     --integrator, --initial, --time-range, --times-per-decade, --integration-tolerance, --csv, \
                     --ncore, --tunneling)"
                ))
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
    let merged = deck.global.solution.overridden_by(&command_line);
    let resolved = merged.resolve()?;
    // Where the number of cores comes from (shown in RUN SETTINGS).
    let cores_source = if command_line.cores.is_some() {
        "--ncore"
    } else if deck.global.solution.cores.is_some() {
        "NCores in the deck"
    } else {
        "default: RAYON_NUM_THREADS, otherwise all logical cores"
    };
    // One solver per run: the steady state to solve (intermediate for SteadyStateAbsorbingBarrier, final for
    // SteadyStateOlzmann, with its thermal eigenpair), the eigen-solver of CSE, or the time-integration plan.
    let (steady_state, thermal, cse, time_integration) = match resolved.solution {
        Solution::SteadyStateAbsorbingBarrier {
            absorbing_barrier_kt,
        } => {
            let barrier = AbsorbingBarrier::BelowLowestThreshold {
                kt_multiple: absorbing_barrier_kt,
            };
            (
                Some((
                    "intermediate steady state",
                    SteadyState::Intermediate { barrier },
                )),
                None,
                None,
                None,
            )
        }
        Solution::SteadyStateOlzmann(thermal) => (
            Some(("final steady state", SteadyState::Final)),
            Some(thermal),
            None,
            None,
        ),
        Solution::ChemicallySignificantEigenvalues { eigen_solver } => {
            (None, None, Some(eigen_solver), None)
        }
        Solution::TimeIntegration(plan) => (None, None, None, Some(plan)),
    };
    let model = chemical_activation_model_from_mess(&deck, &settings)?;
    if model.entrance_channels.is_empty() {
        return Err(format!("{path}: no barrier connects the Reactant of the deck to a well."));
    }

    let network = &model.network;
    let temperatures = model.temperatures_kelvin.clone();
    let pressures = model.pressures_torr.clone();
    let reactant = deck.global.reactant_name.clone();
    // The conditions (T, p) are independent master equations, computed in parallel (`parallel_conditions.rs`).
    // OpenBLAS gets the threads that the conditions leave free, so that workers x BLAS threads stay within
    // the requested number of cores.
    let pool = ConditionPool::new(merged.cores)?;
    let condition_count = temperatures.len() * pressures.len();
    let concurrent = pool.threads().min(condition_count).max(1);
    let blas_threads = set_blas_threads((pool.threads() / concurrent).max(1));
    let k_inf = |t: f64| {
        model
            .entrance_high_pressure_rate
            .as_ref()
            .map_or(f64::NAN, |k| k.rate_cm3_s(t))
    };
    let method = match resolved.solution {
        Solution::SteadyStateAbsorbingBarrier { absorbing_barrier_kt } => format!(
            "SteadyStateAbsorbingBarrier: intermediate steady state, absorbing barrier {absorbing_barrier_kt} k_BT below the lowest threshold"
        ),
        Solution::SteadyStateOlzmann(t) => format!(
            "SteadyStateOlzmann: final steady state, thermal eigenpair by {:?} (sum-rule tolerance {:e})",
            t.eigen_solver, t.sum_rule_tolerance
        ),
        Solution::ChemicallySignificantEigenvalues { eigen_solver } => {
            format!("CSE: chemically significant eigenvalues (MK06; G13), all eigenpairs by {eigen_solver:?}")
        }
        Solution::TimeIntegration(plan) => {
            format!("TimeIntegration: direct time integration, {:?}, {:?}", plan.integrator, plan.initial_state)
        }
    };
    let io = |e: std::io::Error| e.to_string();
    // Human-readable report on stdout; the machine-readable blocks (CSV tables with titles) go to the file of
    // `--csv`.
    let mut report = std::io::stdout().lock();
    let mut machine: Vec<u8> = Vec::new();
    // Every quantity group of the report, with its section in the title, for FILE_tables.csv.
    let mut tables: Vec<QuantityGroup> = Vec::new();

    writeln!(
        machine,
        "# deck {path}: grain {:.3} cm-1",
        network.grain_width_cm1
    )
    .map_err(io)?;
    writeln!(machine, "# solution method: {method}").map_err(io)?;
    for note in &resolved.unused_settings {
        writeln!(machine, "# note: {note}").map_err(io)?;
        eprintln!("note: {note}");
    }
    for well in &network.wells {
        writeln!(
            machine,
            "#   well {:<8} grains {:>6} (absolute {} .. {}), channels: {}",
            well.name,
            well.grain_count(),
            well.bottom_offset_grains,
            well.bottom_offset_grains + well.grain_count() as isize - 1,
            well.channels
                .iter()
                .map(|c| c.name.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        )
        .map_err(io)?;
    }

    // ---- Report header: run settings, chemical network, energetics, units.
    let rule = "=".repeat(100);
    let list = |xs: &[f64]| {
        xs.iter()
            .map(|x| format!("{x}"))
            .collect::<Vec<_>>()
            .join(", ")
    };
    let (method_lines, solver_lines): (Vec<String>, Vec<String>) = match resolved.solution {
        Solution::SteadyStateAbsorbingBarrier { absorbing_barrier_kt } => (
            vec![
                "STEADY STATE, ABSORBING BARRIER (intermediate steady state; steady-state family)".to_string(),
                format!("J N = F with a lower absorbing barrier {absorbing_barrier_kt} k_BT below the lowest threshold of each well"),
                "(Gonzalez-Garcia, Olzmann, PCCP 12, 12290 (2010), eqs. 7, 8 and Sec. 3.2; Olzmann, PCCP 4, 3614 (2002))".to_string(),
            ],
            vec!["linear systems: banded Cholesky factorization of the symmetrized J, relative residual <= 1e-8".to_string()],
        ),
        Solution::SteadyStateOlzmann(th) => (
            vec![
                "STEADY STATE, OLZMANN (final steady state; steady-state family)".to_string(),
                "J N = R F without an absorbing barrier (Gonzalez-Garcia, Olzmann, PCCP 12, 12290 (2010), eqs. 7, 8)".to_string(),
                "with the thermal rate coefficients of the lowest eigenpair of the same J (eq. 12)".to_string(),
                "and the thermal fates of the wells (Boltzmann distribution of one well as source)".to_string(),
            ],
            vec![
                "linear systems: banded Cholesky factorization of the symmetrized J, relative residual <= 1e-8".to_string(),
                format!("thermal eigenpair: {:?}, sum-rule tolerance {:e} (warning only)", th.eigen_solver, th.sum_rule_tolerance),
            ],
        ),
        Solution::ChemicallySignificantEigenvalues { eigen_solver } => (
            vec![
                "CHEMICALLY SIGNIFICANT EIGENVALUES (CSE): phenomenological rate coefficients".to_string(),
                "(Miller, Klippenstein, JPCA 110, 10528 (2006); Georgievskii et al., JPCA 117, 12146 (2013), eqs. 21-30)".to_string(),
            ],
            vec![format!("all eigenpairs of the symmetrized J: {eigen_solver:?}; well-to-well matrix inverted by Gauss-Jordan elimination")],
        ),
        Solution::TimeIntegration(plan) => (
            vec![
                "DIRECT TIME INTEGRATION of the grained populations, dN/dt = R F - J N, and of the yields of every exit".to_string(),
                format!(
                    "initial state: {}",
                    match plan.initial_state {
                        InitialState::Pulse => "pulse N(0) = F (the normalized chemical-activation source)",
                        InitialState::ContinuousFormation => "continuous formation R F, R = 1 s-1, from N(0) = 0",
                    }
                ),
                format!(
                    "output times {:e} .. {:e} s, {} per decade (operator of the final steady state)",
                    plan.time_range_s.0, plan.time_range_s.1, plan.times_per_decade
                ),
            ],
            vec![
                format!(
                    "adaptive L-stable Rosenbrock method {:?} (adapted from KPP; Hairer, Wanner, Solving ODEs II, Sec. IV.7)",
                    plan.integrator
                ),
                format!(
                    "relative tolerance {:e}, absolute {:e}; steps of powers of two, banded Cholesky of s I + S reused",
                    plan.relative_tolerance, plan.absolute_tolerance
                ),
            ],
        ),
    };
    let collisions = match model.collision_model {
        CollisionModel::ExponentialDown { cutoff_in_mean_down } => format!(
            "exponential down, exactly normalized (Robertson (ed.), CCK 43 (2019), eq. 4.16), transitions up to {cutoff_in_mean_down} <dE_down>"
        ),
        CollisionModel::Stepladder => "stepladder, step <dE_down>(T) (Olzmann, Gebhardt, Scherzer, IJCK 23, 825 (1991))".to_string(),
    };
    let tunneling = match settings.eckart_tunneling {
        EckartTunnelingModel::Exact => {
            "exact Eckart transmission (Miller, J. Am. Chem. Soc. 101, 6810 (1979), eq. 8)"
        }
        EckartTunnelingModel::Mess => {
            "MESS semiclassical Eckart model (option for comparison with MESS)"
        }
    };
    let entrances: Vec<String> = model
        .entrance_channels
        .iter()
        .map(|&(w, c)| {
            format!(
                "{} ({} -> {})",
                network.wells[w].channels[c].name,
                reactant.as_deref().unwrap_or("?"),
                network.wells[w].name
            )
        })
        .collect();
    let field =
        |report: &mut std::io::StdoutLock, label: &str, lines: &[String]| -> std::io::Result<()> {
            for (i, line) in lines.iter().enumerate() {
                writeln!(report, " {:<20}{line}", if i == 0 { label } else { "" })?;
            }
            Ok(())
        };
    writeln!(report, "{rule}\n MarXus master equation\n{rule}").map_err(io)?;
    writeln!(report, " RUN SETTINGS\n {}", "-".repeat(12)).map_err(io)?;
    field(&mut report, "Deck:", &[path.clone()]).map_err(io)?;
    field(&mut report, "Method:", &method_lines).map_err(io)?;
    field(&mut report, "Solvers:", &solver_lines).map_err(io)?;
    field(
        &mut report,
        "Parallel:",
        &[
            format!(
                "{} cores ({cores_source}); the {condition_count} conditions (T, p) in batches of {concurrent} at a time (rayon)",
                pool.threads()
            ),
            format!("BLAS threads per LAPACK call: {blas_threads}"),
        ],
    )
    .map_err(io)?;
    field(&mut report, "Temperatures (K):", &[list(&temperatures)]).map_err(io)?;
    field(&mut report, "Pressures (torr):", &[list(&pressures)]).map_err(io)?;
    field(
        &mut report,
        "Energy grid:",
        &[format!(
            "grains of {:.3} cm-1 on a common absolute scale, sums of {} cm-1 counting cells",
            network.grain_width_cm1, settings.cell_width_cm1
        )],
    )
    .map_err(io)?;
    // Low-energy reservoirs: where the normalization of eq. 4.16 fails, the lowest grains of a well form one
    // thermalized state (MESMER manual, Sec. 14.2.1; per temperature, the kernel does not depend on the pressure).
    let mut reservoir_lines = Vec::new();
    for &t in &temperatures {
        let reservoirs = low_energy_reservoirs(network, t, model.collision_model)?;
        if !reservoirs.is_empty() {
            let wells: Vec<String> = reservoirs
                .iter()
                .map(|r| format!("{} {} grains ({:.0} cm-1)", r.well, r.grains, r.grains as f64 * network.grain_width_cm1))
                .collect();
            reservoir_lines.push(format!("  T = {t} K: {}", wells.join(", ")));
        }
    }
    let mut collision_lines = vec![collisions];
    if matches!(model.collision_model, CollisionModel::ExponentialDown { .. }) {
        if reservoir_lines.is_empty() {
            collision_lines.push("eq. 4.16 holds at every grain of every well (no low-energy reservoir)".into());
        } else {
            collision_lines.push(
                "low-energy reservoirs where eq. 4.16 fails (sparse states, Robertson 2019, p. 294): the grains from the"
                    .into(),
            );
            collision_lines.push(
                "  failing one down form one thermalized state (reservoir state, MESMER manual, Sec. 14.2.1):".into(),
            );
            collision_lines.extend(reservoir_lines);
        }
    }
    collision_lines.push("collision frequency: Lennard-Jones (per well, below)".into());
    field(&mut report, "Collisions:", &collision_lines).map_err(io)?;
    field(&mut report, "Tunneling:", &[tunneling.to_string()]).map_err(io)?;
    field(
        &mut report,
        "Source:",
        &[
            format!("thermal reactant {} through {}", reactant.as_deref().unwrap_or("?"), entrances.join(", ")),
            "F(E) ~ rho(E) k(E) exp(-E/kT) of the entrance channels (Pfeifle, Olzmann, IJCK 46, 231 (2014), eqs. 7, 9)".into(),
        ],
    )
    .map_err(io)?;
    for note in &resolved.unused_settings {
        field(&mut report, "Note:", &[note.clone()]).map_err(io)?;
    }
    writeln!(report, "{rule}\n CHEMICAL NETWORK\n {}\n", "-".repeat(16)).map_err(io)?;
    write_network_summary(&mut report, &deck, network, &model.entrance_channels).map_err(io)?;
    write_energetics(&mut report, &deck).map_err(io)?;
    writeln!(
        report,
        "{rule}\nUnits: unimolecular rate coefficients 1/s; bimolecular rate coefficients cm^3/s; yields in %.\n\
         Numbers with six significant digits; *** marks a value that is not available.\n{rule}\n"
    )
    .map_err(io)?;
    let section =
        |report: &mut std::io::StdoutLock, title: &str, text: &str| -> std::io::Result<()> {
            writeln!(
                report,
                "{}\n{title}\n{}\n",
                "_".repeat(100),
                "_".repeat(100)
            )?;
            if !text.is_empty() {
                writeln!(report, "{text}\n")?;
            }
            Ok(())
        };

    // ---- Steady-state method.
    if let Some((label, steady_state)) = steady_state {
        let is_intermediate = matches!(steady_state, SteadyState::Intermediate { .. });
        let run = ChemicalActivationRun {
            temperatures_kelvin: temperatures.clone(),
            pressures_torr: pressures.clone(),
            options: ChemicalActivationOptions {
                collision_model: model.collision_model,
                steady_state: steady_state.clone(),
            },
            solver: LinearSolver::BandedCholesky,
            source: SourceSpecification::ThermalEntrance { channels: model.entrance_channels.clone() },
            tolerance: 1e-8,
        };
        match &steady_state {
            SteadyState::Intermediate {
                barrier: AbsorbingBarrier::BelowLowestThreshold { kt_multiple },
            } => writeln!(
                machine,
                "\n# {label} (absorbing barrier {kt_multiple} k_BT below the lowest threshold)"
            )
            .map_err(io)?,
            _ => writeln!(machine, "\n# {label}").map_err(io)?,
        }
        // Each (T, p) separately, so that a condition without a valid steady state (e.g. the final steady
        // state of a deep well without a sink at low temperature, or an absorbing barrier below the well
        // bottom at high temperature) is reported and the other conditions are still computed.
        let mut results = Vec::new();
        let mut unavailable = Vec::new();
        let per_condition = pool.map_conditions(&temperatures, &pressures, |t, p| {
            run_chemical_activation(
                network,
                &ChemicalActivationRun {
                    temperatures_kelvin: vec![t],
                    pressures_torr: vec![p],
                    ..run.clone()
                },
            )
        });
        for outcome in per_condition {
            match outcome {
                Ok(mut r) => results.append(&mut r),
                Err(e) => {
                    writeln!(machine, "# not available: {e}").map_err(io)?;
                    unavailable.push(e);
                }
            }
        }
        let (title, text) = if is_intermediate {
            (
                "INTERMEDIATE STEADY STATE (absorbing barrier below the lowest threshold of each well)",
                "J N = F with a lower absorbing barrier: the fate of the nascent, chemically activated adducts on their first\n\
                 collisional descent (Gonzalez-Garcia, Olzmann, PCCP 12, 12290 (2010), Sec. 3.2; Olzmann, PCCP 4, 3614 (2002)).\n\
                 A molecule transferred below the barrier counts as stabilized. Bimolecular rate coefficients k(R -> X) = k_inf Phi_X\n\
                 (Pilling, Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003), eq. 44).",
            )
        } else {
            (
                "FINAL STEADY STATE (no absorbing barrier)",
                "J N = R F without a barrier: under continuous formation, after the stabilized population has reached its own\n\
                 steady state, the yields of every exit, chemical activation and thermal reaction together (GO10 eqs. 7-9).",
            )
        };
        section(&mut report, title, text).map_err(io)?;
        for e in &unavailable {
            writeln!(report, "Not available: {e}").map_err(io)?;
        }
        if !results.is_empty() {
            write_results_table(network, &results, &mut machine).map_err(io)?;
            // k(R -> X) = k_inf Phi_X: bimolecular-to-bimolecular and bimolecular-to-well rates (intermediate steady
            // state) or the overall bimolecular-to-bimolecular rates (final steady state).
            let bimolecular = model
                .entrance_high_pressure_rate
                .as_ref()
                .map(|_| &k_inf as &dyn Fn(f64) -> f64);
            let groups = steady_state_groups(
                network,
                &temperatures,
                &pressures,
                &results,
                reactant.as_deref(),
                bimolecular,
            );
            write_groups(&mut report, &temperatures, &pressures, &groups).map_err(io)?;
            tables.extend(prefixed(label, &groups));

            // Bimolecular rate coefficients of the Reactant, k(R -> X) = k_inf Phi_X (cm3/s), machine-readable.
            if let (true, Some(k_inf_rate), Some(reactant)) = (
                is_intermediate,
                &model.entrance_high_pressure_rate,
                deck.global.reactant_name.as_ref(),
            ) {
                writeln!(
                    machine,
                    "\n# bimolecular rate coefficients of {reactant} [cm3/s], {label}"
                )
                .map_err(io)?;
                let mut header = vec![
                    "T[K]".to_string(),
                    "P[Torr]".to_string(),
                    "k_inf".to_string(),
                ];
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
                writeln!(machine, "{}", header.join(",")).map_err(io)?;
                for r in &results {
                    let k = k_inf_rate.rate_cm3_s(r.conditions.temperature_kelvin);
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
                            row.push(format!(
                                "{:.6e}",
                                k * r.result.wells[w].bimolecular_sink_yield
                            ));
                        }
                    }
                    writeln!(machine, "{}", row.join(",")).map_err(io)?;
                }
            }
        }
    }

    // Thermal rate coefficients of the final steady state: lowest eigenpair of its J (GO10 eq. 12).
    if let Some(thermal) = thermal {
        writeln!(
            machine,
            "\n# thermal rate coefficients of the final steady state: lowest eigenpair of J (GO10 eq. 12) by {:?}; \
             sum-rule tolerance {:e}",
            thermal.eigen_solver, thermal.sum_rule_tolerance
        )
        .map_err(io)?;
        let mut results = Vec::new();
        let mut unavailable = Vec::new();
        let per_condition = pool.map_conditions(&temperatures, &pressures, |t, p| {
            run_thermal_rate_coefficients(
                network,
                &[t],
                &[p],
                model.collision_model,
                thermal.eigen_solver,
                thermal.sum_rule_tolerance,
            )
        });
        for outcome in per_condition {
            match outcome {
                Ok(mut r) => {
                    for c in &r {
                        if let Some(warning) = &c.thermal.warning {
                            eprintln!(
                                "warning: T = {} K, p = {} Torr: {warning}",
                                c.conditions.temperature_kelvin, c.conditions.pressure_torr
                            );
                        }
                    }
                    results.append(&mut r)
                }
                Err(e) => {
                    writeln!(machine, "# not available: {e}").map_err(io)?;
                    unavailable.push(e);
                }
            }
        }
        section(
            &mut report,
            "THERMAL RATE COEFFICIENTS OF THE FINAL STEADY STATE (lowest eigenpair of J)",
            &format!(
                "k_uni = sum of the specific rate coefficients averaged over the normalized eigenvector of the lowest eigenvalue\n\
                 lambda_1 of J (Gonzalez-Garcia, Olzmann, PCCP 12, 12290 (2010), text after eq. 12); lambda_1 = k_uni in exact\n\
                 arithmetic (eq. 12), checked as the sum rule (warning above {:e}). Eigen-solver: {:?}. The thermal branching is\n\
                 k_th(channel)/k_uni: the yields of the thermal decay of the network.",
                thermal.sum_rule_tolerance, thermal.eigen_solver
            ),
        )
        .map_err(io)?;
        for e in &unavailable {
            writeln!(report, "Not available: {e}").map_err(io)?;
        }
        for r in &results {
            if let Some(warning) = &r.thermal.warning {
                writeln!(
                    report,
                    "Warning, T = {} K, p = {} torr: {warning}\n",
                    r.conditions.temperature_kelvin, r.conditions.pressure_torr
                )
                .map_err(io)?;
            }
        }
        if !results.is_empty() {
            write_thermal_table(network, &results, &mut machine).map_err(io)?;
            let mut groups = thermal_groups(network, &temperatures, &pressures, &results);
            // Association by detailed balance for a single well with one entrance channel.
            if let (1, [(w, c)], Some(k_inf_rate), Some(reactant)) = (
                network.wells.len(),
                model.entrance_channels.as_slice(),
                &model.entrance_high_pressure_rate,
                deck.global.reactant_name.as_ref(),
            ) {
                writeln!(
                    machine,
                    "\n# bimolecular rate coefficients of {reactant} [cm3/s] by detailed balance from the thermal k_uni of the final steady state"
                )
                .map_err(io)?;
                writeln!(
                    machine,
                    "T[K],P[Torr],k_inf_assoc[cm3/s],k_inf_diss[1/s],k_uni[1/s],k({reactant}->{})",
                    network.wells[*w].name
                )
                .map_err(io)?;
                let association = |r: &ThermalConditionResult| {
                    let k_assoc_inf = k_inf_rate.rate_cm3_s(r.conditions.temperature_kelvin);
                    let entrance = r
                        .thermal
                        .channels
                        .iter()
                        .find(|ch| ch.well == *w && ch.channel == *c)
                        .unwrap();
                    (
                        k_assoc_inf,
                        entrance.high_pressure_rate_s_inv,
                        r.thermal.k_uni_s_inv * k_assoc_inf / entrance.high_pressure_rate_s_inv,
                    )
                };
                for r in &results {
                    let (k_assoc_inf, k_diss_inf, k_assoc) = association(r);
                    writeln!(
                        machine,
                        "{},{},{:.6e},{:.6e},{:.6e},{:.6e}",
                        r.conditions.temperature_kelvin,
                        r.conditions.pressure_torr,
                        k_assoc_inf,
                        k_diss_inf,
                        r.thermal.k_uni_s_inv,
                        k_assoc
                    )
                    .map_err(io)?;
                }
                let lookup = |t: f64, p: f64| {
                    results
                        .iter()
                        .find(|r| {
                            r.conditions.temperature_kelvin == t && r.conditions.pressure_torr == p
                        })
                        .map(|r| association(r).2)
                };
                groups.push(QuantityGroup {
                    title: format!("Association by detailed balance, k({reactant} -> {}) = k_uni k_inf,assoc/k_inf,diss (cm^3/s)", network.wells[*w].name),
                    quantities: vec![Quantity {
                        name: format!("{reactant}->{}", network.wells[*w].name),
                        values: temperatures.iter().map(|&t| pressures.iter().map(|&p| lookup(t, p)).collect()).collect(),
                    }],
                });
            }
            write_groups(&mut report, &temperatures, &pressures, &groups).map_err(io)?;
            tables.extend(prefixed("thermal rate coefficients", &groups));
        }

        // Thermal fate of the molecules thermalized in each well, and chemical activation and thermal reaction
        // separately and together.
        let mut fates = Vec::new();
        let mut unavailable = Vec::new();
        let per_condition = pool.map_conditions(&temperatures, &pressures, |t, p| {
            run_thermal_well_fates(network, &[t], &[p], model.collision_model)
        });
        for outcome in per_condition {
            match outcome {
                Ok(mut f) => fates.append(&mut f),
                Err(e) => unavailable.push(e),
            }
        }
        section(
            &mut report,
            "THERMAL FATES OF THE WELLS (final steady state with the Boltzmann distribution of one well as the source)",
            "The probability that a molecule thermalized in a well ends in each exit of the network, k_x^T J^-1 f0_w.",
        )
        .map_err(io)?;
        for e in &unavailable {
            writeln!(report, "Not available: {e}").map_err(io)?;
        }
        let groups = well_fate_groups(network, &temperatures, &pressures, &fates);
        write_groups(&mut report, &temperatures, &pressures, &groups).map_err(io)?;
        tables.extend(prefixed("thermal fates", &groups));
    }

    // ---- Direct time integration.
    if let Some(plan) = time_integration {
        let times = log_spaced_times(
            plan.time_range_s.0,
            plan.time_range_s.1,
            plan.times_per_decade,
        );
        let integration = TimeIntegrationSettings {
            method: plan.integrator,
            relative_tolerance: plan.relative_tolerance,
            absolute_tolerance: plan.absolute_tolerance,
            times_s: times,
        };
        let options = ChemicalActivationOptions {
            collision_model: model.collision_model,
            steady_state: SteadyState::Final,
        };
        let per_condition = pool.map_conditions(&temperatures, &pressures, |t, p| {
            let source = thermal_entrance_source(network, &model.entrance_channels, KB_CM * t)?;
            let conditions = Conditions {
                temperature_kelvin: t,
                pressure_torr: p,
            };
            integrate_master_equation(
                network,
                &conditions,
                &options,
                &source,
                plan.initial_state,
                &integration,
            )
        });
        let mut evolutions = Vec::new();
        let mut unavailable = Vec::new();
        for outcome in per_condition {
            match outcome {
                Ok(e) => evolutions.push(e),
                Err(e) => {
                    writeln!(machine, "# not available: {e}").map_err(io)?;
                    unavailable.push(e);
                }
            }
        }
        // Machine-readable: one block per condition, time rows, populations and yields (fractions).
        for e in &evolutions {
            writeln!(
                machine,
                "\n# time evolution: T = {} K, p = {} Torr",
                e.conditions.temperature_kelvin, e.conditions.pressure_torr
            )
            .map_err(io)?;
            let mut header = vec!["t[s]".to_string()];
            header.extend(network.wells.iter().map(|w| format!("N({})", w.name)));
            header.extend(e.exits.iter().cloned());
            writeln!(machine, "{}", header.join(",")).map_err(io)?;
            for point in &e.points {
                let mut row = vec![format!("{:.6e}", point.time_s)];
                row.extend(
                    point
                        .well_populations
                        .iter()
                        .chain(&point.exit_yields)
                        .map(|v| format!("{v:.6e}")),
                );
                writeln!(machine, "{}", row.join(",")).map_err(io)?;
            }
        }
        section(
            &mut report,
            "DIRECT TIME INTEGRATION OF THE MASTER EQUATION",
            "dN/dt = R F - J N for the populations of all grains, and dY_x/dt = sum_E k_x(E) N(E) for the yield of every\n\
             exit (product channels, escape sinks), integrated with an adaptive, L-stable Rosenbrock method from the\n\
             initial state above: the early relaxation and the later chemistry in one integration, without a steady-state\n\
             assumption and without a separation of chemical and relaxation modes. For a pulse the total population +\n\
             yields stays 100%, and the yields at long times equal those of the final steady state (k_x^T J^-1 F).",
        )
        .map_err(io)?;
        for e in &unavailable {
            writeln!(report, "Not available: {e}").map_err(io)?;
        }
        writeln!(report, "--- Time evolution (populations of the wells and yields of the exits, % of the formed adducts) ---\n")
            .map_err(io)?;
        write_time_evolution_tables(&mut report, network, &evolutions).map_err(io)?;
        let k_inf_option = model
            .entrance_high_pressure_rate
            .as_ref()
            .map(|_| &k_inf as &dyn Fn(f64) -> f64);
        let groups = time_integration_groups(
            network,
            &temperatures,
            &pressures,
            &evolutions,
            reactant.as_deref(),
            k_inf_option,
        );
        write_groups(&mut report, &temperatures, &pressures, &groups).map_err(io)?;
        tables.extend(prefixed("time integration", &groups));
    }

    // ---- CSE method.
    if let Some(solver) = cse {
        writeln!(machine, "\n# CSE analysis ({solver:?})").map_err(io)?;
        let capture = |t: f64| {
            model
                .entrance_high_pressure_rate
                .as_ref()
                .map_or(f64::NAN, |k| k.rate_cm3_s(t))
        };
        let reactant_name = deck
            .global
            .reactant_name
            .as_deref()
            .filter(|_| model.entrance_high_pressure_rate.is_some());
        let mut all = Vec::new();
        let mut unavailable = Vec::new();
        let per_condition = pool.map_conditions(&temperatures, &pressures, |t, p| {
            run_phenomenological_rates(
                network,
                &[t],
                &[p],
                model.collision_model,
                solver,
                reactant_name.map(|r| (r, &capture as &dyn Fn(f64) -> f64)),
            )
        });
        for outcome in per_condition {
            match outcome {
                Ok(mut results) => {
                    for r in &results {
                        for warning in &r.rates.warnings {
                            eprintln!(
                                "warning: T = {} K, p = {} Torr: {warning}",
                                r.conditions.temperature_kelvin, r.conditions.pressure_torr
                            );
                        }
                    }
                    write_phenomenological_tables(&results, &mut machine).map_err(io)?;
                    all.append(&mut results);
                }
                Err(e) => {
                    writeln!(machine, "# not available: {e}").map_err(io)?;
                    unavailable.push(e);
                }
            }
        }
        section(
            &mut report,
            "PHENOMENOLOGICAL RATE COEFFICIENTS FROM THE CHEMICALLY SIGNIFICANT EIGENVALUES (CSE)",
            &format!(
                "Miller, Klippenstein, J. Phys. Chem. A 110, 10528 (2006); formulation of Georgievskii, Miller, Burke, Klippenstein,\n\
                 J. Phys. Chem. A 117, 12146 (2013), eqs. 21-30. All eigenpairs by {solver:?}. Rows: from, columns: to; wells in 1/s,\n\
                 the reactant row in cm^3/s. Diagonal: total loss of a well; for the reactant its net reaction (capture - return).\n\
                 Yields: the prompt branching of the reactant, the thermal fate of each well (absorbing chain of the well rate\n\
                 coefficients), and the long-time yields, direct + through the wells, which equal the final steady state."
            ),
        )
        .map_err(io)?;
        for e in &unavailable {
            writeln!(report, "Not available: {e}").map_err(io)?;
        }
        writeln!(report, "--- Species-to-species tables ---\n").map_err(io)?;
        write_cse_species_tables(&mut report, &all).map_err(io)?;
        let groups = cse_groups(&temperatures, &pressures, &all);
        write_groups(&mut report, &temperatures, &pressures, &groups).map_err(io)?;
        tables.extend(prefixed("CSE", &groups));
    }

    if let Some(file) = csv_path {
        std::fs::write(&file, &machine).map_err(|e| format!("--csv {file}: {e}"))?;
        // Every quantity group of the report, machine-readable, next to FILE: FILE_tables.csv.
        let tables_file = file.strip_suffix(".csv").map_or_else(
            || format!("{file}_tables.csv"),
            |stem| format!("{stem}_tables.csv"),
        );
        let mut buffer = Vec::new();
        write_groups_csv(&mut buffer, &temperatures, &pressures, &tables).map_err(io)?;
        std::fs::write(&tables_file, &buffer).map_err(|e| format!("{tables_file}: {e}"))?;
    }
    Ok(())
}
