//! Chemical-activation calculation over a list of temperatures and pressures.
//!
//! Workflow per (T, p), temperatures in the outer and pressures in the inner loop:
//!   1. collision frequency, <dE_down>(T), collision kernels and absorbing barriers per well, assembly
//!      of J (`chemical_activation_operator.rs`; PO14 eq. 2);
//!   2. source F (`chemical_activation_sources.rs`; PO14 eqs. 7-13);
//!   3. steady state J N = F (`chemical_activation_steady_state.rs`; PO14 eq. 5);
//!   4. observables (`chemical_activation_observables.rs`; O02 eq. 10, GO10 eqs. 8-10);
//!   5. one row of the results table.
//! Every row is checked: the relative residual of J N = F and the deviation of the mass balance from
//! one must not exceed the given tolerance.
//!
//! Table columns (comma separated):
//!   T[K], P[Torr];
//!   Phi(well:channel) of every product channel (yield, O02 eq. 10);
//!   Phi_stab(well) and Phi_sink(well) of every well;
//!   k_ca(well:channel)[1/s] of every channel, isomerizations included (GO10 eq. 9);
//!   k_tot[1/s] = 1/sum_w sum_i N_i, the total loss rate coefficient of the intermediates (for a normalized
//!   source the total loss flux is one);
//!   fpop(well), the population fractions, and <E>(well)[cm-1], the mean energies above the well bottoms;
//!   mass_balance and residual.
//!
//! Thermal route (`run_thermal_rate_coefficients`, Olzmann's eigenvalue analysis, `chemical_activation_eigen.rs`):
//! per (T, p) the operator without absorbing barrier is assembled and the pressure-dependent thermal rate
//! coefficient k_uni = lambda_1 (GO10 eq. 12) with its channel rate coefficients and the high-pressure
//! limits is tabulated (`write_thermal_table`): T[K], P[Torr], k_uni = sum_j k_j^th + k_c[D] (the average
//! over the thermal eigenvector, GO10 after eq. 12), lambda_1 (GO10 eq. 12), the relative sum-rule deviation
//! |lambda_1 - k_uni|/k_uni, lambda_2/k_uni, k_th(well:channel), k_inf(well:channel), k_sink(well),
//! fpop(well). Comment lines before the header explain k_uni with its reference and list the sum-rule
//! warnings (deviation above the tolerance; the results are kept).
//!
//! References: see `chemical_activation_network.rs`.

use std::io::Write;

use crate::constants::KB_CM;

use super::chemical_activation_eigen::{thermal_rate_coefficients, EigenSolver, ThermalRateCoefficients};
use super::chemical_activation_network::{
    ChannelDestination, ChemicalActivationNetwork, ChemicalActivationOptions, CollisionModel, Conditions, SteadyState,
};
use super::chemical_activation_observables::{evaluate_observables, ChemicalActivationResult};
use super::chemical_activation_operator::{assemble_operator, WellCollisionData};
use super::chemical_activation_sources::thermal_entrance_source;
use super::chemical_activation_steady_state::{project_source, solve_steady_state, LinearSolver};

/// How the intermediates are formed.
#[derive(Debug, Clone)]
pub enum SourceSpecification {
    /// Thermal reactants entering through the reverse of the channels (well, channel)
    /// (F ∝ rho k exp(-E/kT) on the absolute energy scale, PO14 eqs. 7 and 9), recomputed at every
    /// temperature.
    ThermalEntrance { channels: Vec<(usize, usize)> },
    /// A fixed distribution on the grid of every well (e.g. from `chemical_activation_sources.rs` or a
    /// previous master equation), used at every condition.
    Fixed(Vec<Vec<f64>>),
}

/// Specification of a chemical-activation calculation.
#[derive(Debug, Clone)]
pub struct ChemicalActivationRun {
    pub temperatures_kelvin: Vec<f64>,
    pub pressures_torr: Vec<f64>,
    pub options: ChemicalActivationOptions,
    pub solver: LinearSolver,
    pub source: SourceSpecification,
    /// Largest accepted relative residual ||F - J N||/||F|| and deviation of the mass balance from one.
    pub tolerance: f64,
}

/// Results at one temperature and pressure.
#[derive(Debug, Clone)]
pub struct ConditionResult {
    pub conditions: Conditions,
    pub result: ChemicalActivationResult,
    pub relative_residual: f64,
    pub max_relative_asymmetry: f64,
    pub wells: Vec<WellCollisionData>,
}

/// Run the calculation for every (T, p), temperatures outer, pressures inner.
pub fn run_chemical_activation(
    network: &ChemicalActivationNetwork,
    run: &ChemicalActivationRun,
) -> Result<Vec<ConditionResult>, String> {
    network.validate()?;
    match &run.source {
        SourceSpecification::ThermalEntrance { channels } => {
            for &(well, channel) in channels {
                let ok = network.wells.get(well).map_or(false, |w| channel < w.channels.len());
                if !ok {
                    return Err(format!("Entrance channel {channel} of well {well} does not exist."));
                }
            }
        }
        SourceSpecification::Fixed(source) => {
            if source.len() != network.wells.len() {
                return Err(format!("Source given for {} wells, the network has {}.", source.len(), network.wells.len()));
            }
        }
    }

    let mut results = Vec::with_capacity(run.temperatures_kelvin.len() * run.pressures_torr.len());
    for &temperature_kelvin in &run.temperatures_kelvin {
        // The source of thermal reactants depends on T only.
        let source: Vec<Vec<f64>> = match &run.source {
            SourceSpecification::ThermalEntrance { channels } => {
                thermal_entrance_source(network, channels, KB_CM * temperature_kelvin)?
            }
            SourceSpecification::Fixed(source) => source.clone(),
        };

        for &pressure_torr in &run.pressures_torr {
            let conditions = Conditions { temperature_kelvin, pressure_torr };
            let context = |e: String| format!("T = {temperature_kelvin} K, p = {pressure_torr} Torr: {e}");
            let op = assemble_operator(network, &conditions, &run.options).map_err(context)?;
            let projected = project_source(&op, &source).map_err(context)?;
            let solution = solve_steady_state(&op, &projected.on_states, &run.solver).map_err(context)?;
            let result = evaluate_observables(network, &op, &projected, &solution).map_err(context)?;
            if !(solution.relative_residual <= run.tolerance) {
                return Err(context(format!(
                    "relative residual {:e} of J N = F exceeds the tolerance {:e}. In the final steady state \
                     this means that J is numerically singular: the thermal rate coefficient of a well without \
                     a bimolecular sink is negligible compared with the collision frequency (the final steady \
                     state is reached only after about 0.1/k_uni, O02 text after eq. 13); use the intermediate \
                     steady state for such conditions.",
                    solution.relative_residual, run.tolerance
                )));
            }
            if !((result.mass_balance - 1.0).abs() <= run.tolerance) {
                return Err(context(format!(
                    "the yields sum to {} instead of 1 (tolerance {:e}).",
                    result.mass_balance, run.tolerance
                )));
            }
            results.push(ConditionResult {
                conditions,
                result,
                relative_residual: solution.relative_residual,
                max_relative_asymmetry: solution.max_relative_asymmetry,
                wells: op.wells,
            });
        }
    }
    Ok(results)
}

/// Write the results table (see the module documentation for the columns).
pub fn write_results_table<W: Write>(
    network: &ChemicalActivationNetwork,
    results: &[ConditionResult],
    out: &mut W,
) -> std::io::Result<()> {
    let label = |text: &str| text.replace(',', ";");
    let mut header = vec!["T[K]".to_string(), "P[Torr]".to_string()];
    for well in &network.wells {
        for channel in &well.channels {
            if matches!(channel.destination, ChannelDestination::Products { .. }) {
                header.push(format!("Phi({}:{})", label(&well.name), label(&channel.name)));
            }
        }
    }
    for well in &network.wells {
        header.push(format!("Phi_stab({})", label(&well.name)));
        header.push(format!("Phi_sink({})", label(&well.name)));
    }
    for well in &network.wells {
        for channel in &well.channels {
            header.push(format!("k_ca({}:{})[1/s]", label(&well.name), label(&channel.name)));
        }
    }
    header.push("k_tot[1/s]".into());
    for well in &network.wells {
        header.push(format!("fpop({})", label(&well.name)));
    }
    for well in &network.wells {
        header.push(format!("<E>({})[cm-1]", label(&well.name)));
    }
    header.push("mass_balance".into());
    header.push("residual".into());
    writeln!(out, "{}", header.join(","))?;

    for r in results {
        let res = &r.result;
        let mut row = vec![format!("{}", r.conditions.temperature_kelvin), format!("{}", r.conditions.pressure_torr)];
        for c in &res.channels {
            if matches!(c.destination, ChannelDestination::Products { .. }) {
                row.push(format!("{:.6e}", c.flux));
            }
        }
        for w in &res.wells {
            row.push(format!("{:.6e}", w.stabilization_yield));
            row.push(format!("{:.6e}", w.bimolecular_sink_yield));
        }
        for c in &res.channels {
            row.push(format!("{:.6e}", c.ca_rate_constant_s_inv));
        }
        let total_population: f64 = res.wells.iter().map(|w| w.population).sum();
        row.push(format!("{:.6e}", 1.0 / total_population));
        for w in &res.wells {
            row.push(format!("{:.6e}", w.population_fraction));
        }
        for w in &res.wells {
            row.push(format!("{:.6e}", w.mean_energy_cm1));
        }
        row.push(format!("{:.12}", res.mass_balance));
        row.push(format!("{:.3e}", r.relative_residual));
        writeln!(out, "{}", row.join(","))?;
    }
    Ok(())
}

/// Thermal rate coefficients at one temperature and pressure.
#[derive(Debug, Clone)]
pub struct ThermalConditionResult {
    pub conditions: Conditions,
    pub thermal: ThermalRateCoefficients,
    pub wells: Vec<WellCollisionData>,
}

/// Olzmann's eigenvalue analysis at every (T, p), temperatures outer, pressures inner: k_uni from the
/// thermal eigenvector of the operator without absorbing barrier (GO10 eq. 12 and after it), lambda_1,
/// channel and high-pressure rate coefficients. Above `sum_rule_tolerance` (relative deviation of lambda_1
/// from k_uni) the result carries a warning (see `thermal_rate_coefficients`).
pub fn run_thermal_rate_coefficients(
    network: &ChemicalActivationNetwork,
    temperatures_kelvin: &[f64],
    pressures_torr: &[f64],
    collision_model: CollisionModel,
    solver: EigenSolver,
    sum_rule_tolerance: f64,
) -> Result<Vec<ThermalConditionResult>, String> {
    network.validate()?;
    let options = ChemicalActivationOptions { collision_model, steady_state: SteadyState::Final };
    let mut results = Vec::with_capacity(temperatures_kelvin.len() * pressures_torr.len());
    for &temperature_kelvin in temperatures_kelvin {
        for &pressure_torr in pressures_torr {
            let conditions = Conditions { temperature_kelvin, pressure_torr };
            let context = |e: String| format!("T = {temperature_kelvin} K, p = {pressure_torr} Torr: {e}");
            let op = assemble_operator(network, &conditions, &options).map_err(context)?;
            let thermal = thermal_rate_coefficients(network, &op, solver, sum_rule_tolerance).map_err(context)?;
            results.push(ThermalConditionResult { conditions, thermal, wells: op.wells });
        }
    }
    Ok(results)
}

/// Explanation written before the thermal table.
const THERMAL_TABLE_EXPLANATION: &str = "\
# k_uni = sum_j k_th(j) + k_sink: the specific rate coefficients k(E) of the product channels (and the sink)
#   averaged over the normalized thermal eigenvector of J, the eigenvector of the lowest eigenvalue lambda_1,
#   \"analogous to eqn (9) but with Ns = Ns_th being the normalized eigenvector associated with the lowest
#   eigenvalue lambda_1\" (Gonzalez-Garcia, Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010), text after eq. 12).
# lambda_1 = lowest eigenvalue of J (same paper, eq. 12). The column sums of J are the loss rates, so
#   lambda_1 = k_uni exactly; numerically lambda_1 carries an absolute error of order eps*||J|| and is lost when
#   it lies many orders of magnitude below the collision frequency (deep wells, low T), while the eigenvector,
#   and with it k_uni, is much less sensitive. k_uni is therefore the reported rate coefficient.
# precision_floor = eps*max(S_ii), S the symmetrized J: the order of the absolute rounding error of lambda_1.
#   A lambda_1 near or below it, even a negative one, is rounding noise below the double-precision floor; this
#   is not a merging of eigenvalues, which would show as a small lambda_2/k_uni (thermal decay no longer
#   separated from relaxation).
# sum_rule_deviation = |lambda_1 - k_uni| / k_uni; above the tolerance a warning is listed below.
";

/// Table of the thermal route (see the module documentation for the columns).
pub fn write_thermal_table<W: Write>(
    network: &ChemicalActivationNetwork,
    results: &[ThermalConditionResult],
    out: &mut W,
) -> std::io::Result<()> {
    write!(out, "{THERMAL_TABLE_EXPLANATION}")?;
    for r in results {
        if let Some(warning) = &r.thermal.warning {
            writeln!(
                out,
                "# warning: T = {} K, p = {} Torr: {warning}",
                r.conditions.temperature_kelvin, r.conditions.pressure_torr
            )?;
        }
    }
    let label = |text: &str| text.replace(',', ";");
    let mut header: Vec<String> = [
        "T[K]",
        "P[Torr]",
        "k_uni[1/s]",
        "lambda_1[1/s]",
        "precision_floor[1/s]",
        "sum_rule_deviation",
        "lambda_2/k_uni",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    for prefix in ["k_th", "k_inf"] {
        for well in &network.wells {
            for channel in &well.channels {
                header.push(format!("{prefix}({}:{})[1/s]", label(&well.name), label(&channel.name)));
            }
        }
    }
    for well in &network.wells {
        header.push(format!("k_sink({})[1/s]", label(&well.name)));
    }
    for well in &network.wells {
        header.push(format!("fpop({})", label(&well.name)));
    }
    writeln!(out, "{}", header.join(","))?;
    for r in results {
        let th = &r.thermal;
        let mut row = vec![
            format!("{}", r.conditions.temperature_kelvin),
            format!("{}", r.conditions.pressure_torr),
            format!("{:.6e}", th.k_uni_s_inv),
            format!("{:.6e}", th.lambda_1_s_inv),
            format!("{:.3e}", th.precision_floor_s_inv),
            format!("{:.3e}", th.sum_rule_relative_deviation),
            format!("{:.6e}", th.lambda_2_s_inv / th.k_uni_s_inv),
        ];
        row.extend(th.channels.iter().map(|c| format!("{:.6e}", c.thermal_rate_s_inv)));
        row.extend(th.channels.iter().map(|c| format!("{:.6e}", c.high_pressure_rate_s_inv)));
        row.extend(th.sink_rates_s_inv.iter().map(|k| format!("{k:.6e}")));
        row.extend(th.population_fractions.iter().map(|x| format!("{x:.6e}")));
        writeln!(out, "{}", row.join(","))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::masterequation::chemical_activation_eigen::{EigenSolver, DEFAULT_SUM_RULE_TOLERANCE};
    use crate::masterequation::chemical_activation_network::{AbsorbingBarrier, CollisionModel, SteadyState};
    use crate::masterequation::chemical_activation_operator::tests::two_well_network;
    use crate::masterequation::chemical_activation_sources::thermal_source_from_rate;
    use crate::constants::KB_CM;

    /// Two-well network with an entrance channel A -> reactants opening at grain 320.
    fn network_with_entrance() -> ChemicalActivationNetwork {
        let mut network = two_well_network();
        let a = &mut network.wells[0];
        a.channels.push(crate::masterequation::chemical_activation_network::Channel {
            name: "A->reactants".into(),
            destination: ChannelDestination::Products { name: "R".into() },
            threshold_grain: None,
            rate_constant_s_inv: (0..400).map(|i| if i >= 320 { 2.0e7 * ((i - 320) as f64 + 1.0) } else { 0.0 }).collect(),
        });
        network
    }

    fn run_spec(pressures: Vec<f64>) -> ChemicalActivationRun {
        ChemicalActivationRun {
            // At both temperatures 10 k_BT stays below the lowest thresholds of both wells (the
            // intermediate steady state is defined).
            temperatures_kelvin: vec![250.0, 300.0],
            pressures_torr: pressures,
            options: ChemicalActivationOptions {
                collision_model: CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 },
                steady_state: SteadyState::Intermediate { barrier: AbsorbingBarrier::default() },
            },
            solver: LinearSolver::BandedCholesky,
            source: SourceSpecification::ThermalEntrance { channels: vec![(0, 2)] },
            tolerance: 1e-8,
        }
    }

    #[test]
    fn every_temperature_and_pressure_is_solved_in_order() {
        let network = network_with_entrance();
        let results = run_chemical_activation(&network, &run_spec(vec![1.0, 100.0, 10000.0])).unwrap();
        let order: Vec<(f64, f64)> = results.iter().map(|r| (r.conditions.temperature_kelvin, r.conditions.pressure_torr)).collect();
        assert_eq!(order, vec![(250.0, 1.0), (250.0, 100.0), (250.0, 10000.0), (300.0, 1.0), (300.0, 100.0), (300.0, 10000.0)]);
        for r in &results {
            assert!((r.result.mass_balance - 1.0).abs() < 1e-8);
        }
    }

    #[test]
    fn stabilization_grows_with_pressure() {
        let network = network_with_entrance();
        let results = run_chemical_activation(&network, &run_spec(vec![0.1, 10.0, 1000.0, 100000.0])).unwrap();
        for t in results.chunks(4) {
            let stab: Vec<f64> = t.iter().map(|r| r.result.total_stabilization_yield()).collect();
            assert!(stab.windows(2).all(|p| p[1] > p[0]), "Phi_stab not increasing with pressure: {stab:?}");
        }
    }

    #[test]
    fn thermal_entrance_source_is_rebuilt_at_every_temperature() {
        // Same result as a fixed source built by hand at that temperature.
        let network = network_with_entrance();
        let spec = run_spec(vec![100.0]);
        let from_driver = run_chemical_activation(&network, &spec).unwrap();
        for r in &from_driver {
            let t = r.conditions.temperature_kelvin;
            let a = &network.wells[0];
            let f_a = thermal_source_from_rate(&a.density_of_states, &a.channels[2].rate_constant_s_inv, 10.0, KB_CM * t).unwrap();
            let fixed = ChemicalActivationRun {
                temperatures_kelvin: vec![t],
                source: SourceSpecification::Fixed(vec![f_a, vec![0.0; 460]]),
                ..spec.clone()
            };
            let by_hand = run_chemical_activation(&network, &fixed).unwrap();
            for (x, y) in r.result.channels.iter().zip(&by_hand[0].result.channels) {
                assert!((x.flux - y.flux).abs() <= 1e-13 * x.flux.abs().max(1e-300));
            }
        }
    }

    #[test]
    fn results_table_has_a_header_and_one_row_per_condition() {
        let network = network_with_entrance();
        let results = run_chemical_activation(&network, &run_spec(vec![1.0, 100.0])).unwrap();
        let mut out = Vec::new();
        write_results_table(&network, &results, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 1 + 4);
        let columns = lines[0].split(',').count();
        // T, P; 3 product channels; stab and sink of 2 wells; k_ca of 5 channels; k_tot; fpop and <E> of
        // 2 wells; mass balance and residual.
        assert_eq!(columns, 2 + 3 + 4 + 5 + 1 + 4 + 2);
        for line in &lines[1..] {
            assert_eq!(line.split(',').count(), columns);
        }
        assert!(lines[0].starts_with("T[K],P[Torr],Phi(A:A-products)"));
    }

    #[test]
    fn the_thermal_route_gives_falloff_rate_coefficients_at_every_condition() {
        // k_uni (GO10 eq. 12 and the eigenvector average after it) rises with pressure towards the
        // high-pressure limit.
        let network = network_with_entrance();
        let results = run_thermal_rate_coefficients(
            &network,
            &[250.0, 300.0],
            &[10.0, 1000.0, 100000.0],
            CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 },
            EigenSolver::default(),
            DEFAULT_SUM_RULE_TOLERANCE,
        )
        .unwrap();
        assert_eq!(results.len(), 6);
        for t in results.chunks(3) {
            let k: Vec<f64> = t.iter().map(|r| r.thermal.k_uni_s_inv).collect();
            assert!(k.windows(2).all(|p| p[1] > p[0]), "{k:?}");
        }
        let mut out = Vec::new();
        write_thermal_table(&network, &results, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        // The table explains k_uni with its reference in comment lines before the header.
        let comments: Vec<&str> = text.lines().take_while(|l| l.starts_with('#')).collect();
        assert!(comments.iter().any(|l| l.contains("Phys. Chem. Chem. Phys. 12, 12290 (2010)")), "{comments:?}");
        assert!(comments.iter().any(|l| l.contains("rounding noise below the double-precision floor")), "{comments:?}");
        let lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(lines.len(), 7);
        let columns = lines[0].split(',').count();
        assert!(lines[0].starts_with(
            "T[K],P[Torr],k_uni[1/s],lambda_1[1/s],precision_floor[1/s],sum_rule_deviation,lambda_2/k_uni"
        ));
        // The deviation column is the reported relative sum-rule deviation of each condition (printed with
        // four significant digits).
        for (line, r) in lines[1..].iter().zip(&results) {
            let k_uni: f64 = line.split(',').nth(2).unwrap().parse().unwrap();
            assert!((k_uni / r.thermal.k_uni_s_inv - 1.0).abs() < 1e-6);
            let deviation: f64 = line.split(',').nth(5).unwrap().parse().unwrap();
            let reported = r.thermal.sum_rule_relative_deviation;
            assert!((deviation - reported).abs() <= 1e-3 * reported, "{deviation} vs {reported}");
        }
        assert!(lines[1..].iter().all(|l| l.split(',').count() == columns));
    }

    #[test]
    fn sum_rule_warnings_are_written_before_the_thermal_table() {
        let network = network_with_entrance();
        let mut results = run_thermal_rate_coefficients(
            &network,
            &[300.0],
            &[1000.0],
            CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 },
            EigenSolver::default(),
            DEFAULT_SUM_RULE_TOLERANCE,
        )
        .unwrap();
        results[0].thermal.warning = Some("test warning".into());
        let mut out = Vec::new();
        write_thermal_table(&network, &results, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.lines().any(|l| l == "# warning: T = 300 K, p = 1000 Torr: test warning"), "{text}");
    }

    #[test]
    fn a_source_with_the_wrong_number_of_wells_is_an_error() {
        let network = network_with_entrance();
        let mut spec = run_spec(vec![1.0]);
        spec.source = SourceSpecification::Fixed(vec![vec![1.0; 400]]);
        assert!(run_chemical_activation(&network, &spec).is_err());
        spec.source = SourceSpecification::ThermalEntrance { channels: vec![(0, 7)] };
        assert!(run_chemical_activation(&network, &spec).is_err());
    }
}
