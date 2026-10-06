//! Phenomenological rate coefficients from the chemically significant eigenvalues (CSE) of the master
//! equation: the eigenvalue method of Miller and Klippenstein (J. A. Miller, S. J. Klippenstein, J. Phys.
//! Chem. A 110, 10528 (2006), MK06 eqs. 12-29; Bartis, Widom, J. Chem. Phys. 60, 3474 (1974)) in the
//! formulation of Georgievskii, Miller, Burke, Klippenstein, J. Phys. Chem. A 117, 12146 (2013) (G13), in
//! which bimolecular reactants act as decoupled thermal sources and all bimolecular products (and
//! pseudo-first-order sinks) as infinite sinks.
//!
//! With G the relaxation operator J of the final steady state (all wells, no absorbing barrier), its
//! symmetrized form S = D^-1 J D, D = diag(sqrt(f0)), f0 the Boltzmann weights on the absolute energy scale,
//! eigenvalues Lambda_lambda (ascending) and orthonormal eigenvectors u_lambda, the eigenvector of G is
//! f^(lambda) = D u_lambda (normalized with the scalar product of G13 eq. 12). For N wells the N lowest
//! eigenpairs are the chemical eigenstates (G13 Sec. III; MK06 eq. 19), the others the internal-energy
//! relaxation eigenstates. Then (G13 eqs. 15, 23, 25, 27-30, 21, 22):
//!   M_(i,lambda) = (1/sqrt(Q_i)) sum_(grains of i) f^(lambda),   Q_i = sum_(grains of i) f0
//!   p_lambda^(nu) = sum_(all grains) f^(lambda)(E) k_(->nu)(E)        (k = N#/(h rho), eq. 15)
//!   k_(j->i)  = -sqrt(Q_i/Q_j) (M Lambda M^-1)_(i,j),  k_i = (M Lambda M^-1)_(i,i)   (eq. 27; eq. 29)
//!   k_(i->nu) = (1/sqrt(Q_i)) sum_(lambda chem) (M^-1)_(lambda,i) p_lambda^(nu)      (eq. 30)
//!   k_(nu->i) = (sqrt(Q_i)/Q_nu) sum_(lambda chem) M_(i,lambda) p_lambda^(nu)         (eq. 28)
//!   k_(nu->mu) = (1/Q_nu) sum_(lambda relax) p_lambda^(mu) p_lambda^(nu) / Lambda_lambda  (eq. 21)
//!   k_nu^(c) = (1/Q_nu) sum_(all grains) k_(->nu) f0 = k_(nu->nu) + sum_mu k_(nu->mu) + sum_i k_(nu->i)
//!                                                                                  (eqs. 22-23)
//! where nu labels a bimolecular channel (or the escape sink of a well, an energy-independent k). For the
//! bimolecular reactant the capture rate coefficient k^(c) (the high-pressure association rate coefficient
//! of the entrance channels) gives 1/Q_nu = k^(c) / sum k_(->nu) f0, so that the rates (cm3 s-1) need no
//! partition function of the reactants.
//!
//! Diagnostics: the separation Lambda_N / Lambda_(N+1) (the method requires |lambda_N| << |lambda_(N+1)|,
//! MK06 after eq. 24), the projection of every chemical eigenvector on the relaxational subspace,
//! 1 - sum_i M_(i,lambda)^2 (G13 Fig. 2), the loss balance of eq. 29, detailed balance k_(i->j) Q_i =
//! k_(j->i) Q_j, and the double-precision floor eps max S_ii of the eigenvalues.

use super::chemical_activation_eigen::{full_decomposition, require_symmetry, EigenSolver};
use super::chemical_activation_network::{ChannelDestination, ChemicalActivationNetwork};
use super::chemical_activation_operator::ChemicalActivationOperator;
use super::chemical_activation_steady_state::symmetrize;
use crate::numeric::dense_inverse::invert_dense;

/// Ratio Lambda_N / Lambda_(N+1) above which a warning says that the chemical eigenvalues are not well
/// separated from the relaxation ones.
pub const CSE_SEPARATION_WARNING: f64 = 0.1;

/// Rate coefficients from the bimolecular reactant (cm3 s-1).
#[derive(Debug, Clone)]
pub struct ReactantRates {
    pub name: String,
    /// Capture (high-pressure association) rate coefficient k^(c).
    pub capture_cm3_s: f64,
    /// k_(R->i) for every well (G13 eq. 28).
    pub to_well_cm3_s: Vec<f64>,
    /// k_(R->nu) for every bimolecular channel (G13 eq. 21); the entry of the reactant itself is the return
    /// k_(R->R) from the capture balance (eq. 22).
    pub to_bimolecular_cm3_s: Vec<f64>,
}

/// Phenomenological rate coefficients of a network at one temperature and pressure.
#[derive(Debug, Clone)]
pub struct PhenomenologicalRates {
    pub wells: Vec<String>,
    /// Bimolecular channels: product names (the reactant included when a channel leads to it), then the
    /// escape sinks "escape(W)" of the wells with a pseudo-first-order sink.
    pub bimolecular: Vec<String>,
    /// Q_i of the wells (relative, common absolute energy scale).
    pub partition_functions: Vec<f64>,
    /// The N chemical eigenvalues (s-1), ascending.
    pub chemical_eigenvalues_s_inv: Vec<f64>,
    /// The lowest relaxation eigenvalue (s-1); NaN if there is none.
    pub relaxation_eigenvalue_s_inv: f64,
    /// 1 - sum_i M_(i,lambda)^2 of every chemical eigenvector.
    pub relaxational_projection: Vec<f64>,
    /// [i][j] = k_(i->j) (s-1); the diagonal is the total loss k_i.
    pub well_to_well_s_inv: Vec<Vec<f64>>,
    /// [i][nu] = k_(i->nu) (s-1).
    pub well_to_bimolecular_s_inv: Vec<Vec<f64>>,
    pub reactant: Option<ReactantRates>,
    /// eps max S_ii (s-1): absolute rounding error of the eigenvalues.
    pub precision_floor_s_inv: f64,
    /// max_i |k_i - sum_j k_(i->j) - sum_nu k_(i->nu)| / k_i (eq. 29).
    pub loss_balance_max_deviation: f64,
    /// max over well pairs of |k_(i->j) Q_i - k_(j->i) Q_j| / max(k_(i->j) Q_i, k_(j->i) Q_j).
    pub detailed_balance_max_deviation: f64,
    pub warnings: Vec<String>,
}

/// Phenomenological rate coefficients from the chemically significant eigenpairs of the operator `op`
/// (assembled without absorbing barrier). `reactant`: name of the bimolecular channel that is the reactant
/// and its capture rate coefficient k^(c) (cm3 s-1). `solver`: a full decomposition (all eigenpairs).
pub fn phenomenological_rate_coefficients(
    network: &ChemicalActivationNetwork,
    op: &ChemicalActivationOperator,
    reactant: Option<(&str, f64)>,
    solver: EigenSolver,
) -> Result<PhenomenologicalRates, String> {
    let symmetrized = symmetrize(op);
    require_symmetry(symmetrized.max_relative_asymmetry)?;
    let (values, vectors) = full_decomposition(&symmetrized.dense(), solver)?;
    let n_wells = network.wells.len();
    let n_states = op.dimension();
    if n_states <= n_wells {
        return Err("CSE analysis: fewer states than wells.".into());
    }
    let d = &symmetrized.d;
    // Q_i and the states of every well.
    let mut q = vec![0.0; n_wells];
    for (k, &(w, _)) in op.states.iter().enumerate() {
        q[w] += d[k] * d[k];
    }
    if let Some(w) = q.iter().position(|x| !(*x > 0.0)) {
        return Err(format!("CSE analysis: well '{}' has no states.", network.wells[w].name));
    }

    // Bimolecular channels and the rate of every state into each of them.
    let mut bimolecular: Vec<String> = Vec::new();
    for well in &network.wells {
        for channel in &well.channels {
            if let ChannelDestination::Products { name } = &channel.destination {
                if !bimolecular.contains(name) {
                    bimolecular.push(name.clone());
                }
            }
        }
    }
    for well in &network.wells {
        if well.bimolecular_sink_s_inv > 0.0 {
            bimolecular.push(format!("escape({})", well.name));
        }
    }
    let mut k_into = vec![vec![0.0; n_states]; bimolecular.len()];
    for (k, &(w, _)) in op.states.iter().enumerate() {
        let well = &network.wells[w];
        for channel in &well.channels {
            if let ChannelDestination::Products { name } = &channel.destination {
                let nu = bimolecular.iter().position(|b| b == name).unwrap();
                // A low-energy reservoir reacts with the Boltzmann average over its grains.
                k_into[nu][k] += op.state_rate(k, &channel.rate_constant_s_inv);
            }
        }
        if well.bimolecular_sink_s_inv > 0.0 {
            let nu = bimolecular.iter().position(|b| *b == format!("escape({})", well.name)).unwrap();
            k_into[nu][k] += well.bimolecular_sink_s_inv;
        }
    }

    // M (G13 eq. 25) over the chemical eigenstates, its inverse, and K = M Lambda M^-1.
    let n = n_wells;
    let mut m = vec![vec![0.0; n]; n];
    for (lambda, u) in vectors.iter().take(n).enumerate() {
        for (k, &(w, _)) in op.states.iter().enumerate() {
            m[w][lambda] += d[k] * u[k];
        }
    }
    for (w, row) in m.iter_mut().enumerate() {
        row.iter_mut().for_each(|x| *x /= q[w].sqrt());
    }
    let m_inv = invert_dense(&m).map_err(|e| format!("CSE analysis: the projection matrix M is singular ({e})."))?;
    let k_matrix: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| (0..n).map(|l| m[i][l] * values[l] * m_inv[l][j]).sum()).collect())
        .collect();
    // well_to_well[i][j] = k_(i->j) = -sqrt(Q_j/Q_i) K_(j,i); diagonal k_i = K_(i,i) (eqs. 27, 29).
    let well_to_well: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { k_matrix[i][i] } else { -(q[j] / q[i]).sqrt() * k_matrix[j][i] }).collect())
        .collect();

    // p_lambda^(nu) for all eigenstates (eq. 15).
    let p: Vec<Vec<f64>> = k_into
        .iter()
        .map(|k_nu| vectors.iter().map(|u| (0..n_states).map(|k| d[k] * u[k] * k_nu[k]).sum()).collect())
        .collect();
    // k_(i->nu) (eq. 30).
    let well_to_bimolecular: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..bimolecular.len()).map(|nu| (0..n).map(|l| m_inv[l][i] * p[nu][l]).sum::<f64>() / q[i].sqrt()).collect())
        .collect();

    // Diagnostics.
    let loss_balance_max_deviation = (0..n)
        .map(|i| {
            let out: f64 = (0..n).filter(|&j| j != i).map(|j| well_to_well[i][j]).sum::<f64>()
                + well_to_bimolecular[i].iter().sum::<f64>();
            (well_to_well[i][i] - out).abs() / well_to_well[i][i].abs()
        })
        .fold(0.0, f64::max);
    let mut detailed_balance_max_deviation: f64 = 0.0;
    for i in 0..n {
        for j in i + 1..n {
            let (a, b) = (well_to_well[i][j] * q[i], well_to_well[j][i] * q[j]);
            let scale = a.abs().max(b.abs());
            if scale > 0.0 {
                detailed_balance_max_deviation = detailed_balance_max_deviation.max((a - b).abs() / scale);
            }
        }
    }
    let relaxational_projection: Vec<f64> = (0..n).map(|l| 1.0 - (0..n).map(|i| m[i][l] * m[i][l]).sum::<f64>()).collect();
    let relaxation_eigenvalue = values.get(n).copied().unwrap_or(f64::NAN);
    let precision_floor = f64::EPSILON * symmetrized.band_matrix().max_abs_diagonal();
    let mut warnings = Vec::new();
    let separation = values[n - 1] / relaxation_eigenvalue;
    if separation > CSE_SEPARATION_WARNING {
        warnings.push(format!(
            "the largest chemical eigenvalue {:e} s-1 is {separation:.3} of the lowest relaxation eigenvalue {relaxation_eigenvalue:e} \
             s-1: the chemical eigenstates are not well separated (the method requires |lambda_N| << |lambda_(N+1)|, \
             Miller, Klippenstein 2006); wells in fast equilibrium should be merged (Georgievskii et al. 2013, Sec. IV).",
            values[n - 1]
        ));
    }
    for (l, &value) in values.iter().take(n).enumerate() {
        if value < 100.0 * precision_floor {
            warnings.push(format!(
                "chemical eigenvalue {l} = {value:e} s-1 is within a factor 100 of the double-precision floor \
                 eps*max(S_ii) = {precision_floor:e} s-1: the rate coefficients of this mode are not resolved."
            ));
        }
    }

    // Reactant rates (eqs. 28, 21, 22).
    let reactant = match reactant {
        None => None,
        Some((name, capture)) => {
            let r = bimolecular
                .iter()
                .position(|b| b == name)
                .ok_or_else(|| format!("CSE analysis: no channel leads to the reactant '{name}'."))?;
            let denominator: f64 = (0..n_states).map(|k| k_into[r][k] * d[k] * d[k]).sum();
            if !(denominator > 0.0) {
                return Err(format!("CSE analysis: the reactant '{name}' has no association flux."));
            }
            let to_well: Vec<f64> = (0..n)
                .map(|i| capture * q[i].sqrt() * (0..n).map(|l| m[i][l] * p[r][l]).sum::<f64>() / denominator)
                .collect();
            let mut to_bimolecular: Vec<f64> = (0..bimolecular.len())
                .map(|mu| {
                    if mu == r {
                        0.0
                    } else {
                        capture * (n..values.len()).map(|l| p[mu][l] * p[r][l] / values[l]).sum::<f64>() / denominator
                    }
                })
                .collect();
            to_bimolecular[r] = capture - to_well.iter().sum::<f64>() - to_bimolecular.iter().sum::<f64>();
            Some(ReactantRates { name: name.to_string(), capture_cm3_s: capture, to_well_cm3_s: to_well, to_bimolecular_cm3_s: to_bimolecular })
        }
    };

    Ok(PhenomenologicalRates {
        wells: network.wells.iter().map(|w| w.name.clone()).collect(),
        bimolecular,
        partition_functions: q,
        chemical_eigenvalues_s_inv: values[..n].to_vec(),
        relaxation_eigenvalue_s_inv: relaxation_eigenvalue,
        relaxational_projection,
        well_to_well_s_inv: well_to_well,
        well_to_bimolecular_s_inv: well_to_bimolecular,
        reactant,
        precision_floor_s_inv: precision_floor,
        loss_balance_max_deviation,
        detailed_balance_max_deviation,
        warnings,
    })
}

/// Yields that follow from the phenomenological rate coefficients of the reactant and of the wells.
/// "Net reaction" excludes the return to the reactant: every quantity below is a fraction.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactantYields {
    /// Prompt branching of the net reaction of the reactant over the wells (stabilization) and the
    /// bimolecular channels other than the reactant (direct, chemically activated products):
    /// k_(R->X) / sum_(X != R) k_(R->X). Order: wells, then `bimolecular` without the reactant.
    pub prompt_branching: Vec<f64>,
    /// Thermal fate of a molecule in each well: the probabilities to end in each bimolecular channel
    /// (`bimolecular` order, the reactant included), from the absorbing chain of the well rate coefficients.
    pub well_fates: Vec<Vec<f64>>,
    /// Long-time yields as fractions of the eventual net reaction (no return to the reactant), for the
    /// bimolecular channels other than the reactant: the direct part, the part through the wells, and
    /// their sum. They equal k_x^T J^-1 F of the final steady state, normalized without the reactant.
    pub direct: Vec<f64>,
    pub via_wells: Vec<f64>,
    pub total: Vec<f64>,
    /// Names of the columns of `direct`, `via_wells` and `total` (the bimolecular channels without the reactant).
    pub channels: Vec<String>,
}

/// Yields of the reactant and thermal fates of the wells from the phenomenological rate coefficients.
/// None without a reactant.
pub fn reactant_yields(rates: &PhenomenologicalRates) -> Option<Result<ReactantYields, String>> {
    let reactant = rates.reactant.as_ref()?;
    Some((|| {
        let r_index = rates
            .bimolecular
            .iter()
            .position(|b| *b == reactant.name)
            .ok_or_else(|| {
                format!(
                    "Reactant yields: '{}' is not a bimolecular channel.",
                    reactant.name
                )
            })?;
        let others: Vec<usize> = (0..rates.bimolecular.len())
            .filter(|&nu| nu != r_index)
            .collect();

        // Prompt branching of the net reaction of the reactant (wells, then the other bimolecular channels).
        let mut prompt_branching = reactant.to_well_cm3_s.clone();
        prompt_branching.extend(others.iter().map(|&nu| reactant.to_bimolecular_cm3_s[nu]));
        let net: f64 = prompt_branching.iter().sum();
        prompt_branching.iter_mut().for_each(|x| *x /= net);

        // Thermal fates of the wells: absorbing chain B = (I - Q)^-1 A, with Q_ij = k_(i->j)/k_i between the
        // wells, A_i,nu = k_(i->nu)/k_i into the bimolecular channels, k_i the sum of all rates out of well i.
        let n = rates.wells.len();
        let loss: Vec<f64> = (0..n)
            .map(|i| {
                (0..n)
                    .filter(|&j| j != i)
                    .map(|j| rates.well_to_well_s_inv[i][j])
                    .sum::<f64>()
                    + rates.well_to_bimolecular_s_inv[i].iter().sum::<f64>()
            })
            .collect();
        let i_minus_q: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            1.0
                        } else {
                            -rates.well_to_well_s_inv[i][j] / loss[i]
                        }
                    })
                    .collect()
            })
            .collect();
        let fundamental =
            invert_dense(&i_minus_q).map_err(|e| format!("Reactant yields, well chain: {e}"))?;
        let well_fates: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                (0..rates.bimolecular.len())
                    .map(|nu| {
                        (0..n)
                            .map(|j| {
                                fundamental[i][j] * rates.well_to_bimolecular_s_inv[j][nu] / loss[j]
                            })
                            .sum()
                    })
                    .collect()
            })
            .collect();

        // Long-time yields: direct formation plus formation through the wells, normalized to the eventual net
        // reaction (the final return to the reactant excluded).
        let direct_raw: Vec<f64> = others
            .iter()
            .map(|&nu| reactant.to_bimolecular_cm3_s[nu])
            .collect();
        let via_raw: Vec<f64> = others
            .iter()
            .map(|&nu| {
                (0..n)
                    .map(|i| reactant.to_well_cm3_s[i] * well_fates[i][nu])
                    .sum()
            })
            .collect();
        let eventual: f64 = direct_raw.iter().chain(&via_raw).sum();
        let direct: Vec<f64> = direct_raw.iter().map(|x| x / eventual).collect();
        let via_wells: Vec<f64> = via_raw.iter().map(|x| x / eventual).collect();
        let total = direct.iter().zip(&via_wells).map(|(a, b)| a + b).collect();
        let channels = others
            .iter()
            .map(|&nu| rates.bimolecular[nu].clone())
            .collect();
        Ok(ReactantYields {
            prompt_branching,
            well_fates,
            direct,
            via_wells,
            total,
            channels,
        })
    })())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::KB_CM;
    use crate::masterequation::chemical_activation_eigen::{thermal_rate_coefficients, DEFAULT_SUM_RULE_TOLERANCE};
    use crate::masterequation::chemical_activation_network::tests::test_well;
    use crate::masterequation::chemical_activation_network::{
        ChannelDestination, ChemicalActivationNetwork, ChemicalActivationOptions, CollisionModel, Conditions, SteadyState,
    };
    use crate::masterequation::chemical_activation_operator::assemble_operator;
    use crate::masterequation::chemical_activation_operator::tests::{conditions, two_well_network};

    fn final_options() -> ChemicalActivationOptions {
        ChemicalActivationOptions { collision_model: CollisionModel::default(), steady_state: SteadyState::Final }
    }

    /// Boltzmann average of k(E) over the grains of well `w` (the high-pressure rate coefficient).
    /// Canonical average of k(E) of a channel over all grains of well w (a low-energy reservoir keeps its
    /// Boltzmann weight).
    fn high_pressure_rate(network: &ChemicalActivationNetwork, w: usize, channel: usize, kt: f64) -> f64 {
        let well = &network.wells[w];
        let f: Vec<f64> = (0..well.grain_count())
            .map(|i| well.density_of_states[i] * (-network.absolute_energy_cm1(w, i) / kt).exp())
            .collect();
        well.channels[channel].rate_constant_s_inv.iter().zip(&f).map(|(k, f)| k * f).sum::<f64>() / f.iter().sum::<f64>()
    }

    #[test]
    fn a_single_well_gives_the_eigenvector_average_as_its_rate_coefficient() {
        // One chemical eigenstate: k_(W->P) = (1/sqrt(Q)) M^-1 p_1 is the average of k(E) over the thermal
        // eigenvector (Gonzalez-Garcia, Olzmann 2010, after eq. 12), and the total loss k_W = Lambda_1.
        let network = ChemicalActivationNetwork { grain_width_cm1: 10.0, wells: vec![test_well("A", 400, 0, 300)] };
        let cond = Conditions { temperature_kelvin: 600.0, pressure_torr: 760.0 };
        let op = assemble_operator(&network, &cond, &final_options()).unwrap();
        let rates = phenomenological_rate_coefficients(&network, &op, None, EigenSolver::FullDecomposition).unwrap();
        let olzmann = thermal_rate_coefficients(&network, &op, EigenSolver::FullDecomposition, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        assert_eq!(rates.chemical_eigenvalues_s_inv.len(), 1);
        assert!((rates.well_to_bimolecular_s_inv[0][0] / olzmann.k_uni_s_inv - 1.0).abs() < 1e-8);
        assert!((rates.well_to_well_s_inv[0][0] / rates.chemical_eigenvalues_s_inv[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn well_rates_satisfy_the_loss_balance_and_detailed_balance() {
        // Two wells A <-> B, products P from both, an escape sink in B (300 K, 760 Torr).
        let network = two_well_network();
        let op = assemble_operator(&network, &conditions(), &final_options()).unwrap();
        let rates = phenomenological_rate_coefficients(&network, &op, None, EigenSolver::FullDecomposition).unwrap();
        assert_eq!(rates.wells, vec!["A".to_string(), "B".to_string()]);
        assert!(rates.bimolecular.contains(&"P".to_string()) && rates.bimolecular.contains(&"escape(B)".to_string()));
        // Georgievskii et al. (2013) eq. 29: k_i = sum_j k_(i->j) + sum_nu k_(i->nu) (an identity, because the
        // column sums of J are the loss rates).
        assert!(rates.loss_balance_max_deviation < 1e-8, "{}", rates.loss_balance_max_deviation);
        // Detailed balance of the isomerization: k_(A->B) Q_A = k_(B->A) Q_B.
        assert!(rates.detailed_balance_max_deviation < 1e-4, "{}", rates.detailed_balance_max_deviation);
        assert!(rates.well_to_well_s_inv[0][1] > 0.0 && rates.well_to_well_s_inv[1][0] > 0.0);
        // The chemical eigenvalues are well separated from the relaxation ones.
        assert!(rates.chemical_eigenvalues_s_inv[1] < 1e-2 * rates.relaxation_eigenvalue_s_inv);
    }

    #[test]
    fn at_high_pressure_the_rates_approach_the_high_pressure_limits() {
        let network = two_well_network();
        let high = Conditions { temperature_kelvin: 300.0, pressure_torr: 1.0e9 };
        let op = assemble_operator(&network, &high, &final_options()).unwrap();
        let rates = phenomenological_rate_coefficients(&network, &op, None, EigenSolver::FullDecomposition).unwrap();
        let kt = KB_CM * 300.0;
        let isomerization = network.wells[0].channels.iter().position(|c| c.name == "A->B").unwrap();
        let products = network.wells[0].channels.iter().position(|c| matches!(c.destination, ChannelDestination::Products { .. })).unwrap();
        let k_ab = high_pressure_rate(&network, 0, isomerization, kt);
        let k_ap = high_pressure_rate(&network, 0, products, kt);
        let p = rates.bimolecular.iter().position(|n| n == "P").unwrap();
        assert!((rates.well_to_well_s_inv[0][1] / k_ab - 1.0).abs() < 1e-3, "{} vs {k_ab}", rates.well_to_well_s_inv[0][1]);
        assert!((rates.well_to_bimolecular_s_inv[0][p] / k_ap - 1.0).abs() < 1e-3, "{} vs {k_ap}", rates.well_to_bimolecular_s_inv[0][p]);
    }

    #[test]
    fn reactant_rates_satisfy_detailed_balance_with_the_reverse_dissociation() {
        // The products of A are relabelled as the bimolecular reactant R with a capture rate coefficient k_c:
        // k_(R->i) / k_(i->R) = Q_i / Q_R for every well, i.e. k_c / k_inf(A->R) for A (Q_A/Q_R = k_c/k_inf,diss)
        // and (Q_B/Q_A) k_c / k_inf(A->R) for B.
        let mut network = two_well_network();
        let entrance = network.wells[0].channels.iter().position(|c| matches!(c.destination, ChannelDestination::Products { .. })).unwrap();
        network.wells[0].channels[entrance].destination = ChannelDestination::Products { name: "R".into() };
        let op = assemble_operator(&network, &conditions(), &final_options()).unwrap();
        let k_capture = 3.0e-11;
        let rates = phenomenological_rate_coefficients(&network, &op, Some(("R", k_capture)), EigenSolver::FullDecomposition).unwrap();
        let reactant = rates.reactant.as_ref().unwrap();
        let r = rates.bimolecular.iter().position(|n| n == "R").unwrap();
        let kt = KB_CM * conditions().temperature_kelvin;
        let k_diss = high_pressure_rate(&network, 0, entrance, kt);
        let q_ratio_ab = rates.partition_functions[1] / rates.partition_functions[0];
        let expected_a = k_capture / k_diss;
        assert!((reactant.to_well_cm3_s[0] / rates.well_to_bimolecular_s_inv[0][r] / expected_a - 1.0).abs() < 1e-3);
        // B is reached from R only through A: k_(R->B) is a small difference of the two chemical modes, which
        // amplifies the non-orthogonality of M (here ~3e-5, the relaxational projections), which limits detailed
        // balance (Georgievskii et al. 2013, Sec. IV): 0.2% here.
        assert!((reactant.to_well_cm3_s[1] / rates.well_to_bimolecular_s_inv[1][r] / (expected_a * q_ratio_ab) - 1.0).abs() < 1e-2);
        // Capture balance (eq. 22): the return to R closes it and is not negative.
        let total: f64 = reactant.to_well_cm3_s.iter().sum::<f64>() + reactant.to_bimolecular_cm3_s.iter().sum::<f64>();
        assert!((total / k_capture - 1.0).abs() < 1e-12);
        assert!(reactant.to_bimolecular_cm3_s[r] >= 0.0);
    }
}
