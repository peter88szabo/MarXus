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

/// How `ChemicalEigenvalueMax` (0 < value < 1) decides which of the N lowest eigenvectors are chemical.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ChemicalSubspaceCriterion {
    /// From the lowest, the eigenvectors are chemical while their projection on the relaxational subspace,
    /// 1 - F_ne = 1 - sum_i M_(i,lambda)^2, is at most the value: the "relaxation projection threshold" of
    /// MESS's direct method (`MasterEquation::direct_diagonalization_method`, `CalculationMethod direct`).
    #[default]
    RelaxationProjection,
    /// From the lowest, the eigenvalues are chemical while Lambda <= value x Lambda_(N+1): the rule of MESS's
    /// reaction-complex code (`ReactiveComplex::there_are_bound_groups`, `well_reduction_method`).
    EigenvalueRatio,
}

/// Species merging (G13 Sec. IV), with the criteria of MESS (`MasterEquation::direct_diagonalization_method`, and
/// `ReactiveComplex::threshold_well_partition` for the partition).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CseMerging {
    /// `ChemicalEigenvalueMax` (0 < value < 1): the threshold of `criterion`. If fewer than N eigenvectors are
    /// chemical, the wells are partitioned into that many species. MESS has no default (the keyword is required
    /// there); the MarXus default 0.2 is the value of the MESS decks of the validations. MESS's direct method reads a
    /// value above 1 as the threshold Lambda_(N+1)/Lambda >= value; MarXus refuses it (use EigenvalueRatio with
    /// 1/value).
    pub chemical_eigenvalue_max: f64,
    /// `WellProjectionThreshold`: wells whose projection on the chemical subspace is at least this are
    /// primary wells of the partition (MESS default 0.2).
    pub well_projection_threshold: f64,
    /// `ChemicalSubspaceCriterion` (MarXus block): relaxational projection (MESS direct method, default) or
    /// eigenvalue ratio.
    pub criterion: ChemicalSubspaceCriterion,
}

impl Default for CseMerging {
    fn default() -> Self {
        Self { chemical_eigenvalue_max: 0.2, well_projection_threshold: 0.2, criterion: ChemicalSubspaceCriterion::default() }
    }
}

/// Number of chemical eigenvectors by the relaxation projection threshold of MESS's direct method: counted from the
/// lowest while 1 - F_ne <= `threshold` (mess.cc, `direct_diagonalization_method`: `relaxation_projection[itemp] >
/// chemical_threshold` ends the count).
pub fn chemical_projection_count(relaxational_projections: &[f64], threshold: f64) -> usize {
    relaxational_projections.iter().take_while(|&&p| p <= threshold).count()
}

/// Number of chemical eigenvalues: the eigenvalues (ascending) Lambda_0 .. Lambda_(N-1) that are at most
/// `chemical_eigenvalue_max` x Lambda_N, counted from the lowest (MESS: `eigenval[itemp] > chemical_threshold *
/// relax_eval_min` ends the count, with relax_eval_min = eigenval[well_size()]).
pub fn chemical_eigenvalue_count(values: &[f64], n_wells: usize, chemical_eigenvalue_max: f64) -> usize {
    let relaxation = values[n_wells];
    values[..n_wells].iter().take_while(|&&v| v <= chemical_eigenvalue_max * relaxation).count()
}

/// The wells grouped into chemical species, and the wells in equilibrium with the bimolecular species.
#[derive(Debug, Clone, PartialEq)]
pub struct WellPartition {
    /// Every group: well indices, ascending; the groups are the merged species.
    pub groups: Vec<Vec<usize>>,
    /// Wells not in any group: in equilibrium with the bimolecular species (no rate coefficients of their own).
    pub bimolecular_group: Vec<usize>,
    /// Number of chemical eigenvalues minus the total projection of the groups on the chemical subspace.
    pub projection_error: f64,
}

/// Projection of a group of wells on the chemical subspace: |P_chem |g>|^2 with the group vector
/// |g> = sum_(w in g) sqrt(Q_w/Q_g) |w> (G13 eqs. 31-32): sum_lambda (sum_w sqrt(Q_w) M_(w,lambda))^2 / Q_g.
fn group_projection(group: &[usize], pop_chem: &[Vec<f64>], q: &[f64]) -> f64 {
    let q_group: f64 = group.iter().map(|&w| q[w]).sum();
    let chem_size = pop_chem.first().map_or(0, |row| row.len());
    (0..chem_size)
        .map(|l| {
            let x: f64 = group.iter().map(|&w| q[w].sqrt() * pop_chem[w][l]).sum();
            x * x
        })
        .sum::<f64>()
        / q_group
}

/// Partition of the wells into `chem_size` species, as MESS (`ReactiveComplex::threshold_well_partition`):
/// 1. the wells ordered by their own projection on the chemical subspace; the primary wells are taken in that
///    order while their projection is at least `well_projection_threshold`, and at least `chem_size` of them;
/// 2. of all partitions of the primary wells into `chem_size` groups, the one with the largest total
///    projection;
/// 3. the other wells, one at a time, into the group whose projection they raise most, while that raise is
///    not negative; the wells left over form the bimolecular group.
///
/// `pop_chem[w][lambda]` = M_(w,lambda) of the chemical eigenvectors (G13 eq. 25), `q` the Q_w of the wells.
pub fn partition_wells(pop_chem: &[Vec<f64>], q: &[f64], chem_size: usize, well_projection_threshold: f64) -> WellPartition {
    let n = pop_chem.len();
    let mut order: Vec<usize> = (0..n).collect();
    let own: Vec<f64> = (0..n).map(|w| group_projection(&[w], pop_chem, q)).collect();
    order.sort_by(|&a, &b| own[b].total_cmp(&own[a]).then(a.cmp(&b)));
    let mut primary = Vec::new();
    for &w in &order {
        if own[w] < well_projection_threshold && primary.len() >= chem_size {
            break;
        }
        primary.push(w);
    }

    // All partitions of the primary wells into exactly chem_size non-empty groups (restricted growth strings).
    let mut best: Option<(f64, Vec<Vec<usize>>)> = None;
    let mut labels = vec![0usize; primary.len()];
    fn visit(
        k: usize,
        used: usize,
        labels: &mut Vec<usize>,
        primary: &[usize],
        chem_size: usize,
        score: &dyn Fn(&[Vec<usize>]) -> f64,
        best: &mut Option<(f64, Vec<Vec<usize>>)>,
    ) {
        if k == primary.len() {
            if used == chem_size {
                let mut groups = vec![Vec::new(); chem_size];
                for (i, &label) in labels.iter().enumerate() {
                    groups[label].push(primary[i]);
                }
                groups.iter_mut().for_each(|g| g.sort());
                let projection = score(&groups);
                if best.as_ref().map_or(true, |(p, _)| projection > *p) {
                    *best = Some((projection, groups));
                }
            }
            return;
        }
        // Not enough elements left to open the missing groups.
        if chem_size - used > primary.len() - k {
            return;
        }
        for label in 0..(used + 1).min(chem_size) {
            labels[k] = label;
            visit(k + 1, used.max(label + 1), labels, primary, chem_size, score, best);
        }
    }
    let score = |groups: &[Vec<usize>]| groups.iter().map(|g| group_projection(g, pop_chem, q)).sum::<f64>();
    visit(0, 0, &mut labels, &primary, chem_size, &score, &mut best);
    let mut groups = best.map(|(_, g)| g).unwrap_or_default();

    // The other wells, the best (well, group) pair first, while the projection does not decrease.
    let mut rest: Vec<usize> = (0..n).filter(|w| !primary.contains(w)).collect();
    while !rest.is_empty() && !groups.is_empty() {
        let mut choice: Option<(f64, usize, usize)> = None;
        for (r, &w) in rest.iter().enumerate() {
            for (g, group) in groups.iter().enumerate() {
                let mut with = group.clone();
                with.push(w);
                let raise = group_projection(&with, pop_chem, q) - group_projection(group, pop_chem, q);
                if choice.map_or(true, |(best, _, _)| raise > best) {
                    choice = Some((raise, r, g));
                }
            }
        }
        match choice {
            Some((raise, r, g)) if raise >= 0.0 => {
                groups[g].push(rest.remove(r));
                groups[g].sort();
            }
            _ => break,
        }
    }
    let projection_error = chem_size as f64 - score(&groups);
    WellPartition { groups, bimolecular_group: rest, projection_error }
}

/// A prepared distribution F (a pulse) on the CSE description (G13 eqs. 13, 14, 24, 37-41): per unit of F, the populations
/// of the species after the internal-energy relaxation and the yields formed promptly, during the relaxation, in every
/// bimolecular channel and escape sink. The species then follow the phenomenological rate coefficients, so the long-time
/// yield of a channel x is prompt_x + sum_g n_g (fate of g in x), which equals k_x^T J^-1 F (exact with all eigenpairs).
#[derive(Debug, Clone)]
pub struct SourceProjection {
    pub name: String,
    /// n_g after relaxation, per unit of F, for every species (`PhenomenologicalRates::wells` order).
    pub species_populations: Vec<f64>,
    /// Prompt yields per unit of F (`PhenomenologicalRates::bimolecular` order).
    pub prompt_bimolecular: Vec<f64>,
}

/// Rate coefficients from a bimolecular species, the reactant or a product (cm3 s-1).
#[derive(Debug, Clone)]
pub struct ReactantRates {
    pub name: String,
    /// Capture (high-pressure association) rate coefficient k^(c).
    pub capture_cm3_s: f64,
    /// k_(R->i) for every species (G13 eq. 28).
    pub to_well_cm3_s: Vec<f64>,
    /// k_(R->nu) for every bimolecular channel (G13 eq. 21); the entry of the species itself is the return
    /// k_(R->R) from the capture balance (eq. 22).
    pub to_bimolecular_cm3_s: Vec<f64>,
}

/// Phenomenological rate coefficients of a network at one temperature and pressure.
#[derive(Debug, Clone)]
pub struct PhenomenologicalRates {
    /// The kinetic species: the wells, or, where wells are merged (G13 Sec. IV), the merged species, named by
    /// their wells joined with "+". All rate coefficients below refer to these species.
    pub wells: Vec<String>,
    /// The wells of every species (one well each without merging).
    pub well_groups: Vec<Vec<String>>,
    /// Wells in equilibrium with the bimolecular species (MESS's "bimolecular group"): no rate coefficients.
    pub bimolecular_group: Vec<String>,
    /// Number of chemical eigenvalues minus the projection of the species on the chemical subspace (0 without
    /// merging).
    pub partition_projection_error: f64,
    /// Bimolecular channels: product names (the reactant included when a channel leads to it), then the
    /// escape sinks "escape(W)" of the wells with a pseudo-first-order sink.
    pub bimolecular: Vec<String>,
    /// Q_i of the species (relative, common absolute energy scale); the sum over its wells (G13 eq. 32).
    pub partition_functions: Vec<f64>,
    /// The chemical eigenvalues (s-1), ascending: one per species.
    pub chemical_eigenvalues_s_inv: Vec<f64>,
    /// The lowest eigenvalue that is not chemical (s-1); NaN if there is none. Without merging this is
    /// Lambda_(N+1), the lowest relaxation eigenvalue.
    pub relaxation_eigenvalue_s_inv: f64,
    /// 1 - sum_i M_(i,lambda)^2 of every chemical eigenvector.
    pub relaxational_projection: Vec<f64>,
    /// [i][j] = k_(i->j) (s-1); the diagonal is the total loss k_i.
    pub well_to_well_s_inv: Vec<Vec<f64>>,
    /// [i][nu] = k_(i->nu) (s-1).
    pub well_to_bimolecular_s_inv: Vec<Vec<f64>>,
    pub reactant: Option<ReactantRates>,
    /// The rows of the bimolecular products given with a capture rate coefficient, in the order given.
    pub products: Vec<ReactantRates>,
    /// Projections of the prepared distributions given to `phenomenological_rate_coefficients_with_sources`.
    pub source_projections: Vec<SourceProjection>,
    /// The wells of the network, in its order: the rows of `kappa`.
    pub network_wells: Vec<String>,
    /// [w][nu] = kappa_(w,nu) of every well (not merged) and bimolecular channel (G13 eqs. 33-34): close to 1 if the well
    /// is in equilibrium with the bimolecular species nu, close to 0 otherwise; the sum runs over the eigenstates that
    /// are not chemical.
    pub kappa: Vec<Vec<f64>>,
    /// max_w |(1/sqrt(Q_w)) sum_(all lambda) M_(w,lambda) sum_nu p_lambda^(nu) / Lambda_lambda - 1|: the sum rule of eq. 34
    /// extended over every eigenstate and loss channel (the collisions conserve the Boltzmann distribution, so
    /// G f0 = -sum_nu K_nu f0).
    pub kappa_sum_rule_max_deviation: f64,
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
/// and its capture rate coefficient k^(c) (cm3 s-1); `products`: the same for bimolecular products whose rows
/// (k_(P->i), k_(P->nu)) are wanted, with the capture rate coefficient of the reverse association (1/Q_P =
/// k^(c)/sum k_(->P) f0, eq. 23). `solver`: a full decomposition (all eigenpairs).
/// `merging`: when fewer eigenvalues than wells are chemical, the wells are merged into species (G13 Sec. IV,
/// with the criteria of MESS; `CseMerging`).
pub fn phenomenological_rate_coefficients(
    network: &ChemicalActivationNetwork,
    op: &ChemicalActivationOperator,
    reactant: Option<(&str, f64)>,
    products: &[(&str, f64)],
    solver: EigenSolver,
    merging: &CseMerging,
) -> Result<PhenomenologicalRates, String> {
    phenomenological_rate_coefficients_with_sources(network, op, reactant, products, &[], solver, merging)
}

/// `phenomenological_rate_coefficients`, and the projection of prepared distributions (`sources`: name and F on the
/// states of `op`, e.g. from `project_source`) onto the species and the prompt products (`SourceProjection`).
pub fn phenomenological_rate_coefficients_with_sources(
    network: &ChemicalActivationNetwork,
    op: &ChemicalActivationOperator,
    reactant: Option<(&str, f64)>,
    products: &[(&str, f64)],
    sources: &[(&str, Vec<f64>)],
    solver: EigenSolver,
    merging: &CseMerging,
) -> Result<PhenomenologicalRates, String> {
    if !(merging.chemical_eigenvalue_max > 0.0 && merging.chemical_eigenvalue_max < 1.0) {
        return Err(format!(
            "CSE analysis: ChemicalEigenvalueMax = {} is outside 0 < value < 1. MESS's direct method reads a value \
             above 1 as the threshold Lambda_(N+1)/Lambda >= value; for that use ChemicalSubspaceCriterion \
             EigenvalueRatio with 1/value.",
            merging.chemical_eigenvalue_max
        ));
    }
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

    // M (G13 eq. 25) of the wells over all eigenstates (kappa, eq. 34, needs the relaxational ones), and over the
    // N lowest.
    let n_wells_all = n_wells;
    let mut m_all = vec![vec![0.0; vectors.len()]; n_wells_all];
    for (lambda, u) in vectors.iter().enumerate() {
        for (k, &(w, _)) in op.states.iter().enumerate() {
            m_all[w][lambda] += d[k] * u[k];
        }
    }
    for (w, row) in m_all.iter_mut().enumerate() {
        row.iter_mut().for_each(|x| *x /= q[w].sqrt());
    }
    let m_wells: Vec<Vec<f64>> = m_all.iter().map(|row| row[..n_wells_all].to_vec()).collect();

    // Chemical eigenvalues and the species (G13 Sec. IV; MESS direct method and threshold_well_partition).
    // 1 - F_ne of the N lowest eigenvectors (MESS: relaxation_projection[l] = 1 - |eigen_pop row l|^2)
    let projections: Vec<f64> =
        (0..n_wells_all).map(|l| 1.0 - (0..n_wells_all).map(|w| m_wells[w][l] * m_wells[w][l]).sum::<f64>()).collect();
    let chem_size = match merging.criterion {
        ChemicalSubspaceCriterion::RelaxationProjection => chemical_projection_count(&projections, merging.chemical_eigenvalue_max),
        ChemicalSubspaceCriterion::EigenvalueRatio => {
            chemical_eigenvalue_count(&values, n_wells_all, merging.chemical_eigenvalue_max)
        }
    };
    let partition = if chem_size == n_wells_all {
        WellPartition { groups: (0..n_wells_all).map(|w| vec![w]).collect(), bimolecular_group: Vec::new(), projection_error: 0.0 }
    } else if chem_size == 0 {
        WellPartition { groups: Vec::new(), bimolecular_group: (0..n_wells_all).collect(), projection_error: 0.0 }
    } else {
        let pop_chem: Vec<Vec<f64>> = m_wells.iter().map(|row| row[..chem_size].to_vec()).collect();
        partition_wells(&pop_chem, &q, chem_size, merging.well_projection_threshold)
    };
    let groups = &partition.groups;
    let n = groups.len();
    // Q of the species (eq. 32) and M of the species: the basis vector |g> = sum_(w in g) sqrt(Q_w/Q_g) |w>
    // (eq. 31), so M_(g,lambda) = sum_(w in g) sqrt(Q_w/Q_g) M_(w,lambda).
    let q_species: Vec<f64> = groups.iter().map(|g| g.iter().map(|&w| q[w]).sum()).collect();
    let m: Vec<Vec<f64>> = groups
        .iter()
        .zip(&q_species)
        .map(|(g, &qg)| (0..n).map(|l| g.iter().map(|&w| (q[w] / qg).sqrt() * m_wells[w][l]).sum()).collect())
        .collect();
    let q_wells = q;
    let q = q_species;

    // Its inverse and K = M Lambda M^-1.
    let m_inv = if n > 0 {
        invert_dense(&m).map_err(|e| format!("CSE analysis: the projection matrix M is singular ({e})."))?
    } else {
        Vec::new()
    };
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
    // kappa_(w,nu) (eq. 34) over the eigenstates that are not chemical, and its sum rule over all eigenstates.
    let kappa_sum = |w: usize, nu: usize, modes: std::ops::Range<usize>| -> f64 {
        modes.map(|l| m_all[w][l] * p[nu][l] / values[l]).sum::<f64>() / q_wells[w].sqrt()
    };
    let kappa: Vec<Vec<f64>> =
        (0..n_wells_all).map(|w| (0..bimolecular.len()).map(|nu| kappa_sum(w, nu, n..values.len())).collect()).collect();
    let kappa_sum_rule_max_deviation = (0..n_wells_all)
        .map(|w| ((0..bimolecular.len()).map(|nu| kappa_sum(w, nu, 0..values.len())).sum::<f64>() - 1.0).abs())
        .fold(0.0, f64::max);

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
    // Projection of every chemical eigenvector on the relaxational subspace (outside the thermal well vectors).
    let relaxational_projection: Vec<f64> =
        (0..n).map(|l| 1.0 - (0..n_wells_all).map(|w| m_wells[w][l] * m_wells[w][l]).sum::<f64>()).collect();
    let relaxation_eigenvalue = values.get(n).copied().unwrap_or(f64::NAN);
    let precision_floor = f64::EPSILON * symmetrized.band_matrix().max_abs_diagonal();
    let well_names: Vec<String> = network.wells.iter().map(|w| w.name.clone()).collect();
    let species_names: Vec<String> =
        groups.iter().map(|g| g.iter().map(|&w| well_names[w].clone()).collect::<Vec<_>>().join("+")).collect();
    let mut warnings = Vec::new();
    if n < n_wells_all {
        let eigenvalues: Vec<String> = values[..n_wells_all].iter().map(|v| format!("{v:.4e}")).collect();
        let free: Vec<String> = partition.bimolecular_group.iter().map(|&w| well_names[w].clone()).collect();
        let criterion = match merging.criterion {
            ChemicalSubspaceCriterion::RelaxationProjection => format!(
                "only {n} of the {n_wells_all} lowest eigenvectors (eigenvalues {} s-1) have a relaxational projection \
                 1 - F_ne ({}) of at most ChemicalEigenvalueMax = {} (relaxation projection threshold of the MESS direct \
                 method)",
                eigenvalues.join(", "),
                projections.iter().map(|p| format!("{p:.4}")).collect::<Vec<_>>().join(", "),
                merging.chemical_eigenvalue_max
            ),
            ChemicalSubspaceCriterion::EigenvalueRatio => format!(
                "only {n} of the {n_wells_all} lowest eigenvalues ({}) s-1 are at most ChemicalEigenvalueMax = {} x the \
                 lowest relaxation eigenvalue {:.4e} s-1",
                eigenvalues.join(", "),
                merging.chemical_eigenvalue_max,
                values[n_wells_all]
            ),
        };
        warnings.push(format!(
            "species merged (Georgievskii et al. 2013, Sec. IV; MESS partition): {criterion}, so the wells are not all \
             kinetically distinct. Species: {}{}. Partition projection error {:.3e}. Rate coefficients are given for \
             the merged species only; their wells cannot be distinguished at this condition.",
            if species_names.is_empty() { "none".to_string() } else { species_names.join(", ") },
            if free.is_empty() {
                String::new()
            } else {
                format!("; in equilibrium with the bimolecular species (no rate coefficients): {}", free.join(", "))
            },
            partition.projection_error
        ));
    }
    if n > 0 {
        let separation = values[n - 1] / relaxation_eigenvalue;
        if separation > CSE_SEPARATION_WARNING {
            warnings.push(format!(
                "the largest chemical eigenvalue {:e} s-1 is {separation:.3} of the lowest relaxation eigenvalue {relaxation_eigenvalue:e} \
                 s-1: the chemical eigenstates are not well separated (the method requires |lambda_N| << |lambda_(N+1)|, \
                 Miller, Klippenstein 2006). Species are merged only above ChemicalEigenvalueMax = {} (Georgievskii et \
                 al. 2013, Sec. IV).",
                values[n - 1],
                merging.chemical_eigenvalue_max
            ));
        }
    }
    for (l, &value) in values.iter().take(n).enumerate() {
        if value < 100.0 * precision_floor {
            warnings.push(format!(
                "chemical eigenvalue {l} = {value:e} s-1 is within a factor 100 of the double-precision floor \
                 eps*max(S_ii) = {precision_floor:e} s-1: the rate coefficients of this mode are not resolved."
            ));
        }
    }

    // Rows of a bimolecular species with capture rate coefficient k^(c) (eqs. 28, 21, 22); 1/Q_nu = k^(c)/sum k_(->nu) f0 (eq. 23).
    // Eq. 21 runs over every eigenstate that is not chemical.
    let source_rates = |name: &str, capture: f64, role: &str| -> Result<ReactantRates, String> {
        let r = bimolecular
            .iter()
            .position(|b| b == name)
            .ok_or_else(|| format!("CSE analysis: no channel leads to the {role} '{name}'."))?;
        let denominator: f64 = (0..n_states).map(|k| k_into[r][k] * d[k] * d[k]).sum();
        if !(denominator > 0.0) {
            return Err(format!("CSE analysis: the {role} '{name}' has no association flux."));
        }
        let to_well: Vec<f64> =
            (0..n).map(|i| capture * q[i].sqrt() * (0..n).map(|l| m[i][l] * p[r][l]).sum::<f64>() / denominator).collect();
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
        Ok(ReactantRates { name: name.to_string(), capture_cm3_s: capture, to_well_cm3_s: to_well, to_bimolecular_cm3_s: to_bimolecular })
    };
    let reactant = reactant.map(|(name, capture)| source_rates(name, capture, "reactant")).transpose()?;
    let products = products
        .iter()
        .map(|&(name, capture)| source_rates(name, capture, "bimolecular product"))
        .collect::<Result<Vec<_>, _>>()?;

    // Prepared distributions (pulses) on the CSE description: the amplitudes c_l = <f^l|F> = sum_k u_lk F_k/d_k
    // (eqs. 12-13, in logarithms since d spans many orders of magnitude), the species populations after the relaxation
    // n_g = sqrt(Q_g) sum_(l chem) M_(g,l) c_l (eqs. 24, 37) and the prompt yields sum_(l relax) p_l^(nu) c_l/Lambda_l
    // (eqs. 38, 41).
    let half_log_d = &symmetrized.half_log_d;
    let mut source_projections = Vec::new();
    for (name, on_states) in sources {
        if on_states.len() != n_states {
            return Err(format!("CSE analysis: the source '{name}' has {} states, the operator {n_states}.", on_states.len()));
        }
        let c: Vec<f64> = vectors
            .iter()
            .map(|u| {
                (0..n_states)
                    .filter(|&k| on_states[k] != 0.0 && u[k] != 0.0)
                    .map(|k| u[k].signum() * on_states[k].signum() * (u[k].abs().ln() + on_states[k].abs().ln() - half_log_d[k]).exp())
                    .sum()
            })
            .collect();
        let species_populations = (0..n).map(|g| q[g].sqrt() * (0..n).map(|l| m[g][l] * c[l]).sum::<f64>()).collect();
        let prompt_bimolecular =
            (0..bimolecular.len()).map(|nu| (n..values.len()).map(|l| p[nu][l] * c[l] / values[l]).sum()).collect();
        source_projections.push(SourceProjection { name: name.to_string(), species_populations, prompt_bimolecular });
    }

    Ok(PhenomenologicalRates {
        wells: species_names,
        well_groups: groups.iter().map(|g| g.iter().map(|&w| well_names[w].clone()).collect()).collect(),
        bimolecular_group: partition.bimolecular_group.iter().map(|&w| well_names[w].clone()).collect(),
        partition_projection_error: partition.projection_error,
        bimolecular,
        partition_functions: q,
        chemical_eigenvalues_s_inv: values[..n].to_vec(),
        relaxation_eigenvalue_s_inv: relaxation_eigenvalue,
        relaxational_projection,
        well_to_well_s_inv: well_to_well,
        well_to_bimolecular_s_inv: well_to_bimolecular,
        reactant,
        products,
        source_projections,
        network_wells: well_names,
        kappa,
        kappa_sum_rule_max_deviation,
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
pub(crate) mod tests {
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

    // ---- Species merging (G13 Sec. IV; MESS `threshold_well_partition`) ----

    #[test]
    fn the_number_of_chemical_eigenvalues_follows_the_mess_threshold() {
        // Chemical eigenvalues are those with Lambda <= ChemicalEigenvalueMax x Lambda_(N+1), counted from the
        // lowest; Lambda_(N+1) is the lowest relaxation eigenvalue for N wells (MESS, direct method).
        let values = [1.0, 50.0, 300.0, 1000.0];
        assert_eq!(chemical_eigenvalue_count(&values, 3, 0.2), 2);
        assert_eq!(chemical_eigenvalue_count(&values, 3, 0.01), 1);
        assert_eq!(chemical_eigenvalue_count(&values, 3, 0.5), 3);
        assert_eq!(chemical_eigenvalue_count(&values, 3, 0.0005), 0);
    }

    #[test]
    fn wells_in_equilibrium_form_one_group_and_a_decoupled_well_joins_the_bimolecular_group() {
        // One chemical eigenvector, in equilibrium between wells 0 and 1 (Q = 3 and 1), absent from well 2.
        let q = [3.0, 1.0, 2.0];
        let pop = vec![vec![(3.0f64 / 4.0).sqrt()], vec![(1.0f64 / 4.0).sqrt()], vec![0.0]];
        let p = partition_wells(&pop, &q, 1, 0.2);
        assert_eq!(p.groups, vec![vec![0, 1]]);
        assert_eq!(p.bimolecular_group, vec![2]);
        assert!(p.projection_error.abs() < 1e-12, "{}", p.projection_error);
    }

    #[test]
    fn the_partition_with_the_largest_projection_is_chosen() {
        // Two chemical eigenvectors: well 2 alone, and wells 0 + 1 in equilibrium.
        let q = [3.0, 1.0, 2.0];
        let pop = vec![vec![0.0, (3.0f64 / 4.0).sqrt()], vec![0.0, (1.0f64 / 4.0).sqrt()], vec![1.0, 0.0]];
        let p = partition_wells(&pop, &q, 2, 0.2);
        let mut groups = p.groups.clone();
        groups.sort();
        assert_eq!(groups, vec![vec![0, 1], vec![2]]);
        assert!(p.bimolecular_group.is_empty());
        assert!(p.projection_error.abs() < 1e-12, "{}", p.projection_error);
    }

    #[test]
    fn a_well_with_a_small_projection_joins_the_group_whose_projection_it_raises() {
        // One eigenvector in equilibrium among all three wells. Well 2 has a small weight: its own projection
        // 0.01 lies below the threshold (a secondary well), and adding it raises the group projection to 1.
        let q = [3.0, 0.96, 0.04];
        let total: f64 = q.iter().sum();
        let pop: Vec<Vec<f64>> = q.iter().map(|x| vec![(x / total).sqrt()]).collect();
        let p = partition_wells(&pop, &q, 1, 0.2);
        assert_eq!(p.groups, vec![vec![0, 1, 2]]);
        assert!(p.bimolecular_group.is_empty());
        assert!(p.projection_error.abs() < 1e-12, "{}", p.projection_error);
    }

    /// Two wells A, B with a low isomerization barrier (absolute grain 100, k = W#/(h rho) both ways): they
    /// equilibrate much faster than they react (product thresholds at grains 300 and 340).
    pub(crate) fn fast_equilibrium_network() -> ChemicalActivationNetwork {
        let mut network = two_well_network();
        let planck_cm1_s = 3.3356e-11;
        let ts: isize = 100;
        let w_ts = |absolute: isize| if absolute < ts { 0.0 } else { (1.0 + 0.05 * (absolute - ts) as f64).powi(6) };
        for (w, other) in [(0usize, 1usize), (1, 0)] {
            let offset = network.wells[w].bottom_offset_grains;
            let rho = network.wells[w].density_of_states.clone();
            let channel = network.wells[w]
                .channels
                .iter_mut()
                .find(|c| c.destination == ChannelDestination::Well { index: other })
                .unwrap();
            channel.rate_constant_s_inv = (0..rho.len()).map(|i| w_ts(i as isize + offset) / (planck_cm1_s * rho[i])).collect();
        }
        network
    }

    #[test]
    fn by_default_an_eigenvector_is_chemical_while_its_relaxational_projection_is_below_the_threshold() {
        // MESS direct method (mess.cc, direct_diagonalization_method, "relaxation projection threshold"): for
        // 0 < ChemicalEigenvalueMax < 1 the eigenvectors, from the lowest, are chemical while 1 - F_ne <= the value.
        assert_eq!(ChemicalSubspaceCriterion::default(), ChemicalSubspaceCriterion::RelaxationProjection);
        assert_eq!(chemical_projection_count(&[0.05, 0.15, 0.3], 0.2), 2);
        assert_eq!(chemical_projection_count(&[0.25, 0.1], 0.2), 0);
        let network = fast_equilibrium_network();
        let op = assemble_operator(&network, &conditions(), &final_options()).unwrap();
        let rates = |threshold: f64| {
            let merging = CseMerging { chemical_eigenvalue_max: threshold, ..CseMerging::default() };
            phenomenological_rate_coefficients(&network, &op, None, &[], EigenSolver::FullDecomposition, &merging).unwrap()
        };
        let separate = rates(0.999);
        let p = separate.relaxational_projection.clone();
        assert!(p[0] > 0.0 && p[1] > 10.0 * p[0], "the test needs 1 - F_ne of the second eigenvector well above the first: {p:?}");
        assert_eq!(rates(0.5 * p[0]).wells, Vec::<String>::new());
        assert_eq!(rates((p[0] * p[1]).sqrt()).wells, vec!["A+B".to_string()]);
        assert_eq!(separate.wells.len(), 2);
        let merged = rates((p[0] * p[1]).sqrt());
        assert!(merged.warnings.iter().any(|w| w.contains("1 - F_ne")), "{:?}", merged.warnings);
    }

    #[test]
    fn wells_in_fast_equilibrium_are_merged_and_the_group_decays_with_k_uni() {
        let network = fast_equilibrium_network();
        let op = assemble_operator(&network, &conditions(), &final_options()).unwrap();
        let no_merging = CseMerging { chemical_eigenvalue_max: 0.999, ..CseMerging::default() };
        let separate = phenomenological_rate_coefficients(&network, &op, None, &[], EigenSolver::FullDecomposition, &no_merging).unwrap();
        assert_eq!(separate.wells, vec!["A".to_string(), "B".to_string()]);
        let (l1, l2, l3) = (
            separate.chemical_eigenvalues_s_inv[0],
            separate.chemical_eigenvalues_s_inv[1],
            separate.relaxation_eigenvalue_s_inv,
        );
        assert!(l2 > 10.0 * l1, "the test needs Lambda_2 well above Lambda_1: {l1:e} {l2:e} {l3:e}");
        // A threshold between the two leaves one chemical eigenvalue: A and B are merged.
        let merging = CseMerging { chemical_eigenvalue_max: (l1 * l2).sqrt() / l3, criterion: ChemicalSubspaceCriterion::EigenvalueRatio, ..CseMerging::default() };
        let merged = phenomenological_rate_coefficients(&network, &op, None, &[], EigenSolver::FullDecomposition, &merging).unwrap();
        assert_eq!(merged.wells, vec!["A+B".to_string()]);
        assert_eq!(merged.well_groups, vec![vec!["A".to_string(), "B".to_string()]]);
        assert!(merged.bimolecular_group.is_empty());
        assert_eq!(merged.chemical_eigenvalues_s_inv.len(), 1);
        assert!(merged.warnings.iter().any(|w| w.contains("merged")), "{:?}", merged.warnings);
        // Q of the merged species is the sum (G13 eq. 32).
        let q_sum = separate.partition_functions[0] + separate.partition_functions[1];
        assert!((merged.partition_functions[0] / q_sum - 1.0).abs() < 1e-12);
        // One chemical eigenstate: the rate of the merged species into the bimolecular channels is the thermal
        // k_uni, the average of k(E) over the thermal eigenvector of J (GO10 after eq. 12).
        let th = thermal_rate_coefficients(&network, &op, EigenSolver::FullDecomposition, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        let out: f64 = merged.well_to_bimolecular_s_inv[0].iter().sum();
        assert!((out / th.k_uni_s_inv - 1.0).abs() < 1e-8, "{out:e} vs {:e}", th.k_uni_s_inv);
        assert!((merged.well_to_well_s_inv[0][0] / out - 1.0).abs() < 1e-8);
    }

    #[test]
    fn separated_wells_are_not_merged_with_the_default_threshold() {
        let network = two_well_network();
        let op = assemble_operator(&network, &conditions(), &final_options()).unwrap();
        let rates = phenomenological_rate_coefficients(&network, &op, None, &[], EigenSolver::FullDecomposition, &CseMerging::default()).unwrap();
        assert_eq!(rates.wells, vec!["A".to_string(), "B".to_string()]);
        assert!(rates.bimolecular_group.is_empty());
        assert!(!rates.warnings.iter().any(|w| w.contains("merged")));
    }

    #[test]
    fn a_single_well_gives_the_eigenvector_average_as_its_rate_coefficient() {
        // One chemical eigenstate: k_(W->P) = (1/sqrt(Q)) M^-1 p_1 is the average of k(E) over the thermal
        // eigenvector (Gonzalez-Garcia, Olzmann 2010, after eq. 12), and the total loss k_W = Lambda_1.
        let network = ChemicalActivationNetwork { grain_width_cm1: 10.0, wells: vec![test_well("A", 400, 0, 300)] };
        let cond = Conditions { temperature_kelvin: 600.0, pressure_torr: 760.0 };
        let op = assemble_operator(&network, &cond, &final_options()).unwrap();
        let rates = phenomenological_rate_coefficients(&network, &op, None, &[], EigenSolver::FullDecomposition, &CseMerging::default()).unwrap();
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
        let rates = phenomenological_rate_coefficients(&network, &op, None, &[], EigenSolver::FullDecomposition, &CseMerging::default()).unwrap();
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
        let rates = phenomenological_rate_coefficients(&network, &op, None, &[], EigenSolver::FullDecomposition, &CseMerging::default()).unwrap();
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
        let rates = phenomenological_rate_coefficients(&network, &op, Some(("R", k_capture)), &[], EigenSolver::FullDecomposition, &CseMerging::default()).unwrap();
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


    // ---- Product rows (G13 eqs. 28, 21, 22 for every bimolecular species) and kappa (eqs. 33-34) ----

    /// The two-well network with A's products relabelled as the reactant R; B's products P are then formed from B only.
    fn reactant_and_product_network() -> ChemicalActivationNetwork {
        let mut network = two_well_network();
        network.wells[0].channels[0].destination = ChannelDestination::Products { name: "R".into() };
        network
    }

    #[test]
    fn product_rows_satisfy_detailed_balance_and_the_capture_balance() {
        let network = reactant_and_product_network();
        let op = assemble_operator(&network, &conditions(), &final_options()).unwrap();
        let (k_r, k_p) = (3.0e-11, 2.0e-12);
        let rates = phenomenological_rate_coefficients(
            &network,
            &op,
            Some(("R", k_r)),
            &[("P", k_p)],
            EigenSolver::FullDecomposition,
            &CseMerging::default(),
        )
        .unwrap();
        assert_eq!(rates.products.len(), 1);
        let product = &rates.products[0];
        assert_eq!(product.name, "P");
        assert_eq!(product.capture_cm3_s, k_p);
        let (r, p) = (
            rates.bimolecular.iter().position(|n| n == "R").unwrap(),
            rates.bimolecular.iter().position(|n| n == "P").unwrap(),
        );
        let kt = KB_CM * conditions().temperature_kelvin;
        let (k_diss_a, k_diss_b) = (high_pressure_rate(&network, 0, 0, kt), high_pressure_rate(&network, 1, 0, kt));
        // k_(P->B) / k_(B->P) = Q_B/Q_P = k_p / k_inf(B->P), B forming P directly.
        assert!((product.to_well_cm3_s[1] / rates.well_to_bimolecular_s_inv[1][p] / (k_p / k_diss_b) - 1.0).abs() < 1e-3);
        // Capture balance (eq. 22): the return to P closes it and is not negative.
        let total: f64 = product.to_well_cm3_s.iter().sum::<f64>() + product.to_bimolecular_cm3_s.iter().sum::<f64>();
        assert!((total / k_p - 1.0).abs() < 1e-12);
        assert!(product.to_bimolecular_cm3_s[p] >= 0.0);
        // Bimolecular-to-bimolecular (eq. 21, symmetric in nu and mu): k_(R->P) Q_R = k_(P->R) Q_P with
        // Q_R = Q_A k_inf(A->R)/k_r and Q_P = Q_B k_inf(B->P)/k_p.
        let reactant = rates.reactant.as_ref().unwrap();
        let (q_a, q_b) = (rates.partition_functions[0], rates.partition_functions[1]);
        let (q_r, q_p) = (q_a * k_diss_a / k_r, q_b * k_diss_b / k_p);
        assert!((reactant.to_bimolecular_cm3_s[p] * q_r / (product.to_bimolecular_cm3_s[r] * q_p) - 1.0).abs() < 1e-9);
        // The reactant row is the same as without products.
        let alone = phenomenological_rate_coefficients(&network, &op, Some(("R", k_r)), &[], EigenSolver::FullDecomposition, &CseMerging::default())
            .unwrap();
        assert!(alone.products.is_empty());
        assert_eq!(alone.reactant.as_ref().unwrap().to_well_cm3_s, reactant.to_well_cm3_s);
    }

    #[test]
    fn a_product_without_a_channel_is_an_error() {
        let network = reactant_and_product_network();
        let op = assemble_operator(&network, &conditions(), &final_options()).unwrap();
        let e = phenomenological_rate_coefficients(&network, &op, None, &[("Q", 1.0e-12)], EigenSolver::FullDecomposition, &CseMerging::default())
            .unwrap_err();
        assert!(e.contains("'Q'"), "{e}");
    }

    #[test]
    fn kappa_satisfies_the_sum_rule_and_vanishes_for_chemically_distinct_wells() {
        // G13 eq. 34 over the relaxational eigenstates; with the chemical ones the sum over every loss channel is 1
        // for every well (the collisions conserve the Boltzmann distribution).
        let network = reactant_and_product_network();
        let op = assemble_operator(&network, &conditions(), &final_options()).unwrap();
        let rates = phenomenological_rate_coefficients(&network, &op, None, &[], EigenSolver::FullDecomposition, &CseMerging::default()).unwrap();
        assert_eq!(rates.wells, vec!["A", "B"]);
        assert!(rates.kappa_sum_rule_max_deviation < 1e-8, "{}", rates.kappa_sum_rule_max_deviation);
        assert_eq!((rates.kappa.len(), rates.kappa[0].len()), (2, rates.bimolecular.len()));
        for row in &rates.kappa {
            for &kappa in row {
                assert!(kappa.abs() < 0.05, "{:?}", rates.kappa);
            }
        }
    }

    #[test]
    fn wells_in_equilibrium_with_the_bimolecular_species_have_kappas_summing_to_one() {
        // No chemical eigenvalue (both wells in the bimolecular group): every eigenstate is relaxational, so the
        // sum rule is the sum of the kappas themselves.
        let network = reactant_and_product_network();
        let op = assemble_operator(&network, &conditions(), &final_options()).unwrap();
        let separate = phenomenological_rate_coefficients(&network, &op, None, &[], EigenSolver::FullDecomposition, &CseMerging::default()).unwrap();
        let none = CseMerging {
            chemical_eigenvalue_max: 0.5 * separate.chemical_eigenvalues_s_inv[0] / separate.relaxation_eigenvalue_s_inv,
            criterion: ChemicalSubspaceCriterion::EigenvalueRatio,
            ..CseMerging::default()
        };
        let rates = phenomenological_rate_coefficients(&network, &op, None, &[], EigenSolver::FullDecomposition, &none).unwrap();
        assert_eq!(rates.bimolecular_group, vec!["A", "B"]);
        for row in &rates.kappa {
            assert!((row.iter().sum::<f64>() - 1.0).abs() < 1e-8, "{:?}", rates.kappa);
        }
    }


    // ---- Projection of a prepared distribution (G13 eqs. 13, 14, 24, 37-41) ----

    #[test]
    fn the_projection_of_the_reactant_shape_gives_the_reactant_rows_per_capture() {
        // F = k_(->R) f0 / sum k_(->R) f0 has amplitudes c = p^(R)/sum k f0, so the species populations are k_(R->i)/k_c
        // (eq. 28) and the prompt yields k_(R->mu)/k_c (eq. 21).
        use crate::masterequation::chemical_activation_sources::thermal_entrance_source;
        use crate::masterequation::chemical_activation_steady_state::project_source;
        let network = reactant_and_product_network();
        let op = assemble_operator(&network, &conditions(), &final_options()).unwrap();
        let kt = KB_CM * conditions().temperature_kelvin;
        let f = thermal_entrance_source(&network, &[(0, 0)], kt).unwrap();
        let on_states = project_source(&op, &f).unwrap().on_states;
        let k_c = 3.0e-11;
        let rates = phenomenological_rate_coefficients_with_sources(
            &network,
            &op,
            Some(("R", k_c)),
            &[],
            &[("thermal R", on_states)],
            EigenSolver::FullDecomposition,
            &CseMerging::default(),
        )
        .unwrap();
        let reactant = rates.reactant.as_ref().unwrap();
        let projection = &rates.source_projections[0];
        assert_eq!(projection.name, "thermal R");
        for i in 0..rates.wells.len() {
            assert!((projection.species_populations[i] / (reactant.to_well_cm3_s[i] / k_c) - 1.0).abs() < 1e-9, "{i}");
        }
        let r = rates.bimolecular.iter().position(|b| b == "R").unwrap();
        for mu in (0..rates.bimolecular.len()).filter(|&mu| mu != r) {
            let expected = reactant.to_bimolecular_cm3_s[mu] / k_c;
            assert!((projection.prompt_bimolecular[mu] - expected).abs() <= 1e-9 * expected.abs() + 1e-300, "{mu}");
        }
    }

    #[test]
    fn projected_species_and_prompt_yields_reproduce_the_steady_state_yields_of_a_hot_source() {
        // The pulse identity Y_x = k_x^T J^-1 F (exact with all eigenpairs): the prompt yields plus the species
        // populations times their thermal fates (absorbing chain of the rate coefficients).
        use crate::masterequation::chemical_activation_driver::{run_chemical_activation, ChemicalActivationRun, SourceSpecification};
        use crate::masterequation::chemical_activation_steady_state::{project_source, LinearSolver};
        use crate::masterequation::prepared_distributions::{gaussian, EnergyReference, Representation};
        let network = reactant_and_product_network();
        let op = assemble_operator(&network, &conditions(), &final_options()).unwrap();
        let hot = gaussian(&network, 0, 3300.0, 150.0, Representation::Density, EnergyReference::AboveWellGround).unwrap();
        let on_states = project_source(&op, &hot.mass).unwrap().on_states;
        let rates = phenomenological_rate_coefficients_with_sources(
            &network,
            &op,
            Some(("R", 3.0e-11)),
            &[],
            &[("hot", on_states)],
            EigenSolver::FullDecomposition,
            &CseMerging::default(),
        )
        .unwrap();
        let projection = &rates.source_projections[0];
        let fates = reactant_yields(&rates).unwrap().unwrap().well_fates;
        let run = ChemicalActivationRun {
            temperatures_kelvin: vec![conditions().temperature_kelvin],
            pressures_torr: vec![conditions().pressure_torr],
            options: final_options(),
            solver: LinearSolver::BandedCholesky,
            source: SourceSpecification::Fixed(hot.mass.clone()),
            tolerance: 1e-8,
        };
        let steady = run_chemical_activation(&network, &run).unwrap().remove(0).result;
        for (x, name) in rates.bimolecular.iter().enumerate() {
            let expected: f64 = if let Some(well) = name.strip_prefix("escape(").and_then(|n| n.strip_suffix(')')) {
                steady.wells.iter().find(|w| w.name == well).unwrap().bimolecular_sink_yield
            } else {
                steady.channels.iter().filter(|c| matches!(&c.destination, ChannelDestination::Products { name: n } if n == name)).map(|c| c.flux).sum()
            };
            let cse = projection.prompt_bimolecular[x]
                + (0..rates.wells.len()).map(|i| projection.species_populations[i] * fates[i][x]).sum::<f64>();
            assert!((cse - expected).abs() < 1e-7 * expected.max(1e-3), "{name}: CSE {cse:e} vs steady state {expected:e}");
        }
    }
}
