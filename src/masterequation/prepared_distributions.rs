//! Prepared energy distributions on the master-equation grid: the distribution layer of non-thermal
//! preparations and sources (design note "Preparation, Sources, and Experimental Observables in MarXus",
//! papers/Reactant_flux_initiation/MarXus_Nonthermal_Sources.tex, 7 October 2026, cited as N; Sections 5-6).
//!
//! A `GrainDistribution` is a probability mass per grain of every well, normalized over all wells together
//! (N eq. wells: normalizing every well separately and concatenating would give every well the same
//! population). Inputs declare their representation (N Sec. 5.1):
//!   - bin mass: probability per input bin;
//!   - density f(E) (per cm-1): F_i = integral of f over the grain (N eq. binmass);
//!   - per-state weight h(E): f(E) = rho(E) h(E) / normalization (N eq. perstate);
//! and their energy reference (N Sec. 5.2): above the well bottom (grain 0 of the well, centred at 0) or
//! absolute (the common energy scale of the network, grain g of well w centred at (g + offset_w) dE).
//! Grain i of a well covers [(i - 1/2) dE, (i + 1/2) dE] on its own scale. Input bins are mapped onto the
//! grains by their overlaps (conservative rebinning, N eq. rebin); probability outside the grid of the well
//! is reported as the lost fraction before renormalization (N eq. after rebin) and is an error unless
//! truncation is allowed (an intentional conditioning on the represented interval).

use super::chemical_activation_network::ChemicalActivationNetwork;
use super::chemical_activation_sources::thermal_distribution;
use crate::constants::KB_CM;
pub use crate::numeric::special_functions::{erf, normal_interval_probability};

/// What the values of an input distribution are (N Sec. 5.1).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Representation {
    /// Probability mass per input bin.
    BinMass,
    /// Probability density f(E) per cm-1, constant inside each input bin.
    Density,
    /// Per-state weight h(E): the distribution is rho(E) h(E), normalized (N eq. perstate).
    PerStateWeight,
}

/// Energy origin of an input distribution (N Sec. 5.2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyReference {
    /// Above the bottom of the receiving well (the centre of its grain 0).
    AboveWellGround,
    /// The common energy scale of the network: grain i of well w is centred at (i + offset_w) dE.
    Absolute,
}

/// Whether probability outside the grid of the well is allowed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Support {
    /// The input lies on the grid; probability outside it (beyond `MAX_LOST_FRACTION`) is an error.
    Complete,
    /// Probability outside the grid is dropped (an intentional conditioning) and reported as lost.
    Truncate,
}

/// Lost fraction tolerated with `Support::Complete` (rounding at the edges of the grid).
pub const MAX_LOST_FRACTION: f64 = 1.0e-10;

/// Probability mass per grain of every well, normalized over all wells together, with the fraction of the input
/// that lay outside the grid (before renormalization) and a description of its construction.
#[derive(Debug, Clone, PartialEq)]
pub struct GrainDistribution {
    pub mass: Vec<Vec<f64>>,
    pub lost_fraction: f64,
    pub description: String,
}

impl GrainDistribution {
    /// Normalize non-negative masses over all wells.
    fn normalized(mut mass: Vec<Vec<f64>>, lost_fraction: f64, description: String) -> Result<Self, String> {
        if mass.iter().flatten().any(|x| !(*x >= 0.0) || !x.is_finite()) {
            return Err(format!("{description}: the distribution must be finite and >= 0."));
        }
        let total: f64 = mass.iter().flatten().sum();
        if !(total > 0.0) {
            return Err(format!("{description}: the distribution is zero on the grid."));
        }
        mass.iter_mut().flatten().for_each(|x| *x /= total);
        Ok(Self { mass, lost_fraction, description })
    }

    /// Fraction of the distribution in every well.
    pub fn well_fractions(&self) -> Vec<f64> {
        self.mass.iter().map(|m| m.iter().sum()).collect()
    }

    /// Mean energy above the bottom of every well (cm-1); NaN for a well without mass.
    pub fn mean_energy_above_bottom_cm1(&self, network: &ChemicalActivationNetwork) -> Vec<f64> {
        let de = network.grain_width_cm1;
        self.mass
            .iter()
            .map(|m| {
                let total: f64 = m.iter().sum();
                if total > 0.0 {
                    m.iter().enumerate().map(|(i, x)| i as f64 * de * x).sum::<f64>() / total
                } else {
                    f64::NAN
                }
            })
            .collect()
    }

    /// Mean energy on the common absolute scale of the network (cm-1).
    pub fn mean_absolute_energy_cm1(&self, network: &ChemicalActivationNetwork) -> f64 {
        self.mass
            .iter()
            .enumerate()
            .map(|(w, m)| m.iter().enumerate().map(|(i, x)| network.absolute_energy_cm1(w, i) * x).sum::<f64>())
            .sum()
    }
}

/// Conservative rebinning (N eq. rebin): every input bin [lo, hi] with mass m gives m |overlap|/(hi - lo) to
/// each output bin [start + j w, start + (j + 1) w), j < n; a zero-width bin (a stick) goes to the output bin
/// that contains it. Returns the output masses and the mass outside the output range.
pub fn rebin(input: &[(f64, f64)], masses: &[f64], start: f64, width: f64, n: usize) -> (Vec<f64>, f64) {
    let mut out = vec![0.0; n];
    let mut assigned = 0.0;
    let end = start + n as f64 * width;
    for (&(lo, hi), &m) in input.iter().zip(masses) {
        if m == 0.0 {
            continue;
        }
        if hi == lo {
            if lo >= start && lo < end {
                let j = (((lo - start) / width).floor() as usize).min(n - 1);
                out[j] += m;
                assigned += m;
            }
            continue;
        }
        let first = (((lo.max(start) - start) / width).floor().max(0.0)) as usize;
        let last = ((((hi.min(end) - start) / width).ceil()).max(0.0) as usize).min(n);
        for (j, slot) in out.iter_mut().enumerate().take(last).skip(first) {
            let (a, b) = (start + j as f64 * width, start + (j + 1) as f64 * width);
            let overlap = hi.min(b) - lo.max(a);
            if overlap > 0.0 {
                let part = m * (overlap / (hi - lo));
                *slot += part;
                assigned += part;
            }
        }
    }
    let total: f64 = masses.iter().sum();
    (out, (total - assigned).max(0.0))
}

/// Start of the grain grid of well `w` on the scale of `reference`: grain i covers [start + i dE, start + (i+1) dE).
fn grid_start(network: &ChemicalActivationNetwork, w: usize, reference: EnergyReference) -> f64 {
    let de = network.grain_width_cm1;
    match reference {
        EnergyReference::AboveWellGround => -0.5 * de,
        EnergyReference::Absolute => (network.wells[w].bottom_offset_grains as f64 - 0.5) * de,
    }
}

fn check_well(network: &ChemicalActivationNetwork, w: usize, what: &str) -> Result<(), String> {
    if w >= network.wells.len() {
        return Err(format!("{what}: well index {w} is not in the network ({} wells).", network.wells.len()));
    }
    Ok(())
}

fn empty_mass(network: &ChemicalActivationNetwork) -> Vec<Vec<f64>> {
    network.wells.iter().map(|well| vec![0.0; well.grain_count()]).collect()
}

/// Thermal population of well `w` at its own preparation temperature (N eq. boltz), independent of the bath.
pub fn thermal(network: &ChemicalActivationNetwork, w: usize, preparation_temperature_kelvin: f64) -> Result<GrainDistribution, String> {
    let what = format!("Thermal preparation at {preparation_temperature_kelvin} K");
    check_well(network, w, &what)?;
    if !(preparation_temperature_kelvin > 0.0) || !preparation_temperature_kelvin.is_finite() {
        return Err(format!("{what}: the temperature must be positive."));
    }
    let mut mass = empty_mass(network);
    mass[w] = thermal_distribution(&network.wells[w].density_of_states, network.grain_width_cm1, KB_CM * preparation_temperature_kelvin)?;
    GrainDistribution::normalized(mass, 0.0, format!("{what} in well '{}'", network.wells[w].name))
}

/// Gaussian of mean `centre_cm1` and width `sigma_cm1` in well `w` (N eq. gaussian). As a density it is integrated
/// over the grain edges and conditioned on the grid of the well (the fraction outside is reported as lost); as a
/// per-state weight it is multiplied by the density of states (N eq. perstate).
pub fn gaussian(
    network: &ChemicalActivationNetwork,
    w: usize,
    centre_cm1: f64,
    sigma_cm1: f64,
    representation: Representation,
    reference: EnergyReference,
) -> Result<GrainDistribution, String> {
    let scale = match reference {
        EnergyReference::AboveWellGround => "above the well ground",
        EnergyReference::Absolute => "on the common energy scale",
    };
    let what = format!("Gaussian centred at {centre_cm1:.1} cm-1 {scale}, standard deviation {sigma_cm1:.1} cm-1");
    check_well(network, w, &what)?;
    if !(sigma_cm1 > 0.0) || !sigma_cm1.is_finite() || !centre_cm1.is_finite() {
        return Err(format!("{what}: the centre must be finite and the width positive."));
    }
    let de = network.grain_width_cm1;
    let start = grid_start(network, w, reference);
    let well = &network.wells[w];
    let bins: Vec<f64> = (0..well.grain_count())
        .map(|i| normal_interval_probability(start + i as f64 * de, start + (i + 1) as f64 * de, centre_cm1, sigma_cm1))
        .collect();
    let mut mass = empty_mass(network);
    let (description, lost) = match representation {
        Representation::Density | Representation::BinMass => {
            let lost = (1.0 - bins.iter().sum::<f64>()).max(0.0);
            mass[w] = bins;
            (format!("{what}, density, in well '{}'", well.name), lost)
        }
        Representation::PerStateWeight => {
            mass[w] = bins.iter().zip(&well.density_of_states).map(|(b, r)| b * r).collect();
            (format!("{what}, per-state weight, in well '{}'", well.name), 0.0)
        }
    };
    GrainDistribution::normalized(mass, lost, description)
}

/// A tabulated input distribution: bins [lo, hi] (cm-1, on the scale of `reference`) with their values.
#[derive(Debug, Clone, PartialEq)]
pub struct Tabulated {
    pub bins: Vec<(f64, f64)>,
    pub values: Vec<f64>,
    pub representation: Representation,
    pub reference: EnergyReference,
}

/// A tabulated distribution in well `w`, rebinned conservatively onto its grains (N eq. rebin).
pub fn tabulated(network: &ChemicalActivationNetwork, w: usize, table: &Tabulated, support: Support) -> Result<GrainDistribution, String> {
    let what = "Tabulated distribution";
    check_well(network, w, what)?;
    if table.bins.is_empty() || table.bins.len() != table.values.len() {
        return Err(format!("{what}: {} bins and {} values; one value per bin is needed.", table.bins.len(), table.values.len()));
    }
    if table.bins.iter().any(|&(lo, hi)| !lo.is_finite() || !hi.is_finite() || hi < lo) {
        return Err(format!("{what}: every bin needs finite edges with lower <= upper."));
    }
    if table.values.iter().any(|v| !(*v >= 0.0) || !v.is_finite()) {
        return Err(format!("{what}: the values must be finite and >= 0."));
    }
    let masses: Vec<f64> = match table.representation {
        Representation::BinMass => table.values.clone(),
        Representation::Density | Representation::PerStateWeight => {
            if table.bins.iter().any(|&(lo, hi)| hi == lo) {
                return Err(format!("{what}: a density or per-state weight needs bins of nonzero width."));
            }
            table.bins.iter().zip(&table.values).map(|(&(lo, hi), v)| v * (hi - lo)).collect()
        }
    };
    let total: f64 = masses.iter().sum();
    if !(total > 0.0) {
        return Err(format!("{what}: the values are zero everywhere."));
    }
    let well = &network.wells[w];
    let (on_grains, outside) = rebin(&table.bins, &masses, grid_start(network, w, table.reference), network.grain_width_cm1, well.grain_count());
    let lost = outside / total;
    if lost > MAX_LOST_FRACTION && support == Support::Complete {
        return Err(format!(
            "{what}: {:.3e} of the input lies outside the grid of well '{}' ({} grains of {} cm-1 from its bottom); \
             extend the grid, or allow truncation to condition the distribution on the grid.",
            lost,
            well.name,
            well.grain_count(),
            network.grain_width_cm1
        ));
    }
    let mut mass = empty_mass(network);
    let description = match table.representation {
        Representation::PerStateWeight => {
            mass[w] = on_grains.iter().zip(&well.density_of_states).map(|(h, r)| h * r).collect();
            format!("{what} (per-state weight) in well '{}'", well.name)
        }
        _ => {
            mass[w] = on_grains;
            format!("{what} in well '{}'", well.name)
        }
    };
    // The lost probability of a per-state weight would need rho beyond the grid; it is reported as 0.
    let lost = if table.representation == Representation::PerStateWeight { 0.0 } else { lost };
    GrainDistribution::normalized(mass, lost, description)
}

/// Mixture sum_m alpha_m F_m (N eq. mixture), alpha_m >= 0, sum alpha_m = 1; with components confined to single
/// wells this sets the global well fractions (N eq. wells).
pub fn mixture(components: &[(f64, GrainDistribution)]) -> Result<GrainDistribution, String> {
    let what = "Mixture";
    let first = components.first().ok_or_else(|| format!("{what}: no components."))?;
    if components.iter().any(|(a, _)| !(*a >= 0.0) || !a.is_finite()) {
        return Err(format!("{what}: the weights must be finite and >= 0."));
    }
    let sum: f64 = components.iter().map(|(a, _)| a).sum();
    if (sum - 1.0).abs() > 1e-12 {
        return Err(format!("{what}: the weights sum to {sum}, not 1."));
    }
    let shape: Vec<usize> = first.1.mass.iter().map(|m| m.len()).collect();
    if components.iter().any(|(_, d)| d.mass.iter().map(|m| m.len()).collect::<Vec<_>>() != shape) {
        return Err(format!("{what}: the components are on different grids."));
    }
    let mut mass: Vec<Vec<f64>> = shape.iter().map(|&n| vec![0.0; n]).collect();
    for (a, d) in components {
        for (target, source) in mass.iter_mut().zip(&d.mass) {
            for (t, s) in target.iter_mut().zip(source) {
                *t += a * s;
            }
        }
    }
    let lost = components.iter().map(|(a, d)| a * d.lost_fraction).sum();
    let description =
        format!("{what} of {}", components.iter().map(|(a, d)| format!("{a} x [{}]", d.description)).collect::<Vec<_>>().join(", "));
    GrainDistribution::normalized(mass, lost, description)
}

/// The distribution shifted by `delta_cm1` in energy (every well on its own grid), rebinned conservatively.
pub fn shifted(network: &ChemicalActivationNetwork, d: &GrainDistribution, delta_cm1: f64, support: Support) -> Result<GrainDistribution, String> {
    let what = format!("Shift by {delta_cm1:.1} cm-1");
    if d.mass.len() != network.wells.len() {
        return Err(format!("{what}: the distribution has {} wells, the network {}.", d.mass.len(), network.wells.len()));
    }
    let de = network.grain_width_cm1;
    let mut mass = empty_mass(network);
    let mut outside = 0.0;
    for (w, m) in d.mass.iter().enumerate() {
        let bins: Vec<(f64, f64)> = (0..m.len()).map(|i| ((i as f64 - 0.5) * de + delta_cm1, (i as f64 + 0.5) * de + delta_cm1)).collect();
        let (on_grains, lost) = rebin(&bins, m, -0.5 * de, de, network.wells[w].grain_count());
        mass[w] = on_grains;
        outside += lost;
    }
    if outside > MAX_LOST_FRACTION && support == Support::Complete {
        return Err(format!("{what}: {outside:.3e} of the distribution leaves the grid."));
    }
    let lost = 1.0 - (1.0 - d.lost_fraction) * (1.0 - outside);
    GrainDistribution::normalized(mass, lost, format!("{} shifted by {delta_cm1:.1} cm-1", d.description))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::masterequation::chemical_activation_operator::tests::two_well_network;
    use crate::masterequation::chemical_activation_sources::thermal_distribution;
    use crate::constants::KB_CM;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-300)
    }

    #[test]
    fn overlap_rebinning_conserves_mass_and_splits_bins_by_their_overlap() {
        // Input bins [0, 10] and [10, 20] with masses 0.3 and 0.7 onto output bins of width 5 from 0.
        let (out, lost) = rebin(&[(0.0, 10.0), (10.0, 20.0)], &[0.3, 0.7], 0.0, 5.0, 4);
        assert_eq!(out, vec![0.15, 0.15, 0.35, 0.35]);
        assert!(lost < 1e-15, "{lost}");
        // Output covering only [0, 15]: the last half of the second bin is lost.
        let (out, lost) = rebin(&[(0.0, 10.0), (10.0, 20.0)], &[0.3, 0.7], 0.0, 5.0, 3);
        assert_eq!(out, vec![0.15, 0.15, 0.35]);
        assert!(close(lost, 0.35, 1e-15));
        // A zero-width input bin (a stick) goes to the output bin that contains it.
        let (out, lost) = rebin(&[(7.0, 7.0)], &[1.0], 0.0, 5.0, 4);
        assert_eq!((out, lost), (vec![0.0, 1.0, 0.0, 0.0], 0.0));
    }

    #[test]
    fn a_thermal_preparation_uses_its_own_temperature_on_its_well_only() {
        let network = two_well_network();
        let t_prep = 450.0;
        let d = thermal(&network, 1, t_prep).unwrap();
        let reference = thermal_distribution(&network.wells[1].density_of_states, network.grain_width_cm1, KB_CM * t_prep).unwrap();
        assert!(d.mass[1].iter().zip(&reference).all(|(a, b)| close(*a, *b, 1e-14)));
        assert!(d.mass[0].iter().all(|x| *x == 0.0));
        let f = d.well_fractions();
        assert!(f[0] == 0.0 && close(f[1], 1.0, 1e-14), "{f:?}");
        assert_eq!(d.lost_fraction, 0.0);
    }

    #[test]
    fn a_gaussian_density_is_integrated_over_the_grain_edges_and_conditioned_on_the_grid() {
        // N eq. gaussian: f(E) = C exp(-(E - Ec)^2/(2 s^2)) on the allowed domain, integrated over the grains.
        let network = two_well_network();
        let de = network.grain_width_cm1;
        let (centre, sigma) = (1000.0, 50.0);
        let d = gaussian(&network, 0, centre, sigma, Representation::Density, EnergyReference::AboveWellGround).unwrap();
        let cdf = |e: f64| 0.5 * (1.0 + libm_erf((e - centre) / (sigma * 2f64.sqrt())));
        let i = (centre / de).round() as usize;
        let expected = cdf((i as f64 + 0.5) * de) - cdf((i as f64 - 0.5) * de);
        assert!(close(d.mass[0][i], expected, 1e-12), "{} vs {expected}", d.mass[0][i]);
        assert!(close(d.mean_energy_above_bottom_cm1(&network)[0], centre, 1e-6));
        assert!(d.lost_fraction < 1e-15);
        // Centred at the bottom of the well: half of it lies below the grid and is reported as lost.
        let clipped = gaussian(&network, 0, 0.0, sigma, Representation::Density, EnergyReference::AboveWellGround).unwrap();
        // The lost part is the probability below the bottom edge of grain 0, Phi(-dE/2) (the upper tail is negligible).
        assert!(close(clipped.lost_fraction, cdf_at(-0.5 * de, 0.0, sigma), 1e-12), "{}", clipped.lost_fraction);
        assert!(close(clipped.mass[0].iter().sum::<f64>(), 1.0, 1e-14));
    }

    fn cdf_at(e: f64, centre: f64, sigma: f64) -> f64 {
        0.5 * (1.0 + libm_erf((e - centre) / (sigma * 2f64.sqrt())))
    }

    #[test]
    fn a_per_state_gaussian_weight_is_multiplied_by_the_density_of_states() {
        // N eq. perstate: f = rho h / normalization, so F_i/F_j = (rho_i int h)/(rho_j int h).
        let network = two_well_network();
        let rho = &network.wells[0].density_of_states;
        let density = gaussian(&network, 0, 2000.0, 300.0, Representation::Density, EnergyReference::AboveWellGround).unwrap();
        let weight = gaussian(&network, 0, 2000.0, 300.0, Representation::PerStateWeight, EnergyReference::AboveWellGround).unwrap();
        let (i, j) = (150, 250);
        let ratio = (weight.mass[0][i] / weight.mass[0][j]) / (density.mass[0][i] / density.mass[0][j]);
        assert!(close(ratio, rho[i] / rho[j], 1e-10), "{ratio} vs {}", rho[i] / rho[j]);
        // rho rises with E: the per-state weight moves the mean up.
        assert!(weight.mean_energy_above_bottom_cm1(&network)[0] > density.mean_energy_above_bottom_cm1(&network)[0] + 100.0);
    }

    #[test]
    fn a_tabulated_bin_mass_is_rebinned_on_the_scale_of_its_reference() {
        let network = two_well_network();
        let de = network.grain_width_cm1;
        // Bin [995, 1005] above the bottom of B lies exactly on grain 100 of B.
        let table = Tabulated {
            bins: vec![(995.0, 1005.0)],
            values: vec![2.0],
            representation: Representation::BinMass,
            reference: EnergyReference::AboveWellGround,
        };
        let d = tabulated(&network, 1, &table, Support::Complete).unwrap();
        assert_eq!(d.mass[1][100], 1.0);
        // The same bin on the absolute scale: B's grain 0 lies 60 grains below A's, so absolute 1000 cm-1 is
        // grain 100 + 60 of B.
        let absolute = Tabulated { reference: EnergyReference::Absolute, ..table.clone() };
        let d = tabulated(&network, 1, &absolute, Support::Complete).unwrap();
        assert_eq!(d.mass[1][(1000.0 / de) as usize + 60], 1.0);
    }

    #[test]
    fn tabulated_densities_and_per_state_weights_follow_their_definitions() {
        let network = two_well_network();
        let de = network.grain_width_cm1;
        let rho = &network.wells[0].density_of_states;
        // Two bins of one grain width each: densities 1 and 3 per cm-1 give masses 1 : 3.
        let bins = vec![(495.0, 505.0), (505.0, 515.0)];
        let density = Tabulated { bins: bins.clone(), values: vec![1.0, 3.0], representation: Representation::Density, reference: EnergyReference::AboveWellGround };
        let d = tabulated(&network, 0, &density, Support::Complete).unwrap();
        assert!(close(d.mass[0][51] / d.mass[0][50], 3.0, 1e-14));
        assert!(close(d.mass[0][50], 0.25, 1e-14));
        // Per-state weights 1 and 3 are multiplied by rho of the grains.
        let weight = Tabulated { values: vec![1.0, 3.0], representation: Representation::PerStateWeight, ..density.clone() };
        let d = tabulated(&network, 0, &weight, Support::Complete).unwrap();
        assert!(close(d.mass[0][51] / d.mass[0][50], 3.0 * rho[51] / rho[50], 1e-12));
        let _ = de;
    }

    #[test]
    fn mass_beyond_the_grid_is_an_error_unless_truncation_is_allowed() {
        let network = two_well_network();
        let top = network.wells[0].grain_count() as f64 * network.grain_width_cm1;
        let table = Tabulated {
            bins: vec![(1000.0, 1010.0), (top + 100.0, top + 110.0)],
            values: vec![0.9, 0.1],
            representation: Representation::BinMass,
            reference: EnergyReference::AboveWellGround,
        };
        let e = tabulated(&network, 0, &table, Support::Complete).unwrap_err();
        assert!(e.contains("outside"), "{e}");
        let d = tabulated(&network, 0, &table, Support::Truncate).unwrap();
        assert!(close(d.lost_fraction, 0.1, 1e-14));
        assert!(close(d.mass[0].iter().sum::<f64>(), 1.0, 1e-14));
    }

    #[test]
    fn invalid_tables_are_refused() {
        let network = two_well_network();
        let base = Tabulated { bins: vec![(0.0, 10.0)], values: vec![1.0], representation: Representation::BinMass, reference: EnergyReference::AboveWellGround };
        for bad in [
            Tabulated { values: vec![-1.0], ..base.clone() },
            Tabulated { values: vec![f64::NAN], ..base.clone() },
            Tabulated { values: vec![0.0], ..base.clone() },
            Tabulated { bins: vec![(10.0, 0.0)], ..base.clone() },
            Tabulated { values: vec![1.0, 2.0], ..base.clone() },
        ] {
            assert!(tabulated(&network, 0, &bad, Support::Complete).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn a_mixture_keeps_the_given_global_well_fractions() {
        // N eqs. wells and mixture: F = sum_m alpha_m F_m with sum alpha = 1, here pi_A = 0.3 and pi_B = 0.7.
        let network = two_well_network();
        let a = thermal(&network, 0, 300.0).unwrap();
        let b = thermal(&network, 1, 800.0).unwrap();
        let m = mixture(&[(0.3, a.clone()), (0.7, b.clone())]).unwrap();
        let f = m.well_fractions();
        assert!(close(f[0], 0.3, 1e-14) && close(f[1], 0.7, 1e-14), "{f:?}");
        assert!(close(m.mass[1][200], 0.7 * b.mass[1][200], 1e-14));
        // A mixture of one component is that component.
        let one = mixture(&[(1.0, a.clone())]).unwrap();
        assert!(one.mass.iter().flatten().zip(a.mass.iter().flatten()).all(|(x, y)| close(*x, *y, 1e-14)));
        // Weights must be >= 0 and sum to 1.
        assert!(mixture(&[(0.3, a.clone()), (0.6, b.clone())]).is_err());
        assert!(mixture(&[(-0.1, a), (1.1, b)]).is_err());
    }

    #[test]
    fn a_shift_moves_mass_up_and_splits_a_fractional_grain() {
        let network = two_well_network();
        let de = network.grain_width_cm1;
        let mut mass = vec![vec![0.0; network.wells[0].grain_count()], vec![0.0; network.wells[1].grain_count()]];
        mass[0][100] = 1.0;
        let d = GrainDistribution { mass, lost_fraction: 0.0, description: "delta".into() };
        let up = shifted(&network, &d, de, Support::Complete).unwrap();
        assert_eq!(up.mass[0][101], 1.0);
        let half = shifted(&network, &d, -0.5 * de, Support::Complete).unwrap();
        assert!(close(half.mass[0][99], 0.5, 1e-14) && close(half.mass[0][100], 0.5, 1e-14));
    }

    #[test]
    fn descriptions_state_the_energy_reference_and_round_the_energies() {
        let network = two_well_network();
        let name = &network.wells[0].name;
        let d = gaussian(&network, 0, 1234.5678, 50.0, Representation::Density, EnergyReference::AboveWellGround).unwrap();
        assert_eq!(d.description, format!("Gaussian centred at 1234.6 cm-1 above the well ground, standard deviation 50.0 cm-1, density, in well '{name}'"));
        let offset = network.wells[0].bottom_offset_grains as f64 * network.grain_width_cm1;
        let a = gaussian(&network, 0, offset + 1000.04, 50.0, Representation::PerStateWeight, EnergyReference::Absolute).unwrap();
        assert_eq!(
            a.description,
            format!("Gaussian centred at {:.1} cm-1 on the common energy scale, standard deviation 50.0 cm-1, per-state weight, in well '{name}'", offset + 1000.04)
        );
        let s = shifted(&network, &d, 12.3456, Support::Truncate).unwrap();
        assert_eq!(s.description, format!("{} shifted by 12.3 cm-1", d.description));
    }

    /// erf from the standard library is not stable; the test uses the same implementation as the module.
    fn libm_erf(x: f64) -> f64 {
        erf(x)
    }
}
