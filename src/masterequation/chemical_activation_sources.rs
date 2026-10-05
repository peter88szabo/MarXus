//! Nascent (source) distributions F of the chemical-activation master equation.
//!
//! All distributions are probability masses per grain on the grid of the receiving well (grain 0 =
//! well bottom, the energy zero of the master equation, PO14 p. 234), normalized to sum_i F_i = 1.
//! Continuous densities of the papers, integral f(E) dE = 1, correspond to F_i = f(E_i) dE.
//!
//! - Thermal reactants through the entrance channel (PO14 eq. 7; O02 eq. 11):
//!     f(E) = W(E - E0) exp[-(E - E0)/kT] / integral_0^inf W(e) exp(-e/kT) de,   E >= E0,
//!   with W the sum of states of the transition state of the reverse (dissociation) reaction and E0
//!   its threshold. Since k(E) = W(E - E0)/(h rho(E)) (PO14 eq. 9), this equals
//!     f(E) ∝ rho(E) k(E) exp(-E/kT),
//!   i.e. the source follows from the microcanonical rate coefficient of the reverse reaction alone.
//! - Non-thermal reactants with normalized distributions n_A, n_B (PO14 eq. 8):
//!     f(E) = integral_0^{E-E0} n_A(e) n_B(E - E0 - e) de,   E >= E0,
//!   assuming "that the reactive cross section of reaction (a) is independent of the internal energies
//!   of the reactants A and B, which is equivalent to the assumption that the rate of the reverse
//!   reaction (-a) can be described by (vibrational) phase space theory neglecting angular momentum
//!   restrictions" (PO14 p. 235). For consecutive activation, A is the chemically activated
//!   intermediate of a previous master equation and B the thermal partner (PO14 eqs. 10-11, with
//!   E0 = -RE, RE the 0 K reaction energy).
//! - Shift approximation (PO14 eqs. 12-13): neglecting the width of the partner distribution,
//!   n_B(e) = delta(e - <E_B>), gives
//!     f(E) = n_A(E + RE - <E_B>),
//!   "the distribution of the chemically activated reactant shifted to higher energies by the sum of
//!   the reaction energy at 0 K and the average thermal energy" of the partner (PO14 p. 237).
//!
//! References: see `chemical_activation_network.rs`.

use super::chemical_activation_network::ChemicalActivationNetwork;

/// Normalize a non-negative distribution to unit sum.
fn normalized(mut dist: Vec<f64>, what: &str) -> Result<Vec<f64>, String> {
    if dist.iter().any(|x| !(*x >= 0.0) || !x.is_finite()) {
        return Err(format!("{what}: the distribution must be finite and >= 0."));
    }
    let total: f64 = dist.iter().sum();
    if !(total > 0.0) {
        return Err(format!("{what}: the distribution is zero everywhere."));
    }
    dist.iter_mut().for_each(|x| *x /= total);
    Ok(dist)
}

/// Fraction of a normalized distribution that may be lost beyond the top of the receiving grid.
pub const MAX_TRUNCATED_FRACTION: f64 = 1.0e-10;

/// Thermal (Boltzmann) distribution rho(E_i) exp(-E_i/kT) on a grid, normalized.
pub fn thermal_distribution(rho: &[f64], grain_width_cm1: f64, kt_cm1: f64) -> Result<Vec<f64>, String> {
    if rho.iter().any(|r| !(*r >= 0.0) || !r.is_finite()) {
        return Err("Thermal distribution: the density of states must be finite and >= 0.".into());
    }
    let log_weight: Vec<f64> =
        rho.iter().enumerate().map(|(i, r)| r.ln() - i as f64 * grain_width_cm1 / kt_cm1).collect();
    from_log_weights(&log_weight, "Thermal distribution")
}

/// exp(log_weight) normalized, evaluated relative to the largest weight (no under- or overflow).
fn from_log_weights(log_weight: &[f64], what: &str) -> Result<Vec<f64>, String> {
    let max = log_weight.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return Err(format!("{what}: the distribution is zero everywhere."));
    }
    normalized(log_weight.iter().map(|l| (l - max).exp()).collect(), what)
}

/// Mean energy sum_i E_i n_i of a normalized distribution, cm-1.
pub fn mean_energy_cm1(dist: &[f64], grain_width_cm1: f64) -> f64 {
    let total: f64 = dist.iter().sum();
    dist.iter().enumerate().map(|(i, x)| i as f64 * grain_width_cm1 * x).sum::<f64>() / total
}

/// Thermal source through an entrance channel from the microcanonical rate coefficient of its reverse
/// (dissociation) reaction: F_i ∝ rho(E_i) k(E_i) exp(-E_i/kT) (PO14 eqs. 7 and 9).
pub fn thermal_source_from_rate(
    rho: &[f64],
    reverse_rate_s_inv: &[f64],
    grain_width_cm1: f64,
    kt_cm1: f64,
) -> Result<Vec<f64>, String> {
    if rho.len() != reverse_rate_s_inv.len() {
        return Err(format!(
            "Thermal source: {} densities of states for {} rate coefficients.",
            rho.len(),
            reverse_rate_s_inv.len()
        ));
    }
    if reverse_rate_s_inv.iter().any(|k| !(*k >= 0.0) || !k.is_finite()) {
        return Err("Thermal source: the rate coefficients must be finite and >= 0.".into());
    }
    if rho.iter().any(|r| !(*r >= 0.0) || !r.is_finite()) {
        return Err("Thermal source: the density of states must be finite and >= 0.".into());
    }
    // ln(rho k exp(-E/kT)); ln 0 = -inf gives an exact zero after exponentiation.
    let log_weight: Vec<f64> = rho
        .iter()
        .zip(reverse_rate_s_inv)
        .enumerate()
        .map(|(i, (r, k))| r.ln() + k.ln() - i as f64 * grain_width_cm1 / kt_cm1)
        .collect();
    from_log_weights(&log_weight, "Thermal source")
}

/// Thermal source through one or several entrance channels (well, channel) of a network, from the
/// reverse rate coefficients: F_{w,i} ∝ rho_w(E_i) k(E_i) exp(-E_i/kT) with E_i on the ABSOLUTE energy
/// scale of the network, so that channels into different wells are weighted by their thermal fluxes
/// sum W‡ exp(-E/kT) (PO14 eqs. 7 and 9). Normalized over all wells.
pub fn thermal_entrance_source(
    network: &ChemicalActivationNetwork,
    channels: &[(usize, usize)],
    kt_cm1: f64,
) -> Result<Vec<Vec<f64>>, String> {
    if channels.is_empty() {
        return Err("Thermal entrance source: no entrance channel given.".into());
    }
    // ln(rho k exp(-E/kT)) summed over the entrance channels of each well; -inf where no channel is open.
    let mut log_weight: Vec<Vec<f64>> =
        network.wells.iter().map(|w| vec![f64::NEG_INFINITY; w.grain_count()]).collect();
    for &(w, c) in channels {
        let well = network
            .wells
            .get(w)
            .ok_or_else(|| format!("Thermal entrance source: well {w} does not exist."))?;
        let channel = well
            .channels
            .get(c)
            .ok_or_else(|| format!("Thermal entrance source: well '{}' has no channel {c}.", well.name))?;
        for i in 0..well.grain_count() {
            let k = channel.rate_constant_s_inv[i];
            if k > 0.0 {
                let l = well.density_of_states[i].ln() + k.ln() - network.absolute_energy_cm1(w, i) / kt_cm1;
                // ln(exp(a) + exp(b)) for two channels into the same well.
                let a = log_weight[w][i];
                log_weight[w][i] = if a == f64::NEG_INFINITY { l } else { a.max(l) + (-(a - l).abs()).exp().ln_1p() };
            }
        }
    }
    let max = log_weight.iter().flatten().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return Err("Thermal entrance source: the entrance channels are closed on the whole grid.".into());
    }
    let mut f: Vec<Vec<f64>> = log_weight.iter().map(|lw| lw.iter().map(|l| (l - max).exp()).collect()).collect();
    let total: f64 = f.iter().flatten().sum();
    f.iter_mut().flatten().for_each(|x| *x /= total);
    Ok(f)
}

/// Thermal source from the sum of states W(e) of the entrance transition state, e = 0, dE, 2dE, ...
/// above the threshold grain `threshold_grain` of the well (PO14 eq. 7; O02 eq. 11).
pub fn thermal_source_from_sum_of_states(
    grains: usize,
    threshold_grain: usize,
    sum_of_states: &[f64],
    grain_width_cm1: f64,
    kt_cm1: f64,
) -> Result<Vec<f64>, String> {
    if threshold_grain >= grains {
        return Err(format!("Thermal source: threshold grain {threshold_grain} beyond the grid ({grains} grains)."));
    }
    if sum_of_states.len() < grains - threshold_grain {
        return Err(format!(
            "Thermal source: the sum of states covers {} grains above the threshold, the grid needs {}.",
            sum_of_states.len(),
            grains - threshold_grain
        ));
    }
    if sum_of_states.iter().any(|w| !(*w >= 0.0) || !w.is_finite()) {
        return Err("Thermal source: the sum of states must be finite and >= 0.".into());
    }
    let mut log_weight = vec![f64::NEG_INFINITY; grains];
    for i in threshold_grain..grains {
        let e = (i - threshold_grain) as f64 * grain_width_cm1;
        log_weight[i] = sum_of_states[i - threshold_grain].ln() - e / kt_cm1;
    }
    from_log_weights(&log_weight, "Thermal source")
}

/// All population formed in one grain.
pub fn single_grain_source(grains: usize, grain: usize) -> Result<Vec<f64>, String> {
    if grain >= grains {
        return Err(format!("Single-grain source: grain {grain} beyond the grid ({grains} grains)."));
    }
    let mut f = vec![0.0; grains];
    f[grain] = 1.0;
    Ok(f)
}

/// Place the masses `mass[m]` at grains `first + m` of a grid of `grains` grains (`first` may be
/// negative); masses that fall outside the grid may not exceed MAX_TRUNCATED_FRACTION of the total.
fn place_on_grid(grains: usize, first: isize, mass: &[f64], what: &str) -> Result<Vec<f64>, String> {
    let total: f64 = mass.iter().sum();
    let mut f = vec![0.0; grains];
    let mut lost = 0.0;
    for (m, &x) in mass.iter().enumerate() {
        let i = first + m as isize;
        if i >= 0 && (i as usize) < grains {
            f[i as usize] = x;
        } else {
            lost += x;
        }
    }
    if lost > MAX_TRUNCATED_FRACTION * total {
        return Err(format!(
            "{what}: a fraction {:e} of the distribution lies outside the energy grid of the receiving well \
             ({grains} grains); extend the grid.",
            lost / total
        ));
    }
    normalized(f, what)
}

/// Convolution source of PO14 eq. 8 (eqs. 10-11 for consecutive activation): F_i ∝
/// sum_k a_k b_{i - t - k} for i >= t, with t = `threshold_grain` = E0/dE of the reverse reaction
/// (E0 = -RE). `a` and `b` are normalized internal-energy distributions of the two reactants on grids
/// of the same grain width.
pub fn convolution_source(grains: usize, threshold_grain: usize, a: &[f64], b: &[f64]) -> Result<Vec<f64>, String> {
    let a = normalized(a.to_vec(), "Convolution source, first reactant")?;
    let b = normalized(b.to_vec(), "Convolution source, second reactant")?;
    // Discrete form of PO14 eq. 8 for probability masses: c_m = sum_k a_k b_{m-k}, placed at
    // E = E0 + m dE.
    let mut c = vec![0.0; a.len() + b.len() - 1];
    for (k, &x) in a.iter().enumerate() {
        if x == 0.0 {
            continue;
        }
        for (l, &y) in b.iter().enumerate() {
            c[k + l] += x * y;
        }
    }
    place_on_grid(grains, threshold_grain as isize, &c, "Convolution source")
}

/// Shift approximation of PO14 eqs. 12-13: F(E) = n_A(E + RE - <E_B>), with the 0 K reaction energy
/// `reaction_energy_cm1` (RE < 0 for an exothermic step) and the mean thermal energy of the partner.
/// The shift (-RE + <E_B>) is rounded to whole grains.
pub fn shifted_source(
    grains: usize,
    activated_reactant: &[f64],
    reaction_energy_cm1: f64,
    partner_mean_energy_cm1: f64,
    grain_width_cm1: f64,
) -> Result<Vec<f64>, String> {
    let n_a = normalized(activated_reactant.to_vec(), "Shifted source")?;
    // f(E) = n_A(E + RE - <E_B>): grain j of A appears at grain j + shift, shift = (-RE + <E_B>)/dE.
    let shift = ((-reaction_energy_cm1 + partner_mean_energy_cm1) / grain_width_cm1).round() as isize;
    place_on_grid(grains, shift, &n_a, "Shifted source")
}

#[cfg(test)]
mod tests {
    use super::*;

    const D_E: f64 = 10.0;
    const KT: f64 = 0.695_034_76 * 298.0;

    fn assert_close(a: &[f64], b: &[f64], tol: f64) {
        assert_eq!(a.len(), b.len());
        for (i, (x, y)) in a.iter().zip(b).enumerate() {
            assert!((x - y).abs() <= tol * x.abs().max(y.abs()).max(1e-300), "element {i}: {x:e} vs {y:e}");
        }
    }

    #[test]
    fn source_from_the_reverse_rate_equals_the_transition_state_sum_of_states_form() {
        // k(E) = W(E - E0)/(h rho(E)) (PO14 eq. 9) makes rho k exp(-E/kT) ∝ W(E - E0) exp(-(E - E0)/kT).
        let n = 500;
        let t = 130;
        let rho: Vec<f64> = (0..n).map(|i| (1.0 + 0.03 * i as f64).powi(10)).collect();
        let w: Vec<f64> = (0..n - t).map(|e| (1.0 + 0.04 * e as f64).powi(7)).collect();
        let k: Vec<f64> = (0..n).map(|i| if i >= t { w[i - t] / (3.3356e-11 * rho[i]) } else { 0.0 }).collect();
        let from_rate = thermal_source_from_rate(&rho, &k, D_E, KT).unwrap();
        let from_w = thermal_source_from_sum_of_states(n, t, &w, D_E, KT).unwrap();
        assert_close(&from_rate, &from_w, 1e-12);
        assert!(from_w[..t].iter().all(|&x| x == 0.0));
        assert!((from_w.iter().sum::<f64>() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn entrance_channels_into_several_wells_are_weighted_by_their_thermal_fluxes() {
        use crate::masterequation::chemical_activation_operator::tests::two_well_network;
        // Use the products channels of A (index 0) and B (index 0) as two entrance channels.
        let network = two_well_network();
        let f = thermal_entrance_source(&network, &[(0, 0), (1, 0)], KT).unwrap();
        let total: f64 = f.iter().flatten().sum();
        assert!((total - 1.0).abs() < 1e-14);
        let flux = |w: usize| -> f64 {
            let well = &network.wells[w];
            (0..well.grain_count())
                .map(|i| well.density_of_states[i] * well.channels[0].rate_constant_s_inv[i] * (-network.absolute_energy_cm1(w, i) / KT).exp())
                .sum()
        };
        let ratio = f[0].iter().sum::<f64>() / f[1].iter().sum::<f64>();
        assert!((ratio / (flux(0) / flux(1)) - 1.0).abs() < 1e-12, "ratio {ratio}");
        // A single channel reproduces the one-well form.
        let single = thermal_entrance_source(&network, &[(0, 0)], KT).unwrap();
        let a = &network.wells[0];
        let expected = thermal_source_from_rate(&a.density_of_states, &a.channels[0].rate_constant_s_inv, D_E, KT).unwrap();
        assert_close(&single[0], &expected, 1e-12);
        assert!(single[1].iter().all(|&x| x == 0.0));
        assert!(thermal_entrance_source(&network, &[(0, 9)], KT).is_err());
    }

    #[test]
    fn convolution_of_two_single_grains_is_a_single_grain_above_the_threshold() {
        let mut a = vec![0.0; 50];
        a[7] = 1.0;
        let mut b = vec![0.0; 80];
        b[11] = 1.0;
        let f = convolution_source(300, 120, &a, &b).unwrap();
        assert_eq!(f, single_grain_source(300, 120 + 7 + 11).unwrap());
    }

    #[test]
    fn convolution_with_a_sharp_partner_reduces_to_the_shift_approximation() {
        // PO14 eqs. 10 -> 12: n_B(e) = delta(e - <E_B>).
        let n1: Vec<f64> = normalized((0..400).map(|i| (-(i as f64 - 150.0).powi(2) / 800.0).exp()).collect(), "n1").unwrap();
        let partner_grain = 21;
        let mut partner = vec![0.0; 60];
        partner[partner_grain] = 1.0;
        let reaction_energy_cm1 = -6630.0; // E0 = 663 grains
        let convolved = convolution_source(1200, 663, &partner, &n1).unwrap();
        let shifted = shifted_source(1200, &n1, reaction_energy_cm1, partner_grain as f64 * D_E, D_E).unwrap();
        assert_close(&convolved, &shifted, 1e-14);
        // The mean energy rises by -RE + <E_B>.
        let shift = mean_energy_cm1(&shifted, D_E) - mean_energy_cm1(&n1, D_E);
        assert!((shift - (6630.0 + 210.0)).abs() < 1e-9, "shift {shift}");
    }

    #[test]
    fn thermal_distribution_is_normalized_boltzmann() {
        let rho: Vec<f64> = (0..300).map(|i| 1.0 + i as f64).collect();
        let n = thermal_distribution(&rho, D_E, KT).unwrap();
        assert!((n.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        let ratio = n[100] / n[50];
        let expected = rho[100] / rho[50] * (-(50.0 * D_E) / KT).exp();
        assert!((ratio / expected - 1.0).abs() < 1e-12);
    }

    #[test]
    fn sources_that_do_not_fit_on_the_receiving_grid_are_rejected() {
        let n1 = vec![0.5, 0.5];
        assert!(shifted_source(100, &n1, -990.0, 0.0, D_E).is_err()); // shifted to grains 99 and 100
        assert!(shifted_source(101, &n1, -990.0, 0.0, D_E).is_ok());
        assert!(single_grain_source(10, 10).is_err());
        assert!(convolution_source(100, 99, &[0.0, 1.0], &[1.0]).is_err());
    }
}
