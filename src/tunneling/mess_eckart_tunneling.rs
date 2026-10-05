//! Eckart tunneling as implemented in MESS (`mess_eckart_tunneling`), as an alternative to the exact Eckart
//! transmission of `tunneling::eckart_transmission_probability` (Miller, J. Am. Chem. Soc. 101, 6810 (1979),
//! eq. 8), so that MESS decks can be reproduced.
//!
//! Mirrors MESS (Georgievskii, Jasper, Zador, Miller, Burke, Goldsmith, Klippenstein; source
//! src/libmess/model.cc, classes `Model::Tunnel` and `Model::EckartTunnel`; Apache License 2.0), the
//! model that MESS uses for `Tunneling Eckart` (ImaginaryFrequency, two WellDepth values):
//!   - semiclassical transmission probability P(E) = 1 / (1 + exp(-S(E))), E measured from the barrier top,
//!     with P = 1 for S > S_max and P = 0 for S < -S_max (S_max = 100, MESS `Tunnel::_action_max`);
//!   - action of the Eckart barrier, with d_w = V_w/omega (V_0 <= V_1 the two well depths, omega the
//!     magnitude of the imaginary frequency):
//!       S(E) = 4 pi / (d_0^(-1/2) + d_1^(-1/2)) * sum_w [ sqrt(max(E/omega + d_w, 0)) - sqrt(d_w) ],
//!     which tends to 2 pi E/omega (the parabolic barrier) near the top;
//!   - energy range: from -E_c, with E_c = V_0 (the smaller depth) unless -S(-V_0) > S_max, in which case
//!     E_c is lowered by bisection (to 1 cm-1) until -S(-E_c) <= S_max (`Tunnel::_adjust_cutoff`);
//!   - number of states: N_QM(E_j) = N(E_j + E_c) P(-E_c) + sum_i N(E_j + E_c - e_i) [P(e_i) - P(e_(i-1))]
//!     + N(E_j + E_c - e_len) [1 - P(e_(len-1))], e_i = -E_c + i dE up to 2 omega above the top
//!     (`Tunnel::convolute`);
//!   - canonical factor relative to the barrier top (the "tunneling partition function correction factor"
//!     of the MESS log): kappa(T) = (dE/kT) sum P(e) exp(-e/kT), e = -E_c, -E_c + dE, ... < 10 kT,
//!     dE = 0.01 kT (`Tunnel::weight`, without its factor exp(-E_c/kT), which refers the weight to -E_c).
//! Reproduces the factors printed in a MESS log to 1e-4, except for deep tunneling (kappa >> 100) at
//! T <= 300 K, where MESS's ground-state bookkeeping of the barrier lowers the cutoff further (by about
//! 60 cm-1 for the barriers of validation/ZZAllyl+O2_Gamma_Case2; not reproduced): 0.06-0.4% at 300 K.
//! For deep tunneling (kappa >> 1) it is smaller than the exact Eckart factor: by 17-23% for the
//! H-transfer barriers of validation/ZZAllyl+O2_Gamma_Case2, 2-6% for the others.

use crate::constants::PI;

/// Limit of the semiclassical action (MESS `Tunnel::_action_max`, default 100).
pub const MESS_ACTION_MAX: f64 = 100.0;
/// Energy step of the canonical factor in units of kT (MESS `Tunnel::_wtol`, default 0.01).
pub const MESS_WEIGHT_STEP_OVER_TEMPERATURE: f64 = 0.01;
/// Upper limit of the canonical factor in units of kT (MESS `Tunnel::weight`).
pub const MESS_WEIGHT_UPPER_LIMIT_OVER_TEMPERATURE: f64 = 10.0;
/// Range of the convolution above the barrier top in units of omega (MESS `upper_cutoff_factor`).
pub const MESS_CONVOLUTION_RANGE_OVER_FREQUENCY: f64 = 2.0;

/// The MESS Eckart tunneling model of one barrier (energies in cm-1).
#[derive(Debug, Clone)]
pub struct MessEckartTunneling {
    omega_cm1: f64,
    /// Well depths over omega, ascending.
    depths_over_omega: [f64; 2],
    /// 4 pi / (d_0^(-1/2) + d_1^(-1/2)).
    action_factor: f64,
    cutoff_cm1: f64,
}

impl MessEckartTunneling {
    /// Model for the imaginary frequency (magnitude) and the two well depths (cm-1).
    pub fn new(imaginary_frequency_cm1: f64, well_depths_cm1: [f64; 2]) -> Result<Self, String> {
        if !(imaginary_frequency_cm1 > 0.0) || well_depths_cm1.iter().any(|d| !(*d > 0.0)) {
            return Err(format!(
                "MESS Eckart tunneling: the imaginary frequency ({imaginary_frequency_cm1}) and the well depths \
                 ({well_depths_cm1:?}) must be positive."
            ));
        }
        let mut depths = well_depths_cm1;
        if depths[0] > depths[1] {
            depths.swap(0, 1);
        }
        let depths_over_omega = [depths[0] / imaginary_frequency_cm1, depths[1] / imaginary_frequency_cm1];
        let action_factor = 4.0 * PI / depths_over_omega.iter().map(|d| 1.0 / d.sqrt()).sum::<f64>();
        let mut model = Self { omega_cm1: imaginary_frequency_cm1, depths_over_omega, action_factor, cutoff_cm1: depths[0] };
        // Tunnel::_adjust_cutoff: bisection to 1 cm-1 for the largest E_c with -S(-E_c) < S_max.
        if -model.action(-model.cutoff_cm1) > MESS_ACTION_MAX {
            let (mut e_min, mut e_max) = (0.0, model.cutoff_cm1);
            while e_max - e_min > 1.0 {
                let e_test = 0.5 * (e_min + e_max);
                if -model.action(-e_test) < MESS_ACTION_MAX {
                    e_min = e_test;
                } else {
                    e_max = e_test;
                }
            }
            model.cutoff_cm1 = e_min;
        }
        Ok(model)
    }

    /// Lowest energy below the barrier top at which tunneling is counted, E_c (cm-1).
    pub fn cutoff_cm1(&self) -> f64 {
        self.cutoff_cm1
    }

    /// Semiclassical action S(E), E (cm-1) from the barrier top.
    pub fn action(&self, energy_cm1: f64) -> f64 {
        let e = energy_cm1 / self.omega_cm1;
        self.action_factor
            * self.depths_over_omega.iter().map(|d| (e + d).max(0.0).sqrt() - d.sqrt()).sum::<f64>()
    }

    /// Transmission probability P(E) = 1/(1 + exp(-S(E))), clamped at |S| > S_max.
    pub fn transmission(&self, energy_cm1: f64) -> f64 {
        let s = self.action(energy_cm1);
        if s > MESS_ACTION_MAX {
            1.0
        } else if s < -MESS_ACTION_MAX {
            0.0
        } else {
            1.0 / (1.0 + (-s).exp())
        }
    }

    /// Canonical tunneling factor relative to the barrier top at k_B T = `kt_cm1` (cm-1).
    pub fn canonical_factor(&self, kt_cm1: f64) -> f64 {
        let step = kt_cm1 * MESS_WEIGHT_STEP_OVER_TEMPERATURE;
        let e_max = MESS_WEIGHT_UPPER_LIMIT_OVER_TEMPERATURE * kt_cm1;
        // ln P(e) - e/kT on the grid of Tunnel::weight; summed relative to its maximum.
        let mut terms = Vec::new();
        let mut e = -self.cutoff_cm1;
        while e < e_max {
            let p = self.transmission(e);
            if p > 0.0 {
                terms.push(p.ln() - e / kt_cm1);
            }
            e += step;
        }
        let max = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let sum: f64 = terms.iter().map(|t| (t - max).exp()).sum();
        MESS_WEIGHT_STEP_OVER_TEMPERATURE * sum * max.exp()
    }
}

/// Tunneling-corrected number of states of a transition state with the MESS model (`Tunnel::convolute`).
/// `ts_sum_of_states[i]` = N(i dE) above the barrier top. Returns (m, w) with m = round(E_c/dE) and
/// w[k + m] = N_QM(k dE) for -m <= k <= L - 1 - m, the interface of
/// `tunneling::eckart_tunneling_sum_of_states` (N must be given m cells beyond the highest energy at which
/// N_QM is wanted). The grid of the transmission probabilities is aligned with the cells (E_c rounded to a
/// whole number of cells).
pub fn mess_eckart_tunneling_sum_of_states(
    ts_sum_of_states: &[f64],
    grain_width_cm1: f64,
    well_depths_cm1: [f64; 2],
    imaginary_frequency_cm1: f64,
) -> Result<(usize, Vec<f64>), String> {
    let model = MessEckartTunneling::new(imaginary_frequency_cm1, well_depths_cm1)?;
    let m = (model.cutoff_cm1() / grain_width_cm1).round() as usize;
    let len_td = ((m as f64 * grain_width_cm1 + MESS_CONVOLUTION_RANGE_OVER_FREQUENCY * imaginary_frequency_cm1)
        / grain_width_cm1)
        .ceil() as usize;
    let td: Vec<f64> =
        (0..len_td).map(|i| model.transmission((i as f64 - m as f64) * grain_width_cm1)).collect();
    let stat = ts_sum_of_states;
    let w = (0..stat.len())
        .map(|j| {
            let mut value = stat[j] * td[0];
            for i in 1..td.len().min(j + 1) {
                value += stat[j - i] * (td[i] - td[i - 1]);
            }
            if j >= td.len() {
                value += stat[j - td.len()] * (1.0 - td[td.len() - 1]);
            }
            value
        })
        .collect();
    Ok((m, w))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{CM1_TO_KCAL, KB_CM};
    use crate::constants::PI;

    fn kcal(x: f64) -> f64 {
        x / CM1_TO_KCAL
    }

    #[test]
    fn near_the_top_the_transmission_is_that_of_a_parabolic_barrier() {
        // S(E) -> 2 pi E / omega for |E| << depths: P = 1/(1 + exp(-2 pi E/omega)).
        let model = MessEckartTunneling::new(1500.0, [kcal(15.0), kcal(20.0)]).unwrap();
        assert!((model.transmission(0.0) - 0.5).abs() < 1e-15);
        for e in [-20.0, -5.0, 3.0, 10.0] {
            let parabolic = 1.0 / (1.0 + (-2.0 * PI * e / 1500.0).exp());
            assert!((model.transmission(e) / parabolic - 1.0).abs() < 2e-3, "{e}: {} vs {parabolic}", model.transmission(e));
        }
    }

    #[test]
    fn the_transmission_rises_monotonically_and_is_clamped_at_the_action_limit() {
        let model = MessEckartTunneling::new(2742.1, [kcal(16.3), kcal(18.7)]).unwrap();
        let mut last = 0.0;
        for i in 0..400 {
            let e = -model.cutoff_cm1() + i as f64 * 25.0;
            let p = model.transmission(e);
            assert!(p >= last && (0.0..=1.0).contains(&p));
            last = p;
        }
        // |S| above the limit: P = 1 above the top, 0 below.
        assert_eq!(model.transmission(1.0e6), 1.0);
        let below = MessEckartTunneling::new(300.0, [kcal(40.0), kcal(45.0)]).unwrap();
        assert_eq!(below.transmission(-below.cutoff_cm1() - 1.0e3), 0.0);
    }

    #[test]
    fn the_cutoff_is_the_smaller_depth_unless_the_action_there_exceeds_the_limit() {
        // Wide barrier (low frequency, deep): -S(-V0) > action_max, the cutoff is lowered (bisection to
        // 1 cm-1) so that -S(-cutoff) <= action_max.
        let deep = MessEckartTunneling::new(300.0, [kcal(40.0), kcal(45.0)]).unwrap();
        assert!(deep.cutoff_cm1() < kcal(40.0));
        assert!(-deep.action(-deep.cutoff_cm1()) <= MESS_ACTION_MAX);
        assert!(-deep.action(-deep.cutoff_cm1() - 2.0) > MESS_ACTION_MAX);
        // Narrow barrier: the cutoff is the smaller well depth.
        let narrow = MessEckartTunneling::new(2742.1, [kcal(18.7), kcal(16.3)]).unwrap();
        assert!((narrow.cutoff_cm1() - kcal(16.3)).abs() < 1e-9);
    }

    #[test]
    fn canonical_factors_reproduce_the_reference_mess_log_of_the_zz_allyl_o2_case_2_deck() {
        // "tunneling partition function correction factors" of the MESS log (2025-09-29 run) for the
        // barriers of validation/ZZAllyl+O2_Gamma_Case2: (imaginary frequency, well depths in kcal/mol,
        // [(T, factor)], relative tolerance).
        // Reproduced to 1e-4 except for the deep H-transfer barriers (kappa >> 100) at T <= 300 K: there the
        // MESS values correspond to a cutoff about 60 cm-1 below the smaller well depth (from MESS's
        // ground-state bookkeeping of the barrier, not reproduced here); the mimic is higher by 0.06-0.4% at
        // 300 K and 0.7-2.5% at 200 K.
        let cases: [(f64, [f64; 2], [(f64, f64, f64); 3]); 6] = [
            (2742.1, [16.3, 18.7], [(200.0, 3.20644e9, 3e-2), (300.0, 23568.9, 5e-3), (1000.0, 1.84159, 1e-4)]),
            (2191.7, [11.7, 14.8], [(200.0, 705872.0, 3e-2), (300.0, 304.858, 5e-3), (1000.0, 1.45202, 1e-4)]),
            (2666.8, [19.1, 19.8], [(200.0, 3.13499e10, 1e-2), (300.0, 45667.9, 1e-3), (1000.0, 1.80414, 1e-4)]),
            (602.6, [18.7, 18.3], [(200.0, 2.47364, 1e-4), (300.0, 1.43623, 1e-4), (1000.0, 1.03011, 1e-4)]),
            (586.3, [21.0, 7.0], [(200.0, 2.23716, 1e-4), (300.0, 1.38985, 1e-4), (1000.0, 1.02669, 1e-4)]),
            (988.9, [15.8, 23.8], [(200.0, 33.5499, 1e-4), (300.0, 2.96203, 1e-4), (1000.0, 1.08335, 1e-4)]),
        ];
        for (omega, depths, reference) in cases {
            let model = MessEckartTunneling::new(omega, [kcal(depths[0]), kcal(depths[1])]).unwrap();
            for (t, factor, tolerance) in reference {
                let kappa = model.canonical_factor(KB_CM * t);
                assert!((kappa / factor - 1.0).abs() < tolerance, "omega {omega}, T {t}: {kappa:e} vs {factor:e}");
            }
        }
    }

    #[test]
    fn the_convolution_of_a_constant_number_of_states_is_the_transmission_probability() {
        // Stieltjes convolution N_QM(E) = sum N(E - e_i) [P(e_i) - P(e_{i-1})] with N = 1 everywhere gives
        // P(E) on the grid below 2 omega above the top, and 1 beyond.
        let (omega, depths) = (1500.0, [kcal(10.0), kcal(12.0)]);
        let cell = 1.0;
        let len = 12000;
        let (below, w) = mess_eckart_tunneling_sum_of_states(&vec![1.0; len], cell, depths, omega).unwrap();
        let model = MessEckartTunneling::new(omega, depths).unwrap();
        assert_eq!(below, (model.cutoff_cm1() / cell).round() as usize);
        for index in [0usize, 100, below, below + 50, below + 2000, len - 1] {
            let e = (index as f64 - below as f64) * cell;
            let expected = if e < 2.0 * omega { model.transmission(e) } else { 1.0 };
            assert!((w[index] - expected).abs() < 1e-12, "index {index}: {} vs {expected}", w[index]);
        }
    }
}

