use crate::constants::PI;

pub fn wigner(beta: f64, omega: f64) -> f64 {
    let u = omega * beta;

    return 1.0 + u * u / 24.0;
}

pub fn bell(beta: f64, omega: f64) -> f64 {
    let u = omega * beta;

    return (0.5 * u / f64::sin(0.5 * u)).abs();
}

pub fn skodje_truhlar(beta: f64, omega: f64, v0: f64) -> f64 {
    let alpha = 2.0 * PI / omega;

    if beta > alpha {
        // When omega is large and low T
        let kappa = (((beta - alpha) * v0).exp() - 1.0) * beta / (beta - alpha);

        //println!("Skodje, beta*omega: {}", beta * omega);
        //println!("beta*omega > twopi");

        return kappa;
    } else if beta < alpha {
        // When omega is small and high T
        let dum0 = PI * beta / alpha;
        let dum1 = dum0 / dum0.sin();
        let dum2 = ((beta - alpha) * v0).exp() * beta / (beta - alpha);
        let kappa = dum1 + dum2;

        //println!("Skodje, beta*omega: {}", beta * omega);
        //println!("beta*omega < twopi");

        return kappa;
    } else {
        panic!("Wrong mode-argument in Skodje_Truhlar routine");
    }
}

pub fn skodje_truhlar_exact(beta: f64, omega: f64, v0: f64) -> f64 {
    let alpha = 2.0 * PI / omega;

    const NMAX: usize = 100;

    let mut res = 0.0;

    for n in 0..=NMAX {
        let numerator = 1.0 - ((beta - (n + 1) as f64 * alpha) * v0).exp();
        let denom = (n + 1) as f64 * alpha - beta;
        let dum = 1.0 / (n as f64 * alpha + beta);
        let mut alter = 1.0;

        if n % 2 == 1 {
            alter = -1.0;
        }

        res += alter * beta * (numerator / denom + dum);
    }

    return res;
}

/// Canonical Eckart tunneling correction (energies in one consistent unit, beta = 1/kT):
///   kappa(T) = beta exp(beta V_f) integral_0^emax P(E - V_f) exp(-beta E) dE,
/// with E measured from the asymptote on the forward side and P the Eckart transmission probability
/// (`eckart_transmission_probability`, Miller 1979 eq. 8); this is Miller's thermally averaged factor
/// Gamma = integral dE1 P(E1) exp(-E1/kT)/kT (Miller, J. Am. Chem. Soc. 101, 6810 (1979), p. 6811).
/// Because P vanishes below the higher asymptote, kappa is the same in both directions.
/// `vf`, `vb`: barrier heights from the forward and backward sides; `omega`: magnitude of the imaginary
/// frequency; `de`, `emax`: integration step and upper limit (Simpson rule).
pub fn eckart(beta: f64, omega: f64, vf: f64, vb: f64, de: f64, emax: f64) -> f64 {
    let n_emax = (emax / de) as usize;
    let integrand: Vec<f64> = (0..n_emax)
        .map(|i| {
            let e = i as f64 * de;
            (-e * beta).exp() * eckart_transmission_probability(e - vf, vf, vb, omega)
        })
        .collect();
    simpson_integrate(&integrand, 0, n_emax, de) * (beta * vf).exp() * beta
}

fn simpson_integrate(func: &[f64], nmin: usize, nmax: usize, step: f64) -> f64 {
    let mut s0 = 0.0;
    let mut s1 = 0.0;
    let mut s2 = 0.0;

    if nmin > nmax {
        panic!("Wrong boundaries in Simpson Integration: nmin or nmax");
    }

    let ndata: usize = nmax - nmin;

    for i in (nmin..nmax - 2).step_by(2) {
        s1 += func[i];
        s0 += func[i + 1];
        s2 += func[i + 2];
    }

    //println!("sm1 =  {}, s0 = {}, sp1 = {}", s1, s0, s2);
    let mut res = step * (s1 + 4.0 * s0 + s2) / 3.0;
    //println!("step =  {}, res = {}", step, res);

    // If n is even, add the last slice separately
    if ndata % 2 == 0 {
        res += step * (5.0 * func[nmax - 1] + 8.0 * func[nmax - 2] - func[nmax - 3]) / 12.0;
    }

    return res;
}


/// One-dimensional transmission probability P(E1) through an Eckart barrier (energies in cm-1).
///
/// Miller, J. Am. Chem. Soc. 101, 6810 (1979), eq. 8:
///   P(E1) = sinh(a) sinh(b) / [sinh^2((a + b)/2) + cosh^2(c)],
///   a = (4 pi/hw) sqrt(E1 + V0) (V0^-1/2 + V1^-1/2)^-1,
///   b = (4 pi/hw) sqrt(E1 + V1) (V0^-1/2 + V1^-1/2)^-1,
///   c = 2 pi sqrt(V0 V1/hw^2 - 1/16),
/// with E1 the energy in the reaction coordinate relative to the barrier top, V0 and V1 the barrier
/// heights relative to the two sides and hw the magnitude of the imaginary frequency. When
/// V0 V1/hw^2 < 1/16, c is imaginary and cosh(c) becomes cos(|c|) (Johnston, Heicklen, J. Phys. Chem.
/// 66, 532 (1962), text after eq. 13). P is symmetric in V0 and V1 and vanishes for
/// E1 <= -min(V0, V1), below the higher of the two asymptotes.
pub fn eckart_transmission_probability(e1: f64, v0: f64, v1: f64, imaginary_frequency_cm1: f64) -> f64 {
    if !(e1 + v0.min(v1) > 0.0) {
        return 0.0;
    }
    let hw = imaginary_frequency_cm1;
    let s = 1.0 / (v0.powf(-0.5) + v1.powf(-0.5));
    let a = 4.0 * PI / hw * (e1 + v0).sqrt() * s;
    let b = 4.0 * PI / hw * (e1 + v1).sqrt() * s;
    let u = a + b;
    // Overflow-free evaluation with all terms scaled by exp(-M):
    //   sinh(a) sinh(b) = exp(a + b) (1 - exp(-2a)) (1 - exp(-2b)) / 4,
    //   sinh^2(u/2)     = exp(u) (1 - exp(-u))^2 / 4,
    //   cosh^2(c)       = exp(2c) (1 + exp(-2c))^2 / 4   (c real).
    let factor_a = -(-2.0 * a).exp_m1();
    let factor_b = -(-2.0 * b).exp_m1();
    let c_squared = v0 * v1 / (hw * hw) - 1.0 / 16.0;
    let (numerator, denominator) = if c_squared >= 0.0 {
        let c = 2.0 * PI * c_squared.sqrt();
        let m = u.max(2.0 * c);
        (
            (u - m).exp() * factor_a * factor_b / 4.0,
            (u - m).exp() * (-u).exp_m1().powi(2) / 4.0 + (2.0 * c - m).exp() * (1.0 + (-2.0 * c).exp()).powi(2) / 4.0,
        )
    } else {
        // Imaginary c: cosh^2(c) = cos^2(|c|) (Johnston, Heicklen 1962).
        let c = 2.0 * PI * (-c_squared).sqrt();
        (
            factor_a * factor_b / 4.0,
            (-u).exp_m1().powi(2) / 4.0 + c.cos().powi(2) * (-u).exp(),
        )
    };
    numerator / denominator
}

/// Tunneling-corrected number of states of a transition state (Miller 1979, eqs. 6 and 9):
///   N_QM(e) = integral dE1 P'(E1) N(e - E1),
/// e the energy above the zero-point level of the transition state, N(e) its number of states without
/// the reaction coordinate on grains e = i dE (`ts_sum_of_states[i]`), P the Eckart probability.
/// Discretization: the reaction-coordinate energy is divided into cells of width dE centred at
/// E1 = j dE, each carrying the probability increment P((j + 1/2) dE) - P((j - 1/2) dE), so that
///   N_QM(k dE) = sum_j [P((j + 1/2) dE) - P((j - 1/2) dE)] N((k - j) dE);
/// for a step P = h(E1) (no tunneling) this gives N_QM = N exactly.
/// Returns (m, w) with m = floor(min(V0, V1)/dE + 1/2) and w[k + m] = N_QM(k dE) for
/// -m <= k <= L - 1 - m, L = ts_sum_of_states.len(): the states extend m grains below the top of the
/// barrier (down to the higher of the two asymptotes), and N_QM(k dE) needs N up to (k + m) dE, so N must
/// be given m grains beyond the highest energy at which N_QM is wanted.
pub fn eckart_tunneling_sum_of_states(
    ts_sum_of_states: &[f64],
    grain_width_cm1: f64,
    v0: f64,
    v1: f64,
    imaginary_frequency_cm1: f64,
) -> (usize, Vec<f64>) {
    let len = ts_sum_of_states.len();
    let m = (v0.min(v1) / grain_width_cm1 + 0.5).floor() as usize;
    let probability = |x: f64| eckart_transmission_probability(x * grain_width_cm1, v0, v1, imaginary_frequency_cm1);
    // Probability increments of the cells j = -m ..= len-1-m (stored at j + m); cells below -m carry none.
    let increments: Vec<f64> = (0..len)
        .map(|index| {
            let j = index as f64 - m as f64;
            probability(j + 0.5) - probability(j - 0.5)
        })
        .collect();
    // N_QM(k dE), k = index - m: sum over the cells j = -m ..= k of increment(j) N((k - j) dE), where
    // k - j = index - (j + m) <= index <= len - 1.
    let w = (0..len)
        .map(|index| (0..=index).map(|jj| increments[jj] * ts_sum_of_states[index - jj]).sum::<f64>())
        .collect();
    (m, w)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::KB_CM;

    /// Direct evaluation of Miller 1979 eq. 8 (no overflow protection), real c.
    fn eckart_direct(e1: f64, v0: f64, v1: f64, hw: f64) -> f64 {
        let s = 1.0 / (v0.powf(-0.5) + v1.powf(-0.5));
        let a = 4.0 * PI / hw * (e1 + v0).sqrt() * s;
        let b = 4.0 * PI / hw * (e1 + v1).sqrt() * s;
        let c = 2.0 * PI * (v0 * v1 / (hw * hw) - 1.0 / 16.0).sqrt();
        a.sinh() * b.sinh() / ((0.5 * (a + b)).sinh().powi(2) + c.cosh().powi(2))
    }

    #[test]
    fn eckart_probability_matches_the_hyperbolic_form_of_miller_eq_8() {
        for &(v0, v1, hw) in &[(3000.0_f64, 5000.0_f64, 1500.0_f64), (6800.0, 7450.0, 2658.84), (2000.0, 2000.0, 800.0)] {
            for k in -20..40 {
                let e1 = k as f64 * 0.04 * v0.min(v1);
                if e1 <= -v0.min(v1) {
                    continue;
                }
                let p = eckart_transmission_probability(e1, v0, v1, hw);
                let q = eckart_direct(e1, v0, v1, hw);
                assert!((p - q).abs() <= 1e-12 * q.max(1e-300) + 1e-300, "E1 = {e1}: {p:e} vs {q:e}");
            }
        }
    }

    #[test]
    fn eckart_probability_is_symmetric_and_bounded() {
        let (v0, v1, hw) = (4000.0, 9000.0, 1800.0);
        assert_eq!(eckart_transmission_probability(-4000.0, v0, v1, hw), 0.0);
        assert_eq!(eckart_transmission_probability(-5000.0, v0, v1, hw), 0.0);
        assert!((eckart_transmission_probability(20.0 * v1, v0, v1, hw) - 1.0).abs() < 1e-9);
        for k in -39..50 {
            let e1 = k as f64 * 100.0;
            let p = eckart_transmission_probability(e1, v0, v1, hw);
            assert!((0.0..=1.0).contains(&p));
            assert!((p - eckart_transmission_probability(e1, v1, v0, hw)).abs() < 1e-14);
        }
    }

    #[test]
    fn eckart_probability_becomes_the_parabolic_barrier_near_the_top_of_a_wide_barrier() {
        // Near the top of a barrier much higher than hw, P -> 1/(1 + exp(-2 pi E1/hw)) (Miller 1979).
        let hw = 1000.0;
        let v = 2.0e6;
        for k in -5..=5 {
            let e1 = k as f64 * 100.0;
            let parabolic = 1.0 / (1.0 + (-2.0 * PI * e1 / hw).exp());
            let p = eckart_transmission_probability(e1, v, v, hw);
            assert!((p - parabolic).abs() < 2e-3, "E1 = {e1}: {p} vs {parabolic}");
        }
    }

    #[test]
    fn eckart_probability_stays_finite_for_narrow_high_and_for_wide_low_barriers() {
        // Narrow and high: a, b, c of several hundred (sinh, cosh overflow without scaling).
        for k in -29..60 {
            let p = eckart_transmission_probability(k as f64 * 1000.0, 30000.0, 40000.0, 100.0);
            assert!(p.is_finite() && (0.0..=1.0).contains(&p), "{p}");
        }
        // Wide and low: V0 V1/hw^2 < 1/16, c imaginary, cosh^2(c) -> cos^2(|c|).
        let (v0, v1, hw) = (100.0_f64, 120.0_f64, 1000.0_f64);
        let s = 1.0 / (v0.powf(-0.5) + v1.powf(-0.5));
        for k in -9..30 {
            let e1 = k as f64 * 10.0;
            let a = 4.0 * PI / hw * (e1 + v0).sqrt() * s;
            let b = 4.0 * PI / hw * (e1 + v1).sqrt() * s;
            let c = 2.0 * PI * (1.0 / 16.0 - v0 * v1 / (hw * hw)).sqrt();
            let expected = a.sinh() * b.sinh() / ((0.5 * (a + b)).sinh().powi(2) + c.cos().powi(2));
            let p = eckart_transmission_probability(e1, v0, v1, hw);
            assert!((p - expected).abs() < 1e-12, "E1 = {e1}: {p} vs {expected}");
        }
    }

    #[test]
    fn canonical_eckart_correction_is_the_same_in_both_directions() {
        // kappa = beta exp(beta V_from) integral_0^inf P(E - V_from) exp(-beta E) dE
        //       = beta integral P(E1) exp(-beta E1) dE1 over E1 > -min(V_f, V_b): independent of the
        // direction, because P vanishes below the higher asymptote (energies in hartree here).
        let hartree_per_cm1 = 4.556_335e-6;
        let omega = 1500.0 * hartree_per_cm1;
        let (vf, vb) = (5000.0 * hartree_per_cm1, 3000.0 * hartree_per_cm1);
        let beta = 1.0 / (3.166_811_563e-6 * 300.0);
        let de = 1.0 * hartree_per_cm1;
        let emax = 30000.0 * hartree_per_cm1;
        let forward = eckart(beta, omega, vf, vb, de, emax);
        let reverse = eckart(beta, omega, vb, vf, de, emax);
        assert!(forward > 1.0);
        assert!(((forward - reverse) / reverse).abs() < 1e-6, "forward {forward}, reverse {reverse}");
    }

    #[test]
    fn tunneling_sum_of_states_reduces_to_the_classical_count_for_a_step_probability() {
        // A very wide barrier (hw -> 0) makes P a step at E1 = 0: N_QM = N on the grains above the top.
        let d_e = 10.0;
        let n: Vec<f64> = (0..800).map(|i| (1.0 + 0.1 * i as f64).powi(3)).collect();
        let (m, w) = eckart_tunneling_sum_of_states(&n, d_e, 5000.0, 6000.0, 1.0e-3);
        assert_eq!(m, 500);
        assert_eq!(w.len(), 800);
        for k in 0..m {
            assert!(w[k] < 1e-12, "below the top: {}", w[k]);
        }
        for i in 0..300 {
            assert!((w[m + i] - n[i]).abs() < 1e-9 * n[i], "grain {i}: {} vs {}", w[m + i], n[i]);
        }
    }

    #[test]
    fn tunneling_sum_of_states_is_the_stieltjes_convolution_of_miller_eq_9() {
        let d_e = 20.0;
        let (v0, v1, hw) = (3000.0, 4000.0, 1500.0);
        let n: Vec<f64> = (0..300).map(|i| 1.0 + 0.5 * i as f64).collect();
        let (m, w) = eckart_tunneling_sum_of_states(&n, d_e, v0, v1, hw);
        assert_eq!(m, 150);
        let p = |x: f64| eckart_transmission_probability(x * d_e, v0, v1, hw);
        for k in [-150isize, -100, -1, 0, 1, 50, 149] {
            // sum_{j <= k} [P((j + 1/2) dE) - P((j - 1/2) dE)] N(k - j)
            let mut expected = 0.0;
            for j in -(m as isize)..=k {
                expected += (p(j as f64 + 0.5) - p(j as f64 - 0.5)) * n[(k - j) as usize];
            }
            let got = w[(k + m as isize) as usize];
            assert!((got - expected).abs() <= 1e-12 * expected.abs().max(1e-300), "k = {k}: {got} vs {expected}");
        }
        // Tunneling below the top: N_QM > 0 there; above the top N_QM < N + P-weighted tail.
        assert!(w[m - 10] > 0.0);
    }

    #[test]
    fn skodje_truhlar_exact_matches_truncated_parabolic_integral() {
        // The series is the term-by-term integral of
        //   kappa = beta * exp(beta*V0) * int_0^inf exp(-beta*E) P(E) dE,
        //   P(E)  = 1 / (1 + exp(alpha*(V0 - E))),  alpha = 2*pi/omega,
        // so it must agree with a direct quadrature of that integral
        // (up to the truncation of the alternating series at NMAX = 100, ~0.2% here).
        let omega = 1000.0; // cm-1, imaginary barrier frequency
        let temp = 300.0;
        let v0 = 3000.0; // cm-1
        let beta = 1.0 / (KB_CM * temp);
        let alpha = 2.0 * PI / omega;

        let n = 200_000;
        let emax = v0 + 60.0 / beta + 60.0 / alpha;
        let h = emax / n as f64;
        let mut sum = 0.0;
        for i in 0..=n {
            let e = i as f64 * h;
            let weight = if i == 0 || i == n { 0.5 } else { 1.0 };
            let p = 1.0 / (1.0 + (alpha * (v0 - e)).exp());
            sum += weight * (beta * (v0 - e)).exp() * p;
        }
        let kappa_ref = beta * sum * h;

        let kappa = skodje_truhlar_exact(beta, omega, v0);

        assert!(
            ((kappa - kappa_ref) / kappa_ref).abs() < 5.0e-3,
            "kappa = {kappa}, direct integral = {kappa_ref}"
        );
    }
}

