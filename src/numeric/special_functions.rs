//! Error function and complementary error function in double precision.
//!
//! - |x| < 1: erf(x) = (2/sqrt(pi)) exp(-x^2) sum_{n>=0} 2^n x^(2n+1)/(2n+1)!!, a series of positive terms
//!   without cancellation (Abramowitz, Stegun, Handbook of Mathematical Functions (1964), eq. 7.1.6).
//! - x >= 1: erfc(x) = exp(-x^2)/sqrt(pi) * 1/(x + (1/2)/(x + 1/(x + (3/2)/(x + 2/(x + ...))))), the continued
//!   fraction of eq. 7.1.14, evaluated backwards from depth 300 (converged to 2e-16 relative for x >= 0.8).
//! erf(-x) = -erf(x), erfc(-x) = 2 - erfc(x).

const FRAC_2_SQRT_PI: f64 = std::f64::consts::FRAC_2_SQRT_PI;
const SERIES_LIMIT: f64 = 1.0;

fn erf_series(x: f64) -> f64 {
    let x2 = x * x;
    let (mut term, mut sum, mut n) = (x, x, 0.0);
    loop {
        n += 1.0;
        term *= 2.0 * x2 / (2.0 * n + 1.0);
        sum += term;
        if term <= f64::EPSILON * sum {
            break;
        }
    }
    FRAC_2_SQRT_PI * (-x2).exp() * sum
}

fn erfc_continued_fraction(x: f64) -> f64 {
    let mut t = x;
    for k in (1..=300).rev() {
        t = x + 0.5 * k as f64 / t;
    }
    0.5 * FRAC_2_SQRT_PI * (-x * x).exp() / t
}

/// erf(x).
pub fn erf(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let a = x.abs();
    let value = if a < SERIES_LIMIT { erf_series(a) } else { 1.0 - erfc_continued_fraction(a) };
    value.copysign(x)
}

/// erfc(x) = 1 - erf(x), accurate in the upper tail.
pub fn erfc(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x >= SERIES_LIMIT {
        erfc_continued_fraction(x)
    } else if x > -SERIES_LIMIT {
        1.0 - erf(x)
    } else {
        2.0 - erfc_continued_fraction(-x)
    }
}

/// Probability that a normal variable with mean `centre` and standard deviation `sigma` lies in [a, b],
/// accurate in both tails (differences of erfc on the far side of the mean).
pub fn normal_interval_probability(a: f64, b: f64, centre: f64, sigma: f64) -> f64 {
    let s = sigma * std::f64::consts::SQRT_2;
    let (z1, z2) = ((a - centre) / s, (b - centre) / s);
    if z1 >= 0.0 {
        0.5 * (erfc(z1) - erfc(z2))
    } else if z2 <= 0.0 {
        0.5 * (erfc(-z2) - erfc(-z1))
    } else {
        0.5 * (erf(z2) - erf(z1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erf_and_erfc_match_reference_values() {
        // Abramowitz, Stegun, Table 7.1 and standard tabulations.
        for (x, e) in [(0.5, 0.5204998778130465), (1.0, 0.8427007929497149), (2.0, 0.9953222650189527), (2.5, 0.9995930479825550)] {
            assert!((erf(x) - e).abs() < 5e-16, "erf({x}) = {} vs {e}", erf(x));
            assert!((erf(-x) + e).abs() < 5e-16);
        }
        for (x, e) in [(3.0, 2.209049699858544e-5), (5.0, 1.537459794428035e-12), (2.0, 4.677734981047266e-3)] {
            assert!((erfc(x) / e - 1.0).abs() < 1e-14, "erfc({x}) = {:e} vs {e:e}", erfc(x));
        }
        assert_eq!(erf(0.0), 0.0);
        assert!((erfc(-3.0) - (2.0 - 2.209049699858544e-5)).abs() < 1e-15);
    }

    #[test]
    fn normal_interval_probabilities_add_up_and_are_accurate_in_the_tails() {
        let (c, s) = (10.0, 2.0);
        let whole = normal_interval_probability(-1e3, 1e3, c, s);
        assert!((whole - 1.0).abs() < 1e-15);
        let parts: f64 = (0..40).map(|k| normal_interval_probability(k as f64 * 0.5, (k + 1) as f64 * 0.5, c, s)).sum();
        assert!((parts - normal_interval_probability(0.0, 20.0, c, s)).abs() < 1e-15);
        // Far upper tail: erfc(5)/2 - erfc(6)/2 with s*sqrt(2) = 1 at centre 0.
        let tail = normal_interval_probability(5.0, 6.0, 0.0, std::f64::consts::FRAC_1_SQRT_2);
        assert!((tail / (0.5 * (1.537459794428035e-12 - 2.151973671249892e-17)) - 1.0).abs() < 1e-13, "{tail:e}");
    }
}
