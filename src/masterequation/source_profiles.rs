//! Time profiles of source amplitudes R(t) (population per second) for prepared sources (design note N,
//! papers/Reactant_flux_initiation/MarXus_Nonthermal_Sources.tex, Table "Useful time profiles", eqs.
//! effectivesource and injected).
//!
//! Every profile states its amplitude as an amount (integrated) or as a rate, never ambiguously. An impulse is
//! not a rate: it is an exact population jump n(t_p+) = n(t_p-) + N_p F, applied by the integrator as an event
//! (N, Table "Useful time profiles"; Sec. 10.3: "Implement an ideal delta pulse as a population jump"). The
//! continuous part `rate(t)` excludes impulses. `breakpoints` are the times where the rate jumps (rectangular
//! edges, the start and end of a feed); the integrator stops there. The injected amount
//! N_in(0, t) = int_0^t R dt + the impulses at t_p <= t (N eq. injected) is given in closed form.

use crate::numeric::special_functions::normal_interval_probability;

/// Time profile of a source amplitude.
#[derive(Debug, Clone, PartialEq)]
pub enum TimeProfile {
    /// N_p delta(t - t_p): an exact population jump at t_p.
    Impulse { time_s: f64, amount: f64 },
    /// R = N_p/(end - start) on [start, end).
    Rectangular { start_s: f64, end_s: f64, amount: f64 },
    /// R = N_p exp(-(t - t_p)^2/(2 sigma^2))/(sqrt(2 pi) sigma) for t >= 0, renormalized for the part before t = 0
    /// (`truncated_fraction`).
    Gaussian { centre_s: f64, sigma_s: f64, amount: f64 },
    /// R = R_0 on [start, end); open-ended without `end_s`.
    Feed { start_s: f64, end_s: Option<f64>, rate: f64 },
    /// R = k_f N_prec exp(-k_tot (t - start)) for t >= start, a first-order precursor with competing total loss
    /// k_tot >= k_f.
    PrecursorDecay { start_s: f64, formation_rate_s_inv: f64, total_loss_rate_s_inv: f64, precursor_amount: f64 },
    /// Measured (t, R), piecewise linear between the points and zero outside them.
    Tabulated { times_s: Vec<f64>, rates: Vec<f64> },
    /// The sum of several profiles (a pulse train, or a combination).
    Train(Vec<TimeProfile>),
}

fn nonnegative(x: f64) -> bool {
    x >= 0.0 && x.is_finite()
}

impl TimeProfile {
    pub fn validate(&self) -> Result<(), String> {
        let ok = match self {
            TimeProfile::Impulse { time_s, amount } => nonnegative(*time_s) && nonnegative(*amount),
            TimeProfile::Rectangular { start_s, end_s, amount } => {
                nonnegative(*start_s) && end_s.is_finite() && end_s > start_s && nonnegative(*amount)
            }
            TimeProfile::Gaussian { centre_s, sigma_s, amount } => {
                centre_s.is_finite() && *sigma_s > 0.0 && sigma_s.is_finite() && nonnegative(*amount)
            }
            TimeProfile::Feed { start_s, end_s, rate } => {
                nonnegative(*start_s) && nonnegative(*rate) && end_s.map_or(true, |e| e.is_finite() && e > *start_s)
            }
            TimeProfile::PrecursorDecay { start_s, formation_rate_s_inv, total_loss_rate_s_inv, precursor_amount } => {
                nonnegative(*start_s)
                    && nonnegative(*formation_rate_s_inv)
                    && *total_loss_rate_s_inv > 0.0
                    && total_loss_rate_s_inv.is_finite()
                    && total_loss_rate_s_inv >= formation_rate_s_inv
                    && nonnegative(*precursor_amount)
            }
            TimeProfile::Tabulated { times_s, rates } => {
                times_s.len() >= 2
                    && times_s.len() == rates.len()
                    && times_s.iter().all(|t| nonnegative(*t))
                    && times_s.windows(2).all(|w| w[1] > w[0])
                    && rates.iter().all(|r| nonnegative(*r))
            }
            TimeProfile::Train(members) => {
                if members.is_empty() {
                    false
                } else {
                    for m in members {
                        m.validate()?;
                    }
                    true
                }
            }
        };
        if ok {
            Ok(())
        } else {
            Err(format!(
                "Source profile {self:?}: times and amounts must be finite and >= 0, an interval must have end > start, \
                 a width must be positive, a precursor needs k_tot >= k_f >= 0, and a table needs >= 2 increasing times."
            ))
        }
    }

    /// Fraction of a Gaussian pulse that lies before t = 0 and is redistributed by its renormalization (0 otherwise).
    pub fn truncated_fraction(&self) -> f64 {
        match self {
            TimeProfile::Gaussian { centre_s, sigma_s, .. } => normal_interval_probability(f64::NEG_INFINITY, 0.0, *centre_s, *sigma_s),
            _ => 0.0,
        }
    }

    /// The continuous rate R(t) (impulses excluded).
    pub fn rate(&self, t: f64) -> f64 {
        match self {
            TimeProfile::Impulse { .. } => 0.0,
            TimeProfile::Rectangular { start_s, end_s, amount } => {
                if t >= *start_s && t < *end_s {
                    amount / (end_s - start_s)
                } else {
                    0.0
                }
            }
            TimeProfile::Gaussian { centre_s, sigma_s, amount } => {
                if t < 0.0 {
                    return 0.0;
                }
                let z = (t - centre_s) / sigma_s;
                let kept = 1.0 - self.truncated_fraction();
                amount * (-0.5 * z * z).exp() / ((2.0 * std::f64::consts::PI).sqrt() * sigma_s) / kept
            }
            TimeProfile::Feed { start_s, end_s, rate } => {
                if t >= *start_s && end_s.map_or(true, |e| t < e) {
                    *rate
                } else {
                    0.0
                }
            }
            TimeProfile::PrecursorDecay { start_s, formation_rate_s_inv, total_loss_rate_s_inv, precursor_amount } => {
                if t >= *start_s {
                    formation_rate_s_inv * precursor_amount * (-total_loss_rate_s_inv * (t - start_s)).exp()
                } else {
                    0.0
                }
            }
            TimeProfile::Tabulated { times_s, rates } => {
                let n = times_s.len();
                if t < times_s[0] || t > times_s[n - 1] {
                    return 0.0;
                }
                let k = times_s.partition_point(|&x| x <= t).clamp(1, n - 1);
                let (t0, t1, r0, r1) = (times_s[k - 1], times_s[k], rates[k - 1], rates[k]);
                r0 + (r1 - r0) * (t - t0) / (t1 - t0)
            }
            TimeProfile::Train(members) => members.iter().map(|m| m.rate(t)).sum(),
        }
    }

    /// The impulses (time, amount), sorted by time.
    pub fn impulses(&self) -> Vec<(f64, f64)> {
        let mut out = match self {
            TimeProfile::Impulse { time_s, amount } => vec![(*time_s, *amount)],
            TimeProfile::Train(members) => members.iter().flat_map(|m| m.impulses()).collect(),
            _ => Vec::new(),
        };
        out.sort_by(|a, b| a.0.total_cmp(&b.0));
        out
    }

    /// Times where the continuous rate jumps (sorted, without duplicates).
    pub fn breakpoints(&self) -> Vec<f64> {
        let mut out = match self {
            TimeProfile::Rectangular { start_s, end_s, .. } => vec![*start_s, *end_s],
            TimeProfile::Feed { start_s, end_s, .. } => std::iter::once(*start_s).chain(*end_s).collect(),
            TimeProfile::PrecursorDecay { start_s, .. } => vec![*start_s],
            TimeProfile::Tabulated { times_s, .. } => vec![times_s[0], times_s[times_s.len() - 1]],
            TimeProfile::Train(members) => members.iter().flat_map(|m| m.breakpoints()).collect(),
            TimeProfile::Impulse { .. } | TimeProfile::Gaussian { .. } => Vec::new(),
        };
        out.sort_by(|a, b| a.total_cmp(b));
        out.dedup();
        out
    }

    /// Whether the continuous rate is constant on the open interval (a, b).
    pub fn is_constant_on(&self, a: f64, b: f64) -> bool {
        match self {
            TimeProfile::Impulse { .. } => true,
            TimeProfile::Rectangular { .. } | TimeProfile::Feed { .. } => {
                !self.breakpoints().iter().any(|&t| t > a && t < b)
            }
            TimeProfile::Gaussian { amount, .. } => *amount == 0.0,
            TimeProfile::PrecursorDecay { start_s, formation_rate_s_inv, precursor_amount, .. } => {
                b <= *start_s || *formation_rate_s_inv == 0.0 || *precursor_amount == 0.0
            }
            TimeProfile::Tabulated { times_s, .. } => b <= times_s[0] || a >= times_s[times_s.len() - 1],
            TimeProfile::Train(members) => members.iter().all(|m| m.is_constant_on(a, b)),
        }
    }

    /// Injected amount from 0 to t: the integral of the rate plus the impulses at t_p <= t (N eq. injected).
    pub fn injected_until(&self, t: f64) -> f64 {
        match self {
            TimeProfile::Impulse { time_s, amount } => {
                if t >= *time_s {
                    *amount
                } else {
                    0.0
                }
            }
            TimeProfile::Rectangular { start_s, end_s, amount } => amount * ((t - start_s) / (end_s - start_s)).clamp(0.0, 1.0),
            TimeProfile::Gaussian { centre_s, sigma_s, amount } => {
                if t <= 0.0 {
                    return 0.0;
                }
                amount * normal_interval_probability(0.0, t, *centre_s, *sigma_s) / (1.0 - self.truncated_fraction())
            }
            TimeProfile::Feed { start_s, end_s, rate } => {
                let upper = end_s.map_or(t, |e| t.min(e));
                rate * (upper - start_s).max(0.0)
            }
            TimeProfile::PrecursorDecay { start_s, formation_rate_s_inv, total_loss_rate_s_inv, precursor_amount } => {
                if t <= *start_s {
                    return 0.0;
                }
                formation_rate_s_inv / total_loss_rate_s_inv
                    * precursor_amount
                    * -(-total_loss_rate_s_inv * (t - start_s)).exp_m1()
            }
            TimeProfile::Tabulated { times_s, rates } => {
                let mut total = 0.0;
                for k in 1..times_s.len() {
                    let (t0, t1) = (times_s[k - 1], times_s[k]);
                    if t <= t0 {
                        break;
                    }
                    let upper = t.min(t1);
                    total += 0.5 * (rates[k - 1] + self.rate(upper)) * (upper - t0);
                }
                total
            }
            TimeProfile::Train(members) => members.iter().map(|m| m.injected_until(t)).sum(),
        }
    }

    /// Total amount over all times; None for an open-ended feed.
    pub fn total_amount(&self) -> Option<f64> {
        match self {
            TimeProfile::Impulse { amount, .. } | TimeProfile::Rectangular { amount, .. } | TimeProfile::Gaussian { amount, .. } => {
                Some(*amount)
            }
            TimeProfile::Feed { start_s, end_s, rate } => end_s.map(|e| rate * (e - start_s)),
            TimeProfile::PrecursorDecay { formation_rate_s_inv, total_loss_rate_s_inv, precursor_amount, .. } => {
                Some(formation_rate_s_inv / total_loss_rate_s_inv * precursor_amount)
            }
            TimeProfile::Tabulated { times_s, .. } => Some(self.injected_until(times_s[times_s.len() - 1])),
            TimeProfile::Train(members) => members.iter().map(|m| m.total_amount()).sum(),
        }
    }
}

/// A number as written in a deck: plain between 1e-3 and 1e5, otherwise in exponent form.
fn number(x: f64) -> String {
    if x == 0.0 || (1e-3..=1e5).contains(&x.abs()) {
        format!("{x}")
    } else {
        format!("{x:e}")
    }
}

/// The profile in words, with units (report and run-settings lines).
impl std::fmt::Display for TimeProfile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TimeProfile::Impulse { time_s, amount } => write!(f, "impulse at {} s, amount {}", number(*time_s), number(*amount)),
            TimeProfile::Rectangular { start_s, end_s, amount } => {
                write!(f, "rectangular pulse {} .. {} s, amount {}", number(*start_s), number(*end_s), number(*amount))
            }
            TimeProfile::Gaussian { centre_s, sigma_s, amount } => write!(
                f,
                "Gaussian pulse centred at {} s, standard deviation {} s, amount {}",
                number(*centre_s),
                number(*sigma_s),
                number(*amount)
            ),
            TimeProfile::Feed { start_s, end_s, rate } => {
                write!(f, "constant feed {} 1/s from {} s", number(*rate), number(*start_s))?;
                match end_s {
                    Some(e) => write!(f, " to {} s", number(*e)),
                    None => write!(f, ", open-ended"),
                }
            }
            TimeProfile::PrecursorDecay { start_s, formation_rate_s_inv, total_loss_rate_s_inv, precursor_amount } => write!(
                f,
                "precursor decay from {} s: k_f = {} 1/s, k_tot = {} 1/s, precursor amount {}",
                number(*start_s),
                number(*formation_rate_s_inv),
                number(*total_loss_rate_s_inv),
                number(*precursor_amount)
            ),
            TimeProfile::Tabulated { times_s, .. } => write!(
                f,
                "tabulated rate, {} points {} .. {} s",
                times_s.len(),
                number(times_s[0]),
                number(times_s[times_s.len() - 1])
            ),
            TimeProfile::Train(members) => {
                write!(f, "train of {}: ", members.len())?;
                for (k, m) in members.iter().enumerate() {
                    write!(f, "{}{m}", if k > 0 { "; " } else { "" })?;
                }
                Ok(())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-300)
    }

    /// Trapezoid integral of the continuous rate on [a, b] with n intervals.
    fn integral(p: &TimeProfile, a: f64, b: f64, n: usize) -> f64 {
        let h = (b - a) / n as f64;
        (0..=n).map(|k| p.rate(a + k as f64 * h) * if k == 0 || k == n { 0.5 } else { 1.0 }).sum::<f64>() * h
    }

    #[test]
    fn an_impulse_is_a_jump_and_not_a_rate() {
        let p = TimeProfile::Impulse { time_s: 1e-6, amount: 2.0 };
        p.validate().unwrap();
        assert_eq!(p.rate(1e-6), 0.0);
        assert_eq!(p.impulses(), vec![(1e-6, 2.0)]);
        assert_eq!(p.injected_until(0.5e-6), 0.0);
        assert_eq!(p.injected_until(1e-6), 2.0);
        assert_eq!(p.total_amount(), Some(2.0));
        assert!(p.breakpoints().is_empty());
    }

    #[test]
    fn a_rectangular_pulse_has_its_amount_spread_over_its_width() {
        let p = TimeProfile::Rectangular { start_s: 1.0, end_s: 3.0, amount: 4.0 };
        p.validate().unwrap();
        assert_eq!((p.rate(0.5), p.rate(1.5), p.rate(3.5)), (0.0, 2.0, 0.0));
        assert_eq!(p.injected_until(2.0), 2.0);
        assert_eq!(p.injected_until(10.0), 4.0);
        assert_eq!(p.breakpoints(), vec![1.0, 3.0]);
        assert!(p.is_constant_on(1.0, 3.0) && p.is_constant_on(3.0, 9.0) && !p.is_constant_on(0.0, 2.0));
    }

    #[test]
    fn a_gaussian_pulse_is_renormalized_for_its_part_before_the_start() {
        // Centred at 1 sigma: 15.87% of the full Gaussian lies before t = 0 and is redistributed.
        let p = TimeProfile::Gaussian { centre_s: 2e-8, sigma_s: 2e-8, amount: 1.0 };
        p.validate().unwrap();
        assert!(close(p.injected_until(1.0), 1.0, 1e-14));
        assert!(close(integral(&p, 0.0, 2e-7, 20000), 1.0, 1e-7));
        assert!(close(p.injected_until(2e-8), integral(&p, 0.0, 2e-8, 20000), 1e-8));
        assert!(!p.is_constant_on(0.0, 1e-8));
        assert!(close(p.truncated_fraction(), 0.15865525393145707, 1e-12));
    }

    #[test]
    fn a_feed_has_a_constant_rate_on_its_interval() {
        let p = TimeProfile::Feed { start_s: 1.0, end_s: Some(4.0), rate: 0.5 };
        p.validate().unwrap();
        assert_eq!((p.rate(0.5), p.rate(2.0), p.rate(4.0)), (0.0, 0.5, 0.0));
        assert_eq!(p.injected_until(3.0), 1.0);
        assert_eq!(p.total_amount(), Some(1.5));
        let open = TimeProfile::Feed { start_s: 0.0, end_s: None, rate: 2.0 };
        assert_eq!(open.total_amount(), None);
        assert_eq!(open.injected_until(3.0), 6.0);
        assert_eq!(open.breakpoints(), vec![0.0]);
    }

    #[test]
    fn a_precursor_decay_forms_its_share_of_the_precursor() {
        // R = k_f N_prec exp(-k_tot (t - t0)): in total (k_f/k_tot) N_prec.
        let p = TimeProfile::PrecursorDecay { start_s: 1e-6, formation_rate_s_inv: 2e5, total_loss_rate_s_inv: 5e5, precursor_amount: 3.0 };
        p.validate().unwrap();
        assert_eq!(p.rate(0.5e-6), 0.0);
        assert!(close(p.rate(1e-6), 6e5, 1e-14));
        assert!(close(p.total_amount().unwrap(), 1.2, 1e-14));
        let t = 3e-6;
        assert!(close(p.injected_until(t), 1.2 * (1.0 - (-5e5f64 * 2e-6).exp()), 1e-13));
        assert!(close(integral(&p, 1e-6, t, 200000), p.injected_until(t), 1e-8));
        let bad = TimeProfile::PrecursorDecay { start_s: 0.0, formation_rate_s_inv: 5.0, total_loss_rate_s_inv: 2.0, precursor_amount: 1.0 };
        assert!(bad.validate().is_err());
    }

    #[test]
    fn a_tabulated_rate_is_piecewise_linear_and_zero_outside() {
        let p = TimeProfile::Tabulated { times_s: vec![1.0, 2.0, 4.0], rates: vec![0.0, 2.0, 0.0] };
        p.validate().unwrap();
        assert_eq!((p.rate(0.5), p.rate(1.5), p.rate(3.0), p.rate(5.0)), (0.0, 1.0, 1.0, 0.0));
        assert_eq!(p.injected_until(2.0), 1.0);
        assert_eq!(p.injected_until(3.0), 2.5);
        assert_eq!(p.total_amount(), Some(3.0));
        assert!(TimeProfile::Tabulated { times_s: vec![1.0, 1.0], rates: vec![0.0, 1.0] }.validate().is_err());
        assert!(TimeProfile::Tabulated { times_s: vec![1.0, 2.0], rates: vec![-1.0, 1.0] }.validate().is_err());
    }

    #[test]
    fn a_pulse_train_is_the_sum_of_its_members() {
        let p = TimeProfile::Train(vec![
            TimeProfile::Impulse { time_s: 0.0, amount: 1.0 },
            TimeProfile::Rectangular { start_s: 1.0, end_s: 2.0, amount: 0.5 },
            TimeProfile::Impulse { time_s: 3.0, amount: 0.25 },
        ]);
        p.validate().unwrap();
        assert_eq!(p.impulses(), vec![(0.0, 1.0), (3.0, 0.25)]);
        assert_eq!(p.rate(1.5), 0.5);
        assert_eq!(p.injected_until(2.5), 1.5);
        assert_eq!(p.injected_until(3.0), 1.75);
        assert_eq!(p.breakpoints(), vec![1.0, 2.0]);
    }

    #[test]
    fn invalid_profiles_are_refused() {
        for bad in [
            TimeProfile::Impulse { time_s: -1.0, amount: 1.0 },
            TimeProfile::Impulse { time_s: 0.0, amount: -1.0 },
            TimeProfile::Rectangular { start_s: 2.0, end_s: 2.0, amount: 1.0 },
            TimeProfile::Gaussian { centre_s: 0.0, sigma_s: 0.0, amount: 1.0 },
            TimeProfile::Feed { start_s: 1.0, end_s: Some(0.5), rate: 1.0 },
            TimeProfile::Feed { start_s: 0.0, end_s: None, rate: f64::NAN },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
    }

    #[test]
    fn profiles_are_described_in_words_with_units() {
        let cases = [
            (TimeProfile::Impulse { time_s: 1e-9, amount: 0.5 }, "impulse at 1e-9 s, amount 0.5"),
            (TimeProfile::Rectangular { start_s: 0.0, end_s: 2e-6, amount: 1.0 }, "rectangular pulse 0 .. 2e-6 s, amount 1"),
            (
                TimeProfile::Gaussian { centre_s: 1e-8, sigma_s: 1e-9, amount: 2.0 },
                "Gaussian pulse centred at 1e-8 s, standard deviation 1e-9 s, amount 2",
            ),
            (TimeProfile::Feed { start_s: 0.0, end_s: None, rate: 1e3 }, "constant feed 1000 1/s from 0 s, open-ended"),
            (TimeProfile::Feed { start_s: 1e-6, end_s: Some(1e-3), rate: 0.25 }, "constant feed 0.25 1/s from 1e-6 s to 0.001 s"),
            (
                TimeProfile::PrecursorDecay { start_s: 0.0, formation_rate_s_inv: 1e5, total_loss_rate_s_inv: 2e5, precursor_amount: 1.0 },
                "precursor decay from 0 s: k_f = 100000 1/s, k_tot = 2e5 1/s, precursor amount 1",
            ),
            (
                TimeProfile::Tabulated { times_s: vec![0.0, 1e-7, 1e-6], rates: vec![0.0, 1.0, 0.0] },
                "tabulated rate, 3 points 0 .. 1e-6 s",
            ),
        ];
        for (p, text) in &cases {
            assert_eq!(p.to_string(), *text);
        }
        let train = TimeProfile::Train(vec![cases[0].0.clone(), cases[3].0.clone()]);
        assert_eq!(train.to_string(), "train of 2: impulse at 1e-9 s, amount 0.5; constant feed 1000 1/s from 0 s, open-ended");
    }
}
