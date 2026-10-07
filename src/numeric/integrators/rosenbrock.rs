//! Adaptive Rosenbrock integrator for stiff ordinary differential equations dy/dt = f(t, y), adapted from
//! the Rosenbrock integrator of KPP, the Kinetic PreProcessor (int/rosenbrock.f90, subroutine
//! ros_Integrator: (C) Adrian Sandu, August 2004, Virginia Polytechnic Institute and State University;
//! revised by Philipp Miehe and Adrian Sandu, May 2006; GNU General Public License v3, as MarXus). The
//! implementation follows Hairer, Wanner, Solving Ordinary Differential Equations II, Springer (1991,
//! 1996), Section IV.7; the methods are in `rosenbrock_methods.rs`.
//!
//! One step from t to t + h (KPP ros_Integrator):
//! 1. f(t, y), df/dt (by finite differences, unless the system is autonomous) and the matrix
//!    G = 1/(h gamma_1) - df/dy, prepared (factorized) by the system (`StiffSystem::prepare`); if the
//!    factorization fails, h is halved, at most five times.
//! 2. The stages G K_i = f(T_i, Y_i) + sum_(j<i) (C_ij/h) K_j + h gamma_i df/dt.
//! 3. y_new = y + sum M_i K_i and the error estimate sum E_i K_i, measured in the scaled norm
//!    err = sqrt(mean_k [e_k / (atol_k + rtol_k max(|y_k|, |y_new_k|))]^2) (KPP ros_ErrorNorm, at least 1e-10).
//! 4. Step-size control: h_new = h min(fac_max, max(fac_min, fac_safety / err^(1/ELO))). The step is
//!    accepted if err <= 1 (or h <= h_min); after an accepted step that followed a rejection h does not
//!    grow; after two successive rejections h is multiplied by fac_rejection.

use super::rosenbrock_methods::RosenbrockMethod;

/// The system dy/dt = f(t, y) and its linear algebra.
pub trait StiffSystem {
    /// Number of unknowns.
    fn dimension(&self) -> usize;
    /// f(t, y) into `dydt`.
    fn rhs(&self, t: f64, y: &[f64], dydt: &mut [f64]);
    /// Prepares (factorizes) G = shift I - df/dy(t, y) for the following `solve` calls, with
    /// shift = 1/(h gamma_1). An error if G is singular.
    fn prepare(&mut self, t: f64, y: &[f64], shift: f64) -> Result<(), String>;
    /// b := G^-1 b with the prepared G.
    fn solve(&self, b: &mut [f64]) -> Result<(), String>;
}

/// Parameters of the integration (KPP ICNTRL and RCNTRL, with the KPP defaults).
#[derive(Debug, Clone, PartialEq)]
pub struct RosenbrockOptions {
    pub method: RosenbrockMethod,
    /// Absolute tolerance: one value for all components, or one per component.
    pub absolute_tolerance: Vec<f64>,
    /// Relative tolerance: one value for all components, or one per component.
    pub relative_tolerance: Vec<f64>,
    /// f independent of t: df/dt = 0 is not computed (KPP ICNTRL(1) = 1).
    pub autonomous: bool,
    /// Lower bound of the step (0: none; KPP recommends 0).
    pub h_min: f64,
    /// Upper bound of the step (0: the integration interval).
    pub h_max: f64,
    /// First step (0: max(h_min, 1e-5)).
    pub h_start: f64,
    pub fac_min: f64,
    pub fac_max: f64,
    pub fac_rejection: f64,
    pub fac_safety: f64,
    pub max_steps: usize,
    /// Set negative components of the solution to zero after every accepted step (KPP ICNTRL(16) = 1).
    pub clip_negative: bool,
    /// MarXus extension (default off, KPP behaviour): round every new step down to a power of two. The steps
    /// are never larger than KPP's, and step sizes repeat, so that a system with a constant Jacobian (the
    /// linear master equation) can reuse the factorization of 1/(h gamma) - J.
    pub power_of_two_steps: bool,
}

impl RosenbrockOptions {
    /// KPP defaults with the given method and scalar tolerances.
    pub fn new(method: RosenbrockMethod, relative_tolerance: f64, absolute_tolerance: f64) -> Self {
        Self {
            method,
            absolute_tolerance: vec![absolute_tolerance],
            relative_tolerance: vec![relative_tolerance],
            autonomous: false,
            h_min: 0.0,
            h_max: 0.0,
            h_start: 0.0,
            fac_min: 0.2,
            fac_max: 6.0,
            fac_rejection: 0.1,
            fac_safety: 0.9,
            max_steps: 200_000,
            clip_negative: false,
            power_of_two_steps: false,
        }
    }
}

/// Work done by an integration (KPP ISTATUS and RSTATUS).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct IntegrationStatistics {
    pub function_evaluations: usize,
    pub factorizations: usize,
    pub solves: usize,
    pub steps: usize,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    /// Last accepted step.
    pub last_step: f64,
    /// Next predicted step: the first step of a following integration from the end point.
    pub next_step: f64,
}

/// Integrates dy/dt = f(t, y) from `t_start` to `t_end`; `y` holds the initial values and receives the
/// solution at `t_end`.
pub fn integrate<S: StiffSystem>(
    system: &mut S,
    y: &mut [f64],
    t_start: f64,
    t_end: f64,
    options: &RosenbrockOptions,
) -> Result<IntegrationStatistics, String> {
    let n = system.dimension();
    if y.len() != n {
        return Err(format!(
            "Rosenbrock: {} initial values for a system of dimension {n}.",
            y.len()
        ));
    }
    let tableau = options.method.tableau();
    let roundoff = f64::EPSILON;
    check_tolerances(options, n, roundoff)?;
    // MarXus extension: the round-off of time values relative to the time scale of the interval (KPP compares
    // with the absolute machine epsilon, which assumes times of order one; master-equation times reach 1e-12 s).
    let time_roundoff = roundoff * t_start.abs().max(t_end.abs()).max(f64::MIN_POSITIVE);
    let atol = |k: usize| {
        options.absolute_tolerance[if options.absolute_tolerance.len() == 1 {
            0
        } else {
            k
        }]
    };
    let rtol = |k: usize| {
        options.relative_tolerance[if options.relative_tolerance.len() == 1 {
            0
        } else {
            k
        }]
    };

    // Step bounds (KPP Rosenbrock: Hmin, Hmax, Hstart with DeltaMin = 1e-5).
    const DELTA_MIN: f64 = 1.0e-5;
    let interval = (t_end - t_start).abs();
    let h_min = options.h_min.max(0.0);
    let h_max = if options.h_max > 0.0 {
        options.h_max.min(interval)
    } else {
        interval
    };
    let h_start = if options.h_start > 0.0 {
        options.h_start.min(interval)
    } else {
        h_min.max(DELTA_MIN)
    };

    let mut stats = IntegrationStatistics::default();
    let direction = if t_end >= t_start { 1.0 } else { -1.0 };
    let mut t = t_start;
    // MarXus extension: steps rounded down to powers of two (`RosenbrockOptions::power_of_two_steps`).
    let grid = |h: f64| {
        if options.power_of_two_steps && h > 0.0 {
            2f64.powf(h.log2().floor()).max(h_min)
        } else {
            h
        }
    };
    let mut h = h_start.max(h_min).min(h_max);
    if h <= 10.0 * time_roundoff {
        h = DELTA_MIN;
    }
    h = grid(h);
    let mut reject_last = false;
    let mut reject_more = false;

    let s = tableau.stages;
    let mut k = vec![vec![0.0; n]; s];
    let mut f0 = vec![0.0; n];
    let mut f = vec![0.0; n];
    let mut df_dt = vec![0.0; n];
    let mut y_new = vec![0.0; n];
    let mut y_stage = vec![0.0; n];

    while (direction > 0.0 && (t - t_end) + time_roundoff <= 0.0)
        || (direction < 0.0 && (t_end - t) + time_roundoff <= 0.0)
    {
        if stats.steps > options.max_steps {
            return Err(format!(
                "Rosenbrock: more than {} steps (t = {t:e}, h = {h:e}).",
                options.max_steps
            ));
        }
        if t + 0.1 * h * direction == t || h <= time_roundoff {
            return Err(format!(
                "Rosenbrock: step size too small (t = {t:e}, h = {h:e})."
            ));
        }
        // Do not step beyond t_end.
        h = h.min((t_end - t).abs());

        system.rhs(t, y, &mut f0);
        stats.function_evaluations += 1;
        if !options.autonomous {
            // df/dt by finite differences (KPP ros_FunTimeDerivative).
            let delta = roundoff.sqrt() * 1.0e-6_f64.max(t.abs());
            system.rhs(t + delta, y, &mut df_dt);
            stats.function_evaluations += 1;
            for (d, &f0k) in df_dt.iter_mut().zip(&f0) {
                *d = (*d - f0k) / delta;
            }
        }

        // Repeat the step until it is accepted.
        loop {
            prepare_matrix(
                system,
                t,
                y,
                &mut h,
                direction,
                tableau.gamma[0],
                &mut stats,
            )?;
            for i in 0..s {
                if i == 0 {
                    f.copy_from_slice(&f0);
                } else if tableau.new_function[i] {
                    y_stage.copy_from_slice(y);
                    for j in 0..i {
                        let a = tableau.a(i, j);
                        if a != 0.0 {
                            for (ys, kj) in y_stage.iter_mut().zip(&k[j]) {
                                *ys += a * kj;
                            }
                        }
                    }
                    system.rhs(t + tableau.alpha[i] * direction * h, &y_stage, &mut f);
                    stats.function_evaluations += 1;
                }
                let mut stage = f.clone();
                for j in 0..i {
                    let hc = tableau.c(i, j) / (direction * h);
                    for (st, kj) in stage.iter_mut().zip(&k[j]) {
                        *st += hc * kj;
                    }
                }
                if !options.autonomous && tableau.gamma[i] != 0.0 {
                    let hg = direction * h * tableau.gamma[i];
                    for (st, d) in stage.iter_mut().zip(&df_dt) {
                        *st += hg * d;
                    }
                }
                system.solve(&mut stage)?;
                stats.solves += 1;
                k[i] = stage;
            }

            // New solution and error estimate.
            y_new.copy_from_slice(y);
            let mut error_sum = 0.0;
            for (c, yn) in y_new.iter_mut().enumerate() {
                let mut e = 0.0;
                for i in 0..s {
                    *yn += tableau.m[i] * k[i][c];
                    e += tableau.e[i] * k[i][c];
                }
                let scale = atol(c) + rtol(c) * y[c].abs().max(yn.abs());
                error_sum += (e / scale).powi(2);
            }
            let error = (error_sum / n as f64).sqrt().max(1.0e-10);

            // New step bounded by fac_min <= h_new/h <= fac_max.
            let factor = options.fac_max.min(
                options
                    .fac_min
                    .max(options.fac_safety / error.powf(1.0 / tableau.error_order)),
            );
            let mut h_new = h * factor;
            stats.steps += 1;
            if error <= 1.0 || h <= h_min {
                stats.accepted_steps += 1;
                for (yk, &yn) in y.iter_mut().zip(&y_new) {
                    *yk = if options.clip_negative {
                        yn.max(0.0)
                    } else {
                        yn
                    };
                }
                t += direction * h;
                h_new = h_min.max(h_new.min(h_max));
                if reject_last {
                    // No step increase after a rejected step.
                    h_new = h_new.min(h);
                }
                h_new = grid(h_new);
                stats.last_step = h;
                stats.next_step = h_new;
                reject_last = false;
                reject_more = false;
                h = h_new;
                break;
            }
            if reject_more {
                h_new = h * options.fac_rejection;
            }
            reject_more = reject_last;
            reject_last = true;
            h = grid(h_new);
            if stats.accepted_steps >= 1 {
                stats.rejected_steps += 1;
            }
        }
    }
    Ok(stats)
}

/// Tolerances must be positive and the relative ones between 10 eps and 1 (KPP Rosenbrock).
fn check_tolerances(options: &RosenbrockOptions, n: usize, roundoff: f64) -> Result<(), String> {
    for (name, values) in [
        ("absolute", &options.absolute_tolerance),
        ("relative", &options.relative_tolerance),
    ] {
        if values.len() != 1 && values.len() != n {
            return Err(format!(
                "Rosenbrock: {} {name} tolerances for {n} unknowns (1 or {n}).",
                values.len()
            ));
        }
    }
    if options.absolute_tolerance.iter().any(|&a| !(a > 0.0)) {
        return Err("Rosenbrock: the absolute tolerance must be positive.".into());
    }
    if options
        .relative_tolerance
        .iter()
        .any(|&r| !(r > 10.0 * roundoff && r < 1.0))
    {
        return Err(format!(
            "Rosenbrock: the relative tolerance must lie between {:e} and 1.",
            10.0 * roundoff
        ));
    }
    Ok(())
}

/// G = 1/(h gamma_1) - df/dy prepared by the system; if it fails, h is halved, at most five times
/// (KPP ros_PrepareMatrix).
fn prepare_matrix<S: StiffSystem>(
    system: &mut S,
    t: f64,
    y: &[f64],
    h: &mut f64,
    direction: f64,
    gamma_1: f64,
    stats: &mut IntegrationStatistics,
) -> Result<(), String> {
    let mut failures = 0;
    loop {
        stats.factorizations += 1;
        match system.prepare(t, y, 1.0 / (direction * *h * gamma_1)) {
            Ok(()) => return Ok(()),
            Err(e) => {
                failures += 1;
                if failures > 5 {
                    return Err(format!("Rosenbrock: the matrix 1/(h gamma) - J could not be factorized ({e}); t = {t:e}, h = {h:e}."));
                }
                *h *= 0.5;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numeric::dense_inverse::invert_dense;

    /// A small dense system with an explicit Jacobian; G is inverted densely.
    struct Dense<F, J>
    where
        F: Fn(f64, &[f64], &mut [f64]),
        J: Fn(f64, &[f64]) -> Vec<Vec<f64>>,
    {
        n: usize,
        f: F,
        jacobian: J,
        inverse: Vec<Vec<f64>>,
    }

    impl<F, J> Dense<F, J>
    where
        F: Fn(f64, &[f64], &mut [f64]),
        J: Fn(f64, &[f64]) -> Vec<Vec<f64>>,
    {
        fn new(n: usize, f: F, jacobian: J) -> Self {
            Self {
                n,
                f,
                jacobian,
                inverse: Vec::new(),
            }
        }
    }

    impl<F, J> StiffSystem for Dense<F, J>
    where
        F: Fn(f64, &[f64], &mut [f64]),
        J: Fn(f64, &[f64]) -> Vec<Vec<f64>>,
    {
        fn dimension(&self) -> usize {
            self.n
        }
        fn rhs(&self, t: f64, y: &[f64], dydt: &mut [f64]) {
            (self.f)(t, y, dydt)
        }
        fn prepare(&mut self, t: f64, y: &[f64], shift: f64) -> Result<(), String> {
            let jac = (self.jacobian)(t, y);
            let g: Vec<Vec<f64>> = (0..self.n)
                .map(|i| {
                    (0..self.n)
                        .map(|j| if i == j { shift } else { 0.0 } - jac[i][j])
                        .collect()
                })
                .collect();
            self.inverse = invert_dense(&g)?;
            Ok(())
        }
        fn solve(&self, b: &mut [f64]) -> Result<(), String> {
            let x: Vec<f64> = self
                .inverse
                .iter()
                .map(|row| row.iter().zip(b.iter()).map(|(a, v)| a * v).sum())
                .collect();
            b.copy_from_slice(&x);
            Ok(())
        }
    }

    #[test]
    fn exponential_decay_is_reproduced_by_every_method() {
        for method in RosenbrockMethod::ALL {
            let mut system = Dense::new(
                1,
                |_t, y: &[f64], d: &mut [f64]| d[0] = -2.0 * y[0],
                |_t, _y: &[f64]| vec![vec![-2.0]],
            );
            let mut y = [1.0];
            let mut options = RosenbrockOptions::new(method, 1e-8, 1e-14);
            options.autonomous = true;
            let stats = integrate(&mut system, &mut y, 0.0, 3.0, &options).unwrap();
            let exact = (-6.0_f64).exp();
            assert!(
                (y[0] / exact - 1.0).abs() < 1e-5,
                "{method:?}: {} vs {exact}",
                y[0]
            );
            assert!(
                stats.accepted_steps > 3 && stats.steps >= stats.accepted_steps,
                "{method:?}: {stats:?}"
            );
        }
    }

    #[test]
    fn a_stiff_non_autonomous_problem_takes_large_steps() {
        // Prothero-Robinson: y' = -lambda (y - sin t) + cos t, y(0) = 0, solution sin t; lambda = 1e6
        // (stiffness ratio 1e6 over t in [0, 10]). df/dt = lambda cos t - sin t by finite differences.
        let lambda = 1.0e6;
        for method in RosenbrockMethod::ALL {
            let mut system = Dense::new(
                1,
                move |t, y: &[f64], d: &mut [f64]| d[0] = -lambda * (y[0] - t.sin()) + t.cos(),
                move |_t, _y: &[f64]| vec![vec![-lambda]],
            );
            let mut y = [0.0];
            let options = RosenbrockOptions::new(method, 1e-6, 1e-10);
            let stats = integrate(&mut system, &mut y, 0.0, 10.0, &options).unwrap();
            assert!(
                (y[0] - 10.0_f64.sin()).abs() < 1e-5,
                "{method:?}: {} vs {}",
                y[0],
                10.0_f64.sin()
            );
            // Ros2 (order 2) needs about 6000 steps at this tolerance; an explicit method would need ~1e7.
            assert!(stats.accepted_steps < 20_000, "{method:?}: {stats:?}");
        }
    }

    #[test]
    fn the_robertson_problem_reaches_its_reference_values() {
        // Robertson's chemical kinetics (Hairer, Wanner, Section IV.1): y1' = -0.04 y1 + 1e4 y2 y3,
        // y2' = 0.04 y1 - 1e4 y2 y3 - 3e7 y2^2, y3' = 3e7 y2^2; y(0) = (1, 0, 0). Reference at t = 40:
        // (0.7158270687, 9.185534764e-6, 0.2841637457).
        let f = |_t: f64, y: &[f64], d: &mut [f64]| {
            d[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
            d[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1];
            d[2] = 3.0e7 * y[1] * y[1];
        };
        let jac = |_t: f64, y: &[f64]| {
            vec![
                vec![-0.04, 1.0e4 * y[2], 1.0e4 * y[1]],
                vec![0.04, -1.0e4 * y[2] - 6.0e7 * y[1], -1.0e4 * y[1]],
                vec![0.0, 6.0e7 * y[1], 0.0],
            ]
        };
        for method in [
            RosenbrockMethod::Ros3,
            RosenbrockMethod::Ros4,
            RosenbrockMethod::Rodas3,
            RosenbrockMethod::Rodas4,
        ] {
            let mut system = Dense::new(3, f, jac);
            let mut y = [1.0, 0.0, 0.0];
            let mut options = RosenbrockOptions::new(method, 1e-8, 1e-14);
            options.autonomous = true;
            integrate(&mut system, &mut y, 0.0, 40.0, &options).unwrap();
            let reference = [0.7158270687, 9.185534764e-6, 0.2841637457];
            for k in 0..3 {
                assert!(
                    (y[k] / reference[k] - 1.0).abs() < 1e-5,
                    "{method:?} y{}: {} vs {}",
                    k + 1,
                    y[k],
                    reference[k]
                );
            }
            // Linear invariant y1 + y2 + y3 = 1 (preserved by Rosenbrock methods with the exact Jacobian).
            assert!((y.iter().sum::<f64>() - 1.0).abs() < 1e-12, "{method:?}");
        }
    }

    #[test]
    fn every_method_has_its_nominal_order() {
        // y' = -y^3, y(0) = 1: y(1) = 1/sqrt(3). Fixed steps (h_min = h_max = h_start = h, every step
        // accepted); the observed order log2(err(h)/err(h/2)). (y' = -y^2 is unsuitable: Rodas3 integrates
        // it exactly, to rounding.)
        for method in RosenbrockMethod::ALL {
            let error = |h: f64| {
                let mut system = Dense::new(
                    1,
                    |_t, y: &[f64], d: &mut [f64]| d[0] = -y[0].powi(3),
                    |_t, y: &[f64]| vec![vec![-3.0 * y[0] * y[0]]],
                );
                let mut y = [1.0];
                let mut options = RosenbrockOptions::new(method, 1e-3, 1e-3);
                options.autonomous = true;
                options.h_min = h;
                options.h_max = h;
                options.h_start = h;
                integrate(&mut system, &mut y, 0.0, 1.0, &options).unwrap();
                (y[0] - 1.0 / 3.0_f64.sqrt()).abs()
            };
            let observed = (error(0.05) / error(0.025)).log2();
            let nominal = method.tableau().order as f64;
            eprintln!("{method:?}: observed order {observed:.2}, nominal {nominal}");
            // At least the nominal order (a coding error lowers it; Rodas4 shows 4.5 at these steps, where
            // higher-order terms still contribute).
            assert!(
                observed > nominal - 0.35 && observed < nominal + 1.0,
                "{method:?}: observed order {observed:.2}, nominal {nominal}"
            );
        }
    }

    #[test]
    fn power_of_two_steps_keep_the_accuracy_and_repeat_the_step_sizes() {
        // Robertson as above, with steps rounded down to powers of two: the same reference values, and the
        // matrix G = 1/(h gamma) - J is prepared for few distinct step sizes.
        let f = |_t: f64, y: &[f64], d: &mut [f64]| {
            d[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
            d[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1];
            d[2] = 3.0e7 * y[1] * y[1];
        };
        let jac = |_t: f64, y: &[f64]| {
            vec![
                vec![-0.04, 1.0e4 * y[2], 1.0e4 * y[1]],
                vec![0.04, -1.0e4 * y[2] - 6.0e7 * y[1], -1.0e4 * y[1]],
                vec![0.0, 6.0e7 * y[1], 0.0],
            ]
        };
        struct Counting<S> {
            inner: S,
            shifts: Vec<f64>,
        }
        impl<S: StiffSystem> StiffSystem for Counting<S> {
            fn dimension(&self) -> usize {
                self.inner.dimension()
            }
            fn rhs(&self, t: f64, y: &[f64], d: &mut [f64]) {
                self.inner.rhs(t, y, d)
            }
            fn prepare(&mut self, t: f64, y: &[f64], shift: f64) -> Result<(), String> {
                self.shifts.push(shift);
                self.inner.prepare(t, y, shift)
            }
            fn solve(&self, b: &mut [f64]) -> Result<(), String> {
                self.inner.solve(b)
            }
        }
        let mut system = Counting {
            inner: Dense::new(3, f, jac),
            shifts: Vec::new(),
        };
        let mut y = [1.0, 0.0, 0.0];
        let mut options = RosenbrockOptions::new(RosenbrockMethod::Rodas4, 1e-8, 1e-14);
        options.autonomous = true;
        options.power_of_two_steps = true;
        let stats = integrate(&mut system, &mut y, 0.0, 40.0, &options).unwrap();
        let reference = [0.7158270687, 9.185534764e-6, 0.2841637457];
        for k in 0..3 {
            assert!(
                (y[k] / reference[k] - 1.0).abs() < 1e-5,
                "y{}: {} vs {}",
                k + 1,
                y[k],
                reference[k]
            );
        }
        // Every step except the last one (truncated at t_end) has h = 2^k: shift = 2^-k / gamma_1.
        let gamma_1 = RosenbrockMethod::Rodas4.tableau().gamma[0];
        let off_grid = system
            .shifts
            .iter()
            .filter(|&&s| ((1.0 / (s * gamma_1)).log2().fract()).abs() > 1e-12)
            .count();
        assert!(off_grid <= 1, "{off_grid} steps not a power of two");
        let mut distinct: Vec<f64> = system.shifts.clone();
        distinct.sort_by(|a, b| a.partial_cmp(b).unwrap());
        distinct.dedup();
        assert!(
            distinct.len() * 2 < stats.steps,
            "{} distinct step sizes in {} steps",
            distinct.len(),
            stats.steps
        );
    }

    #[test]
    fn unreasonable_tolerances_are_refused() {
        let mut system = Dense::new(
            1,
            |_t, y: &[f64], d: &mut [f64]| d[0] = -y[0],
            |_t, _y: &[f64]| vec![vec![-1.0]],
        );
        let mut y = [1.0];
        for (rtol, atol) in [(1e-20, 1e-10), (1.0, 1e-10), (1e-6, 0.0)] {
            let options = RosenbrockOptions::new(RosenbrockMethod::Ros4, rtol, atol);
            assert!(integrate(&mut system, &mut y, 0.0, 1.0, &options)
                .unwrap_err()
                .contains("tolerance"));
        }
    }
}
