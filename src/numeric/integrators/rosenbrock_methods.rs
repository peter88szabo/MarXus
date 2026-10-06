//! Coefficients of the Rosenbrock methods, adapted from the Rosenbrock integrator of KPP, the Kinetic
//! PreProcessor (int/rosenbrock.f90: (C) Adrian Sandu, August 2004, Virginia Polytechnic Institute and State
//! University; revised by Philipp Miehe and Adrian Sandu, May 2006; GNU General Public License v3, as MarXus).
//!
//! The methods are written in the form used by KPP (Hairer, Wanner, Solving Ordinary Differential
//! Equations II, Springer (1991, 1996), Section IV.7):
//!   G = 1/(h gamma_1) - df/dy(t0, y0),
//!   T_i = t0 + alpha_i h,  Y_i = y0 + sum_(j<i) A_ij K_j,
//!   G K_i = f(T_i, Y_i) + sum_(j<i) (C_ij / h) K_j + h gamma_i df/dt(t0, y0),
//!   y1 = y0 + sum_i M_i K_i,  error estimate sum_i E_i K_i.
//! A and C are strictly lower triangular, stored row-wise: A(i, j) = a[(i-1)(i-2)/2 + j - 1] for i > j
//! (1-based i, j). `new_function[i]` is false when stage i reuses the function value of stage i - 1.
//! `error_order` (KPP ros_ELO) is the smaller order of the main and embedded method plus one; the new step
//! is h min(fac_max, max(fac_min, fac_safety / err^(1/error_order))).
//!
//! References of the methods (from the KPP documentation):
//! - Ros2: J. G. Verwer, E. J. Spee, J. G. Blom, W. Hundsdorfer, SIAM J. Sci. Comput. 20, 1456 (1999).
//! - Ros3, Rodas3: A. Sandu, J. G. Verwer, J. G. Blom, E. J. Spee, G. R. Carmichael, F. A. Potra, Atmos.
//!   Environ. 31, 3459 (1997).
//! - Ros4, Rodas4: E. Hairer, G. Wanner, Solving Ordinary Differential Equations II, Springer (1991, 1996).
//!
//! Not adopted: KPP's Rang3 (W-method of J. Rang, L. Angermann, BIT Numer. Math. 45, 761 (2005)). With its
//! coefficients the method itself is of order 3, but its embedded error estimate vanishes for a linear
//! problem (y' = -2y, h = 0.1 ... 2: estimate below 1e-15, true local error 3e-5 ... 9e-2), so that the
//! step-size control accepts any step (checked 2026-10-06, `reports/direct_time_integration.md`).

/// A Rosenbrock method of KPP.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RosenbrockMethod {
    /// L-stable, 2 stages, order 2 (Verwer et al. 1999).
    Ros2,
    /// L-stable, 3 stages, order 3, 2 function evaluations (Sandu et al. 1997).
    Ros3,
    /// L-stable, 4 stages, order 4, embedded order 3 (Hairer, Wanner).
    #[default]
    Ros4,
    /// Stiffly accurate, 4 stages, order 3 (Sandu et al. 1997); KPP's default.
    Rodas3,
    /// Stiffly accurate, 6 stages, order 4 (Hairer, Wanner).
    Rodas4,
}

/// Coefficients of a Rosenbrock method (see the module documentation for their meaning).
#[derive(Debug, Clone, PartialEq)]
pub struct RosenbrockTableau {
    pub name: &'static str,
    pub stages: usize,
    /// Order of the main method.
    pub order: usize,
    pub a: Vec<f64>,
    pub c: Vec<f64>,
    pub m: Vec<f64>,
    pub e: Vec<f64>,
    pub alpha: Vec<f64>,
    pub gamma: Vec<f64>,
    pub new_function: Vec<bool>,
    pub error_order: f64,
}

impl RosenbrockTableau {
    /// A(i, j), 0-based stage indices, i > j.
    pub fn a(&self, i: usize, j: usize) -> f64 {
        self.a[i * (i - 1) / 2 + j]
    }

    /// C(i, j), 0-based stage indices, i > j.
    pub fn c(&self, i: usize, j: usize) -> f64 {
        self.c[i * (i - 1) / 2 + j]
    }
}

impl RosenbrockMethod {
    /// All methods.
    pub const ALL: [RosenbrockMethod; 5] = [
        Self::Ros2,
        Self::Ros3,
        Self::Ros4,
        Self::Rodas3,
        Self::Rodas4,
    ];

    /// The coefficients of the method (KPP int/rosenbrock.f90, subroutines Ros2 ... Rodas4).
    pub fn tableau(self) -> RosenbrockTableau {
        match self {
            Self::Ros2 => {
                let g = 1.0 + 1.0 / 2.0_f64.sqrt();
                RosenbrockTableau {
                    name: "ROS-2",
                    stages: 2,
                    order: 2,
                    a: vec![1.0 / g],
                    c: vec![-2.0 / g],
                    m: vec![3.0 / (2.0 * g), 1.0 / (2.0 * g)],
                    e: vec![1.0 / (2.0 * g), 1.0 / (2.0 * g)],
                    alpha: vec![0.0, 1.0],
                    gamma: vec![g, -g],
                    new_function: vec![true, true],
                    error_order: 2.0,
                }
            }
            Self::Ros3 => RosenbrockTableau {
                name: "ROS-3",
                stages: 3,
                order: 3,
                a: vec![1.0, 1.0, 0.0],
                c: vec![
                    -0.10156171083877702091975600115545E+01,
                    0.40759956452537699824805835358067E+01,
                    0.92076794298330791242156818474003E+01,
                ],
                m: vec![
                    0.1E+01,
                    0.61697947043828245592553615689730E+01,
                    -0.42772256543218573326238373806514,
                ],
                e: vec![
                    0.5,
                    -0.29079558716805469821718236208017E+01,
                    0.22354069897811569627360909276199,
                ],
                alpha: vec![
                    0.0,
                    0.43586652150845899941601945119356,
                    0.43586652150845899941601945119356,
                ],
                gamma: vec![
                    0.43586652150845899941601945119356,
                    0.24291996454816804366592249683314,
                    0.21851380027664058511513169485832E+01,
                ],
                new_function: vec![true, true, false],
                error_order: 3.0,
            },
            Self::Ros4 => {
                let a2 = 0.1867943637803922E+01;
                let a3 = 0.2344449711399156;
                let alpha3 = 0.6552168638155900;
                RosenbrockTableau {
                    name: "ROS-4",
                    stages: 4,
                    order: 4,
                    a: vec![0.2000000000000000E+01, a2, a3, a2, a3, 0.0],
                    c: vec![
                        -0.7137615036412310E+01,
                        0.2580708087951457E+01,
                        0.6515950076447975,
                        -0.2137148994382534E+01,
                        -0.3214669691237626,
                        -0.6949742501781779,
                    ],
                    m: vec![
                        0.2255570073418735E+01,
                        0.2870493262186792,
                        0.4353179431840180,
                        0.1093502252409163E+01,
                    ],
                    e: vec![
                        -0.2815431932141155,
                        -0.7276199124938920E-01,
                        -0.1082196201495311,
                        -0.1093502252409163E+01,
                    ],
                    alpha: vec![0.0, 0.1145640000000000E+01, alpha3, alpha3],
                    gamma: vec![
                        0.5728200000000000,
                        -0.1769193891319233E+01,
                        0.7592633437920482,
                        -0.1049021087100450,
                    ],
                    new_function: vec![true, true, true, false],
                    error_order: 4.0,
                }
            }
            Self::Rodas3 => RosenbrockTableau {
                name: "RODAS-3",
                stages: 4,
                order: 3,
                a: vec![0.0, 2.0, 0.0, 2.0, 0.0, 1.0],
                c: vec![4.0, 1.0, -1.0, 1.0, -1.0, -(8.0 / 3.0)],
                m: vec![2.0, 0.0, 1.0, 1.0],
                e: vec![0.0, 0.0, 0.0, 1.0],
                alpha: vec![0.0, 0.0, 1.0, 1.0],
                gamma: vec![0.5, 1.5, 0.0, 0.0],
                new_function: vec![true, false, true, true],
                error_order: 3.0,
            },
            Self::Rodas4 => {
                let a7 = 0.1221224509226641E+01;
                let a8 = 0.6019134481288629E+01;
                let a9 = 0.1253708332932087E+02;
                let a10 = -0.6878860361058950;
                RosenbrockTableau {
                    name: "RODAS-4",
                    stages: 6,
                    order: 4,
                    a: vec![
                        0.1544000000000000E+01,
                        0.9466785280815826,
                        0.2557011698983284,
                        0.3314825187068521E+01,
                        0.2896124015972201E+01,
                        0.9986419139977817,
                        a7,
                        a8,
                        a9,
                        a10,
                        a7,
                        a8,
                        a9,
                        a10,
                        1.0,
                    ],
                    c: vec![
                        -0.5668800000000000E+01,
                        -0.2430093356833875E+01,
                        -0.2063599157091915,
                        -0.1073529058151375,
                        -0.9594562251023355E+01,
                        -0.2047028614809616E+02,
                        0.7496443313967647E+01,
                        -0.1024680431464352E+02,
                        -0.3399990352819905E+02,
                        0.1170890893206160E+02,
                        0.8083246795921522E+01,
                        -0.7981132988064893E+01,
                        -0.3152159432874371E+02,
                        0.1631930543123136E+02,
                        -0.6058818238834054E+01,
                    ],
                    m: vec![a7, a8, a9, a10, 1.0, 1.0],
                    e: vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                    alpha: vec![0.000, 0.386, 0.210, 0.630, 1.000, 1.000],
                    gamma: vec![
                        0.2500000000000000,
                        -0.1043000000000000,
                        0.1035000000000000,
                        -0.3620000000000023E-01,
                        0.0,
                        0.0,
                    ],
                    new_function: vec![true; 6],
                    error_order: 4.0,
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_tableau_has_consistent_dimensions() {
        for method in RosenbrockMethod::ALL {
            let t = method.tableau();
            let s = t.stages;
            let lower = s * (s - 1) / 2;
            assert_eq!((t.a.len(), t.c.len()), (lower, lower), "{}", t.name);
            for v in [&t.m, &t.e, &t.alpha, &t.gamma] {
                assert_eq!(v.len(), s, "{}", t.name);
            }
            assert_eq!(t.new_function.len(), s, "{}", t.name);
            assert!(t.new_function[0], "{}", t.name);
            assert_eq!(t.alpha[0], 0.0, "{}", t.name);
            assert!(
                t.gamma[0] > 0.0,
                "{}: G = 1/(h gamma_1) - J needs gamma_1 > 0",
                t.name
            );
        }
    }

    #[test]
    fn row_wise_storage_maps_stage_indices() {
        // Rodas4: A(5, 1..4) (1-based) = a7..a10, A(6, 5) = 1 (KPP ros_A(7..10), ros_A(15)).
        let t = RosenbrockMethod::Rodas4.tableau();
        assert_eq!(t.a(4, 0), 0.1221224509226641E+01);
        assert_eq!(t.a(4, 3), -0.6878860361058950);
        assert_eq!(t.a(5, 4), 1.0);
        assert_eq!(t.c(1, 0), -0.5668800000000000E+01);
    }
}
