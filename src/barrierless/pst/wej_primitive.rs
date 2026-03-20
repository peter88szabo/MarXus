#![allow(non_snake_case)]

use crate::numeric::lanczos_gamma::gamma_func;
use crate::rrkm::sum_and_density::beyer_swinehart_counting;
use crate::barrierless::pst::collision_type::{ReactantType, CollisionType, classify_collision};

// ALL Wl(E,J) values are meant for the rotational/orbital contribution
// of the collision partners.

/*
===============================================================================
Here we compute 

    dW[i] = W(E_i, J) - W(E_i - dE, J)

instead of storing W(E_i, J) directly
===============================================================================

The functions that evaluate W(E, J) return cumulative sums of states:

    W(E, J) = number of rotational/orbital states with energy <= E

This is an integrated quantity. It does not represent how many states lie
within a specific energy interval.

The Beyer–Swinehart vibrational convolution requires the number of states
in each discrete energy grain, not the cumulative total. For an energy grid

    E_i = i * dE

the number of states in the i-th bin is

    ΔW(E_i, J) = W(E_i, J) - W(E_i - dE, J)

This is what is stored in dW[i].

-------------------------------------------------------------------------------
Interpretation
-------------------------------------------------------------------------------

dW[i] represents the number of states in the interval

    (E_i - dE, E_i]

not the total number of states up to E_i.

In the continuous limit:

    ρ(E, J) = dW(E, J) / dE

so on a discrete grid:

    W(E_i, J) - W(E_i - dE, J) ≈ ρ(E_i, J) * dE

Thus dW[i] can be interpreted as a discrete density-of-states contribution
(density × bin width), and is dimensionless.

-------------------------------------------------------------------------------
Why this form is required
-------------------------------------------------------------------------------

The total PST sum of states is constructed via convolution:

    W_all(E, J) = ∫ ρ_rot(E', J) * W_vib(E - E') dE'

In discrete form:

    W_all(E_i, J) = Σ_k [ΔW_rot(E_k, J)] * W_vib(E_i - E_k)

where

    ΔW_rot(E_k, J) = W_rot(E_k, J) - W_rot(E_k - dE, J)

Therefore the rotational/orbital contribution must be provided as bin counts,
not as cumulative W(E, J).

-------------------------------------------------------------------------------
Why W(E) is evaluated at both E and E - dE
-------------------------------------------------------------------------------

The two evaluations correspond to different cumulative quantities:

    W(E, ...)       = all states up to the upper edge of the bin
    W(E - dE, ...)  = all states up to the lower edge of the bin

Their difference isolates the states belonging to the current bin.
This is not redundant work, but the standard way to convert a cumulative
function into per-bin values.

-------------------------------------------------------------------------------
After convolution
-------------------------------------------------------------------------------

The rotational/orbital dW array is used as the starting vector for the
Beyer–Swinehart vibrational counting. Depending on the starting vector,
Beyer–Swinehart returns either a density-like grain count or a cumulative
sum of states.

Workflow:

    cumulative W(E, J)
        -> convert to ΔW(E, J) per bin
        -> vibrational Beyer–Swinehart counting
        -> internal rotor convolution (if needed)

-------------------------------------------------------------------------------
Interpolation
-------------------------------------------------------------------------------

When interpolation replaces exact counting, the same procedure applies:

    dW[i] = W_interp(E_i, J) - W_interp(E_i - dE, J)

Only the source of W(E, J) changes (exact vs interpolated), not the logic.
===============================================================================
*/

// -----------------------------------------------------------------------------
// NOTE ON CENTRIFUGAL BARRIER TREATMENT
// -----------------------------------------------------------------------------
// This implementation performs pure phase-space counting of rotational/orbital
// states under the constraint of total energy E.
//
// No explicit centrifugal barrier E0(J) is included here. In the formal PST
// framework (e.g. Olzmann), the number of states should depend on the available
// energy
//
//     E* = E - E0(J)
//
// where E0(J) is the maximum of the effective radial potential
// V_eff(r) = V(r) + L(L+1)/(2μr²).
//
// In the present code, no interaction potential V(r) is defined and no barrier
// search is performed. Therefore this implementation implicitly assumes
//
//     E0(J) = 0
//
// and counts all energetically allowed rotational/orbital states directly from E.
//
// As a consequence, the centrifugal barrier is not treated explicitly, and the
// result corresponds to a "free-rotor" phase-space model.
//
// If a physically correct PST treatment is required, the barrier height E0(J)
// must be computed externally (e.g. from an effective potential) and the energy
// should be shifted before calling these routines:
//
//     E → E - E0(J)
//
// States with E ≤ E0(J) should then be excluded.
// -----------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CountingMode {
    Exact,
    Interpol,
}

impl CountingMode {
    pub fn from_str(mode: &str) -> Result<Self, String> {
        match mode {
            "exact" => Ok(Self::Exact),
            "interpol" => Ok(Self::Interpol),
            _ => Err(format!(
                "invalid counting mode '{}'; use \"exact\" or \"interpol\"",
                mode
            )),
        }
    }
}

fn number_of_bins(emax: f64, dE: f64) -> usize {
    (emax / dE + 0.5).floor() as usize
}

fn energy_of_bin(bin_index: usize, dE: f64) -> f64 {
    (bin_index as f64) * dE
}

fn relative_error_percent(exact: f64, approx: f64) -> f64 {
    if exact.abs() < 1.0e-14 {
        0.0
    } else {
        ((exact - approx).abs() / exact.abs()) * 100.0
    }
}

/// Build rovibrational W(E,J) by:
/// 1. generating rotational/orbital grain counts dW(E,J),
/// 2. applying Beyer–Swinehart for vibrational convolution,
/// 3. convolving with free internal rotors if present.
/// Build rovibrational W(E,J) by:
/// 1. generating rotational/orbital grain counts dW(E,J),
/// 2. applying Beyer–Swinehart for vibrational convolution,
/// 3. convolving with free internal rotors if present.
///
/// Rotational constants:
/// - Linear: use brot
/// - Spherical: use brot
/// - Prolate: use arot and brot
/// - Oblate: use brot and crot
pub fn get_pst_rovibrational_WEJ(
    fragA: ReactantType,
    fragB: ReactantType,
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    arot1: f64,
    brot1: f64,
    crot1: f64,
    arot2: f64,
    brot2: f64,
    crot2: f64,
    nvib1: usize,
    nvib2: usize,
    vib_freqs1: &[f64],
    vib_freqs2: &[f64],
    n_internal_rotors1: usize,
    n_internal_rotors2: usize,
    internal_rot_consts1: &[f64],
    internal_rot_consts2: &[f64],
) -> Vec<f64> {
    let collision_type = classify_collision(fragA, fragB);
    let nebin = (emax / dE + 0.5).floor() as usize;

    let mut rotational_bin_counts = vec![0.0; nebin + 1];

    match collision_type {
        CollisionType::AtomAtom => {
            panic!("AtomAtom is not handled in get_pst_rovibrational_WEJ()");
        }

        CollisionType::AtomLinear => {
            let b_linear = match (fragA, fragB) {
                (ReactantType::Atom, ReactantType::Linear) => brot2,
                (ReactantType::Linear, ReactantType::Atom) => brot1,
                _ => panic!("Inconsistent AtomLinear classification"),
            };

            build_dW0_atom_and_linear(
                emax,
                dE,
                total_j,
                b_linear,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::AtomSpherical => {
            let b_spherical = match (fragA, fragB) {
                (ReactantType::Atom, ReactantType::Spherical) => brot2,
                (ReactantType::Spherical, ReactantType::Atom) => brot1,
                _ => panic!("Inconsistent AtomSpherical classification"),
            };

            build_dW0_atom_and_sphericaltop(
                emax,
                dE,
                total_j,
                b_spherical,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::AtomProlate => {
            let (a_prol, b_prol) = match (fragA, fragB) {
                (ReactantType::Atom, ReactantType::Prolate) => (arot2, brot2),
                (ReactantType::Prolate, ReactantType::Atom) => (arot1, brot1),
                _ => panic!("Inconsistent AtomProlate classification"),
            };

            build_dW0_atom_and_prolate(
                counting_mode,
                emax,
                dE,
                total_j,
                a_prol,
                b_prol,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::AtomOblate => {
            let (b_obl, c_obl) = match (fragA, fragB) {
                (ReactantType::Atom, ReactantType::Oblate) => (brot2, crot2),
                (ReactantType::Oblate, ReactantType::Atom) => (brot1, crot1),
                _ => panic!("Inconsistent AtomOblate classification"),
            };

            build_dW0_atom_and_oblate(
                counting_mode,
                emax,
                dE,
                total_j,
                b_obl,
                c_obl,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::LinearLinear => {
            build_dW0_linear_and_linear(
                counting_mode,
                emax,
                dE,
                total_j,
                brot1,
                brot2,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::LinearSpherical => {
            let (b_spherical, b_linear) = match (fragA, fragB) {
                (ReactantType::Linear, ReactantType::Spherical) => (brot2, brot1),
                (ReactantType::Spherical, ReactantType::Linear) => (brot1, brot2),
                _ => panic!("Inconsistent LinearSpherical classification"),
            };

            build_dW0_linear_and_sphericaltop(
                counting_mode,
                emax,
                dE,
                total_j,
                b_spherical,
                b_linear,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::LinearProlate => {
            let (a_prol, b_prol, b_linear) = match (fragA, fragB) {
                (ReactantType::Linear, ReactantType::Prolate) => (arot2, brot2, brot1),
                (ReactantType::Prolate, ReactantType::Linear) => (arot1, brot1, brot2),
                _ => panic!("Inconsistent LinearProlate classification"),
            };

            build_dW0_prolate_and_linear(
                counting_mode,
                emax,
                dE,
                total_j,
                a_prol,
                b_prol,
                b_linear,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::LinearOblate => {
            let (b_obl, c_obl, b_linear) = match (fragA, fragB) {
                (ReactantType::Linear, ReactantType::Oblate) => (brot2, crot2, brot1),
                (ReactantType::Oblate, ReactantType::Linear) => (brot1, crot1, brot2),
                _ => panic!("Inconsistent LinearOblate classification"),
            };

            build_dW0_oblate_and_linear(
                counting_mode,
                emax,
                dE,
                total_j,
                b_obl,
                c_obl,
                b_linear,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::SphericalSpherical => {
            build_dW0_sphericaltop_and_sphericaltop(
                counting_mode,
                emax,
                dE,
                total_j,
                brot1,
                brot2,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::SphericalProlate => {
            let (a_prol, b_prol, b_spherical) = match (fragA, fragB) {
                (ReactantType::Spherical, ReactantType::Prolate) => (arot2, brot2, brot1),
                (ReactantType::Prolate, ReactantType::Spherical) => (arot1, brot1, brot2),
                _ => panic!("Inconsistent SphericalProlate classification"),
            };

            build_dW0_prolate_and_sphericaltop(
                counting_mode,
                emax,
                dE,
                total_j,
                a_prol,
                b_prol,
                b_spherical,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::SphericalOblate => {
            let (b_obl, c_obl, b_spherical) = match (fragA, fragB) {
                (ReactantType::Spherical, ReactantType::Oblate) => (brot2, crot2, brot1),
                (ReactantType::Oblate, ReactantType::Spherical) => (brot1, crot1, brot2),
                _ => panic!("Inconsistent SphericalOblate classification"),
            };

            build_dW0_oblate_and_sphericaltop(
                counting_mode,
                emax,
                dE,
                total_j,
                b_obl,
                c_obl,
                b_spherical,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::ProlateProlate => {
            build_dW0_prolate_and_prolate(
                counting_mode,
                emax,
                dE,
                total_j,
                arot1,
                brot1,
                arot2,
                brot2,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::OblateOblate => {
            build_dW0_oblate_and_oblate(
                counting_mode,
                emax,
                dE,
                total_j,
                brot1,
                crot1,
                brot2,
                crot2,
                &mut rotational_bin_counts,
            );
        }

        CollisionType::ProlateOblate => {
            let (a_prol, b_prol, b_obl, c_obl) = match (fragA, fragB) {
                (ReactantType::Prolate, ReactantType::Oblate) => {
                    (arot1, brot1, brot2, crot2)
                }
                (ReactantType::Oblate, ReactantType::Prolate) => {
                    (arot2, brot2, brot1, crot1)
                }
                _ => panic!("Inconsistent ProlateOblate classification"),
            };

            build_dW0_prolate_and_oblate(
                counting_mode,
                emax,
                dE,
                total_j,
                a_prol,
                b_prol,
                b_obl,
                c_obl,
                &mut rotational_bin_counts,
            );
        }
    }

    let total_vibrations = nvib1 + nvib2;
    let mut vib_freq_bins = vec![0usize; total_vibrations];

    for i in 0..nvib1 {
        vib_freq_bins[i] = (vib_freqs1[i] / dE + 0.5).floor() as usize;
    }
    for i in 0..nvib2 {
        vib_freq_bins[nvib1 + i] = (vib_freqs2[i] / dE + 0.5).floor() as usize;
    }

    let rovibrational_sum_states = beyer_swinehart_counting(
        total_vibrations,
        nebin,
        &vib_freq_bins,
        &rotational_bin_counts,
    );

    if n_internal_rotors1 == 0 && n_internal_rotors2 == 0 {
        return rovibrational_sum_states;
    }

    let mut internal_rotor_density = vec![0.0; nebin + 1];
    internal_rotor(
        nebin,
        dE,
        n_internal_rotors1,
        n_internal_rotors2,
        internal_rot_consts1,
        internal_rot_consts2,
        &mut internal_rotor_density,
    );

    let mut final_sum_states = vec![0.0; nebin + 1];
    final_sum_states[0] = rovibrational_sum_states[0];

    for i in 1..=nebin {
        let mut value = 0.0;
        for l in 0..=i {
            value += rovibrational_sum_states[l] * internal_rotor_density[i - l];
        }
        final_sum_states[i] = value * dE;
    }

    final_sum_states
}

//================================================================================
// Atom + linear
//================================================================================
pub fn build_dW0_atom_and_linear(
    emax: f64,
    dE: f64,
    total_j: usize,
    b: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    for i in 1..=nebin {
        let e = energy_of_bin(i, dE);
        let em = e - dE;
        dW[i] =
            W_atom_linear(e, total_j, b) - W_atom_linear(em, total_j, b);
    }
}

fn W_atom_linear(e: f64, total_j: usize, b: f64) -> f64 {
    let jmax = ((0.25 + e / b).sqrt() - 0.5).floor() as usize;

    if total_j > jmax {
        ((jmax + 1) * (jmax + 1)) as f64
    } else {
        ((2 * total_j + 1) * (jmax + 1) - total_j * (total_j + 1)) as f64
    }
}

//================================================================================
// Atom + spherical top
//================================================================================
pub fn build_dW0_atom_and_sphericaltop(
    emax: f64,
    dE: f64,
    total_j: usize,
    b: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    for i in 1..=nebin {
        let e = energy_of_bin(i, dE);
        let em = e - dE;
        dW[i] = W_atom_sphericaltop(e, total_j, b)
            - W_atom_sphericaltop(em, total_j, b);
    }
}

fn W_atom_sphericaltop(e: f64, total_j: usize, b: f64) -> f64 {
    let jmax = ((0.25 + e / b).sqrt() - 0.5).floor() as usize;

    if total_j > jmax {
        ((jmax + 1) * (2 * jmax + 3) * (2 * jmax + 1)) as f64 / 3.0
    } else {
        ((2 * total_j + 1) * (jmax + 1) * (jmax + 1)) as f64
            - ((2 * total_j + 1) * (total_j + 1) * total_j) as f64 / 3.0
    }
}

//================================================================================
// Linear + linear
//================================================================================
pub fn build_dW0_linear_and_linear(
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    b1: f64,
    b2: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    match counting_mode {
        CountingMode::Exact => {
            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;
                dW[i] = W_linear_linear_exact(e, total_j, b1, b2)
                    - W_linear_linear_exact(em, total_j, b1, b2);
            }
        }

        CountingMode::Interpol => {
            let mut switched_to_interpolation = false;

            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;

                if !switched_to_interpolation {
                    let exact = W_linear_linear_exact(e, total_j, b1, b2);
                    let approx = W_linear_linear_interpol(e, total_j, b1, b2);

                    if relative_error_percent(exact, approx) <= 2.0 {
                        switched_to_interpolation = true;
                    }

                    dW[i] = exact - W_linear_linear_exact(em, total_j, b1, b2);
                } else {
                    dW[i] = W_linear_linear_interpol(e, total_j, b1, b2)
                        - W_linear_linear_interpol(em, total_j, b1, b2);
                }
            }
        }
    }
}


#![allow(non_snake_case)]

// ============================================================================
// Helpers
// ============================================================================

#[inline]
fn number_of_bins(emax: f64, dE: f64) -> usize {
    (emax / dE + 0.5).floor() as usize
}

#[inline]
fn energy_of_bin(i: usize, dE: f64) -> f64 {
    (i as f64) * dE
}

// ============================================================================
// Atom + prolate
// ============================================================================

pub fn build_dW0_atom_and_prolate(
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    a_prol: f64,
    b_prol: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    match counting_mode {
        CountingMode::Exact | CountingMode::Interpol => {
            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;
                dW[i] = W_atom_prolate_exact(e, total_j, a_prol, b_prol)
                    - W_atom_prolate_exact(em, total_j, a_prol, b_prol);
            }
        }
    }
}

// ============================================================================
// Atom + oblate
// ============================================================================

pub fn build_dW0_atom_and_oblate(
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    b_obl: f64,
    c_obl: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    match counting_mode {
        CountingMode::Exact | CountingMode::Interpol => {
            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;
                dW[i] = W_atom_oblate_exact(e, total_j, b_obl, c_obl)
                    - W_atom_oblate_exact(em, total_j, b_obl, c_obl);
            }
        }
    }
}

// ============================================================================
// Prolate + linear
// ============================================================================

pub fn build_dW0_prolate_and_linear(
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    a_prol: f64,
    b_prol: f64,
    b_lin: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    match counting_mode {
        CountingMode::Exact | CountingMode::Interpol => {
            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;
                dW[i] = W_prolate_linear_exact(e, total_j, a_prol, b_prol, b_lin)
                    - W_prolate_linear_exact(em, total_j, a_prol, b_prol, b_lin);
            }
        }
    }
}

// ============================================================================
// Prolate + spherical top
// ============================================================================

pub fn build_dW0_prolate_and_sphericaltop(
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    a_prol: f64,
    b_prol: f64,
    b_sph: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    match counting_mode {
        CountingMode::Exact | CountingMode::Interpol => {
            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;
                dW[i] = W_prolate_sphericaltop_exact(e, total_j, a_prol, b_prol, b_sph)
                    - W_prolate_sphericaltop_exact(em, total_j, a_prol, b_prol, b_sph);
            }
        }
    }
}

// ============================================================================
// Prolate + prolate
// ============================================================================

pub fn build_dW0_prolate_and_prolate(
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    a1: f64,
    b1: f64,
    a2: f64,
    b2: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    match counting_mode {
        CountingMode::Exact | CountingMode::Interpol => {
            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;
                dW[i] = W_prolate_prolate_exact(e, total_j, a1, b1, a2, b2)
                    - W_prolate_prolate_exact(em, total_j, a1, b1, a2, b2);
            }
        }
    }
}

// ============================================================================
// Oblate + linear
// ============================================================================

pub fn build_dW0_oblate_and_linear(
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    b_obl: f64,
    c_obl: f64,
    b_lin: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    match counting_mode {
        CountingMode::Exact | CountingMode::Interpol => {
            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;
                dW[i] = W_oblate_linear_exact(e, total_j, b_obl, c_obl, b_lin)
                    - W_oblate_linear_exact(em, total_j, b_obl, c_obl, b_lin);
            }
        }
    }
}

// ============================================================================
// Oblate + spherical top
// ============================================================================

pub fn build_dW0_oblate_and_sphericaltop(
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    b_obl: f64,
    c_obl: f64,
    b_sph: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    match counting_mode {
        CountingMode::Exact | CountingMode::Interpol => {
            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;
                dW[i] = W_oblate_sphericaltop_exact(e, total_j, b_obl, c_obl, b_sph)
                    - W_oblate_sphericaltop_exact(em, total_j, b_obl, c_obl, b_sph);
            }
        }
    }
}

// ============================================================================
// Oblate + oblate
// ============================================================================

pub fn build_dW0_oblate_and_oblate(
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    b1: f64,
    c1: f64,
    b2: f64,
    c2: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    match counting_mode {
        CountingMode::Exact | CountingMode::Interpol => {
            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;
                dW[i] = W_oblate_oblate_exact(e, total_j, b1, c1, b2, c2)
                    - W_oblate_oblate_exact(em, total_j, b1, c1, b2, c2);
            }
        }
    }
}

// ============================================================================
// Prolate + oblate
// ============================================================================

pub fn build_dW0_prolate_and_oblate(
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    a_prol: f64,
    b_prol: f64,
    b_obl: f64,
    c_obl: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    match counting_mode {
        CountingMode::Exact | CountingMode::Interpol => {
            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;
                dW[i] = W_prolate_oblate_exact(e, total_j, a_prol, b_prol, b_obl, c_obl)
                    - W_prolate_oblate_exact(em, total_j, a_prol, b_prol, b_obl, c_obl);
            }
        }
    }
}




fn W_linear_linear_exact(e: f64, total_j: usize, b1: f64, b2: f64) -> f64 {
    let mut total_states = 0.0;

    let j1max = ((0.25 + e / b1).sqrt() - 0.5).floor() as usize;

    for j1 in 0..=j1max {
        let e_remaining = e - b1 * (j1 * (j1 + 1)) as f64;
        if e_remaining < 0.0 {
            continue;
        }

        let j2max = ((0.25 + e_remaining / b2).sqrt() - 0.5).floor() as usize;

        for j2 in 0..=j2max {
            let j_coupled_min = j1.abs_diff(j2);
            let j_coupled_max = j1 + j2;

            for j_coupled in j_coupled_min..=j_coupled_max {
                let l_min = total_j.abs_diff(j_coupled);
                let l_max = total_j + j_coupled;

                for _l in l_min..=l_max {
                    total_states += 1.0;
                }
            }
        }
    }

    total_states
}

fn W_linear_linear_j0(e: f64, b1: f64, b2: f64) -> f64 {
    2.0 * e * e.sqrt() * (b1.sqrt() + b2.sqrt() - (b1 + b2).sqrt()) / (3.0 * b1 * b2)
}

fn W_linear_linear_high_j(e: f64, b1: f64, b2: f64) -> f64 {
    e * e / (2.0 * b1 * b2)
}

fn W_linear_linear_interpol(e: f64, total_j: usize, b1: f64, b2: f64) -> f64 {
    let high_j_limit = W_linear_linear_high_j(e, b1, b2);
    let reduced_x =
        ((2 * total_j + 1) as f64) * W_linear_linear_j0(e, b1, b2) / high_j_limit;
    high_j_limit * interpolate_reduced_W(reduced_x)
}

//================================================================================
// Linear + spherical top
//================================================================================
pub fn build_dW0_linear_and_sphericaltop(
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    bst: f64,
    bli: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    match counting_mode {
        CountingMode::Exact => {
            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;
                dW[i] = W_linear_sphericaltop_exact(e, total_j, bst, bli)
                    - W_linear_sphericaltop_exact(em, total_j, bst, bli);
            }
        }

        CountingMode::Interpol => {
            let mut switched_to_interpolation = false;

            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;

                if !switched_to_interpolation {
                    let exact = W_linear_sphericaltop_exact(e, total_j, bst, bli);
                    let approx = W_linear_sphericaltop_interpol(e, total_j, bst, bli);

                    if relative_error_percent(exact, approx) <= 2.0 {
                        switched_to_interpolation = true;
                    }

                    dW[i] = exact - W_linear_sphericaltop_exact(em, total_j, bst, bli);
                } else {
                    dW[i] = W_linear_sphericaltop_interpol(e, total_j, bst, bli)
                        - W_linear_sphericaltop_interpol(em, total_j, bst, bli);
                }
            }
        }
    }
}

fn W_linear_sphericaltop_exact(
    e: f64,
    total_j: usize,
    bst: f64,
    bli: f64,
) -> f64 {
    let mut total_states = 0.0;

    let j1max = ((0.25 + e / bst).sqrt() - 0.5).floor() as usize;

    for j1 in 0..=j1max {
        let e_remaining = e - bst * (j1 * (j1 + 1)) as f64;
        if e_remaining < 0.0 {
            continue;
        }

        let j2max = ((0.25 + e_remaining / bli).sqrt() - 0.5).floor() as usize;
        let degeneracy_j1 = (2 * j1 + 1) as f64;

        for j2 in 0..=j2max {
            let j_coupled_min = j1.abs_diff(j2);
            let j_coupled_max = j1 + j2;

            for j_coupled in j_coupled_min..=j_coupled_max {
                let l_min = total_j.abs_diff(j_coupled);
                let l_max = total_j + j_coupled;

                for _l in l_min..=l_max {
                    total_states += degeneracy_j1;
                }
            }
        }
    }

    total_states
}

fn W_linear_sphericaltop_j0(e: f64, bst: f64, bli: f64) -> f64 {
    e * e * (bst / (bst + bli)).sqrt().asin()
        / (2.0 * bst * bst.sqrt() * bli.sqrt())
}

fn W_linear_sphericaltop_high_j(e: f64, bst: f64, bli: f64) -> f64 {
    8.0 * e * e * e.sqrt() / (15.0 * bst * bst.sqrt() * bli)
}

fn W_linear_sphericaltop_interpol(
    e: f64,
    total_j: usize,
    bst: f64,
    bli: f64,
) -> f64 {
    let high_j_limit = W_linear_sphericaltop_high_j(e, bst, bli);
    let reduced_x = ((2 * total_j + 1) as f64)
        * W_linear_sphericaltop_j0(e, bst, bli)
        / high_j_limit;
    high_j_limit * interpolate_reduced_W(reduced_x)
}

//================================================================================
// Spherical top + spherical top
//================================================================================
pub fn build_dW0_sphericaltop_and_sphericaltop(
    counting_mode: CountingMode,
    emax: f64,
    dE: f64,
    total_j: usize,
    bst1: f64,
    bst2: f64,
    dW: &mut [f64],
) {
    let nebin = number_of_bins(emax, dE);
    dW[0] = 1.0;

    match counting_mode {
        CountingMode::Exact => {
            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;
                dW[i] = W_sphericaltop_sphericaltop_exact(e, total_j, bst1, bst2)
                    - W_sphericaltop_sphericaltop_exact(em, total_j, bst1, bst2);
            }
        }

        CountingMode::Interpol => {
            let mut switched_to_interpolation = false;

            for i in 1..=nebin {
                let e = energy_of_bin(i, dE);
                let em = e - dE;

                if !switched_to_interpolation {
                    let exact = W_sphericaltop_sphericaltop_exact(e, total_j, bst1, bst2);
                    let approx =
                        W_sphericaltop_sphericaltop_interpol(e, total_j, bst1, bst2);

                    if relative_error_percent(exact, approx) <= 2.0 {
                        switched_to_interpolation = true;
                    }

                    dW[i] = exact
                        - W_sphericaltop_sphericaltop_exact(em, total_j, bst1, bst2);
                } else {
                    dW[i] = W_sphericaltop_sphericaltop_interpol(e, total_j, bst1, bst2)
                        - W_sphericaltop_sphericaltop_interpol(em, total_j, bst1, bst2);
                }
            }
        }
    }
}

fn W_sphericaltop_sphericaltop_exact(
    e: f64,
    total_j: usize,
    bst1: f64,
    bst2: f64,
) -> f64 {
    let mut total_states = 0.0;

    let j1max = ((0.25 + e / bst1).sqrt() - 0.5).floor() as usize;

    for j1 in 0..=j1max {
        let e_remaining = e - bst1 * (j1 * (j1 + 1)) as f64;
        if e_remaining < 0.0 {
            continue;
        }

        let j2max = ((0.25 + e_remaining / bst2).sqrt() - 0.5).floor() as usize;
        let degeneracy_j1 = (2 * j1 + 1) as f64;

        for j2 in 0..=j2max {
            let degeneracy = (2 * j2 + 1) as f64 * degeneracy_j1;
            let j_coupled_min = j1.abs_diff(j2);
            let j_coupled_max = j1 + j2;

            for j_coupled in j_coupled_min..=j_coupled_max {
                let l_min = total_j.abs_diff(j_coupled);
                let l_max = total_j + j_coupled;

                for _l in l_min..=l_max {
                    total_states += degeneracy;
                }
            }
        }
    }

    total_states
}



#![allow(non_snake_case)]

// ============================================================================
// Helpers
// ============================================================================

#[inline]
fn k_degeneracy(k: usize) -> f64 {
    if k == 0 { 1.0 } else { 2.0 }
}

#[inline]
fn jmax_rigid_rotor(e: f64, b: f64) -> usize {
    if e <= 0.0 || b <= 0.0 {
        0
    } else {
        ((0.25 + e / b).sqrt() - 0.5).floor() as usize
    }
}

// Prolate symmetric top: E(j,K) = B j(j+1) + (A-B) K^2,  A > B
#[inline]
fn prolate_energy(j: usize, k: usize, a: f64, b: f64) -> f64 {
    b * (j * (j + 1)) as f64 + (a - b) * (k * k) as f64
}

// Oblate symmetric top: E(j,K) = B j(j+1) - (B-C) K^2,  B > C
#[inline]
fn oblate_energy(j: usize, k: usize, b: f64, c: f64) -> f64 {
    b * (j * (j + 1)) as f64 - (b - c) * (k * k) as f64
}

// For oblate tops the minimum energy at fixed j occurs at K = j:
// Emin(j) = B j(j+1) - (B-C) j^2 = B j + C j^2
#[inline]
fn jmax_oblate(e: f64, b: f64, c: f64) -> usize {
    if e <= 0.0 || b <= 0.0 || c <= 0.0 {
        0
    } else {
        ((-b + (b * b + 4.0 * c * e).sqrt()) / (2.0 * c)).floor() as usize
    }
}

// Smallest allowed nonnegative K for an oblate top at fixed (j,e):
// B j(j+1) - (B-C) K^2 <= e
// K^2 >= [B j(j+1) - e] / (B-C)
#[inline]
fn oblate_kmin(j: usize, e: f64, b: f64, c: f64) -> usize {
    let numerator = b * (j * (j + 1)) as f64 - e;

    if numerator <= 0.0 {
        return 0;
    }

    let denom = b - c;
    if denom <= 0.0 {
        return j + 1;
    }

    (numerator / denom).sqrt().ceil() as usize
}

// ============================================================================
// Atom + prolate
// ============================================================================

pub fn W_atom_prolate_exact(
    e: f64,
    total_j: usize,
    a_prol: f64,
    b_prol: f64,
) -> f64 {
    let mut total_states = 0.0;

    let jmax = jmax_rigid_rotor(e, b_prol);

    for j in 0..=jmax {
        for k in 0..=j {
            let e_rot = prolate_energy(j, k, a_prol, b_prol);
            if e_rot > e {
                continue;
            }

            let weight = k_degeneracy(k);
            let l_min = total_j.abs_diff(j);
            let l_max = total_j + j;

            for _l in l_min..=l_max {
                total_states += weight;
            }
        }
    }

    total_states
}

// ============================================================================
// Atom + oblate
// ============================================================================

pub fn W_atom_oblate_exact(
    e: f64,
    total_j: usize,
    b_obl: f64,
    c_obl: f64,
) -> f64 {
    let mut total_states = 0.0;

    let jmax = jmax_oblate(e, b_obl, c_obl);

    for j in 0..=jmax {
        let kmin = oblate_kmin(j, e, b_obl, c_obl);
        if kmin > j {
            continue;
        }

        for k in kmin..=j {
            let e_rot = oblate_energy(j, k, b_obl, c_obl);
            if e_rot > e {
                continue;
            }

            let weight = k_degeneracy(k);
            let l_min = total_j.abs_diff(j);
            let l_max = total_j + j;

            for _l in l_min..=l_max {
                total_states += weight;
            }
        }
    }

    total_states
}

// ============================================================================
// Prolate + linear
// ============================================================================

pub fn W_prolate_linear_exact(
    e: f64,
    total_j: usize,
    a_prol: f64,
    b_prol: f64,
    b_lin: f64,
) -> f64 {
    let mut total_states = 0.0;

    let j1max = jmax_rigid_rotor(e, b_prol);

    for j1 in 0..=j1max {
        for k1 in 0..=j1 {
            let e1 = prolate_energy(j1, k1, a_prol, b_prol);
            if e1 > e {
                continue;
            }

            let weight1 = k_degeneracy(k1);
            let e_rem = e - e1;
            let j2max = jmax_rigid_rotor(e_rem, b_lin);

            for j2 in 0..=j2max {
                let j_coupled_min = j1.abs_diff(j2);
                let j_coupled_max = j1 + j2;

                for j_coupled in j_coupled_min..=j_coupled_max {
                    let l_min = total_j.abs_diff(j_coupled);
                    let l_max = total_j + j_coupled;

                    for _l in l_min..=l_max {
                        total_states += weight1;
                    }
                }
            }
        }
    }

    total_states
}

// ============================================================================
// Prolate + spherical top
// ============================================================================

pub fn W_prolate_sphericaltop_exact(
    e: f64,
    total_j: usize,
    a_prol: f64,
    b_prol: f64,
    b_sph: f64,
) -> f64 {
    let mut total_states = 0.0;

    let j1max = jmax_rigid_rotor(e, b_prol);

    for j1 in 0..=j1max {
        for k1 in 0..=j1 {
            let e1 = prolate_energy(j1, k1, a_prol, b_prol);
            if e1 > e {
                continue;
            }

            let weight1 = k_degeneracy(k1);
            let e_rem = e - e1;
            let j2max = jmax_rigid_rotor(e_rem, b_sph);

            for j2 in 0..=j2max {
                let weight2 = (2 * j2 + 1) as f64;
                let j_coupled_min = j1.abs_diff(j2);
                let j_coupled_max = j1 + j2;

                for j_coupled in j_coupled_min..=j_coupled_max {
                    let l_min = total_j.abs_diff(j_coupled);
                    let l_max = total_j + j_coupled;

                    for _l in l_min..=l_max {
                        total_states += weight1 * weight2;
                    }
                }
            }
        }
    }

    total_states
}

// ============================================================================
// Prolate + prolate
// ============================================================================

pub fn W_prolate_prolate_exact(
    e: f64,
    total_j: usize,
    a1: f64,
    b1: f64,
    a2: f64,
    b2: f64,
) -> f64 {
    let mut total_states = 0.0;

    let j1max = jmax_rigid_rotor(e, b1);

    for j1 in 0..=j1max {
        for k1 in 0..=j1 {
            let e1 = prolate_energy(j1, k1, a1, b1);
            if e1 > e {
                continue;
            }

            let weight1 = k_degeneracy(k1);
            let e_rem1 = e - e1;
            let j2max = jmax_rigid_rotor(e_rem1, b2);

            for j2 in 0..=j2max {
                for k2 in 0..=j2 {
                    let e2 = prolate_energy(j2, k2, a2, b2);
                    if e2 > e_rem1 {
                        continue;
                    }

                    let weight2 = k_degeneracy(k2);
                    let j_coupled_min = j1.abs_diff(j2);
                    let j_coupled_max = j1 + j2;

                    for j_coupled in j_coupled_min..=j_coupled_max {
                        let l_min = total_j.abs_diff(j_coupled);
                        let l_max = total_j + j_coupled;

                        for _l in l_min..=l_max {
                            total_states += weight1 * weight2;
                        }
                    }
                }
            }
        }
    }

    total_states
}

// ============================================================================
// Oblate + linear
// ============================================================================

pub fn W_oblate_linear_exact(
    e: f64,
    total_j: usize,
    b_obl: f64,
    c_obl: f64,
    b_lin: f64,
) -> f64 {
    let mut total_states = 0.0;

    let j1max = jmax_oblate(e, b_obl, c_obl);

    for j1 in 0..=j1max {
        let kmin1 = oblate_kmin(j1, e, b_obl, c_obl);
        if kmin1 > j1 {
            continue;
        }

        for k1 in kmin1..=j1 {
            let e1 = oblate_energy(j1, k1, b_obl, c_obl);
            if e1 > e {
                continue;
            }

            let weight1 = k_degeneracy(k1);
            let e_rem = e - e1;
            let j2max = jmax_rigid_rotor(e_rem, b_lin);

            for j2 in 0..=j2max {
                let j_coupled_min = j1.abs_diff(j2);
                let j_coupled_max = j1 + j2;

                for j_coupled in j_coupled_min..=j_coupled_max {
                    let l_min = total_j.abs_diff(j_coupled);
                    let l_max = total_j + j_coupled;

                    for _l in l_min..=l_max {
                        total_states += weight1;
                    }
                }
            }
        }
    }

    total_states
}

// ============================================================================
// Oblate + spherical top
// ============================================================================

pub fn W_oblate_sphericaltop_exact(
    e: f64,
    total_j: usize,
    b_obl: f64,
    c_obl: f64,
    b_sph: f64,
) -> f64 {
    let mut total_states = 0.0;

    let j1max = jmax_oblate(e, b_obl, c_obl);

    for j1 in 0..=j1max {
        let kmin1 = oblate_kmin(j1, e, b_obl, c_obl);
        if kmin1 > j1 {
            continue;
        }

        for k1 in kmin1..=j1 {
            let e1 = oblate_energy(j1, k1, b_obl, c_obl);
            if e1 > e {
                continue;
            }

            let weight1 = k_degeneracy(k1);
            let e_rem = e - e1;
            let j2max = jmax_rigid_rotor(e_rem, b_sph);

            for j2 in 0..=j2max {
                let weight2 = (2 * j2 + 1) as f64;
                let j_coupled_min = j1.abs_diff(j2);
                let j_coupled_max = j1 + j2;

                for j_coupled in j_coupled_min..=j_coupled_max {
                    let l_min = total_j.abs_diff(j_coupled);
                    let l_max = total_j + j_coupled;

                    for _l in l_min..=l_max {
                        total_states += weight1 * weight2;
                    }
                }
            }
        }
    }

    total_states
}

// ============================================================================
// Oblate + oblate
// ============================================================================

pub fn W_oblate_oblate_exact(
    e: f64,
    total_j: usize,
    b1: f64,
    c1: f64,
    b2: f64,
    c2: f64,
) -> f64 {
    let mut total_states = 0.0;

    let j1max = jmax_oblate(e, b1, c1);

    for j1 in 0..=j1max {
        let kmin1 = oblate_kmin(j1, e, b1, c1);
        if kmin1 > j1 {
            continue;
        }

        for k1 in kmin1..=j1 {
            let e1 = oblate_energy(j1, k1, b1, c1);
            if e1 > e {
                continue;
            }

            let weight1 = k_degeneracy(k1);
            let e_rem1 = e - e1;
            let j2max = jmax_oblate(e_rem1, b2, c2);

            for j2 in 0..=j2max {
                let kmin2 = oblate_kmin(j2, e_rem1, b2, c2);
                if kmin2 > j2 {
                    continue;
                }

                for k2 in kmin2..=j2 {
                    let e2 = oblate_energy(j2, k2, b2, c2);
                    if e2 > e_rem1 {
                        continue;
                    }

                    let weight2 = k_degeneracy(k2);
                    let j_coupled_min = j1.abs_diff(j2);
                    let j_coupled_max = j1 + j2;

                    for j_coupled in j_coupled_min..=j_coupled_max {
                        let l_min = total_j.abs_diff(j_coupled);
                        let l_max = total_j + j_coupled;

                        for _l in l_min..=l_max {
                            total_states += weight1 * weight2;
                        }
                    }
                }
            }
        }
    }

    total_states
}

// ============================================================================
// Prolate + oblate
// ============================================================================

pub fn W_prolate_oblate_exact(
    e: f64,
    total_j: usize,
    a_prol: f64,
    b_prol: f64,
    b_obl: f64,
    c_obl: f64,
) -> f64 {
    let mut total_states = 0.0;

    let j1max = jmax_rigid_rotor(e, b_prol);

    for j1 in 0..=j1max {
        for k1 in 0..=j1 {
            let e1 = prolate_energy(j1, k1, a_prol, b_prol);
            if e1 > e {
                continue;
            }

            let weight1 = k_degeneracy(k1);
            let e_rem1 = e - e1;
            let j2max = jmax_oblate(e_rem1, b_obl, c_obl);

            for j2 in 0..=j2max {
                let kmin2 = oblate_kmin(j2, e_rem1, b_obl, c_obl);
                if kmin2 > j2 {
                    continue;
                }

                for k2 in kmin2..=j2 {
                    let e2 = oblate_energy(j2, k2, b_obl, c_obl);
                    if e2 > e_rem1 {
                        continue;
                    }

                    let weight2 = k_degeneracy(k2);
                    let j_coupled_min = j1.abs_diff(j2);
                    let j_coupled_max = j1 + j2;

                    for j_coupled in j_coupled_min..=j_coupled_max {
                        let l_min = total_j.abs_diff(j_coupled);
                        let l_max = total_j + j_coupled;

                        for _l in l_min..=l_max {
                            total_states += weight1 * weight2;
                        }
                    }
                }
            }
        }
    }

    total_states
}

fn W_sphericaltop_sphericaltop_j0(e: f64, bst1: f64, bst2: f64) -> f64 {
    8.0 * e * e * e.sqrt() / (15.0 * bst1 * bst2 * (bst1 + bst2).sqrt())
}

fn W_sphericaltop_sphericaltop_high_j(e: f64, bst1: f64, bst2: f64) -> f64 {
    3.1416 * e * e * e / (6.0 * bst1 * bst2 * (bst1 * bst2).sqrt())
}

fn W_sphericaltop_sphericaltop_interpol(
    e: f64,
    total_j: usize,
    bst1: f64,
    bst2: f64,
) -> f64 {
    let high_j_limit = W_sphericaltop_sphericaltop_high_j(e, bst1, bst2);
    let reduced_x = ((2 * total_j + 1) as f64)
        * W_sphericaltop_sphericaltop_j0(e, bst1, bst2)
        / high_j_limit;
    high_j_limit * interpolate_reduced_W(reduced_x)
}

//================================================================================
// Shared interpolation helper
//================================================================================
pub fn interpolate_reduced_W(x: f64) -> f64 {
    x.tanh() * (0.08 * x.sqrt() * (-1.1 * (x - 1.0) * (x - 1.0)).exp() + 1.0)
}

//================================================================================
// Internal rotor density
//================================================================================
pub fn internal_rotor(
    nebin: usize,
    dE: f64,
    n_internal_rotors1: usize,
    n_internal_rotors2: usize,
    internal_rot_consts1: &[f64],
    internal_rot_consts2: &[f64],
    internal_rotor_density: &mut [f64],
) {
    let total_internal_rotors = (n_internal_rotors1 + n_internal_rotors2) as f64;

    internal_rotor_density[0] = 0.0;

    let mut prod_Brot_sqrt = 1.0;

    for i in 0..n_internal_rotors1 {
        prod_Brot_sqrt *= internal_rot_consts1[i].sqrt();
    }

    for i in 0..n_internal_rotors2 {
        prod_Brot_sqrt *= internal_rot_consts2[i].sqrt();
    }

    let exponent = total_internal_rotors / 2.0 - 1.0;
    let gamma = gamma_func(total_internal_rotors / 2.0);

    for i in 1..=nebin {
        let e = (i as f64) * dE;
        internal_rotor_density[i] =
            1.7725_f64.powf(total_internal_rotors) / (gamma * prod_Brot_sqrt) * e.powf(exponent);
    }
}
