use crate::numeric::lanczos_gamma::gamma_func;

use super::hindered_rotor::fqhind;
use super::threshold_energy::morse_threshold_energy;
use super::types::{ProductRotorCase, ReactantRotorCase};

#[derive(Debug, Clone)]
pub struct SacmThermalInput {
    pub temperatures: Vec<f64>,
    pub max_temperature: f64,
    pub rotor_case_reactant: ReactantRotorCase,
    pub rotor_case_products: ProductRotorCase,
    pub alpha_over_beta: f64,
    pub anisotropy_scale: f64,
    pub transition_symmetry: f64,
    pub reactant_symmetry: f64,
    pub frag1_symmetry: f64,
    pub frag2_symmetry: f64,
    pub dissociation_energy: f64,
    pub product_zpe: f64,
    pub zpe_shift: f64,
    pub enthalpy_0k: f64,
    pub base_rot_const: f64,
    pub centrifugal_a1: f64,
    pub centrifugal_a2: f64,
    pub external_rot_a: f64,
    pub external_rot_b: f64,
    pub external_rot_c: f64,
    pub frag_rot_b1: f64,
    pub frag_rot_b2: f64,
    pub red_mass: f64,
    pub vib_reactant: Vec<f64>,
    pub vib_frag1: Vec<f64>,
    pub vib_frag2: Vec<f64>,
    pub internal_rot_reactant: Vec<f64>,
    pub internal_rot_frag1: Vec<f64>,
    pub internal_rot_frag2: Vec<f64>,
    pub hindered_barrier_reactant: f64,
    pub hindered_barrier_frag1: f64,
    pub hindered_barrier_frag2: f64,
    pub hindered_count_reactant: f64,
    pub hindered_count_frag1: f64,
    pub hindered_count_frag2: f64,
    pub electronic_q_reactant: f64,
    pub electronic_q_frag1: f64,
    pub electronic_q_frag2: f64,
    pub qstar_eps: Vec<f64>,
    pub qstar_x: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct SacmThermalRate {
    pub temperature: f64,
    pub dissociation_rate: f64,
    pub recombination_rate: f64,
    pub equilibrium_constant: f64,
}

#[derive(Debug, Clone)]
pub struct SacmThermalDebugRate {
    pub temperature: f64,
    pub kt: f64,
    pub qcent: f64,
    pub enull_0: f64,
    pub enull: Vec<f64>,
    pub qera: f64,
    pub qerb: f64,
    pub qerc: f64,
    pub qira: f64,
    pub qirb: f64,
    pub qirc: f64,
    pub qviba: f64,
    pub qvibb: f64,
    pub qvibc: f64,
    pub qela: f64,
    pub qelb: f64,
    pub qelc: f64,
    pub qvra: f64,
    pub qvrb: f64,
    pub qvrc: f64,
    pub qstarp: f64,
    pub qprod: f64,
    pub dissociation_rate: f64,
    pub recombination_rate: f64,
    pub equilibrium_constant: f64,
}

/// Compute high-pressure thermal rates using the QCENT/QSTAR pathway.
pub fn compute_thermal_rates(input: &SacmThermalInput) -> Vec<SacmThermalRate> {
    let debug_rates = compute_thermal_rates_debug(input);
    debug_rates
        .into_iter()
        .map(|rate| SacmThermalRate {
            temperature: rate.temperature,
            dissociation_rate: rate.dissociation_rate,
            recombination_rate: rate.recombination_rate,
            equilibrium_constant: rate.equilibrium_constant,
        })
        .collect()
}

/// Compute high-pressure thermal rates and return intermediate terms for debugging.
pub fn compute_thermal_rates_debug(input: &SacmThermalInput) -> Vec<SacmThermalDebugRate> {
    let mut rates = Vec::new();
    let enull = compute_enull_j(
        input.dissociation_energy,
        input.product_zpe,
        input.zpe_shift,
        input.alpha_over_beta,
        input.base_rot_const,
        input.centrifugal_a1,
        input.centrifugal_a2,
        input.max_temperature,
    );
    if enull.is_empty() {
        return rates;
    }

    let enull_0 = enull[0];
    let delta_e0 = enull_0 - input.dissociation_energy - input.product_zpe;
    let famc = compute_famc(
        input.rotor_case_reactant,
        input.rotor_case_products,
        input.anisotropy_scale,
        input.alpha_over_beta,
    );

    for &temp in &input.temperatures {
        if temp <= 0.0 {
            continue;
        }
        let kt = 0.69503 * temp;
        let qcent = compute_qcent(&enull, enull_0, kt);
        let qviba = vibrational_partition(&input.vib_reactant, kt);
        let qvibb = vibrational_partition(&input.vib_frag1, kt);
        let qvibc = vibrational_partition(&input.vib_frag2, kt);
        let qira = internal_rotor_partition(
            temp,
            &input.internal_rot_reactant,
            input.hindered_barrier_reactant,
            input.hindered_count_reactant,
        );
        let qirb = internal_rotor_partition(
            temp,
            &input.internal_rot_frag1,
            input.hindered_barrier_frag1,
            input.hindered_count_frag1,
        );
        let qirc = internal_rotor_partition(
            temp,
            &input.internal_rot_frag2,
            input.hindered_barrier_frag2,
            input.hindered_count_frag2,
        );

        let qera = external_rotor_partition_reactant(
            input.rotor_case_reactant,
            temp,
            input.external_rot_a,
            input.external_rot_b,
            input.external_rot_c,
        );
        let (qerb, qerc) = external_rotor_partition_products(
            input.rotor_case_products,
            temp,
            input.frag_rot_b1,
            input.frag_rot_b2,
        );

        let qvra = qviba * qira * qera / input.reactant_symmetry;
        let qvrb = qvibb * qirb * qerb / input.frag1_symmetry;
        let qvrc = qvibc * qirc * qerc / input.frag2_symmetry;

        let mut qstarp = qstar_product(&input.qstar_eps, &input.qstar_x, kt);
        if input.internal_rot_frag1.len() + input.internal_rot_frag2.len() > 0 {
            qstarp *=
                (input.internal_rot_frag1.len() + input.internal_rot_frag2.len()) as f64 * 2.0;
        }

        let qprod = qcent * famc / input.transition_symmetry * qstarp * (-delta_e0 / kt).exp();
        let k_diss = kt / 3.3356e-11 * qprod / qvra * (-input.enthalpy_0k / kt).exp();
        let k_recomb = kt / 3.3356e-11 * 3.0835e-21
            / (input.red_mass * kt * (input.red_mass * kt).sqrt())
            * input.electronic_q_reactant
            / (input.electronic_q_frag1 * input.electronic_q_frag2)
            * qprod
            / (qvrb * qvrc);
        let k_eq = if k_recomb > 0.0 {
            k_diss / k_recomb
        } else {
            0.0
        };

        rates.push(SacmThermalDebugRate {
            temperature: temp,
            kt,
            qcent,
            enull_0,
            enull: enull.clone(),
            qera,
            qerb,
            qerc,
            qira,
            qirb,
            qirc,
            qviba,
            qvibb,
            qvibc,
            qela: input.electronic_q_reactant,
            qelb: input.electronic_q_frag1,
            qelc: input.electronic_q_frag2,
            qvra,
            qvrb,
            qvrc,
            qstarp,
            qprod,
            dissociation_rate: k_diss,
            recombination_rate: k_recomb,
            equilibrium_constant: k_eq,
        });
    }

    rates
}

/// Precompute E0(J) values to build QCENT.
fn compute_enull_j(
    dmorse: f64,
    ezp: f64,
    zpe_shift: f64,
    alpha_over_beta: f64,
    be: f64,
    a1: f64,
    a2: f64,
    max_temp: f64,
) -> Vec<f64> {
    let mut enull = Vec::new();
    let mut j = 0;
    let enull_0 = morse_threshold_energy(dmorse, ezp, zpe_shift, alpha_over_beta, be, a1, a2, 0);
    if enull_0 < 0.0 {
        return enull;
    }
    enull.push(enull_0);

    let kb = 0.69503;
    loop {
        j += 1;
        let e = morse_threshold_energy(
            dmorse,
            ezp,
            zpe_shift,
            alpha_over_beta,
            be,
            a1,
            a2,
            j as i32,
        );
        if e < 0.0 {
            break;
        }
        let term = (2 * j + 1) as f64 * (-(e - enull_0) / (kb * max_temp)).exp();
        enull.push(e);
        if term < 1.0e-2 {
            break;
        }
        if (e - enull_0) / (kb * max_temp) > 87.0 {
            break;
        }
    }
    enull
}

/// Angular-momentum coupling factor for the thermal partition (FAMC).
fn compute_famc(
    reactant_rotor_case: ReactantRotorCase,
    product_rotor_case: ProductRotorCase,
    c3: f64,
    alpha_over_beta: f64,
) -> f64 {
    let famec = [1.0, 2.0, 2.0, 2.0];
    let famuc = [1.2732, 2.5465, 1.6211, 3.2423, 6.4846];
    let a = famec
        .get(reactant_rotor_case.famc_index())
        .copied()
        .unwrap_or(1.0);
    let b = famuc
        .get(product_rotor_case.famc_index())
        .copied()
        .unwrap_or(1.0);
    b + (a - b) * (-(c3 * alpha_over_beta)).exp()
}

/// Centrifugal partition function QCENT.
fn compute_qcent(enull: &[f64], enull_0: f64, kt: f64) -> f64 {
    let mut qcent = 0.0;
    for (j, &e) in enull.iter().enumerate() {
        if (e - enull_0) / kt > 87.0 {
            break;
        }
        qcent += (2 * j + 1) as f64 * (-(e - enull_0) / kt).exp();
    }
    qcent
}

/// Harmonic vibrational partition function.
fn vibrational_partition(freqs: &[f64], kt: f64) -> f64 {
    let mut q = 1.0;
    for &w in freqs {
        q /= 1.0 - (-w / kt).exp();
    }
    q
}

/// Internal-rotor partition with hindered-rotor correction.
fn internal_rotor_partition(temp: f64, brot: &[f64], v0: f64, n: f64) -> f64 {
    if brot.is_empty() {
        return 1.0;
    }
    let mut q = 1.478_f64.powi(brot.len() as i32) * temp.powf(brot.len() as f64 / 2.0);
    for &b in brot {
        q /= b.sqrt();
    }
    let kt = 0.69503 * temp;
    q * fqhind(kt, v0, q, n)
}

/// External-rotor partition for the reactant (linear/spherical/symmetric top).
fn external_rotor_partition_reactant(
    rotor_case: ReactantRotorCase,
    temp: f64,
    rot_a: f64,
    rot_b: f64,
    rot_c: f64,
) -> f64 {
    match rotor_case {
        ReactantRotorCase::Linear => 0.695 * temp / rot_b,
        ReactantRotorCase::SphericalTop => 1.027 * temp / rot_b * (temp / rot_b).sqrt(),
        _ => 1.027 * (temp * temp * temp / (rot_a * rot_b * rot_c)).sqrt(),
    }
}

/// External-rotor partition for the fragments (linear/spherical/symmetric top cases).
fn external_rotor_partition_products(
    rotor_case: ProductRotorCase,
    temp: f64,
    rot_b1: f64,
    rot_b2: f64,
) -> (f64, f64) {
    match rotor_case {
        ProductRotorCase::LinearAtom => (0.695 * temp / rot_b1, 1.0),
        ProductRotorCase::SphericalAtom => (1.027 * temp / rot_b1 * (temp / rot_b1).sqrt(), 1.0),
        ProductRotorCase::LinearLinear => (0.695 * temp / rot_b1, 0.695 * temp / rot_b2),
        ProductRotorCase::LinearSpherical => (
            1.027 * temp / rot_b1 * (temp / rot_b1).sqrt(),
            0.695 * temp / rot_b2,
        ),
        ProductRotorCase::SphericalSpherical => (
            1.027 * temp / rot_b1 * (temp / rot_b1).sqrt(),
            1.027 * temp / rot_b2 * (temp / rot_b2).sqrt(),
        ),
    }
}

/// Product of QSTAR terms for disappearing degrees of freedom.
fn qstar_product(eps: &[f64], x: &[f64], kt: f64) -> f64 {
    let mut q = 1.0;
    for (e, xval) in eps.iter().zip(x.iter()) {
        let term = gamma_func(1.0 + xval) * (1.0 - (-e / kt).exp()).powf(-xval);
        q *= term;
    }
    q
}
