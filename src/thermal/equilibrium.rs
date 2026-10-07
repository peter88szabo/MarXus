//! Equilibrium constants from the thermochemistry: the standard Gibbs energies of `eval_all_therm_func`
//! (G = gtherm + dh0, translation at the standard pressure p°) give, for sum_r nu_r R_r -> sum_p nu_p P_p,
//!   dG° = sum_p nu_p G_p - sum_r nu_r G_r,   K_p = exp(-dG°/RT)   (referred to p°),
//!   K_c = K_p (p°/k_B T)^dn,   dn = sum_p nu_p - sum_r nu_r,   in (molecule cm-3)^dn
//! (ideal gases; McQuarrie, Statistical Mechanics, 1976, chapter on chemical equilibrium).

use crate::constants::{AU_TO_KCAL, BOLTZMANN_SI, RGAS_AU};
use crate::molecule::MoleculeStruct;

/// Equilibrium constant of one reaction at T and the standard pressure p°.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquilibriumConstant {
    /// dG° (kcal/mol), the standard Gibbs energies at p° (dh0 included).
    pub delta_g_kcal_mol: f64,
    /// dn = sum_p nu_p - sum_r nu_r.
    pub delta_n: f64,
    /// K_p = exp(-dG°/RT), dimensionless, referred to p°.
    pub k_p: f64,
    /// K_c = K_p (p°/k_B T)^dn in (molecule cm-3)^dn.
    pub k_c: f64,
}

/// K_p and K_c of sum_r nu_r R_r -> sum_p nu_p P_p from the thermochemistry of the molecules, evaluated here at
/// `temperature_kelvin`, the standard pressure `standard_pressure_pa` and `freq_cutoff_cm1` (0: harmonic
/// oscillators; above 0: Grimme's quasi-RRHO entropy). Entries are (nu, molecule).
pub fn equilibrium_constant_from_thermochemistry(
    reactants: &mut [(f64, &mut MoleculeStruct)],
    products: &mut [(f64, &mut MoleculeStruct)],
    temperature_kelvin: f64,
    standard_pressure_pa: f64,
    freq_cutoff_cm1: f64,
) -> Result<EquilibriumConstant, String> {
    if !(temperature_kelvin > 0.0) {
        return Err(format!("Equilibrium constant: temperature {temperature_kelvin} K is not positive."));
    }
    if !(standard_pressure_pa > 0.0) {
        return Err(format!("Equilibrium constant: standard pressure {standard_pressure_pa} Pa is not positive."));
    }
    // sum of nu G (Hartree) and of nu
    let side = |molecules: &mut [(f64, &mut MoleculeStruct)]| -> (f64, f64) {
        molecules.iter_mut().fold((0.0, 0.0), |(g, n), (nu, molecule)| {
            molecule.eval_all_therm_func(temperature_kelvin, standard_pressure_pa, freq_cutoff_cm1);
            (g + *nu * molecule.thermo.gtot, n + *nu)
        })
    };
    let (g_reactants, n_reactants) = side(reactants);
    let (g_products, n_products) = side(products);
    let delta_g = g_products - g_reactants;
    let delta_n = n_products - n_reactants;
    let k_p = (-delta_g / (RGAS_AU * temperature_kelvin)).exp();
    // p°/k_B T in molecule cm-3
    let number_density_cm3 = standard_pressure_pa / (BOLTZMANN_SI * temperature_kelvin) * 1.0e-6;
    Ok(EquilibriumConstant {
        delta_g_kcal_mol: delta_g * AU_TO_KCAL,
        delta_n,
        k_p,
        k_c: k_p * number_density_cm3.powf(delta_n),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::{MolType, MoleculeBuilder, MoleculeStruct};

    // Independent closed forms (SI): second radiation constant hc/k (cm K), h, k, amu.
    const C2: f64 = 1.438_776_877;
    const H: f64 = 6.626_070_15e-34;
    const K: f64 = 1.380_649e-23;
    const AMU: f64 = 1.660_539_066_60e-27;

    fn molecule(name: &str, freq: Vec<f64>, brot: Vec<f64>, mass: f64, multi: f64, dh0: f64) -> MoleculeStruct {
        MoleculeBuilder::new(name.to_string(), MolType::mol).freq(freq).brot(brot).mass(mass).multi(multi).dh0(dh0).build()
    }

    fn q_vib(freqs: &[f64], t: f64) -> f64 {
        freqs.iter().map(|nu| 1.0 / (1.0 - (-C2 * nu / t).exp())).product()
    }

    #[test]
    fn an_isomerization_constant_is_the_ratio_of_the_partition_functions() {
        // A <=> B, dn = 0: K = (q_B/q_A) exp(-(E0_B - E0_A)/kT); same mass and rotors, so only the vibrations differ.
        let t = 500.0;
        let mut a = molecule("A", vec![1000.0], vec![2.0, 1.0, 0.5], 30.0, 1.0, 0.0);
        let mut b = molecule("B", vec![500.0], vec![2.0, 1.0, 0.5], 30.0, 1.0, 300.0);
        let k = equilibrium_constant_from_thermochemistry(&mut [(1.0, &mut a)], &mut [(1.0, &mut b)], t, 1.0e5, 0.0).unwrap();
        let expected = q_vib(&[500.0], t) / q_vib(&[1000.0], t) * (-C2 * 300.0 / t).exp();
        assert_eq!(k.delta_n, 0.0);
        assert!((k.k_p / expected - 1.0).abs() < 1e-4, "{} vs {expected}", k.k_p);
        assert_eq!(k.k_c, k.k_p);
    }

    #[test]
    fn an_association_constant_has_the_relative_translation_per_volume() {
        // A + B <=> AB (A nonlinear doublet, B linear triplet, AB nonlinear doublet; m_AB = m_A + m_B):
        //   K_c = q_AB / (q_A q_B) (h^2 / (2 pi mu k T))^(3/2) exp(-dE0/kT)   in cm3,
        // classical rigid rotors sqrt(pi) (T^3/(thA thB thC))^(1/2) and T/th, harmonic vibrations from the ground level.
        let (t, p0): (f64, f64) = (500.0, 1.0e5);
        let mut a = molecule("A", vec![1000.0], vec![2.0, 1.0, 0.5], 30.0, 2.0, 0.0);
        let mut b = molecule("B", vec![1500.0], vec![1.5], 32.0, 3.0, 0.0);
        let mut ab = molecule("AB", vec![800.0, 1200.0], vec![0.5, 0.2, 0.15], 62.0, 2.0, -8000.0);
        let k = equilibrium_constant_from_thermochemistry(&mut [(1.0, &mut a), (1.0, &mut b)], &mut [(1.0, &mut ab)], t, p0, 0.0)
            .unwrap();
        let q_rot_nonlinear = |b: [f64; 3]| std::f64::consts::PI.sqrt() * (t.powi(3) / (C2.powi(3) * b[0] * b[1] * b[2])).sqrt();
        let q_a = q_rot_nonlinear([2.0, 1.0, 0.5]) * q_vib(&[1000.0], t) * 2.0;
        let q_b = t / (C2 * 1.5) * q_vib(&[1500.0], t) * 3.0;
        let q_ab = q_rot_nonlinear([0.5, 0.2, 0.15]) * q_vib(&[800.0, 1200.0], t) * 2.0;
        let mu = 30.0 * 32.0 / 62.0 * AMU;
        let lambda3_m3 = (H * H / (2.0 * std::f64::consts::PI * mu * K * t)).powf(1.5);
        let expected_cm3 = q_ab / (q_a * q_b) * lambda3_m3 * 1.0e6 * (C2 * 8000.0 / t).exp();
        assert_eq!(k.delta_n, -1.0);
        assert!((k.k_c / expected_cm3 - 1.0).abs() < 1e-4, "{} vs {expected_cm3}", k.k_c);
        // K_p = K_c (p°/kT)^dn with p°/kT in molecule cm-3.
        assert!((k.k_p / (k.k_c * p0 / (K * t) * 1.0e-6) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn a_non_positive_temperature_or_pressure_is_refused() {
        let mut a = molecule("A", vec![1000.0], vec![2.0, 1.0, 0.5], 30.0, 1.0, 0.0);
        let mut b = molecule("B", vec![500.0], vec![2.0, 1.0, 0.5], 30.0, 1.0, 300.0);
        assert!(equilibrium_constant_from_thermochemistry(&mut [(1.0, &mut a)], &mut [(1.0, &mut b)], 0.0, 1.0e5, 0.0).is_err());
        assert!(equilibrium_constant_from_thermochemistry(&mut [(1.0, &mut a)], &mut [(1.0, &mut b)], 300.0, 0.0, 0.0).is_err());
    }
}
