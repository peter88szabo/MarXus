use super::anharmcorr::anharmonic_product_w_e;
use super::hindered_rotor::hindered_rotor_sum_factor;
use super::threshold_energy::morse_threshold_energy;
use super::types::{SacmEnergyGrid, SacmJResolved};
use super::system::{prepare_sacm, SacmInput, SacmPrepared};
use super::thermal_rates::{compute_thermal_rates, SacmThermalInput, SacmThermalRate};
use super::fame_angmomcoupling::factor_angmom;

#[derive(Debug, Clone, Copy)]
pub enum SacmDensitySource {
    ReactantA,
    ReactantB,
}

#[derive(Debug, Clone, Copy)]
pub struct SacmThreshold {
    pub dissociation_energy: f64,
    pub product_zpe: f64,
    pub zpe_shift: f64,
    pub alpha_over_beta: f64,
    pub base_rot_const: f64,
    pub centrifugal_a1: Option<f64>,
    pub centrifugal_a2: Option<f64>,
}

#[derive(Debug, Clone, Copy)]
pub struct SacmSymmetry {
    pub reactant_symmetry: f64,
    pub transition_symmetry: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct SacmAngularCouplingInput {
    pub rotor_case: i32,
    pub spin_factor: f64,
    pub rot_const_a: f64,
    pub rot_const_b: f64,
    pub anisotropy_scale: f64,
    pub wr_beta: f64,
    pub reactant_freq_sum: f64,
}

#[derive(Debug, Clone)]
pub struct SacmAnharmonicInput {
    pub anharmonic_mode: i32,
    pub frag1_freq: Vec<f64>,
    pub frag2_freq: Vec<f64>,
    pub frag1_anharm: Vec<f64>,
    pub frag2_anharm: Vec<f64>,
}

#[derive(Debug, Clone, Copy)]
pub struct SacmHinderedInput {
    pub vibrational_count: usize,
    pub barrier_height: f64,
}

#[derive(Debug, Clone)]
pub struct SacmCorrections {
    pub alpha_over_beta: f64,
    pub angular_coupling: Option<SacmAngularCouplingInput>,
    pub anharmonic: Option<SacmAnharmonicInput>,
    pub hindered: Option<SacmHinderedInput>,
    pub faminf_baseline: Option<Vec<f64>>,
}

#[derive(Debug, Clone)]
pub struct SacmRateConfig {
    pub threshold: SacmThreshold,
    pub symmetry: SacmSymmetry,
    pub density_source: SacmDensitySource,
    pub capture_geometry: Option<super::types::SacmCaptureGeometry>,
    pub capture_model: Option<SacmCaptureModel>,
    pub thermal_input: Option<SacmThermalInput>,
}

#[derive(Debug, Clone, Copy)]
pub struct SacmCaptureModel {
    pub beta: f64,
    pub use_triatomic: bool,
}

#[derive(Debug)]
pub struct SacmRateCurve {
    pub j: usize,
    pub k_e: Vec<f64>,
    pub we_corr: Vec<f64>,
    pub energy_offset: f64,
}

#[derive(Debug)]
pub struct SacmEnergyRate {
    pub k_e: Vec<f64>,
    pub rho_e: Vec<f64>,
}


/// Run SACM: PST baseline plus correction factors, returning k(E) per J.
pub fn run_sacm(
    input: &SacmInput,
    config: SacmRateConfig,
    corrections: SacmCorrections,
) -> (SacmPrepared, Vec<SacmRateCurve>, SacmEnergyRate, Vec<SacmThermalRate>) {
    let prepared = prepare_sacm(input);
    let mut config = config;
    if config.capture_geometry.is_none() {
        config.capture_geometry = prepared.capture_geometry.clone();
    }
    if let (Some(geom), Some(model)) = (&config.capture_geometry, config.capture_model) {
        let derived = derive_centrifugal_params(geom, model);
        config.threshold.centrifugal_a1 = Some(derived.0);
        config.threshold.centrifugal_a2 = Some(derived.1);
    }

    let reactant_zpe = match config.density_source {
        SacmDensitySource::ReactantA => input.reactants.reactant_a.molecule.zpe,
        SacmDensitySource::ReactantB => input.reactants.reactant_b.molecule.zpe,
    };
    let thermal_rates = match config.thermal_input.as_ref() {
        Some(thermal_input) => compute_thermal_rates(thermal_input),
        None => Vec::new(),
    };
    let (curves, energy_rate) =
        sacm_rates(&prepared, input.grid, reactant_zpe, config, &corrections);
    (prepared, curves, energy_rate, thermal_rates)
}

/// Build J-resolved rate curves from PST and correction factors.
fn sacm_rates(
    prepared: &SacmPrepared,
    grid: SacmEnergyGrid,
    reactant_zpe: f64,
    config: SacmRateConfig,
    corrections: &SacmCorrections,
) -> (Vec<SacmRateCurve>, SacmEnergyRate) {
    let nbin = (grid.emax / grid.dE + 0.5) as usize;
    let density = match config.density_source {
        SacmDensitySource::ReactantA => &prepared.reactant_a_j,
        SacmDensitySource::ReactantB => &prepared.reactant_b_j,
    };

    let mut numerator = vec![0.0; nbin + 1];
    let mut rho_total = vec![0.0; nbin + 1];

    let curves: Vec<SacmRateCurve> = density
        .iter()
        .map(|j_states| {
            let e0_j = morse_threshold_energy(
                config.threshold.dissociation_energy,
                config.threshold.product_zpe,
                config.threshold.zpe_shift,
                config.threshold.alpha_over_beta,
                config.threshold.base_rot_const,
                config.threshold.centrifugal_a1.unwrap_or(0.0),
                config.threshold.centrifugal_a2.unwrap_or(0.0),
                j_states.j as i32,
            );

            if e0_j < 0.0 {
                return SacmRateCurve {
                    j: j_states.j,
                    k_e: vec![0.0; nbin + 1],
                    we_corr: vec![0.0; nbin + 1],
                    energy_offset: 0.0,
                };
            }

            let energy_offset = e0_j - reactant_zpe;
            let start_bin = ((energy_offset / grid.dE) + 0.5).floor().max(0.0) as usize;
            let base_we = &prepared.phase_space.we_ts;
            let we_corr = apply_corrections(base_we, j_states, grid.dE, corrections);

            let factor =
                config.symmetry.reactant_symmetry / (config.symmetry.transition_symmetry * 3.3356e-11);
            let mut k_e = vec![0.0; nbin + 1];
            for i in 0..=nbin {
                let idx = i + start_bin;
                if idx > nbin || idx >= j_states.states.rho_ej.len() {
                    break;
                }
                let rho = j_states.states.rho_ej[idx];
                if rho > 0.0 {
                    k_e[i] = factor * we_corr[i] / rho;
                    numerator[idx] += k_e[i] * rho;
                }
                rho_total[idx] += rho;
            }

            SacmRateCurve {
                j: j_states.j,
                k_e,
                we_corr,
                energy_offset,
            }
        })
        .collect();

    let mut k_e_integrated = vec![0.0; nbin + 1];
    for i in 0..=nbin {
        if rho_total[i] > 0.0 {
            k_e_integrated[i] = numerator[i] / rho_total[i];
        }
    }

    (
        curves,
        SacmEnergyRate {
            k_e: k_e_integrated,
            rho_e: rho_total,
        },
    )
}

/// Apply Fam/anharmonic/hindered-rotor corrections to PST sum of states.
fn apply_corrections(
    base_we: &[f64],
    j_states: &SacmJResolved,
    dE: f64,
    corrections: &SacmCorrections,
) -> Vec<f64> {
    let mut out = vec![0.0; base_we.len()];
    let use_corrections = corrections.alpha_over_beta > 0.0;
    let faminf = corrections
        .faminf_baseline
        .as_deref()
        .unwrap_or(&[]);

    for (i, &we) in base_we.iter().enumerate() {
        let ei = i as f64 * dE;
        let mut factor = 1.0;

        if use_corrections {
            if let Some(fam) = corrections.angular_coupling {
                let faminf_val = faminf.get(i).copied().unwrap_or(1.0);
                let eif = if fam.reactant_freq_sum > 0.0 {
                    let wr = whitten_rabinovitch_factor(2.0 * ei / fam.reactant_freq_sum, fam.wr_beta);
                    ei + wr * fam.reactant_freq_sum / 2.0
                } else {
                    ei
                };
                let fame_val = factor_angmom(
                    fam.rotor_case,
                    eif,
                    j_states.j as f64,
                    fam.spin_factor,
                    fam.rot_const_a,
                    fam.rot_const_b,
                );
                let blend = (-(fam.anisotropy_scale * corrections.alpha_over_beta)).exp();
                let wpst = faminf_val + (fame_val - faminf_val) * blend;
                factor *= wpst;
            }

            if let Some(anh) = &corrections.anharmonic {
                let nvib1 = anh.frag1_freq.len();
                let nvib2 = anh.frag2_freq.len();
                factor *= anharmonic_product_w_e(
                    anh.anharmonic_mode,
                    ei,
                    nvib1,
                    nvib2,
                    &anh.frag1_freq,
                    &anh.frag2_freq,
                    &anh.frag1_anharm,
                    &anh.frag2_anharm,
                );
            }

            if let Some(hind) = corrections.hindered {
                factor *= hindered_rotor_sum_factor(hind.vibrational_count as i32, hind.barrier_height, ei);
            }
        }

        out[i] = we * factor;
    }

    out
}

/// Whitten-Rabinovitch correction factor used for the energy shift in Fam.
fn whitten_rabinovitch_factor(reduced_energy: f64, beta: f64) -> f64 {
    let w = if reduced_energy > 1.0 {
        10.0_f64.powf(-1.0506 * reduced_energy.powf(0.25))
    } else {
        1.0 / (5.0 * reduced_energy + 2.73 * reduced_energy.sqrt() + 3.51)
    };
    1.0 - beta * w
}

/// Derive centrifugal A1/A2 from capture geometry for pseudo-diatomic/triatomic cases.
fn derive_centrifugal_params(
    geom: &super::types::SacmCaptureGeometry,
    model: SacmCaptureModel,
) -> (f64, f64) {
    let qe = geom.distance_primary;
    if qe <= 0.0 || model.beta <= 0.0 {
        return (0.0, 0.0);
    }

    let base_a1 = 2.0 / (model.beta * qe);
    let base_a2 = 1.0 / (model.beta * qe * model.beta * qe);

    if !model.use_triatomic {
        return (base_a1, base_a2);
    }

    let (q2e, angle_deg, mz) = match (geom.distance_secondary, geom.angle_deg, geom.mass_z) {
        (Some(q2e), Some(angle_deg), Some(mz)) => (q2e, angle_deg, mz),
        _ => return (base_a1, base_a2),
    };

    let mx = geom.mass_x;
    let my = geom.mass_y;
    let theta = angle_deg.to_radians();
    let cos_theta = theta.cos();

    let anen = mx * (my + mz) * q2e * q2e
        - 2.0 * mx * mz * q2e * qe * cos_theta
        + mz * (mx + my) * qe * qe;

    if anen == 0.0 {
        return (base_a1, base_a2);
    }

    let a1 = base_a1 * (mz * (mx + my) * qe * qe - mx * mz * q2e * qe * cos_theta) / anen;
    let a2 = base_a2 * mz * (mx + my) * qe * qe / anen;
    (a1, a2)
}
