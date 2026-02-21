use super::anharmcorr::anharmonic_product_w_e;
use super::fame_angmomcoupling::factor_angmom;
use super::hindered_rotor::hindered_rotor_sum_factor;
use super::interpol_react_prod::ChannelEigenvalues;
use super::pst_channels::PstChannels;
use super::pst_detailed::{build_pst_detailed_with_interpolation, PstDetailedInput};
use super::types::{ReactantRotorCase, SacmEnergyGrid};

#[derive(Debug, Clone, Copy)]
pub struct OlzmanmAngularInput {
    pub rotor_case: ReactantRotorCase,
    pub spin_factor: f64,
    pub rot_const_a: f64,
    pub rot_const_b: f64,
    pub anisotropy_scale: f64,
    pub alpha_over_beta: f64,
    pub wr_beta: f64,
    pub reactant_freq_sum: f64,
}

#[derive(Debug, Clone)]
pub struct OlzmanmAnharmonicInput {
    pub anharmonic_mode: i32,
    pub frag1_freq: Vec<f64>,
    pub frag2_freq: Vec<f64>,
    pub frag1_anharm: Vec<f64>,
    pub frag2_anharm: Vec<f64>,
}

#[derive(Debug, Clone, Copy)]
pub struct OlzmanmHinderedInput {
    pub vibrational_count: usize,
    pub barrier_height: f64,
}

#[derive(Debug, Clone)]
pub struct OlzmanmCorrections {
    pub angular_coupling: Option<OlzmanmAngularInput>,
    pub anharmonic: Option<OlzmanmAnharmonicInput>,
    pub hindered: Option<OlzmanmHinderedInput>,
    pub faminf_baseline: Option<Vec<f64>>,
}

/// Apply Olzmann-style SACM rigidity (FAME/FAMINF blend) to PST W(E,J).
pub fn apply_olzmanm_rigidity(
    pst: &PstChannels,
    grid: SacmEnergyGrid,
    j: f64,
    corrections: &OlzmanmCorrections,
) -> Vec<f64> {
    let mut out = vec![0.0; pst.w_e.len()];
    let faminf = corrections.faminf_baseline.as_deref().unwrap_or(&[]);

    for (i, &we) in pst.w_e.iter().enumerate() {
        let ei = i as f64 * grid.dE;
        let mut factor = 1.0;

        if let Some(ang) = corrections.angular_coupling {
            let faminf_val = faminf.get(i).copied().unwrap_or(1.0);
            let eif = if ang.reactant_freq_sum > 0.0 {
                let wr = whitten_rabinovitch_factor(2.0 * ei / ang.reactant_freq_sum, ang.wr_beta);
                ei + wr * ang.reactant_freq_sum / 2.0
            } else {
                ei
            };
            let fame_val = factor_angmom(
                ang.rotor_case,
                eif,
                j,
                ang.spin_factor,
                ang.rot_const_a,
                ang.rot_const_b,
            );
            let blend = (-(ang.anisotropy_scale * ang.alpha_over_beta)).exp();
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
            factor *=
                hindered_rotor_sum_factor(hind.vibrational_count as i32, hind.barrier_height, ei);
        }

        out[i] = we * factor;
    }

    out
}

#[derive(Debug, Clone)]
pub struct SacmOlzmanmInput<'a> {
    pub pst: PstDetailedInput<'a>,
    pub j: f64,
    pub corrections: OlzmanmCorrections,
}

pub fn build_olzmanm_channels(
    input: SacmOlzmanmInput<'_>,
) -> (PstChannels, Vec<f64>, Option<ChannelEigenvalues>) {
    let grid = input.pst.grid;
    let (pst, channel) = build_pst_detailed_with_interpolation(input.pst);
    let corrected = apply_olzmanm_rigidity(&pst, grid, input.j, &input.corrections);
    (pst, corrected, channel)
}

/// Whitten-Rabinovitch correction factor used for the energy shift in Fame.
fn whitten_rabinovitch_factor(reduced_energy: f64, beta: f64) -> f64 {
    let w = if reduced_energy > 1.0 {
        10.0_f64.powf(-1.0506 * reduced_energy.powf(0.25))
    } else {
        1.0 / (5.0 * reduced_energy + 2.73 * reduced_energy.sqrt() + 3.51)
    };
    1.0 - beta * w
}
