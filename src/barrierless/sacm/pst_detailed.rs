use super::interpol_react_prod::{
    channel_eigval_interpol_from_input, ChannelEigenvalues, SacmInterpolationInput,
};
use super::pst_channels::{convolve_states, PstChannels};
use super::types::SacmEnergyGrid;
use crate::rrkm::sum_and_density::get_rovib_WE_or_rhoE;

/// Build PST W0(E,J) via convolution of conserved and transitional mode counts.
pub fn build_w0_convolution(conserved: &[f64], transitional: &[f64], out_len: usize) -> Vec<f64> {
    convolve_states(conserved, transitional, out_len)
}

/// Apply a FAMINF baseline to W0(E,J)
pub fn apply_faminf(w0: &[f64], faminf: Option<&[f64]>) -> Vec<f64> {
    let mut out = vec![0.0; w0.len()];
    match faminf {
        Some(factors) => {
            for (i, &w) in w0.iter().enumerate() {
                let f = factors.get(i).copied().unwrap_or(1.0);
                out[i] = w * f;
            }
        }
        None => out.copy_from_slice(w0),
    }
    out
}

/// PST-detailed (Fortran-style): convolution plus optional FAMINF baseline.
pub fn build_pst_detailed(
    conserved: &[f64],
    transitional: &[f64],
    out_len: usize,
    faminf: Option<&[f64]>,
    energy_offset: f64,
) -> PstChannels {
    let w0 = build_w0_convolution(conserved, transitional, out_len);
    let w_e = apply_faminf(&w0, faminf);
    PstChannels { w_e, energy_offset }
}

#[derive(Debug, Clone)]
pub struct PstDetailedInput<'a> {
    pub conserved_freq: &'a [f64],
    pub conserved_brot: &'a [f64],
    pub transitional_freq: &'a [f64],
    pub use_interpolation: bool,
    pub interpolation: Option<&'a SacmInterpolationInput>,
    pub grid: SacmEnergyGrid,
    pub faminf: Option<&'a [f64]>,
    pub energy_offset: f64,
}

/// Compute a sum-of-states vector using the existing Beyer-Swinehart routine.
pub fn rovib_sum_from_freqs(freq: &[f64], brot: &[f64], grid: SacmEnergyGrid) -> Vec<f64> {
    let nbin = (grid.emax / grid.dE + 0.5) as usize;
    let freq_bins: Vec<usize> = freq.iter().map(|&w| (w / grid.dE + 0.5) as usize).collect();
    get_rovib_WE_or_rhoE(
        "sum".to_string(),
        freq.len(),
        nbin,
        grid.dE,
        brot.len(),
        &freq_bins,
        brot,
    )
}

/// Build PST detailed channels with optional SPOL interpolation for the transitional modes.
pub fn build_pst_detailed_with_interpolation(
    input: PstDetailedInput<'_>,
) -> (PstChannels, Option<ChannelEigenvalues>) {
    let conserved_sum = rovib_sum_from_freqs(input.conserved_freq, input.conserved_brot, input.grid);
    let channel = if input.use_interpolation {
        input.interpolation.map(channel_eigval_interpol_from_input)
    } else {
        None
    };
    let transitional_freq = match channel.as_ref() {
        Some(channel) => channel.freq_ts.as_slice(),
        None => input.transitional_freq,
    };
    let transitional_sum = rovib_sum_from_freqs(transitional_freq, &[], input.grid);
    let out_len = conserved_sum.len();
    let pst = build_pst_detailed(&conserved_sum, &transitional_sum, out_len, input.faminf, input.energy_offset);
    (pst, channel)
}
