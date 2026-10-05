use super::interpol_react_prod::{
    channel_eigval_interpol_from_input, ChannelEigenvalues, SacmInterpolationInput,
};
use super::pst_channels::{convolve_states, PstChannels};
use super::types::SacmEnergyGrid;
use crate::rrkm::sum_and_density::get_rovib_WE_or_rhoE;

/// Build PST W0(E,J): conserved-mode density convolved with the transitional-mode sum of states.
pub fn build_w0_convolution(
    conserved_density: &[f64],
    transitional_sum: &[f64],
    out_len: usize,
    d_e: f64,
) -> Vec<f64> {
    convolve_states(conserved_density, transitional_sum, out_len, d_e)
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
    conserved_density: &[f64],
    transitional_sum: &[f64],
    out_len: usize,
    faminf: Option<&[f64]>,
    energy_offset: f64,
    d_e: f64,
) -> PstChannels {
    let w0 = build_w0_convolution(conserved_density, transitional_sum, out_len, d_e);
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
    rovib_states_from_freqs("sum", freq, brot, grid)
}

/// Compute a density-of-states vector (states per cm-1) using the Beyer-Swinehart routine.
pub fn rovib_density_from_freqs(freq: &[f64], brot: &[f64], grid: SacmEnergyGrid) -> Vec<f64> {
    rovib_states_from_freqs("den", freq, brot, grid)
}

fn rovib_states_from_freqs(what: &str, freq: &[f64], brot: &[f64], grid: SacmEnergyGrid) -> Vec<f64> {
    let nbin = (grid.emax / grid.dE + 0.5) as usize;
    let freq_bins: Vec<usize> = freq.iter().map(|&w| (w / grid.dE + 0.5) as usize).collect();
    get_rovib_WE_or_rhoE(
        what.to_string(),
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
    let conserved_density =
        rovib_density_from_freqs(input.conserved_freq, input.conserved_brot, input.grid);
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
    let out_len = conserved_density.len();
    let pst = build_pst_detailed(
        &conserved_density,
        &transitional_sum,
        out_len,
        input.faminf,
        input.energy_offset,
        input.grid.dE,
    );
    (pst, channel)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conserved_and_transitional_oscillators_count_exactly() {
        // One conserved and one transitional oscillator of 1000 cm-1: at E = 4000 cm-1 the open
        // states are v1 + v2 <= 4, i.e. 15 (W = int rho_cons(x) W_trans(E - x) dx, Forst,
        // Chem. Rev. 71, 339 (1971), eq. 33).
        let grid = SacmEnergyGrid { dE: 10.0, emax: 5000.0 };
        let input = PstDetailedInput {
            conserved_freq: &[1000.0],
            conserved_brot: &[],
            transitional_freq: &[1000.0],
            use_interpolation: false,
            interpolation: None,
            grid,
            faminf: None,
            energy_offset: 0.0,
        };
        let (pst, _) = build_pst_detailed_with_interpolation(input);
        assert!((pst.w_e[400] - 15.0).abs() < 1e-9, "W0(4000) = {}", pst.w_e[400]);
    }
}
