use super::pst_channels::{capture_probability, PstCaptureModel, PstChannels};
use super::pst_detailed::{build_pst_detailed_with_interpolation, PstDetailedInput};
use super::types::SacmEnergyGrid;

/// PST-simplified (Troe-Ushakov 2006): W(E,J) = W0(E,J) * w(E,J).
pub fn build_pst_simplified(
    w0: &[f64],
    grid: SacmEnergyGrid,
    e0_j: f64,
    model: PstCaptureModel,
) -> PstChannels {
    let mut w_e = vec![0.0; w0.len()];

    for (i, &w) in w0.iter().enumerate() {
        let ei = i as f64 * grid.dE;
        let e_total = e0_j + ei;
        let wcap = capture_probability(e_total, e0_j, model);
        w_e[i] = w * wcap;
    }

    PstChannels {
        w_e,
        energy_offset: e0_j,
    }
}

#[derive(Debug, Clone)]
pub struct PstSimplifiedInput<'a> {
    pub pst: PstDetailedInput<'a>,
    pub e0_j: f64,
    pub model: PstCaptureModel,
}

/// PST-simplified from raw mode data with optional SPOL interpolation.
pub fn build_pst_simplified_from_modes(
    input: PstSimplifiedInput<'_>,
) -> (
    PstChannels,
    Option<super::interpol_react_prod::ChannelEigenvalues>,
) {
    let grid = input.pst.grid;
    let (w0, channel) = build_pst_detailed_with_interpolation(input.pst);
    let pst = build_pst_simplified(&w0.w_e, grid, input.e0_j, input.model);
    (pst, channel)
}
