use super::types::SacmEnergyGrid;

pub const HCM: f64 = 3.3356e-11;
pub const KB_CM: f64 = 0.69503;

#[derive(Debug, Clone, Copy)]
pub enum PstCaptureModel {
    Step,
    PowerLaw { exponent_n: f64 },
}

#[derive(Debug, Clone)]
pub struct PstChannels {
    pub w_e: Vec<f64>,
    pub energy_offset: f64,
}

/// Convolve two state-count arrays to combine conserved and transitional modes.
pub fn convolve_states(a: &[f64], b: &[f64], out_len: usize) -> Vec<f64> {
    let mut out = vec![0.0; out_len];
    for i in 0..out_len {
        let mut sum = 0.0;
        for j in 0..=i {
            if j < a.len() && (i - j) < b.len() {
                sum += a[j] * b[i - j];
            }
        }
        out[i] = sum;
    }
    out
}

/// Build PST open-channel counts by convolving conserved and transitional sums of states.
pub fn build_convolved_channels(
    conserved: &[f64],
    transitional: &[f64],
    out_len: usize,
) -> Vec<f64> {
    convolve_states(conserved, transitional, out_len)
}

/// Capture probability from Troe-Ushakov (2006), either step or power-law form.
pub fn capture_probability(e_total: f64, e0: f64, model: PstCaptureModel) -> f64 {
    if e_total <= 0.0 || e_total <= e0 {
        return 0.0;
    }
    match model {
        PstCaptureModel::Step => 1.0,
        PstCaptureModel::PowerLaw { exponent_n } => {
            let term = 1.0 - e0 / e_total;
            if term <= 0.0 {
                0.0
            } else {
                term.powf(exponent_n)
            }
        }
    }
}

/// Apply the capture probability to a convolved W(E) to build PST W(E,J).
pub fn build_pst_channels(
    w_convolved: &[f64],
    grid: SacmEnergyGrid,
    e0_j: f64,
    model: PstCaptureModel,
) -> PstChannels {
    let mut w_e = vec![0.0; w_convolved.len()];
    for (i, &w) in w_convolved.iter().enumerate() {
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

/// Compute k(E,J) from W(E,J) and rho(E,J) with a J-dependent energy offset.
pub fn microcanonical_rate(
    w_e: &[f64],
    rho_ej: &[f64],
    grid: SacmEnergyGrid,
    energy_offset: f64,
    symmetry_factor: f64,
) -> Vec<f64> {
    let mut k_e = vec![0.0; w_e.len()];
    let start_bin = ((energy_offset / grid.dE) + 0.5).floor().max(0.0) as usize;
    for i in 0..w_e.len() {
        let idx = i + start_bin;
        if idx >= rho_ej.len() {
            break;
        }
        let rho = rho_ej[idx];
        if rho > 0.0 {
            k_e[i] = symmetry_factor * w_e[i] / rho;
        }
    }
    k_e
}
