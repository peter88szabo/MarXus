use super::types::{SacmEnergyGrid, SacmPhaseSpaceStates, SacmReactantPair};

/// Convolve two state-count arrays to combine independent fragments.
fn convolve(a: &[f64], b: &[f64], out_len: usize) -> Vec<f64> {
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

/// Build PST transition-state states by convolving independent reactant states.
pub fn phase_space_states(pair: &SacmReactantPair, grid: SacmEnergyGrid) -> SacmPhaseSpaceStates {
    let states_a = pair.reactant_a.rovib_states(grid);
    let states_b = pair.reactant_b.rovib_states(grid);
    let out_len = (grid.emax / grid.dE + 0.5) as usize + 1;

    let rho_ts = convolve(&states_a.rho_e, &states_b.rho_e, out_len);
    let we_ts = convolve(&states_a.rho_e, &states_b.we, out_len);

    SacmPhaseSpaceStates { rho_ts, we_ts }
}
