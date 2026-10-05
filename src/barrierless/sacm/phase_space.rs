use super::pst_channels::convolve_states;
use super::types::{SacmEnergyGrid, SacmPhaseSpaceStates, SacmReactantPair};

/// Build PST transition-state states by convolving independent reactant states
/// (rho_ts = rho_A * rho_B, W_ts = rho_A * W_B; Forst, Chem. Rev. 71, 339 (1971), eqs. 30, 33).
pub fn phase_space_states(pair: &SacmReactantPair, grid: SacmEnergyGrid) -> SacmPhaseSpaceStates {
    let states_a = pair.reactant_a.rovib_states(grid);
    let states_b = pair.reactant_b.rovib_states(grid);
    let out_len = (grid.emax / grid.dE + 0.5) as usize + 1;

    let rho_ts = convolve_states(&states_a.rho_e, &states_b.rho_e, out_len, grid.dE);
    let we_ts = convolve_states(&states_a.rho_e, &states_b.we, out_len, grid.dE);

    SacmPhaseSpaceStates { rho_ts, we_ts }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::{MolType, MoleculeBuilder};
    use crate::rrkm::sum_and_density::RotorSymmetry;
    use crate::barrierless::sacm::types::SacmReactant;

    // Linear fragment: a 2D rotor; the single oscillator lies above Emax and contributes nothing.
    fn linear_fragment(b: f64) -> SacmReactant {
        let molecule = MoleculeBuilder::new("frag".to_string(), MolType::mol)
            .freq(vec![1.0e6])
            .brot(vec![b, b])
            .mass(30.0)
            .ene(0.0)
            .build();
        SacmReactant { molecule, rotor: RotorSymmetry::SphericalTop, centrifugal: None }
    }

    #[test]
    fn two_linear_fragments_follow_forst_eq43() {
        // Two 2D rotors: r = 4, Q'_4 = 1/(B1 B2), W(E) = E^2/(2 B1 B2), rho(E) = E/(B1 B2)
        // (Forst, Chem. Rev. 71, 339 (1971), eqs. 30, 33, 40, 43); grain error of order dE/E.
        let (b1, b2) = (1.5, 0.3);
        let pair = SacmReactantPair { reactant_a: linear_fragment(b1), reactant_b: linear_fragment(b2) };
        // dE = 2 cm-1 so that a missing dE factor in the convolution cannot go unnoticed.
        let grid = SacmEnergyGrid { dE: 2.0, emax: 4000.0 };
        let states = phase_space_states(&pair, grid);
        let e = 4000.0;
        let w_ref = e * e / (2.0 * b1 * b2);
        let rho_ref = e / (b1 * b2);
        let w = states.we_ts[2000];
        let rho = states.rho_ts[2000];
        assert!(((w - w_ref) / w_ref).abs() < 2e-3, "W = {w} vs {w_ref}");
        assert!(((rho - rho_ref) / rho_ref).abs() < 2e-3, "rho = {rho} vs {rho_ref}");
    }
}
