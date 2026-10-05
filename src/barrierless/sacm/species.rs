use crate::rrkm::sum_and_density::{
    get_Jres_rovib_WEJ_or_rhoEJ, get_rovib_WE_or_rhoE, RotorSymmetry,
};

use super::types::{SacmEnergyGrid, SacmJRange, SacmJResolved, SacmReactant, SacmReactantStates};

/// Map harmonic frequencies to energy-bin indices for Beyer-Swinehart counting.
fn build_frequency_bins(freq: &[f64], dE: f64) -> Vec<usize> {
    freq.iter().map(|&w| (w / dE + 0.5) as usize).collect()
}

impl SacmReactant {
    /// Sum of states for the reactant on the given energy grid.
    pub fn we(&self, grid: SacmEnergyGrid) -> Vec<f64> {
        self.rovib_states(grid).we
    }

    /// Density of states for the reactant on the given energy grid.
    pub fn rho_e(&self, grid: SacmEnergyGrid) -> Vec<f64> {
        self.rovib_states(grid).rho_e
    }

    /// Compute rovibrational sum/density of states using RRKM counting.
    pub fn rovib_states(&self, grid: SacmEnergyGrid) -> SacmReactantStates {
        let nbin = (grid.emax / grid.dE + 0.5) as usize;
        let freq_bins = build_frequency_bins(&self.molecule.freq, grid.dE);

        let rho_e = get_rovib_WE_or_rhoE(
            "den".to_string(),
            self.molecule.freq.len(),
            nbin,
            grid.dE,
            self.molecule.brot.len(),
            &freq_bins,
            &self.molecule.brot,
        );

        let we = get_rovib_WE_or_rhoE(
            "sum".to_string(),
            self.molecule.freq.len(),
            nbin,
            grid.dE,
            self.molecule.brot.len(),
            &freq_bins,
            &self.molecule.brot,
        );

        SacmReactantStates { rho_e, we }
    }

    /// Compute J-resolved states, using Bcent for prolate tops when needed.
    pub fn j_resolved_states(
        &self,
        grid: SacmEnergyGrid,
        j_range: SacmJRange,
    ) -> Vec<SacmJResolved> {
        let nbin = (grid.emax / grid.dE + 0.5) as usize;
        let freq_bins = build_frequency_bins(&self.molecule.freq, grid.dE);

        // The J,K-resolved sums run over the VIBRATIONAL sum/density at E - E_rot(J,K)
        // (Forst, Comput. Chem. 20, 419 (1996), eqs. 5-7): the overall rotation is added by
        // get_Jres_rovib_WEJ_or_rhoEJ and must not already be contained in the base arrays.
        let rho_e = get_rovib_WE_or_rhoE(
            "den".to_string(),
            self.molecule.freq.len(),
            nbin,
            grid.dE,
            0,
            &freq_bins,
            &[],
        );

        let we = get_rovib_WE_or_rhoE(
            "sum".to_string(),
            self.molecule.freq.len(),
            nbin,
            grid.dE,
            0,
            &freq_bins,
            &[],
        );

        let mut results = Vec::new();
        let mut j = j_range.j_start;
        while j <= j_range.j_end {
            let b_effective = match self.rotor {
                RotorSymmetry::ProlateSymmetricTop => Some(
                    self.centrifugal
                        .expect("Prolate SACM requires centrifugal parameters for Bcent")
                        .bcentrifugal(j),
                ),
                _ => None,
            };
            let states = get_Jres_rovib_WEJ_or_rhoEJ(
                self.rotor,
                j,
                grid.dE,
                nbin,
                &self.molecule.brot,
                &rho_e,
                &we,
                b_effective,
            );
            results.push(SacmJResolved { j, states });
            j += j_range.j_step;
        }

        results
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::{MolType, MoleculeBuilder};

    #[test]
    fn j_zero_spherical_top_counts_only_vibrations() {
        // Forst, Comput. Chem. 20, 419 (1996), eqs. 5-7: the J,K-resolved sums run over the
        // VIBRATIONAL G_v, N_v at E - E_r. For a spherical top at J = 0 only the single K = 0
        // rotational state exists, so W(E,0) = G_v(E) and rho(E,0) = N_v(E).
        // One 1000 cm-1 oscillator: at E = 3000 cm-1, G_v = 4 and N_v = 1/dE.
        let molecule = MoleculeBuilder::new("top".to_string(), MolType::mol)
            .freq(vec![1000.0])
            .brot(vec![1.0, 1.0, 1.0])
            .mass(30.0)
            .ene(0.0)
            .build();
        let reactant = SacmReactant { molecule, rotor: RotorSymmetry::SphericalTop, centrifugal: None };
        let grid = SacmEnergyGrid { dE: 10.0, emax: 5000.0 };
        let jres = reactant.j_resolved_states(grid, SacmJRange { j_start: 0, j_end: 0, j_step: 1 });
        let states = &jres[0].states;
        assert!((states.wej[300] - 4.0).abs() < 1e-9, "W(3000, J=0) = {}", states.wej[300]);
        assert!((states.rho_ej[300] - 0.1).abs() < 1e-12, "rho(3000, J=0) = {}", states.rho_ej[300]);
    }
}
