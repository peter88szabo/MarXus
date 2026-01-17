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
