//! Energy graining: fine counting cells -> master-equation grains.
//!
//! States are counted on fine cells (1 cm-1 by default), and every convolution (transition-state sums
//! of states, tunneling, inverse Laplace transforms, fragment pairs) is done on the cells. The master
//! equation then works with grains, contiguous blocks of cells: "the energy axis [is divided] into a
//! set of contiguous intervals or grains for which mean values of energy, microcanonical rate
//! coefficient and density of states are assigned" (Pilling, Robertson, Annu. Rev. Phys. Chem. 54, 245
//! (2003), p. 254). Grain values are cell averages:
//!   rho_g = (1/n) sum_{c in g} rho_c                (states per cm-1, averaged over the grain)
//!   W_g   = (1/n) sum_{c in g} N_c                  (number of states of the transition state)
//!   k_g   = W_g/(h rho_g) = sum_c N_c / (h sum_c rho_c),
//! i.e. the grain rate coefficient is the flux-weighted mean (sum of cell fluxes over sum of cell
//! states). Grains lie on one absolute grid common to all wells: grain g collects the absolute cells
//! c = g n - floor(n/2) ... g n - floor(n/2) + n - 1 and is centred at g dE (exactly for odd n, half a cell
//! above for even n), dE = n x cell width; Boltzmann factors and the collision kernel use these common
//! centres, which keeps detailed balance between wells exact. Cells below the first counted cell of a
//! species hold no states, so the lowest grain of a well may be partially filled.

/// Cells and grains of the master equation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GrainGrid {
    pub cell_width_cm1: f64,
    pub cells_per_grain: usize,
}

impl GrainGrid {
    /// Grid with the requested grain width rounded to a whole number (>= 1) of cells.
    pub fn new(cell_width_cm1: f64, grain_width_cm1: f64) -> Result<Self, String> {
        if !(cell_width_cm1 > 0.0) || !(grain_width_cm1 > 0.0) {
            return Err(format!(
                "Cell width ({cell_width_cm1} cm-1) and grain width ({grain_width_cm1} cm-1) must be positive."
            ));
        }
        let cells_per_grain = ((grain_width_cm1 / cell_width_cm1).round() as usize).max(1);
        Ok(Self { cell_width_cm1, cells_per_grain })
    }

    pub fn grain_width_cm1(&self) -> f64 {
        self.cell_width_cm1 * self.cells_per_grain as f64
    }

    /// Absolute cell of the energy (cm-1) on the common scale.
    pub fn cell_of_energy(&self, energy_cm1: f64) -> isize {
        (energy_cm1 / self.cell_width_cm1).round() as isize
    }

    /// First absolute cell of grain g.
    pub fn first_cell_of_grain(&self, grain: isize) -> isize {
        let n = self.cells_per_grain as isize;
        grain * n - n / 2
    }

    /// Grain that contains the absolute cell.
    pub fn grain_of_cell(&self, cell: isize) -> isize {
        let n = self.cells_per_grain as isize;
        (cell + n / 2).div_euclid(n)
    }
}

/// Averages over the grains first_grain..=last_grain of cell values given for the absolute cells
/// first_cell, first_cell + 1, ...; cells below first_cell count as zero. The values must cover the
/// highest cell of last_grain.
pub fn average_over_grains(
    grid: &GrainGrid,
    values: &[f64],
    first_cell: isize,
    first_grain: isize,
    last_grain: isize,
) -> Result<Vec<f64>, String> {
    let n = grid.cells_per_grain;
    let last_cell = grid.first_cell_of_grain(last_grain) + n as isize - 1;
    if last_cell - first_cell >= values.len() as isize {
        return Err(format!(
            "Cell values cover the absolute cells {first_cell}..{} but grain {last_grain} needs cell {last_cell}.",
            first_cell + values.len() as isize - 1
        ));
    }
    Ok((first_grain..=last_grain)
        .map(|g| {
            let start = grid.first_cell_of_grain(g);
            (start..start + n as isize)
                .filter(|&c| c >= first_cell)
                .map(|c| values[(c - first_cell) as usize])
                .sum::<f64>()
                / n as f64
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grains_are_centred_on_the_common_grid() {
        let grid = GrainGrid::new(1.0, 5.0).unwrap();
        assert_eq!(grid.first_cell_of_grain(0), -2);
        assert_eq!(grid.first_cell_of_grain(1), 3);
        assert_eq!(grid.first_cell_of_grain(-1), -7);
        for cell in -12..13 {
            let g = grid.grain_of_cell(cell);
            assert!(grid.first_cell_of_grain(g) <= cell && cell < grid.first_cell_of_grain(g + 1), "cell {cell}");
        }
        assert_eq!(grid.cell_of_energy(41.6), 42);
        assert_eq!(grid.cell_of_energy(-3.4), -3);
    }

    #[test]
    fn grain_width_is_rounded_to_whole_cells() {
        let grid = GrainGrid::new(1.0, 41.7).unwrap();
        assert_eq!(grid.cells_per_grain, 42);
        assert_eq!(grid.grain_width_cm1(), 42.0);
        assert_eq!(GrainGrid::new(1.0, 0.3).unwrap().cells_per_grain, 1);
        assert!(GrainGrid::new(0.0, 10.0).is_err());
        assert!(GrainGrid::new(1.0, -5.0).is_err());
    }

    #[test]
    fn grain_averages_conserve_the_cell_sum() {
        let grid = GrainGrid::new(1.0, 5.0).unwrap();
        // Cells 3 .. 3+29 hold the values 1..30; grains 0..=6 cover cells -2..=32.
        let values: Vec<f64> = (1..=30).map(|v| v as f64).collect();
        let averages = average_over_grains(&grid, &values, 3, 0, 6).unwrap();
        assert_eq!(averages.len(), 7);
        assert_eq!(averages[0], 0.0); // cells -2..2 lie below the first cell
        assert_eq!(averages[1], (1.0 + 2.0 + 3.0 + 4.0 + 5.0) / 5.0);
        let total: f64 = averages.iter().map(|a| a * 5.0).sum();
        assert_eq!(total, values.iter().sum::<f64>());
    }

    #[test]
    fn a_partially_filled_lowest_grain_averages_over_the_whole_grain() {
        let grid = GrainGrid::new(1.0, 5.0).unwrap();
        // Values start at cell 5: grain 1 (cells 3..7) holds cells 5, 6, 7.
        let averages = average_over_grains(&grid, &[1.0; 10], 5, 1, 2).unwrap();
        assert_eq!(averages, vec![3.0 / 5.0, 1.0]);
    }

    #[test]
    fn values_must_cover_the_highest_grain() {
        let grid = GrainGrid::new(1.0, 5.0).unwrap();
        // Ten values cover the cells 3..=12: grain 2 (cells 8..=12) is covered, grain 3 (13..=17) is not.
        assert!(average_over_grains(&grid, &[1.0; 10], 3, 0, 2).is_ok());
        assert!(average_over_grains(&grid, &[1.0; 10], 3, 0, 3).is_err());
    }
}
