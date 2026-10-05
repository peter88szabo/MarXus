//! Chemical-activation network from an input deck in the MESS input format (`mess_input.rs`).
//!
//! Energy grid (`energy_graining.rs`). All state counting and every convolution is done on fine cells
//! (1 cm-1 by default): densities of states, transition-state numbers of states, Eckart tunneling,
//! inverse Laplace transforms and fragment-pair densities. Every energy of the deck is placed once on
//! the absolute cell grid. The master-equation grains are contiguous blocks of cells on one absolute grid
//! common to all wells, centred at g dE, with grain width dE from the settings or
//! EnergyStepOverTemperature x k_B x (lowest temperature of the deck), rounded to whole cells. Grain values
//! are cell averages: rho_g = <rho>_g and k_g = <W>_g/(h <rho>_g), the sum of the cell fluxes over the
//! sum of the cell states (Pilling, Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003), p. 254). Both
//! directions of an isomerization use the same <W>_g, so rho_a k_ab = rho_b k_ba holds exactly (detailed
//! balance, `chemical_activation_operator.rs`). All grids end at the same top: the settings, otherwise
//! the highest barrier or asymptote + ExcessEnergyOverTemperature x k_B T at the highest temperature of
//! the deck, otherwise ModelEnergyLimit. Each well starts at its lowest grain that contains states.
//!
//! Microcanonical rate coefficients. Tight transition states: k(E) = W‡(E - E0)/(h rho(E)) (RRKM;
//! PO14 eq. 9; Forst 1973 Sec. 4.5 for the symmetry and degeneracy factors carried by W‡ and rho).
//! Eckart tunneling (`Tunneling Eckart` block): W‡ is replaced by the convolution of Miller, J. Am. Chem.
//! Soc. 101, 6810 (1979), eqs. 8-9, N_QM(E) = integral dE1 P'(E1) N‡(E - E1) with the exact Eckart
//! transmission probability of the deck's imaginary frequency and two well depths
//! (`tunneling::eckart_tunneling_sum_of_states`); k(E) is then non-zero down to the higher of the two
//! asymptotes, min(well depths) below the transition state, while the channel's classical threshold
//! (the transition-state grain) remains the reference of the absorbing barrier. Other tunneling models
//! are refused; `ignore_tunneling` leaves tunneling out.
//! Barrierless channels: the inverse Laplace transform of the high-pressure rate coefficient given in
//! the MarXus `InverseLaplaceTransform` block of the barrier (`barrierless::ilt::ilt_barrierless`;
//! Davies, Green, Pilling, Chem. Phys. Lett. 126, 373 (1986)):
//!   association, k_inf in cm3 s-1: W(E - E_th) from the convolved density of the two fragments
//!     (rho_AB(E) = sum_E' rho_A(E') rho_B(E - E') dE), E_th = E(asymptote) + E_inf;
//!   dissociation, k_inf in s-1: W(E - E_th) from the density of the well, E_th = E(well) + E_inf,
//!     which may not lie below the dissociation asymptote.
//! A phase-space-theory core without an ILT block is refused: the barrierless module is not yet
//! connected for chemical activation.
//!
//! Collisions: Lennard-Jones parameters combined as sigma = (sigma_1 + sigma_2)/2,
//! eps = sqrt(eps_1 eps_2) (Troe, J. Chem. Phys. 66, 4758 (1977), Sec. III), reduced mass of the two
//! Masses[amu]; <dE_down>(T) = Factor (T/300 K)^Power (Factor is the value at 300 K in this input
//! format); exponential down with the deck's ExponentCutoff. The `Escape` pseudo-first-order rate
//! constant of a well is the bimolecular sink k_c[D] (PO14 eq. 2).

use std::collections::HashMap;

use crate::barrierless::ilt::ilt_barrierless::{ilt_sum_of_states_association, ilt_sum_of_states_dissociation};
use crate::constants::{CM1_TO_KCAL, H_PLANCK_CM, KB_CM};
use crate::tunneling::tunneling::eckart_tunneling_sum_of_states;
use crate::utils::atomic_masses::mass_vector_from_symbols_amu;

use super::chemical_activation_network::{
    Channel, ChannelDestination, ChemicalActivationNetwork, CollisionModel, EnergyTransferParameters,
    LennardJonesPair, Well,
};
use super::energy_graining::{average_over_grains, GrainGrid};
use super::mess_input::{
    rotational_constants_from_geometry_cm1, IltDirection, MessBarrierCore, MessDeck, MessSpeciesRrho,
    TunnelingSpecification,
};
use super::microcanonical_builder::{rrho_density_of_states, rrho_sum_of_states, SpeciesMicroModel};

/// Conversion factor from Epsilons[1/cm] to K (hc/k_B in cm K).
const CM1_TO_KELVIN: f64 = 1.438_776_877;
/// Reference temperature of the energy-transfer Factor in this input format, K.
const ENERGY_TRANSFER_REFERENCE_TEMPERATURE: f64 = 300.0;

/// Choices that the input deck does not fix.
#[derive(Debug, Clone)]
pub struct MessNetworkSettings {
    /// Grain width of the master equation in cm-1 (rounded to whole cells); None:
    /// EnergyStepOverTemperature x k_B x lowest temperature of the deck.
    pub grain_width_cm1: Option<f64>,
    /// Width of the counting cells in cm-1 (all state counting and convolutions), default 1 cm-1.
    pub cell_width_cm1: f64,
    /// Absolute top of all grids in cm-1 on the energy scale of the deck; None: highest barrier or
    /// asymptote + ExcessEnergyOverTemperature x k_B x highest temperature, otherwise ModelEnergyLimit.
    pub top_energy_cm1: Option<f64>,
    /// Leave tunneling out of k(E) although the deck has Tunneling blocks.
    pub ignore_tunneling: bool,
}

impl Default for MessNetworkSettings {
    fn default() -> Self {
        Self { grain_width_cm1: None, cell_width_cm1: 1.0, top_energy_cm1: None, ignore_tunneling: false }
    }
}

/// Network and conditions built from an input deck.
#[derive(Debug, Clone)]
pub struct MessChemicalActivationModel {
    pub network: ChemicalActivationNetwork,
    pub temperatures_kelvin: Vec<f64>,
    pub pressures_torr: Vec<f64>,
    /// Exponential down with the ExponentCutoff of the deck.
    pub collision_model: CollisionModel,
    /// Channels (well, channel) through which the `Reactant` of the deck forms the wells.
    pub entrance_channels: Vec<(usize, usize)>,
}

/// Build the chemical-activation network of an input deck.
pub fn chemical_activation_model_from_mess(
    deck: &MessDeck,
    settings: &MessNetworkSettings,
) -> Result<MessChemicalActivationModel, String> {
    let global = &deck.global;
    let lowest_temperature = global.temperatures_kelvin.iter().cloned().fold(f64::INFINITY, f64::min);
    let highest_temperature = global.temperatures_kelvin.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Cells and grains.
    let requested_grain = match settings.grain_width_cm1 {
        Some(width) => width,
        None => {
            global
                .energy_step_over_temperature
                .ok_or("Input deck: no EnergyStepOverTemperature and no grain width in the settings.")?
                * KB_CM
                * lowest_temperature
        }
    };
    let grid = GrainGrid::new(settings.cell_width_cm1, requested_grain)?;
    let cell = grid.cell_width_cm1;

    // Top of the grids: the highest barrier or asymptote plus ExcessEnergyOverTemperature kT at the
    // highest temperature (the master-equation top of the input format), otherwise ModelEnergyLimit.
    let top_cm1 = match settings.top_energy_cm1 {
        Some(top) => top,
        None => match global.excess_energy_over_temperature {
            Some(excess) => {
                let highest = deck
                    .barriers
                    .iter()
                    .map(|b| b.rrho.zero_energy_cm1)
                    .chain(deck.bimolecular.values().map(|b| b.ground_energy_cm1))
                    .fold(f64::NEG_INFINITY, f64::max);
                if !highest.is_finite() {
                    return Err("Input deck: no barrier or bimolecular species to place the top of the grid.".into());
                }
                highest + excess * KB_CM * highest_temperature
            }
            None => {
                global
                    .model_energy_limit_kcal_mol
                    .ok_or("Input deck: no ExcessEnergyOverTemperature, ModelEnergyLimit or top energy in the settings.")?
                    / CM1_TO_KCAL
            }
        },
    };
    let top_grain = grid.grain_of_cell(grid.cell_of_energy(top_cm1));
    let last_cell = grid.first_cell_of_grain(top_grain + 1) - 1;
    // Number of cells from absolute cell `from` up to the last cell of the top grain.
    let cells_from = |from: isize, what: &str| -> Result<usize, String> {
        if last_cell - from < 1 {
            return Err(format!("{what} lies at or above the top of the energy grid ({top_cm1} cm-1)."));
        }
        Ok((last_cell - from + 1) as usize)
    };

    // Collision parameters (same for all wells in this input format).
    let (eps_1, eps_2) = global.lj_epsilons_cm1.ok_or("Input deck: missing Epsilons[1/cm].")?;
    let (sigma_1, sigma_2) = global.lj_sigmas_angstrom.ok_or("Input deck: missing Sigmas[angstrom].")?;
    let (m_1, m_2) = global.lj_masses_amu.ok_or("Input deck: missing Masses[amu].")?;
    let lennard_jones = LennardJonesPair {
        sigma_angstrom: 0.5 * (sigma_1 + sigma_2),
        epsilon_kelvin: (eps_1 * eps_2).sqrt() * CM1_TO_KELVIN,
        reduced_mass_amu: m_1 * m_2 / (m_1 + m_2),
    };
    let energy_transfer = EnergyTransferParameters {
        mean_down_at_reference_cm1: global.alpha_factor_cm1.ok_or("Input deck: missing Exponential Factor[1/cm].")?,
        reference_temperature_kelvin: ENERGY_TRANSFER_REFERENCE_TEMPERATURE,
        temperature_exponent: global.alpha_power.ok_or("Input deck: missing Exponential Power.")?,
    };
    let cutoff = global.exponent_cutoff.ok_or("Input deck: missing Exponential ExponentCutoff.")?;

    // Wells: density on the cells from the ground state to the top, averaged over the grains; leading
    // grains without states (the ground-state cell of classical rotors is empty) are dropped.
    struct WellGrid {
        zero_cell: isize,
        rho_cells: Vec<f64>,
        first_grain: isize,
        rho_grains: Vec<f64>,
    }
    let well_index: HashMap<&str, usize> =
        deck.well_order.iter().enumerate().map(|(i, name)| (name.as_str(), i)).collect();
    let mut grids = Vec::with_capacity(deck.well_order.len());
    for name in &deck.well_order {
        let species = &deck.wells[name];
        let zero_cell = grid.cell_of_energy(species.zero_energy_cm1);
        let n = cells_from(zero_cell, &format!("Well '{name}'"))?;
        let rho_cells = rrho_density_of_states(n, cell, &species_model(species)?)?;
        let ground_grain = grid.grain_of_cell(zero_cell);
        let averages = average_over_grains(&grid, &rho_cells, zero_cell, ground_grain, top_grain)?;
        let skip = averages
            .iter()
            .position(|r| *r > 0.0)
            .ok_or_else(|| format!("Well '{name}' has no states below the top of the grid."))?;
        grids.push(WellGrid {
            zero_cell,
            rho_cells,
            first_grain: ground_grain + skip as isize,
            rho_grains: averages[skip..].to_vec(),
        });
    }

    // Grain rate coefficients of well w from the cell numbers of states W (cells from `first_cell`):
    // k_g = <W>_g / (h <rho>_g), the sum of cell fluxes over the sum of cell states of the grain.
    let rates = |w: usize, first_cell: isize, w_cells: &[f64]| -> Result<Vec<f64>, String> {
        let well = &grids[w];
        let w_first_grain = grid.grain_of_cell(first_cell);
        let w_grains = if w_first_grain <= top_grain {
            average_over_grains(&grid, w_cells, first_cell, w_first_grain, top_grain)?
        } else {
            Vec::new()
        };
        Ok((0..well.rho_grains.len())
            .map(|i| {
                let g = well.first_grain + i as isize;
                if g < w_first_grain {
                    0.0
                } else {
                    w_grains[(g - w_first_grain) as usize] / (H_PLANCK_CM * well.rho_grains[i])
                }
            })
            .collect())
    };
    // Classical threshold of a channel opening at the absolute cell `threshold_cell`, as a well grain.
    let threshold_grain = |w: usize, threshold_cell: isize| -> usize {
        (grid.grain_of_cell(threshold_cell) - grids[w].first_grain).max(0) as usize
    };

    let mut channels: Vec<Vec<Channel>> = vec![Vec::new(); grids.len()];
    for barrier in &deck.barriers {
        let name = &barrier.name;
        let phase_space_core = matches!(barrier.core, MessBarrierCore::PhaseSpaceTheory { .. });
        // Tight transition state: (first cell of W, classical threshold cell, W on the cells), with the
        // Eckart tunneling convolution of Miller (1979, eqs. 8-9) when the barrier has one. Tunneling
        // reaches energies only down to `floor_cell`, the higher ground state of the two sides (below it
        // one side has no states).
        let tight_sum_of_states = |floor_cell: isize| -> Result<(isize, isize, Vec<f64>), String> {
            if phase_space_core {
                return Err(format!(
                    "Barrier '{name}' is barrierless (phase-space-theory core): give its high-pressure rate \
                     coefficient in an InverseLaplaceTransform block; the barrierless module is not yet \
                     connected for chemical activation."
                ));
            }
            let threshold = grid.cell_of_energy(barrier.rrho.zero_energy_cm1);
            let n = cells_from(threshold, &format!("Barrier '{name}'"))?;
            let model = species_model(&barrier.rrho)?;
            match (&barrier.tunneling, settings.ignore_tunneling) {
                (Some(TunnelingSpecification::Eckart { imaginary_frequency_cm1, well_depths_cm1 }), false) => {
                    // N‡ is needed m cells beyond the top (contract of eckart_tunneling_sum_of_states).
                    let m = (well_depths_cm1[0].min(well_depths_cm1[1]) / cell + 0.5).floor() as usize;
                    let ts_states = rrho_sum_of_states(n + m, cell, &model)?;
                    let (below, w_cells) = eckart_tunneling_sum_of_states(
                        &ts_states,
                        cell,
                        well_depths_cm1[0],
                        well_depths_cm1[1],
                        *imaginary_frequency_cm1,
                    );
                    debug_assert_eq!(below, m);
                    let first_cell = threshold - below as isize;
                    let drop = (floor_cell - first_cell).clamp(0, below as isize) as usize;
                    Ok((first_cell + drop as isize, threshold, w_cells[drop..].to_vec()))
                }
                (Some(TunnelingSpecification::Unsupported { model }), false) => Err(format!(
                    "Barrier '{name}': tunneling model '{model}' is not implemented (only Eckart); set \
                     ignore_tunneling to run without tunneling."
                )),
                _ => Ok((threshold, threshold, rrho_sum_of_states(n, cell, &model)?)),
            }
        };

        match (well_index.get(barrier.left.as_str()), well_index.get(barrier.right.as_str())) {
            (Some(&a), Some(&b)) => {
                if barrier.inverse_laplace_transform.is_some() {
                    return Err(format!("Barrier '{name}': an ILT channel between two wells is not supported."));
                }
                let floor_cell = grids[a].zero_cell.max(grids[b].zero_cell);
                let (first_cell, threshold, w_cells) = tight_sum_of_states(floor_cell)?;
                for (from, to) in [(a, b), (b, a)] {
                    channels[from].push(Channel {
                        name: name.clone(),
                        destination: ChannelDestination::Well { index: to },
                        threshold_grain: Some(threshold_grain(from, threshold)),
                        rate_constant_s_inv: rates(from, first_cell, &w_cells)?,
                    });
                }
            }
            (Some(&w), None) | (None, Some(&w)) => {
                let other = if well_index.contains_key(barrier.left.as_str()) { &barrier.right } else { &barrier.left };
                let bimolecular = deck.bimolecular.get(other).ok_or_else(|| {
                    format!("Barrier '{name}' connects to '{other}', which is neither a Well nor a Bimolecular species.")
                })?;
                let (first_cell, threshold, w_cells) = match &barrier.inverse_laplace_transform {
                    None => tight_sum_of_states(grids[w].zero_cell.max(grid.cell_of_energy(bimolecular.ground_energy_cm1)))?,
                    Some(ilt) => {
                        let e_inf = ilt.high_pressure_rate.activation_energy_cm1;
                        match ilt.direction {
                            IltDirection::Association => {
                                let threshold = grid.cell_of_energy(bimolecular.ground_energy_cm1 + e_inf);
                                let n = cells_from(threshold, &format!("Barrier '{name}'"))?;
                                let rho_a = rrho_density_of_states(n, cell, &species_model(&bimolecular.fragment_a)?)?;
                                let rho_b = rrho_density_of_states(n, cell, &species_model(&bimolecular.fragment_b)?)?;
                                // rho_AB(E) = sum_E' rho_A(E') rho_B(E - E') dE on the cells.
                                let rho_ab: Vec<f64> =
                                    (0..n).map(|i| (0..=i).map(|j| rho_a[j] * rho_b[i - j]).sum::<f64>() * cell).collect();
                                let mass_a: f64 = mass_vector_from_symbols_amu(&bimolecular.fragment_a.geometry_symbols)?.iter().sum();
                                let mass_b: f64 = mass_vector_from_symbols_amu(&bimolecular.fragment_b.geometry_symbols)?.iter().sum();
                                let w_cells = ilt_sum_of_states_association(
                                    &ilt.high_pressure_rate,
                                    &rho_ab,
                                    mass_a * mass_b / (mass_a + mass_b),
                                    cell,
                                )
                                .map_err(|e| format!("Barrier '{name}': {e}"))?;
                                (threshold, threshold, w_cells)
                            }
                            IltDirection::Dissociation => {
                                let well = &grids[w];
                                let threshold = grid.cell_of_energy(deck.wells[&deck.well_order[w]].zero_energy_cm1 + e_inf);
                                if threshold < grid.cell_of_energy(bimolecular.ground_energy_cm1) {
                                    return Err(format!(
                                        "Barrier '{name}': the ILT threshold E(well) + E_inf lies below the asymptote \
                                         '{other}'; k(E) would be non-zero below the dissociation energy."
                                    ));
                                }
                                let n = cells_from(threshold, &format!("Barrier '{name}'"))?;
                                let w_cells = ilt_sum_of_states_dissociation(&ilt.high_pressure_rate, &well.rho_cells[..n], cell)
                                    .map_err(|e| format!("Barrier '{name}': {e}"))?;
                                (threshold, threshold, w_cells)
                            }
                        }
                    }
                };
                channels[w].push(Channel {
                    name: name.clone(),
                    destination: ChannelDestination::Products { name: other.clone() },
                    threshold_grain: Some(threshold_grain(w, threshold)),
                    rate_constant_s_inv: rates(w, first_cell, &w_cells)?,
                });
            }
            (None, None) => {
                return Err(format!(
                    "Barrier '{name}' connects '{}' and '{}', neither of which is a Well.",
                    barrier.left, barrier.right
                ));
            }
        }
    }

    let wells: Vec<Well> = deck
        .well_order
        .iter()
        .zip(grids)
        .zip(channels)
        .map(|((name, grid), channels)| Well {
            name: name.clone(),
            bottom_offset_grains: grid.first_grain,
            density_of_states: grid.rho_grains,
            channels,
            lennard_jones: lennard_jones.clone(),
            energy_transfer: energy_transfer.clone(),
            bimolecular_sink_s_inv: deck.well_escape_rate_s_inv.get(name).copied().unwrap_or(0.0),
        })
        .collect();
    let network = ChemicalActivationNetwork { grain_width_cm1: grid.grain_width_cm1(), wells };
    network.validate()?;

    let entrance_channels = match &global.reactant_name {
        Some(reactant) => network
            .wells
            .iter()
            .enumerate()
            .flat_map(|(w, well)| {
                well.channels.iter().enumerate().filter_map(move |(c, ch)| match &ch.destination {
                    ChannelDestination::Products { name } if name == reactant => Some((w, c)),
                    _ => None,
                })
            })
            .collect(),
        None => Vec::new(),
    };

    Ok(MessChemicalActivationModel {
        network,
        temperatures_kelvin: global.temperatures_kelvin.clone(),
        pressures_torr: global.pressures_torr.clone(),
        collision_model: CollisionModel::ExponentialDown { cutoff_in_mean_down: cutoff },
        entrance_channels,
    })
}

/// RRHO counting model of a species of the deck (chirality 1; the deck's SymmetryFactor and ground
/// electronic degeneracy enter as g_e/sigma).
fn species_model(species: &MessSpeciesRrho) -> Result<SpeciesMicroModel, String> {
    let rotational_constants_cm1 = if species.geometry_symbols.is_empty() {
        return Err(format!("Species '{}' has no geometry for its rotational constants.", species.name));
    } else {
        rotational_constants_from_geometry_cm1(&species.geometry_symbols, &species.geometry_angstrom)?
    };
    Ok(SpeciesMicroModel {
        name: species.name.clone(),
        vibrational_frequencies_cm1: species.vibrational_frequencies_cm1.clone(),
        rotational_constants_cm1,
        symmetry_number: species.symmetry_factor,
        chirality_number: 1.0,
        electronic_degeneracy: species.electronic_degeneracy_ground,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::masterequation::chemical_activation_driver::{run_chemical_activation, ChemicalActivationRun, SourceSpecification};
    use crate::masterequation::chemical_activation_network::{AbsorbingBarrier, ChemicalActivationOptions, SteadyState};
    use crate::masterequation::chemical_activation_operator::isomerization_detailed_balance;
    use crate::masterequation::chemical_activation_steady_state::LinearSolver;
    use crate::masterequation::mess_input::parse_mess_input;
    use crate::masterequation::energy_graining::GrainGrid;
    use crate::masterequation::chemical_activation_network::Well;

    /// HCO + O2 (R) -> W1 (ILT association, barrierless) <-> W2 (tight) -> OH + CO2 (P, tight);
    /// W2 escapes with 1e5 s-1. Energies in kcal/mol relative to R.
    const DECK: &str = r#"
TemperatureList[K]            300. 500.
PressureList[torr]            10. 760.
EnergyStepOverTemperature     0.2
ModelEnergyLimit[kcal/mol]    40
ExcessEnergyOverTemperature   15
Reactant                      R
Model
  EnergyRelaxation
    Exponential
      Factor[1/cm]            200
      Power                   .85
      ExponentCutoff          15
    End
  CollisionFrequency
    LennardJones
      Epsilons[1/cm]          417.0  33.4
      Sigmas[angstrom]        6.5    3.9
      Masses[amu]             149    28
    End
  Bimolecular R
    Fragment HCO
      RRHO
        Geometry[angstrom] 3
        H 0.0 0.0 0.0
        C 0.0 0.0 1.1
        O 1.0 0.0 1.6
        Core RigidRotor
          SymmetryFactor 1
        End
        Frequencies[1/cm] 3
        1080 1868 2434
        ZeroEnergy[1/cm] 0
        ElectronicLevels[1/cm] 1
          0 2
      End
    Fragment O2
      RRHO
        Geometry[angstrom] 2
        O 0.0 0.0 0.0
        O 0.0 0.0 1.21
        Core RigidRotor
          SymmetryFactor 2
        End
        Frequencies[1/cm] 1
        1580
        ZeroEnergy[1/cm] 0
        ElectronicLevels[1/cm] 1
          0 3
      End
    GroundEnergy[kcal/mol] 0.0
  End
  Well W1
    Species
      RRHO
        Geometry[angstrom] 5
        C  0.0  0.0  0.0
        O  1.2  0.0  0.0
        O -0.7  1.1  0.0
        O -0.3  2.4  0.3
        H -0.6 -0.9  0.2
        Core RigidRotor
          SymmetryFactor 1
        End
        Frequencies[1/cm] 9
        250 400 600 900 1000 1100 1400 1800 2900
        ZeroEnergy[kcal/mol] -30
        ElectronicLevels[1/cm] 1
          0 2
      End
    End
  Well W2
    Escape Constant
      PseudoFirstOrderRateConstant[1/sec]  1.0E5
    End
    Species
      RRHO
        Geometry[angstrom] 5
        C  0.0  0.0  0.0
        O  1.3  0.0  0.0
        O -0.6  1.2  0.0
        O -0.2  2.3  0.5
        H  1.8  0.9  0.1
        Core RigidRotor
          SymmetryFactor 1
        End
        Frequencies[1/cm] 9
        200 350 650 850 950 1150 1300 1700 3500
        ZeroEnergy[kcal/mol] -25
        ElectronicLevels[1/cm] 1
          0 2
      End
    End
  Barrier B0 R W1
    RRHO
      Stoichiometry C1H1O3
      Core PhaseSpaceTheory
        FragmentGeometry[angstrom] 3
        H 0.0 0.0 0.0
        C 0.0 0.0 1.1
        O 1.0 0.0 1.6
        FragmentGeometry[angstrom] 2
        O 0.0 0.0 0.0
        O 0.0 0.0 1.21
        SymmetryFactor 2.0
        PotentialPrefactor[au] 2.4
        PotentialPowerExponent 6.
      End
      InverseLaplaceTransform
        Direction                   Association
        PreExponential[cm^3/s]      5.0e-12
        TemperatureExponent         0.0
        ReferenceTemperature[K]     298.0
        ActivationEnergy[kcal/mol]  0.0
      End
      Frequencies[1/cm] 4
      1080 1868 2434 1580
      ZeroEnergy[kcal/mol] 0.0
      ElectronicLevels[1/cm] 1
        0 6
    End
  Barrier B12 W1 W2
    RRHO
      Geometry[angstrom] 5
      C  0.0  0.0  0.0
      O  1.25 0.0  0.0
      O -0.65 1.15 0.0
      O -0.25 2.35 0.4
      H  0.6  1.2  0.1
      Core RigidRotor
        SymmetryFactor 1
      End
      Frequencies[1/cm] 8
      300 500 700 900 1000 1300 1700 2800
      ZeroEnergy[kcal/mol] -5
      ElectronicLevels[1/cm] 1
        0 2
    End
  Barrier B2P W2 P
    RRHO
      Geometry[angstrom] 5
      C  0.0  0.0  0.0
      O  1.3  0.0  0.0
      O -0.6  1.3  0.0
      O -0.2  2.5  0.5
      H  2.0  1.0  0.1
      Core RigidRotor
        SymmetryFactor 1
      End
      Frequencies[1/cm] 8
      250 450 600 800 1000 1200 1600 3400
      ZeroEnergy[kcal/mol] -2
      ElectronicLevels[1/cm] 1
        0 2
    End
  Bimolecular P
    Fragment OH
      RRHO
        Geometry[angstrom] 2
        O 0.0 0.0 0.0
        H 0.0 0.0 0.97
        Core RigidRotor
          SymmetryFactor 1
        End
        Frequencies[1/cm] 1
        3700
        ZeroEnergy[1/cm] 0
        ElectronicLevels[1/cm] 1
          0 2
      End
    Fragment CO2
      RRHO
        Geometry[angstrom] 3
        O 0.0 0.0 -1.16
        C 0.0 0.0 0.0
        O 0.0 0.0 1.16
        Core RigidRotor
          SymmetryFactor 2
        End
        Frequencies[1/cm] 4
        667 667 1388 2349
        ZeroEnergy[1/cm] 0
        ElectronicLevels[1/cm] 1
          0 1
      End
    GroundEnergy[kcal/mol] -20.0
  End
End
"#;


    fn tunneling_block(depth_1: f64, depth_2: f64) -> String {
        format!(
            "      Tunneling Eckart\n        ImaginaryFrequency[1/cm] 1500\n        WellDepth[kcal/mol] {depth_1}\n        WellDepth[kcal/mol] {depth_2}\n      End\n      ZeroEnergy[kcal/mol] -5"
        )
    }

    /// The deck with an Eckart tunneling block in the isomerization barrier B12 (depths 25 and 20 kcal/mol).
    fn deck_with_tunneling() -> String {
        DECK.replace("      ZeroEnergy[kcal/mol] -5", &tunneling_block(25.0, 20.0))
    }

    fn build(deck: &str, settings: &MessNetworkSettings) -> Result<MessChemicalActivationModel, String> {
        chemical_activation_model_from_mess(&parse_mess_input(deck).unwrap(), settings)
    }

    fn model() -> MessChemicalActivationModel {
        build(DECK, &MessNetworkSettings::default()).unwrap()
    }

    /// 1 cm-1 cells, grain 0.2 kT(300 K) = 41.7 cm-1 rounded to 42 cells.
    fn grid() -> GrainGrid {
        GrainGrid::new(1.0, 0.2 * KB_CM * 300.0).unwrap()
    }

    fn grain_of_energy_kcal(e_kcal: f64) -> isize {
        let g = grid();
        g.grain_of_cell(g.cell_of_energy(e_kcal / CM1_TO_KCAL))
    }

    fn first_open_grain(well: &Well, channel: usize) -> isize {
        well.channels[channel].rate_constant_s_inv.iter().position(|k| *k > 0.0).unwrap() as isize + well.bottom_offset_grains
    }

    #[test]
    fn wells_share_one_grid_up_to_the_excess_energy_above_the_highest_barrier() {
        // Top = highest barrier or asymptote of the deck (R at 0) + ExcessEnergyOverTemperature kT(T_max).
        let m = model();
        let g = grid();
        assert_eq!(m.network.grain_width_cm1, 42.0);
        let top_grain = g.grain_of_cell(g.cell_of_energy(15.0 * KB_CM * 500.0));
        for (w, e_kcal) in [(0, -30.0), (1, -25.0)] {
            let well = &m.network.wells[w];
            assert_eq!(well.bottom_offset_grains + well.grain_count() as isize - 1, top_grain, "well {w}");
            let ground = grain_of_energy_kcal(e_kcal);
            assert!((0..=1).contains(&(well.bottom_offset_grains - ground)), "well {w}");
            assert!(well.density_of_states.iter().all(|r| *r > 0.0));
        }
        assert_eq!(m.temperatures_kelvin, vec![300.0, 500.0]);
        assert_eq!(m.pressures_torr, vec![10.0, 760.0]);
    }

    #[test]
    fn grain_densities_are_cell_averages() {
        // The grain density times the grain width is the number of states counted on the 1 cm-1 cells.
        let m = model();
        let deck = parse_mess_input(DECK).unwrap();
        let well = &m.network.wells[0];
        let g = grid();
        let zero_cell = g.cell_of_energy(deck.wells["W1"].zero_energy_cm1);
        let last_cell = g.first_cell_of_grain(well.bottom_offset_grains + well.grain_count() as isize) - 1;
        let cells = rrho_density_of_states((last_cell - zero_cell + 1) as usize, 1.0, &species_model(&deck.wells["W1"]).unwrap()).unwrap();
        let states_cells: f64 = cells.iter().sum::<f64>() * 1.0;
        let states_grains: f64 = well.density_of_states.iter().sum::<f64>() * 42.0;
        assert!(((states_grains - states_cells) / states_cells).abs() < 1e-12);
    }

    #[test]
    fn isomerization_rates_obey_detailed_balance_exactly() {
        for deck in [DECK.to_string(), deck_with_tunneling()] {
            let balance = isomerization_detailed_balance(&build(&deck, &MessNetworkSettings::default()).unwrap().network);
            assert_eq!(balance.len(), 1);
            assert!(balance[0].max_relative_deviation < 1e-12, "{}", balance[0].max_relative_deviation);
        }
    }

    #[test]
    fn channels_follow_the_barriers_of_the_deck() {
        let m = model();
        let names: Vec<Vec<(String, ChannelDestination)>> = m
            .network
            .wells
            .iter()
            .map(|w| w.channels.iter().map(|c| (c.name.clone(), c.destination.clone())).collect())
            .collect();
        assert_eq!(m.network.wells[0].name, "W1");
        assert_eq!(names[0], vec![
            ("B0".to_string(), ChannelDestination::Products { name: "R".into() }),
            ("B12".to_string(), ChannelDestination::Well { index: 1 }),
        ]);
        assert_eq!(names[1], vec![
            ("B12".to_string(), ChannelDestination::Well { index: 0 }),
            ("B2P".to_string(), ChannelDestination::Products { name: "P".into() }),
        ]);
        // Without tunneling the isomerization opens in the grain of the transition state (or the next one
        // if the TS lies in the last cell of its grain), and that grain is the classical threshold.
        let ts = grain_of_energy_kcal(-5.0);
        for (w, c) in [(0, 1), (1, 0)] {
            let well = &m.network.wells[w];
            assert!((0..=1).contains(&(first_open_grain(well, c) - ts)), "well {w}");
            assert_eq!(well.channels[c].threshold_grain, Some((ts - well.bottom_offset_grains) as usize));
        }
    }

    #[test]
    fn association_ilt_forms_the_entrance_channel_at_the_asymptote() {
        let m = model();
        assert_eq!(m.entrance_channels, vec![(0, 0)]);
        let well = &m.network.wells[0];
        let asymptote = grain_of_energy_kcal(0.0);
        assert!((0..=1).contains(&(first_open_grain(well, 0) - asymptote)));
        assert_eq!(well.channels[0].threshold_grain, Some((asymptote - well.bottom_offset_grains) as usize));
        assert!(well.channels[0].rate_constant_s_inv.iter().all(|k| k.is_finite()));
    }

    #[test]
    fn collision_parameters_and_sink_come_from_the_deck() {
        let m = model();
        let w2 = &m.network.wells[1];
        assert_eq!(m.network.wells[0].bimolecular_sink_s_inv, 0.0);
        assert_eq!(w2.bimolecular_sink_s_inv, 1.0e5);
        assert!((w2.lennard_jones.sigma_angstrom - 5.2).abs() < 1e-12);
        assert!((w2.lennard_jones.epsilon_kelvin - (417.0_f64 * 33.4).sqrt() * 1.438_776_877).abs() < 1e-9);
        assert!((w2.lennard_jones.reduced_mass_amu - 149.0 * 28.0 / 177.0).abs() < 1e-12);
        assert_eq!(w2.energy_transfer.mean_down_at_reference_cm1, 200.0);
        assert_eq!(w2.energy_transfer.reference_temperature_kelvin, 300.0);
        assert_eq!(w2.energy_transfer.temperature_exponent, 0.85);
        assert_eq!(m.collision_model, CollisionModel::ExponentialDown { cutoff_in_mean_down: 15.0 });
    }

    #[test]
    fn a_phase_space_barrier_without_an_ilt_block_is_refused() {
        let start = DECK.find("      InverseLaplaceTransform").unwrap();
        let end = start + DECK[start..].find("      End\n").unwrap() + "      End\n".len();
        let deck = format!("{}{}", &DECK[..start], &DECK[end..]);
        let err = build(&deck, &MessNetworkSettings::default()).unwrap_err();
        assert!(err.contains("B0"), "{err}");
    }

    #[test]
    fn eckart_tunneling_opens_the_isomerization_below_the_transition_state() {
        // Miller (1979) eqs. 8-9: k(E) > 0 below the classical threshold, larger above it, and the
        // classical threshold stays the reference of the absorbing barrier.
        let classical = model();
        let tunneling = build(&deck_with_tunneling(), &MessNetworkSettings::default()).unwrap();
        let ts = grain_of_energy_kcal(-5.0);
        for (w, c) in [(0, 1), (1, 0)] {
            let well = &tunneling.network.wells[w];
            assert!(first_open_grain(well, c) < ts - 10, "well {w}: tunneling must open the channel well below the TS");
            assert_eq!(well.channels[c].threshold_grain, Some((ts - well.bottom_offset_grains) as usize));
            let i = (ts + 2 - well.bottom_offset_grains) as usize;
            let k_tun = well.channels[c].rate_constant_s_inv[i];
            let k_cl = classical.network.wells[w].channels[c].rate_constant_s_inv[i];
            assert!(k_tun > k_cl, "well {w}: {k_tun:e} vs {k_cl:e}");
        }
        // The tunneling region ends at the higher of the two well bottoms (20 kcal/mol below the TS).
        let well = &tunneling.network.wells[0];
        assert!(first_open_grain(well, 1) >= grain_of_energy_kcal(-25.0));
    }

    #[test]
    fn tunneling_stops_at_the_highest_ground_state_of_the_two_sides() {
        // Well depths larger than the energy differences of the deck (30 and 28 kcal/mol, while W1 and W2
        // lie 25 and 20 kcal/mol below the TS): tunneling can only reach energies at which both wells
        // have states, i.e. down to the ground state of W2 at -25 kcal/mol.
        let deck = DECK.replace("      ZeroEnergy[kcal/mol] -5", &tunneling_block(30.0, 28.0));
        let m = build(&deck, &MessNetworkSettings::default()).unwrap();
        let w2_ground = grain_of_energy_kcal(-25.0);
        for (w, c) in [(0, 1), (1, 0)] {
            assert!(first_open_grain(&m.network.wells[w], c) >= w2_ground, "well {w}");
        }
        assert!(isomerization_detailed_balance(&m.network)[0].max_relative_deviation < 1e-12);
    }

    #[test]
    fn tunneling_can_be_switched_off_explicitly() {
        let settings = MessNetworkSettings { ignore_tunneling: true, ..Default::default() };
        let off = build(&deck_with_tunneling(), &settings).unwrap();
        let classical = model();
        for w in 0..2 {
            for c in 0..2 {
                assert_eq!(off.network.wells[w].channels[c].rate_constant_s_inv, classical.network.wells[w].channels[c].rate_constant_s_inv);
            }
        }
    }

    #[test]
    fn unsupported_tunneling_models_are_refused() {
        let deck = deck_with_tunneling().replace("Tunneling Eckart", "Tunneling Read");
        assert!(build(&deck, &MessNetworkSettings::default()).is_err());
        let settings = MessNetworkSettings { ignore_tunneling: true, ..Default::default() };
        assert!(build(&deck, &settings).is_ok());
    }

    #[test]
    fn final_steady_state_of_the_c2h3_deck_reproduces_canonical_transition_state_theory() {
        // Thermal formation through the only channel: the final steady state is equilibrium and k^ca is
        // the canonical average of k(E), i.e. transition-state theory, here evaluated independently from
        // partition functions (quantum harmonic vibrations from the zero-point level, classical rigid
        // rotors, sigma = 1, equal electronic degeneracies):
        //   k_inf = (kT/h) [Q_vib(TS) Q_rot(TS)] / [Q_vib(W1) Q_rot(W1)] exp(-E0/kT).
        use crate::masterequation::chemical_activation_steady_state::LinearSolver;
        let deck = parse_mess_input(include_str!("../../examples/c2h3_chemical_activation.inp")).unwrap();
        let m = chemical_activation_model_from_mess(&deck, &MessNetworkSettings::default()).unwrap();
        let run = ChemicalActivationRun {
            temperatures_kelvin: vec![1000.0],
            pressures_torr: vec![760.0],
            options: ChemicalActivationOptions { collision_model: m.collision_model, steady_state: SteadyState::Final },
            solver: LinearSolver::BandedCholesky,
            source: SourceSpecification::ThermalEntrance { channels: m.entrance_channels.clone() },
            tolerance: 1e-8,
        };
        let k_me = run_chemical_activation(&m.network, &run).unwrap()[0].result.channels[0].ca_rate_constant_s_inv;

        let kt = KB_CM * 1000.0;
        let q_vib = |f: &[f64]| f.iter().map(|w| 1.0 / (1.0 - (-w / kt).exp())).product::<f64>();
        let q_rot = |b: &[f64]| std::f64::consts::PI.sqrt() * kt.powf(1.5) / (b[0] * b[1] * b[2]).sqrt();
        let well = &deck.wells["W1"];
        let ts = &deck.barriers[0].rrho;
        let b_well = rotational_constants_from_geometry_cm1(&well.geometry_symbols, &well.geometry_angstrom).unwrap();
        let b_ts = rotational_constants_from_geometry_cm1(&ts.geometry_symbols, &ts.geometry_angstrom).unwrap();
        let e0 = ts.zero_energy_cm1 - well.zero_energy_cm1;
        let k_tst = kt / H_PLANCK_CM * q_vib(&ts.vibrational_frequencies_cm1) * q_rot(&b_ts)
            / (q_vib(&well.vibrational_frequencies_cm1) * q_rot(&b_well))
            * (-e0 / kt).exp();
        assert!((k_me / k_tst - 1.0).abs() < 5e-3, "master equation {k_me:e} vs transition-state theory {k_tst:e}");
    }

    #[test]
    fn the_deck_runs_through_the_chemical_activation_driver() {
        for deck in [DECK.to_string(), deck_with_tunneling()] {
            let m = build(&deck, &MessNetworkSettings::default()).unwrap();
            let run = ChemicalActivationRun {
                temperatures_kelvin: m.temperatures_kelvin.clone(),
                pressures_torr: m.pressures_torr.clone(),
                options: ChemicalActivationOptions {
                    collision_model: m.collision_model,
                    steady_state: SteadyState::Intermediate { barrier: AbsorbingBarrier::default() },
                },
                solver: LinearSolver::BandedCholesky,
                source: SourceSpecification::ThermalEntrance { channels: m.entrance_channels.clone() },
                tolerance: 1e-8,
            };
            let results = run_chemical_activation(&m.network, &run).unwrap();
            assert_eq!(results.len(), 4);
            for r in &results {
                let products = r.result.channels.iter().find(|c| c.name == "B2P").unwrap();
                assert!(products.flux > 0.0);
                assert!((r.result.mass_balance - 1.0).abs() < 1e-8);
            }
        }
    }
}
