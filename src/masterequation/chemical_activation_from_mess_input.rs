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
//! are refused; `ignore_tunneling` leaves tunneling out. `eckart_tunneling` selects the transmission model of
//! the Eckart blocks: the exact Eckart probability (default) or the MESS semiclassical model
//! (`tunneling::mess_eckart_tunneling`, same convolution with its own P(E) and cutoff energy), to reproduce
//! MESS results.
//! Barrierless channels: the inverse Laplace transform of the high-pressure rate coefficient given in
//! the MarXus `InverseLaplaceTransform` block of the barrier (`barrierless::ilt::ilt_barrierless`;
//! Davies, Green, Pilling, Chem. Phys. Lett. 126, 373 (1986)):
//!   association, k_inf in cm3 s-1: W(E - E_th) from the convolved density of the two fragments
//!     (rho_AB(E) = sum_E' rho_A(E') rho_B(E - E') dE), E_th = E(asymptote) + E_inf;
//!   dissociation, k_inf in s-1: W(E - E_th) from the density of the well, E_th = E(well) + E_inf,
//!     which may not lie below the dissociation asymptote.
//! Without an ILT block, a barrier with a phase-space-theory core (`Core PhaseSpaceTheory`) uses that core
//! (`barrierless::phasespace`, `TSTLevel` T, E or EJ, default EJ as in MESS): W(E - E_barrier) is the
//! core number of states times the electronic degeneracy of the barrier, convolved with its conserved
//! vibrations (`microcanonical_builder::transition_state_sum_of_states`); for an isotropic -V0/R^n
//! potential its canonical capture rate is Georgievskii, Klippenstein, J. Chem. Phys. 122, 194103 (2005),
//! eq. 55. An ILT block, when present, takes precedence.
//!
//! Collisions: Lennard-Jones parameters combined as sigma = (sigma_1 + sigma_2)/2,
//! eps = sqrt(eps_1 eps_2) (Troe, J. Chem. Phys. 66, 4758 (1977), Sec. III), reduced mass of the two
//! Masses[amu]; <dE_down>(T) = Factor (T/300 K)^Power (Factor is the value at 300 K in this input
//! format); exponential down with the deck's ExponentCutoff. The `Escape` pseudo-first-order rate
//! constant of a well is the bimolecular sink k_c[D] (PO14 eq. 2).

use std::collections::HashMap;

use crate::barrierless::ilt::ilt_barrierless::{
    ilt_sum_of_states_association, ilt_sum_of_states_dissociation, translational_partition_constant,
};
use crate::constants::{CM1_TO_KCAL, H_PLANCK_CM, KB_CM};
use crate::tunneling::mess_eckart_tunneling::{mess_eckart_tunneling_sum_of_states, MessEckartTunneling};
use crate::tunneling::tunneling::eckart_tunneling_sum_of_states;
use crate::utils::atomic_masses::mass_vector_from_symbols_amu;

use super::chemical_activation_network::{
    Channel, ChannelDestination, ChemicalActivationNetwork, CollisionModel, EnergyTransferParameters,
    LennardJonesPair, Well,
};
use super::collisional_relaxation::CollisionIntegral;
use super::energy_graining::{average_over_grains, GrainGrid};
use super::mess_input::{
    rotational_constants_from_geometry_cm1, IltDirection, MessBarrierCore, MessDeck, MessSpeciesRrho,
    TunnelingSpecification,
};
use super::microcanonical_builder::{
    rrho_density_of_states, rrho_sum_of_states, transition_state_sum_of_states, SpeciesMicroModel,
    TransitionStateModel,
};
use crate::barrierless::phasespace::phase_space_theory::PhaseSpaceTheoryModel;
use crate::barrierless::phasespace::types::{CaptureFragment, CaptureFragmentRotorModel, PhaseSpaceTheoryInput};
use crate::rrkm::internal_rotor::{rotational_constant_cm1, HinderedRotor, ReducedMomentModel};

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
    /// Transmission model of `Tunneling Eckart` blocks: the exact Eckart probability (default) or the MESS
    /// semiclassical model (`tunneling::mess_eckart_tunneling`).
    pub eckart_tunneling: EckartTunnelingModel,
    /// Form of the Lennard-Jones collision integral Omega(2,2)* of every well.
    pub collision_integral: CollisionIntegral,
    /// Reduced moment of the internal rotors computed from a geometry (Kilpatrick-Pitzer by default).
    pub rotor_reduced_moment: ReducedMomentModel,
}

/// Setup of the internal rotors of the deck species.
#[derive(Debug, Clone, Copy, Default)]
struct RotorContext {
    reduced_moment: ReducedMomentModel,
    /// Energy above the zero of any species up to which its rotor levels must reach (cm-1).
    energy_range_cm1: f64,
}

/// Transmission model of an Eckart barrier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EckartTunnelingModel {
    /// Exact Eckart transmission probability (Miller, J. Am. Chem. Soc. 101, 6810 (1979), eq. 8).
    #[default]
    Exact,
    /// MESS semiclassical Eckart model (`tunneling::mess_eckart_tunneling`), to reproduce MESS decks.
    Mess,
}

impl Default for MessNetworkSettings {
    fn default() -> Self {
        Self {
            grain_width_cm1: None,
            cell_width_cm1: 1.0,
            top_energy_cm1: None,
            ignore_tunneling: false,
            eckart_tunneling: EckartTunnelingModel::default(),
            collision_integral: CollisionIntegral::default(),
            rotor_reduced_moment: ReducedMomentModel::default(),
        }
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
    /// High-pressure rate coefficient of the Reactant forming the wells (None without a bimolecular
    /// Reactant).
    pub entrance_high_pressure_rate: Option<EntranceHighPressureRate>,
    /// High-pressure association rate coefficient of every bimolecular species that is not Dummy and has channels, by
    /// name (the Reactant's included): the capture rate coefficients of the CSE product rows.
    pub bimolecular_high_pressure_rates: Vec<(String, EntranceHighPressureRate)>,
    /// Internal rotors of the deck species as used in the state counts: wells, barriers, then fragments.
    pub internal_rotors: Vec<InternalRotorSummary>,
    /// Partition functions of the wells (deck order), the tight barriers (deck order) and the bimolecular species
    /// (not Dummy; by name), from the cell densities.
    pub species_partition_functions: Vec<DeckSpeciesPartition>,
}

/// An internal rotor of a deck species: B, levels and the range of its potential.
#[derive(Debug, Clone)]
pub struct InternalRotorSummary {
    pub species: String,
    /// 1-based index of the rotor in its species.
    pub rotor: usize,
    pub levels: HinderedRotor,
    pub potential_minimum_cm1: f64,
    pub potential_maximum_cm1: f64,
}

/// Canonical high-pressure rate coefficient of the bimolecular Reactant A + B forming the wells through
/// all its entrance channels, from the same cell numbers of states W(E) that give k(E):
///   k_inf(T) = sum_E W(E) exp(-(E - E_AB)/kT) dE / (h C'(mu) (kT)^(3/2) Q_A(T) Q_B(T)),
/// the transition-state (or ILT) flux over the reactant partition function per unit volume, with
/// C'(mu)(kT)^(3/2) the translational partition function of the relative motion
/// (`ilt_barrierless::translational_partition_constant`) and Q_A, Q_B the internal partition functions
/// of the fragments from their cell densities. With the yields Phi_X of the steady state, the
/// bimolecular rate coefficients are k(A + B -> X) = k_inf Phi_X; for stabilization this is the
/// association rate coefficient obtained "from the rate into the absorbing barrier" (Pilling, Robertson,
/// Annu. Rev. Phys. Chem. 54, 245 (2003), eq. 44 and text).
#[derive(Debug, Clone)]
pub struct EntranceHighPressureRate {
    cell_width_cm1: f64,
    /// Sum of W over the entrance channels on the absolute cells from `first_cell`.
    w_cells: Vec<f64>,
    first_cell: isize,
    /// Absolute cell of the asymptote A + B.
    asymptote_cell: isize,
    /// Internal densities of states of A and B on cells from their ground states (per cm-1).
    density_a: Vec<f64>,
    density_b: Vec<f64>,
    reduced_mass_amu: f64,
}

impl EntranceHighPressureRate {
    /// k_inf(T) in cm3 s-1.
    pub fn rate_cm3_s(&self, temperature_kelvin: f64) -> f64 {
        let kt = KB_CM * temperature_kelvin;
        let d = self.cell_width_cm1;
        let flux: f64 = self
            .w_cells
            .iter()
            .enumerate()
            .map(|(i, w)| w * (-((self.first_cell + i as isize - self.asymptote_cell) as f64 * d) / kt).exp())
            .sum::<f64>()
            * d
            / H_PLANCK_CM;
        let partition_function =
            |rho: &[f64]| rho.iter().enumerate().map(|(i, r)| r * (-(i as f64 * d) / kt).exp()).sum::<f64>() * d;
        flux / (translational_partition_constant(self.reduced_mass_amu)
            * kt.powf(1.5)
            * partition_function(&self.density_a)
            * partition_function(&self.density_b))
    }
}

/// Kind of a deck species with a partition function.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeckSpeciesKind {
    Well,
    /// A tight transition state (RRHO core; not phase-space theory, not an inverse Laplace transform).
    Barrier,
    /// Two fragments with their relative translation.
    Bimolecular,
}

/// Partition function of a deck species from the cell densities of states that also give k(E) (cells of
/// `cell_width_cm1`, counted from the species' ground state up to the top of the grid).
#[derive(Debug, Clone)]
pub struct DeckSpeciesPartition {
    pub name: String,
    pub kind: DeckSpeciesKind,
    /// Ground energy on the deck scale (cm-1): ZeroEnergy of a well or barrier, GroundEnergy of a bimolecular species.
    pub ground_energy_cm1: f64,
    cell_width_cm1: f64,
    /// Densities of states (per cm-1) on the cells from the ground state: the species, or fragments A and B.
    densities: Vec<Vec<f64>>,
    /// Reduced mass of the fragments (amu) of a bimolecular species.
    reduced_mass_amu: Option<f64>,
}

impl DeckSpeciesPartition {
    /// Q(T) = sum rho(E) exp(-E/kT) dE from the ground state; for a bimolecular species
    /// Q_A Q_B C'(mu)(kT)^(3/2), with the relative translation per cm3 (`translational_partition_constant`).
    pub fn partition_function(&self, temperature_kelvin: f64) -> f64 {
        let kt = KB_CM * temperature_kelvin;
        let d = self.cell_width_cm1;
        let internal: f64 = self
            .densities
            .iter()
            .map(|rho| rho.iter().enumerate().map(|(i, r)| r * (-(i as f64 * d) / kt).exp()).sum::<f64>() * d)
            .product();
        match self.reduced_mass_amu {
            Some(mu) => internal * translational_partition_constant(mu) * kt.powf(1.5),
            None => internal,
        }
    }

    /// ln Q - E_ground/kT: the Boltzmann weight on the deck scale.
    fn log_weight(&self, temperature_kelvin: f64) -> f64 {
        self.partition_function(temperature_kelvin).ln() - self.ground_energy_cm1 / (KB_CM * temperature_kelvin)
    }
}

/// K_(X/Y)(T) = Q_X exp(-E_X/kT) / (Q_Y exp(-E_Y/kT)): [X]/[Y] in equilibrium, in cm3 for a well X and a bimolecular
/// species Y (the "Real equilibrium constants" of MESS).
pub fn deck_equilibrium_constant(x: &DeckSpeciesPartition, y: &DeckSpeciesPartition, temperature_kelvin: f64) -> f64 {
    (x.log_weight(temperature_kelvin) - y.log_weight(temperature_kelvin)).exp()
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
    // Internal rotors: levels from every species zero up to the top of the grid, plus the cells of an Eckart
    // tunneling convolution below a barrier.
    let lowest_cm1 = deck
        .wells
        .values()
        .map(|w| w.zero_energy_cm1)
        .chain(deck.barriers.iter().map(|b| b.rrho.zero_energy_cm1))
        .chain(deck.bimolecular.values().map(|b| b.ground_energy_cm1))
        .fold(f64::INFINITY, f64::min);
    let tunneling_cm1 = deck
        .barriers
        .iter()
        .filter_map(|b| match &b.tunneling {
            Some(TunnelingSpecification::Eckart { well_depths_cm1, .. }) => Some(well_depths_cm1[0].max(well_depths_cm1[1])),
            _ => None,
        })
        .fold(0.0, f64::max);
    let rotors = RotorContext {
        reduced_moment: settings.rotor_reduced_moment,
        energy_range_cm1: (last_cell + 1) as f64 * cell - lowest_cm1.min(top_cm1) + tunneling_cm1 + cell,
    };
    let species_model = |species: &MessSpeciesRrho| species_model(species, &rotors);

    // Collision parameters (same for all wells in this input format).
    let (eps_1, eps_2) = global.lj_epsilons_cm1.ok_or("Input deck: missing Epsilons[1/cm].")?;
    let (sigma_1, sigma_2) = global.lj_sigmas_angstrom.ok_or("Input deck: missing Sigmas[angstrom].")?;
    let (m_1, m_2) = global.lj_masses_amu.ok_or("Input deck: missing Masses[amu].")?;
    let lennard_jones = LennardJonesPair {
        sigma_angstrom: 0.5 * (sigma_1 + sigma_2),
        epsilon_kelvin: (eps_1 * eps_2).sqrt() * CM1_TO_KELVIN,
        reduced_mass_amu: m_1 * m_2 / (m_1 + m_2),
        collision_integral: settings.collision_integral,
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
    // Partition functions of the wells from their cells (barriers and bimolecular species are added below).
    let mut species_partition_functions: Vec<DeckSpeciesPartition> = deck
        .well_order
        .iter()
        .zip(&grids)
        .map(|(name, well)| DeckSpeciesPartition {
            name: name.clone(),
            kind: DeckSpeciesKind::Well,
            ground_energy_cm1: deck.wells[name].zero_energy_cm1,
            cell_width_cm1: cell,
            densities: vec![well.rho_cells.clone()],
            reduced_mass_amu: None,
        })
        .collect();

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

    let reactant = global.reactant_name.as_deref();
    // (first cell, W on the cells) of every channel from a well to each bimolecular species (not Dummy).
    let mut bimolecular_states: HashMap<String, Vec<(isize, Vec<f64>)>> = HashMap::new();
    let mut channels: Vec<Vec<Channel>> = vec![Vec::new(); grids.len()];
    for barrier in &deck.barriers {
        let name = &barrier.name;
        // Tight transition state: (first cell of W, classical threshold cell, W on the cells), with the
        // Eckart tunneling convolution of Miller (1979, eqs. 8-9) when the barrier has one. Tunneling
        // reaches energies only down to `floor_cell`, the higher ground state of the two sides (below it
        // one side has no states).
        let tight_sum_of_states = |floor_cell: isize| -> Result<(isize, isize, Vec<f64>), String> {
            // Barrierless transition state with a phase-space-theory core (MESS `Core PhaseSpaceTheory`):
            // the core number of states (`barrierless::phasespace`, TSTLevel as in the deck, default EJ),
            // times the electronic degeneracy, convolved with the conserved vibrations of the barrier
            // (`microcanonical_builder::transition_state_sum_of_states`), from the barrier energy on.
            if let MessBarrierCore::PhaseSpaceTheory {
                fragment_a_geometry_symbols,
                fragment_a_geometry_angstrom,
                fragment_b_geometry_symbols,
                fragment_b_geometry_angstrom,
                symmetry_operations,
                potential_prefactor_au,
                potential_power_exponent,
                tst_level,
            } = &barrier.core
            {
                if !barrier.rrho.internal_rotors.is_empty() {
                    return Err(format!(
                        "Barrier '{name}': internal rotors in a barrier with a PhaseSpaceTheory core are not supported; \
                         give the conserved modes as Frequencies."
                    ));
                }
                let fragment = |symbols: &Vec<String>, coordinates: &Vec<[f64; 3]>| CaptureFragment {
                    mass_amu: None,
                    rotor: CaptureFragmentRotorModel::GeometryAngstrom {
                        symbols: symbols.clone(),
                        coordinates_angstrom: coordinates.clone(),
                    },
                };
                let pst_core = PhaseSpaceTheoryModel::new(PhaseSpaceTheoryInput {
                    fragment_a: fragment(fragment_a_geometry_symbols, fragment_a_geometry_angstrom),
                    fragment_b: fragment(fragment_b_geometry_symbols, fragment_b_geometry_angstrom),
                    symmetry_operations: *symmetry_operations,
                    potential_prefactor_au: *potential_prefactor_au,
                    potential_power_exponent: *potential_power_exponent,
                    tst_level: *tst_level,
                })
                .map_err(|e| format!("Barrier '{name}': {e}"))?;
                let threshold = grid.cell_of_energy(barrier.rrho.zero_energy_cm1);
                let n = cells_from(threshold, &format!("Barrier '{name}'"))?;
                let ts = TransitionStateModel::PhaseSpaceTheoryRRHO {
                    pst_core,
                    vibrational_frequencies_cm1: barrier.rrho.vibrational_frequencies_cm1.clone(),
                    electronic_degeneracy: barrier.rrho.electronic_degeneracy_ground,
                    excited_electronic_levels: excited_levels(&barrier.rrho),
                };
                let w_cells = transition_state_sum_of_states(n, cell, &ts).map_err(|e| format!("Barrier '{name}': {e}"))?;
                return Ok((threshold, threshold, w_cells));
            }
            let threshold = grid.cell_of_energy(barrier.rrho.zero_energy_cm1);
            let n = cells_from(threshold, &format!("Barrier '{name}'"))?;
            let model = species_model(&barrier.rrho)?;
            match (&barrier.tunneling, settings.ignore_tunneling) {
                (Some(TunnelingSpecification::Eckart { imaginary_frequency_cm1, well_depths_cm1 }), false) => {
                    let (below, w_cells) = match settings.eckart_tunneling {
                        EckartTunnelingModel::Exact => {
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
                            (below, w_cells)
                        }
                        EckartTunnelingModel::Mess => {
                            // Same contract, with the MESS cutoff energy E_c (m = round(E_c / cell)).
                            let tunnel = MessEckartTunneling::new(*imaginary_frequency_cm1, *well_depths_cm1)
                                .map_err(|e| format!("Barrier '{name}': {e}"))?;
                            let m = (tunnel.cutoff_cm1() / cell).round() as usize;
                            let ts_states = rrho_sum_of_states(n + m, cell, &model)?;
                            mess_eckart_tunneling_sum_of_states(&ts_states, cell, *well_depths_cm1, *imaginary_frequency_cm1)
                                .map_err(|e| format!("Barrier '{name}': {e}"))?
                        }
                    };
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
                let is_dummy = deck.dummy_bimolecular.iter().any(|d| d == other);
                let bimolecular = deck.bimolecular.get(other);
                if bimolecular.is_none() && !is_dummy {
                    return Err(format!(
                        "Barrier '{name}' connects to '{other}', which is neither a Well nor a Bimolecular species."
                    ));
                }
                let (first_cell, threshold, w_cells) = match &barrier.inverse_laplace_transform {
                    // A Dummy product has no asymptote: k(E) is floored at the well bottom only.
                    None => tight_sum_of_states(match bimolecular {
                        Some(b) => grids[w].zero_cell.max(grid.cell_of_energy(b.ground_energy_cm1)),
                        None => grids[w].zero_cell,
                    })?,
                    Some(ilt) => {
                        let bimolecular = bimolecular.ok_or_else(|| {
                            format!("Barrier '{name}': an inverse Laplace transform needs the fragments of '{other}', not a Dummy species.")
                        })?;
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
                                let mass_a = species_mass_amu(&bimolecular.fragment_a)?;
                                let mass_b = species_mass_amu(&bimolecular.fragment_b)?;
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
                if bimolecular.is_some() {
                    bimolecular_states.entry(other.clone()).or_default().push((first_cell, w_cells.clone()));
                }
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

    // High-pressure association rate coefficient of every bimolecular species with channels (by name): the sum of W over
    // its channels and the fragment densities on the cells from the asymptote to the top.
    let mut bimolecular_high_pressure_rates: Vec<(String, EntranceHighPressureRate)> = Vec::new();
    let mut with_channels: Vec<&String> = bimolecular_states.keys().collect();
    with_channels.sort();
    for name in with_channels {
        let (pair, states) = (&deck.bimolecular[name], &bimolecular_states[name]);
        let first_cell = states.iter().map(|(c, _)| *c).min().unwrap_or(0);
        let mut w_cells = vec![0.0; (last_cell - first_cell + 1).max(0) as usize];
        for (c, w) in states {
            for (i, x) in w.iter().enumerate() {
                w_cells[(c - first_cell) as usize + i] += x;
            }
        }
        let asymptote_cell = grid.cell_of_energy(pair.ground_energy_cm1);
        let n = cells_from(asymptote_cell, &format!("The asymptote of '{name}'"))?;
        let (mass_a, mass_b) = (species_mass_amu(&pair.fragment_a)?, species_mass_amu(&pair.fragment_b)?);
        bimolecular_high_pressure_rates.push((
            name.clone(),
            EntranceHighPressureRate {
                cell_width_cm1: cell,
                w_cells,
                first_cell,
                asymptote_cell,
                density_a: rrho_density_of_states(n, cell, &species_model(&pair.fragment_a)?)?,
                density_b: rrho_density_of_states(n, cell, &species_model(&pair.fragment_b)?)?,
                reduced_mass_amu: mass_a * mass_b / (mass_a + mass_b),
            },
        ));
    }
    // The bimolecular Reactant's.
    let entrance_high_pressure_rate =
        reactant.and_then(|r| bimolecular_high_pressure_rates.iter().find(|(name, _)| name == r)).map(|(_, rate)| rate.clone());

    // Partition functions of the tight barriers and the bimolecular species (not Dummy), from their cells (the wells
    // follow their grids above).
    for barrier in &deck.barriers {
        if barrier.inverse_laplace_transform.is_some() || !matches!(barrier.core, MessBarrierCore::TightRrho) {
            continue;
        }
        let zero = barrier.rrho.zero_energy_cm1;
        let n = cells_from(grid.cell_of_energy(zero), &format!("Barrier '{}'", barrier.name))?;
        species_partition_functions.push(DeckSpeciesPartition {
            name: barrier.name.clone(),
            kind: DeckSpeciesKind::Barrier,
            ground_energy_cm1: zero,
            cell_width_cm1: cell,
            densities: vec![rrho_density_of_states(n, cell, &species_model(&barrier.rrho)?)?],
            reduced_mass_amu: None,
        });
    }
    let mut bimolecular_names: Vec<&String> = deck.bimolecular.keys().collect();
    bimolecular_names.sort();
    for name in bimolecular_names {
        let pair = &deck.bimolecular[name];
        let n = cells_from(grid.cell_of_energy(pair.ground_energy_cm1), &format!("Bimolecular '{name}'"))?;
        let (mass_a, mass_b) = (species_mass_amu(&pair.fragment_a)?, species_mass_amu(&pair.fragment_b)?);
        species_partition_functions.push(DeckSpeciesPartition {
            name: name.clone(),
            kind: DeckSpeciesKind::Bimolecular,
            ground_energy_cm1: pair.ground_energy_cm1,
            cell_width_cm1: cell,
            densities: vec![
                rrho_density_of_states(n, cell, &species_model(&pair.fragment_a)?)?,
                rrho_density_of_states(n, cell, &species_model(&pair.fragment_b)?)?,
            ],
            reduced_mass_amu: Some(mass_a * mass_b / (mass_a + mass_b)),
        });
    }

    let mut fragments: Vec<_> = deck.bimolecular.iter().collect();
    fragments.sort_by(|a, b| a.0.cmp(b.0));
    let species_with_rotors = deck
        .well_order
        .iter()
        .map(|name| &deck.wells[name])
        .chain(deck.barriers.iter().map(|b| &b.rrho))
        .chain(fragments.iter().flat_map(|(_, b)| [&b.fragment_a, &b.fragment_b]))
        .filter(|species| !species.internal_rotors.is_empty());
    let mut internal_rotors = Vec::new();
    for species in species_with_rotors {
        let model = species_model(species)?;
        for (r, (levels, rotor)) in model.internal_rotors.into_iter().zip(&species.internal_rotors).enumerate() {
            let (potential_minimum_cm1, potential_maximum_cm1) = rotor.potential.extrema();
            internal_rotors.push(InternalRotorSummary {
                species: species.name.clone(),
                rotor: r + 1,
                levels,
                potential_minimum_cm1,
                potential_maximum_cm1,
            });
        }
    }

    Ok(MessChemicalActivationModel {
        species_partition_functions,
        internal_rotors,
        entrance_high_pressure_rate,
        bimolecular_high_pressure_rates,
        network,
        temperatures_kelvin: global.temperatures_kelvin.clone(),
        pressures_torr: global.pressures_torr.clone(),
        collision_model: CollisionModel::ExponentialDown { cutoff_in_mean_down: cutoff },
        entrance_channels,
    })
}

/// Mass of a species of the deck (amu): the Atom mass or the sum of the atomic masses of its geometry.
fn species_mass_amu(species: &MessSpeciesRrho) -> Result<f64, String> {
    match species.atom_mass_amu.or(species.mass_amu) {
        Some(mass) => Ok(mass),
        None => Ok(mass_vector_from_symbols_amu(&species.geometry_symbols)?.iter().sum()),
    }
}

/// `species_model` for other users of deck species (`photoion`): internal-rotor levels up to `energy_range_cm1` above
/// the species zero, reduced moments of `reduced_moment`.
pub(crate) fn deck_species_model(
    species: &MessSpeciesRrho,
    reduced_moment: ReducedMomentModel,
    energy_range_cm1: f64,
) -> Result<SpeciesMicroModel, String> {
    species_model(species, &RotorContext { reduced_moment, energy_range_cm1 })
}

/// RRHO counting model of a species of the deck (chirality 1; the deck's SymmetryFactor and ground
/// electronic degeneracy enter as g_e/sigma), with its internal rotors: B given in the deck, or from the reduced
/// moment of the geometry; levels up to the energy range of `rotors`.
fn species_model(species: &MessSpeciesRrho, rotors: &RotorContext) -> Result<SpeciesMicroModel, String> {
    let rotational_constants_cm1 = if species.atom_mass_amu.is_some() {
        Vec::new()
    } else if let Some(constants) = &species.rotational_constants_cm1 {
        constants.clone()
    } else if species.geometry_symbols.is_empty() {
        return Err(format!("Species '{}' has no geometry for its rotational constants.", species.name));
    } else {
        rotational_constants_from_geometry_cm1(&species.geometry_symbols, &species.geometry_angstrom)?
    };
    let internal_rotors = species
        .internal_rotors
        .iter()
        .enumerate()
        .map(|(r, rotor)| {
            let context = format!("Species '{}', rotor {}", species.name, r + 1);
            let b_cm1 = match rotor.rotational_constant_cm1 {
                Some(b) => b,
                None => {
                    let masses = mass_vector_from_symbols_amu(&species.geometry_symbols)?;
                    let moment = rotors
                        .reduced_moment
                        .reduced_moment(&masses, &species.geometry_angstrom, &rotor.group, rotor.axis)
                        .map_err(|e| format!("{context}: {e}"))?;
                    rotational_constant_cm1(moment)
                }
            };
            HinderedRotor::for_energy_range(
                b_cm1,
                rotor.symmetry,
                rotor.potential.clone(),
                rotor.hamilton_size_min,
                rotor.hamilton_size_max,
                rotors.energy_range_cm1,
            )
            .map_err(|e| format!("{context}: {e}"))
        })
        .collect::<Result<Vec<_>, String>>()?;
    Ok(SpeciesMicroModel {
        name: species.name.clone(),
        vibrational_frequencies_cm1: species.vibrational_frequencies_cm1.clone(),
        rotational_constants_cm1,
        symmetry_number: species.symmetry_factor,
        chirality_number: 1.0,
        electronic_degeneracy: species.electronic_degeneracy_ground,
        internal_rotors,
        excited_electronic_levels: excited_levels(species),
    })
}

/// The electronic levels of a deck species above its ground level (energy in cm-1, degeneracy).
fn excited_levels(species: &MessSpeciesRrho) -> Vec<(f64, f64)> {
    species.electronic_levels.iter().copied().filter(|&(e, _)| e > 0.0).collect()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::masterequation::chemical_activation_driver::{run_chemical_activation, ChemicalActivationRun, SourceSpecification};
    use crate::masterequation::chemical_activation_network::{AbsorbingBarrier, ChemicalActivationOptions, SteadyState};
    use crate::masterequation::chemical_activation_operator::isomerization_detailed_balance;
    use crate::masterequation::chemical_activation_steady_state::LinearSolver;
    use crate::masterequation::mess_input::parse_mess_input;
    use crate::masterequation::energy_graining::GrainGrid;
    use crate::masterequation::chemical_activation_network::Well;
    use crate::rrkm::internal_rotor::{
        convolve_rotor_levels, potential_from_equidistant_points, reduced_moment_bond_axis, reduced_moment_pitzer, DEFAULT_BASIS_SIZE,
    };

    /// HCO + O2 (R) -> W1 (ILT association, barrierless) <-> W2 (tight) -> OH + CO2 (P, tight);
    /// W2 escapes with 1e5 s-1. Energies in kcal/mol relative to R.
    pub(crate) const DECK: &str = r#"
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

    const W1_GEOMETRY: &str = "        Geometry[angstrom] 5\n        C  0.0  0.0  0.0\n        O  1.2  0.0  0.0\n        O -0.7  1.1  0.0\n        O -0.3  2.4  0.3\n        H -0.6 -0.9  0.2\n";
    const B12_GEOMETRY: &str = "      Geometry[angstrom] 5\n      C  0.0  0.0  0.0\n      O  1.25 0.0  0.0\n      O -0.65 1.15 0.0\n      O -0.25 2.35 0.4\n      H  0.6  1.2  0.1\n";

    fn assert_same_values(a: &[f64], b: &[f64], what: &str) {
        assert_eq!(a.len(), b.len(), "{what}");
        for (i, (x, y)) in a.iter().zip(b).enumerate() {
            assert!((x - y).abs() <= 1e-12 * x.abs().max(y.abs()), "{what}[{i}]: {x:e} vs {y:e}");
        }
    }

    #[test]
    fn rotational_constants_and_a_mass_give_the_same_network_as_the_geometry() {
        // A species may be given by its rotational constants and its mass instead of a geometry (MarXus extension of
        // the deck format, for species data given that way): the same constants give the same densities of states,
        // sums of states and k(E).
        let reference = model();
        let parsed = parse_mess_input(DECK).unwrap();
        let w1 = &parsed.wells["W1"];
        let ts = &parsed.barriers.iter().find(|b| b.name == "B12").unwrap().rrho;
        let constants = |sp: &MessSpeciesRrho, pad: &str| {
            let b = rotational_constants_from_geometry_cm1(&sp.geometry_symbols, &sp.geometry_angstrom).unwrap();
            let mass: f64 = mass_vector_from_symbols_amu(&sp.geometry_symbols).unwrap().iter().sum();
            format!("{pad}RotationalConstants[1/cm] 3\n{pad}  {} {} {}\n{pad}Mass[amu] {mass}\n", b[0], b[1], b[2])
        };
        assert!(DECK.contains(W1_GEOMETRY) && DECK.contains(B12_GEOMETRY));
        let deck = DECK.replace(W1_GEOMETRY, &constants(w1, "        ")).replace(B12_GEOMETRY, &constants(ts, "      "));
        let replaced = build(&deck, &MessNetworkSettings::default()).unwrap();
        let w = reference.network.wells.iter().position(|w| w.name == "W1").unwrap();
        assert_same_values(
            &reference.network.wells[w].density_of_states,
            &replaced.network.wells[w].density_of_states,
            "rho(W1)",
        );
        for (a, b) in reference.network.wells[w].channels.iter().zip(&replaced.network.wells[w].channels) {
            assert_same_values(&a.rate_constant_s_inv, &b.rate_constant_s_inv, &a.name);
        }
    }

    /// W1 of the test deck with a hindered rotor (group O O about the C-O axis) after its Core block; `extra` is added
    /// inside the Rotor block.
    fn deck_with_w1_rotor(extra: &str) -> String {
        let core = "          SymmetryFactor 1\n        End\n        Frequencies[1/cm] 9\n        250 400";
        assert_eq!(DECK.matches(core).count(), 1);
        DECK.replace(
            core,
            &format!(
                "          SymmetryFactor 1\n        End\n        Rotor Hindered\n          Group 3 4\n          Axis 1 2\n          \
                 Symmetry 1\n          Potential[kcal/mol] 4\n          0. 1.5 3.0 1.5\n{extra}        End\n        \
                 Frequencies[1/cm] 9\n        250 400"
            ),
        )
    }

    #[test]
    fn a_hindered_rotor_of_the_deck_is_convolved_into_the_well_density() {
        // The well density of states is the rotor-free count convolved with the levels of the rotor, whose B comes from
        // the reduced moment of the geometry (Pitzer by default, or the bond-axis definition) or is given in the deck.
        let parsed = parse_mess_input(DECK).unwrap();
        let w1 = &parsed.wells["W1"];
        let masses = mass_vector_from_symbols_amu(&w1.geometry_symbols).unwrap();
        let points: Vec<f64> = [0.0, 1.5, 3.0, 1.5].iter().map(|v| v / CM1_TO_KCAL).collect();
        let potential = potential_from_equidistant_points(&points).unwrap();
        let pitzer = rotational_constant_cm1(reduced_moment_pitzer(&masses, &w1.geometry_angstrom, &[2, 3], (0, 1)).unwrap());
        let bond_axis = rotational_constant_cm1(reduced_moment_bond_axis(&masses, &w1.geometry_angstrom, &[2, 3], (0, 1)).unwrap());
        let reference = model();
        let w = reference.network.wells.iter().position(|w| w.name == "W1").unwrap();
        let g = grid();
        let well = &reference.network.wells[w];
        let zero_cell = g.cell_of_energy(w1.zero_energy_cm1);
        let last_cell = g.first_cell_of_grain(well.bottom_offset_grains + well.grain_count() as isize) - 1;
        let n = (last_cell - zero_cell + 1) as usize;
        let rotor_free = rrho_density_of_states(n, 1.0, &species_model(w1, &RotorContext::default()).unwrap()).unwrap();
        for (extra, reduced_moment, b) in [
            ("", ReducedMomentModel::Pitzer, pitzer),
            ("", ReducedMomentModel::BondAxis, bond_axis),
            ("          RotationalConstant[1/cm] 2.5\n", ReducedMomentModel::Pitzer, 2.5),
        ] {
            let settings = MessNetworkSettings { rotor_reduced_moment: reduced_moment, ..Default::default() };
            let with_rotor = build(&deck_with_w1_rotor(extra), &settings).unwrap();
            let rotor = HinderedRotor::new(b, 1, potential.clone(), DEFAULT_BASIS_SIZE).unwrap();
            let expected: f64 = convolve_rotor_levels(&rotor_free, &rotor.levels_above_ground_cm1, 1.0).iter().sum();
            let well = &with_rotor.network.wells[w];
            let last = g.first_cell_of_grain(well.bottom_offset_grains + well.grain_count() as isize) - 1;
            assert_eq!(last, last_cell);
            let states: f64 = well.density_of_states.iter().sum::<f64>() * 42.0;
            assert!((states / expected - 1.0).abs() < 1e-12, "B = {b}: {states:e} vs {expected:e}");
        }
        assert!((pitzer / bond_axis - 1.0).abs() > 1e-3, "{pitzer} vs {bond_axis}");
    }

    #[test]
    fn the_network_model_lists_the_internal_rotors_as_counted() {
        // The rotors of the deck species, with B, the potential range and the levels used in the state counts.
        assert!(model().internal_rotors.is_empty());
        let parsed = parse_mess_input(DECK).unwrap();
        let w1 = &parsed.wells["W1"];
        let masses = mass_vector_from_symbols_amu(&w1.geometry_symbols).unwrap();
        let b = rotational_constant_cm1(reduced_moment_pitzer(&masses, &w1.geometry_angstrom, &[2, 3], (0, 1)).unwrap());
        let points: Vec<f64> = [0.0, 1.5, 3.0, 1.5].iter().map(|v| v / CM1_TO_KCAL).collect();
        let potential = potential_from_equidistant_points(&points).unwrap();
        let rotor = HinderedRotor::new(b, 1, potential.clone(), DEFAULT_BASIS_SIZE).unwrap();
        let m = build(&deck_with_w1_rotor(""), &MessNetworkSettings::default()).unwrap();
        let [summary] = &m.internal_rotors[..] else { panic!("{:?}", m.internal_rotors.len()) };
        assert_eq!((summary.species.as_str(), summary.rotor), ("W1", 1));
        assert_eq!(summary.levels.rotational_constant_cm1, b);
        assert_eq!(summary.levels.ground_energy_cm1, rotor.ground_energy_cm1);
        assert_eq!(summary.levels.levels_above_ground_cm1[..10], rotor.levels_above_ground_cm1[..10]);
        assert_eq!((summary.potential_minimum_cm1, summary.potential_maximum_cm1), potential.extrema());
    }

    #[test]
    fn the_excited_electronic_levels_of_the_deck_reach_the_species_model() {
        let deck = DECK.replace(
            "        ZeroEnergy[kcal/mol] -30\n        ElectronicLevels[1/cm] 1\n          0 2\n",
            "        ZeroEnergy[kcal/mol] -30\n        ElectronicLevels[1/cm] 2\n          0 2\n          139.7 2\n",
        );
        assert_ne!(deck, DECK);
        let parsed = parse_mess_input(&deck).unwrap();
        let model = species_model(&parsed.wells["W1"], &RotorContext::default()).unwrap();
        assert_eq!((model.electronic_degeneracy, model.excited_electronic_levels.clone()), (2.0, vec![(139.7, 2.0)]));
        assert!(build(&deck, &MessNetworkSettings::default()).is_ok());
    }

    #[test]
    fn a_rotor_in_a_phase_space_theory_barrier_is_an_error() {
        let start = DECK.find("      InverseLaplaceTransform\n").unwrap();
        let end = DECK[start..].find("      End\n").unwrap() + start + "      End\n".len();
        let rotor = "      Rotor Hindered\n        Group 1\n        Axis 2 3\n        RotationalConstant[1/cm] 5\n        \
Potential[kcal/mol] 2\n        0. 1.\n      End\n";
        let deck = format!("{}{rotor}{}", &DECK[..start], &DECK[end..]);
        let err = build(&deck, &MessNetworkSettings::default()).unwrap_err();
        assert!(err.contains("B0") && err.contains("PhaseSpaceTheory"), "{err}");
    }

    #[test]
    fn a_dummy_bimolecular_product_is_a_sink_without_molecular_data() {
        // MESS's `Dummy` bimolecular species: a product reached through a tight barrier, whose molecular data are
        // never needed. The network is the same as with the full product.
        let start = DECK.find("  Bimolecular P\n").unwrap();
        let end = DECK[start..].find("    GroundEnergy[kcal/mol] -20.0\n  End\n").unwrap() + start
            + "    GroundEnergy[kcal/mol] -20.0\n  End\n".len();
        let deck = format!("{}  Bimolecular P\n    Dummy\n{}", &DECK[..start], &DECK[end..]);
        let parsed = parse_mess_input(&deck).unwrap();
        assert_eq!(parsed.dummy_bimolecular, vec!["P".to_string()]);
        assert!(!parsed.bimolecular.contains_key("P"));
        let reference = model();
        let dummy = build(&deck, &MessNetworkSettings::default()).unwrap();
        let w = reference.network.wells.iter().position(|w| w.name == "W2").unwrap();
        let channel = |m: &MessChemicalActivationModel| {
            m.network.wells[w].channels.iter().find(|c| c.name == "B2P").unwrap().clone()
        };
        assert_eq!(channel(&dummy).destination, ChannelDestination::Products { name: "P".into() });
        assert_same_values(&channel(&reference).rate_constant_s_inv, &channel(&dummy).rate_constant_s_inv, "k(B2P)");
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
        let cells = rrho_density_of_states((last_cell - zero_cell + 1) as usize, 1.0, &species_model(&deck.wells["W1"], &RotorContext::default()).unwrap()).unwrap();
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
    fn a_phase_space_barrier_without_an_ilt_block_uses_the_phase_space_core() {
        // B0 (R -> W1) without its ILT block: the phase-space-theory core (default TSTLevel EJ) gives the
        // transition-state sum of states. In the test deck the core fragments are the reactant fragments
        // (same geometries), the TS frequencies are the fragment frequencies, g_TS/(g_A g_B) = 6/(2*3) = 1 and
        // sigma_A sigma_B / sigma_PST = 1*2/2 = 1, so the entrance high-pressure rate coefficient is the
        // isotropic capture rate of Georgievskii, Klippenstein, J. Chem. Phys. 122, 194103 (2005), eq. 55:
        //   k(T) = (8 pi)^(1/2) ((n-2)/2)^(2/n) Gamma(1 - 2/n) mu^(-1/2) V0^(2/n) T^(1/2 - 2/n)  (atomic units).
        use crate::numeric::lanczos_gamma::gamma_func;
        use crate::utils::atomic_masses::mass_vector_from_symbols_amu;
        let start = DECK.find("      InverseLaplaceTransform").unwrap();
        let end = start + DECK[start..].find("      End\n").unwrap() + "      End\n".len();
        let deck = format!("{}{}", &DECK[..start], &DECK[end..]);
        let m = build(&deck, &MessNetworkSettings::default()).unwrap();
        let rate = m.entrance_high_pressure_rate.as_ref().expect("the deck has a bimolecular Reactant");
        let mass = |symbols: &[&str]| -> f64 {
            mass_vector_from_symbols_amu(&symbols.iter().map(|s| s.to_string()).collect::<Vec<_>>()).unwrap().iter().sum()
        };
        let (m_a, m_b) = (mass(&["H", "C", "O"]), mass(&["O", "O"]));
        let mu_au = m_a * m_b / (m_a + m_b) * crate::constants::AMU_TO_ELECTRON_MASS;
        let (v0, n) = (2.4_f64, 6.0_f64);
        // Atomic units (CODATA 2018): k_B = 3.166811563e-6 Hartree/K, a0 = 0.529177210903e-8 cm,
        // t_au = 2.4188843265857e-17 s.
        let au_to_cm3_s = 0.529177210903e-8_f64.powi(3) / 2.4188843265857e-17;
        for t in [300.0_f64, 500.0] {
            let kt_au = 3.166811563e-6 * t;
            let k_capture = (8.0 * std::f64::consts::PI).sqrt() * ((n - 2.0) / 2.0).powf(2.0 / n) * gamma_func(1.0 - 2.0 / n)
                * mu_au.powf(-0.5) * v0.powf(2.0 / n) * kt_au.powf(0.5 - 2.0 / n) * au_to_cm3_s;
            let k = rate.rate_cm3_s(t);
            assert!((k / k_capture - 1.0).abs() < 5e-3, "T = {t}: {k:e} vs capture {k_capture:e}");
        }
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
    fn the_mess_eckart_model_scales_the_channel_by_its_canonical_tunneling_factor() {
        // B12 (W1 <-> W2) with Eckart tunneling (1500 cm-1, depths 25 and 20 kcal/mol): the high-pressure rate
        // coefficient, the Boltzmann average of k(E) over W1, is kappa(T) times the classical one, so the ratio
        // of the two models is kappa_MESS/kappa_exact (tunneling::mess_eckart_tunneling, tunneling::eckart).
        use crate::tunneling::mess_eckart_tunneling::MessEckartTunneling;
        use crate::tunneling::tunneling::eckart;
        let exact = build(&deck_with_tunneling(), &MessNetworkSettings::default()).unwrap();
        let settings = MessNetworkSettings { eckart_tunneling: EckartTunnelingModel::Mess, ..Default::default() };
        let mess = build(&deck_with_tunneling(), &settings).unwrap();
        let k_inf = |m: &MessChemicalActivationModel, t: f64| -> f64 {
            let well = &m.network.wells[0];
            let channel = well.channels.iter().find(|c| c.name == "B12").unwrap();
            let kt = KB_CM * t;
            let f: Vec<f64> = (0..well.grain_count())
                .map(|i| well.density_of_states[i] * (-(i as f64) * m.network.grain_width_cm1 / kt).exp())
                .collect();
            channel.rate_constant_s_inv.iter().zip(&f).map(|(k, f)| k * f).sum::<f64>() / f.iter().sum::<f64>()
        };
        let depths = [25.0 / CM1_TO_KCAL, 20.0 / CM1_TO_KCAL];
        for t in [300.0, 500.0] {
            let kt = KB_CM * t;
            let kappa_exact = eckart(1.0 / kt, 1500.0, depths[0], depths[1], 0.1, depths[0] + 40.0 * kt);
            let kappa_mess = MessEckartTunneling::new(1500.0, depths).unwrap().canonical_factor(kt);
            let ratio = k_inf(&mess, t) / k_inf(&exact, t);
            assert!((ratio / (kappa_mess / kappa_exact) - 1.0).abs() < 1e-2, "T {t}: {ratio} vs {}", kappa_mess / kappa_exact);
            assert!(ratio < 1.0);
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
    fn entrance_high_pressure_rate_reproduces_the_ilt_input() {
        // The ILT inverts k_inf(T) = A (T/T_ref)^n exp(-E_inf/kT); the canonical average of the same cell
        // numbers of states gives it back (test deck: A = 5e-12 cm3/s, n = 0, E_inf = 0).
        let m = model();
        let rate = m.entrance_high_pressure_rate.as_ref().expect("the deck has a bimolecular Reactant");
        for t in [300.0, 500.0] {
            let k = rate.rate_cm3_s(t);
            assert!((k / 5.0e-12 - 1.0).abs() < 1e-2, "T = {t}: {k:e}");
        }
    }

    #[test]
    fn entrance_high_pressure_rate_of_a_tight_transition_state_is_transition_state_theory() {
        // C2H3 deck: k_inf(H + C2H2 -> C2H3) = (kT/h) Q‡ exp(-(E‡ - E_asym)/kT) / (C'(mu) (kT)^(3/2) Q_C2H2 Q_H),
        // quantum harmonic vibrations, classical rigid rotors (C2H2 linear, Q = kT/(sigma B), sigma = 2),
        // electronic degeneracies TS 2, C2H2 1, H 2.
        use crate::barrierless::ilt::ilt_barrierless::translational_partition_constant;
        let deck = parse_mess_input(include_str!("../../examples/c2h3_chemical_activation.inp")).unwrap();
        let m = chemical_activation_model_from_mess(&deck, &MessNetworkSettings::default()).unwrap();
        let t = 1000.0;
        let kt = KB_CM * t;
        let q_vib = |f: &[f64]| f.iter().map(|w| 1.0 / (1.0 - (-w / kt).exp())).product::<f64>();
        let ts = &deck.barriers[0].rrho;
        let b_ts = rotational_constants_from_geometry_cm1(&ts.geometry_symbols, &ts.geometry_angstrom).unwrap();
        let q_ts = 2.0 * q_vib(&ts.vibrational_frequencies_cm1) * std::f64::consts::PI.sqrt() * kt.powf(1.5)
            / (b_ts[0] * b_ts[1] * b_ts[2]).sqrt();
        let pair = &deck.bimolecular["P1"];
        let c2h2 = &pair.fragment_a;
        let b_c2h2 = rotational_constants_from_geometry_cm1(&c2h2.geometry_symbols, &c2h2.geometry_angstrom).unwrap();
        let q_c2h2 = q_vib(&c2h2.vibrational_frequencies_cm1) * kt / (2.0 * b_c2h2[0]);
        let q_h = 2.0;
        let mass = |symbol: &str| crate::utils::atomic_masses::atomic_mass_amu(symbol).unwrap();
        let (m_c2h2, m_h) = (2.0 * mass("C") + 2.0 * mass("H"), mass("H"));
        let mu = m_c2h2 * m_h / (m_c2h2 + m_h);
        let k_tst = kt / H_PLANCK_CM * q_ts * (-(ts.zero_energy_cm1 - pair.ground_energy_cm1) / kt).exp()
            / (translational_partition_constant(mu) * kt.powf(1.5) * q_c2h2 * q_h);
        let k = m.entrance_high_pressure_rate.as_ref().unwrap().rate_cm3_s(t);
        assert!((k / k_tst - 1.0).abs() < 5e-3, "{k:e} vs transition-state theory {k_tst:e}");
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

    /// A + B <=> AB given by rotational constants and masses (m_AB = m_A + m_B), with excited electronic levels in AB and
    /// B: the deck species of the agreement test of the two equilibrium-constant routes (reports/equilibrium_constants.md).
    const ASSOCIATION_DECK: &str = "
TemperatureList[K]                  500.
PressureList[atm]                   1.0
EnergyStepOverTemperature           0.1
ModelEnergyLimit[kcal/mol]          150
Model
  EnergyRelaxation
    Exponential
      Factor[1/cm]                  200
      Power                         .85
      ExponentCutoff                15
    End
  CollisionFrequency
    LennardJones
      Epsilons[1/cm]                6.95   292.0
      Sigmas[angstrom]              2.55   4.36
      Masses[amu]                   4.0    62.0
    End
  Well AB
    Species
      RRHO
        RotationalConstants[1/cm]   3
          0.5 0.2 0.15
        Mass[amu]                   62.0
        Core RigidRotor
          SymmetryFactor            1
        End
        Frequencies[1/cm]           2
          800.0 1200.0
        ZeroEnergy[1/cm]            -8000
        ElectronicLevels[1/cm]      2
          0    2
          300  2
      End
    End
  Barrier TS AB P
    RRHO
      RotationalConstants[1/cm]     3
        0.45 0.18 0.14
      Mass[amu]                     62.0
      Core RigidRotor
        SymmetryFactor              1
      End
      Frequencies[1/cm]             1
        700.0
      ZeroEnergy[1/cm]              500
      ElectronicLevels[1/cm]        1
        0  2
    End
  Bimolecular P
    Fragment A
      RRHO
        RotationalConstants[1/cm]   3
          2.0 1.0 0.5
        Mass[amu]                   30.0
        Core RigidRotor
          SymmetryFactor            1
        End
        Frequencies[1/cm]           1
          1000.0
        ZeroEnergy[1/cm]            0
        ElectronicLevels[1/cm]      1
          0  2
      End
    Fragment B
      RRHO
        RotationalConstants[1/cm]   1
          1.5
        Mass[amu]                   32.0
        Core RigidRotor
          SymmetryFactor            1
        End
        Frequencies[1/cm]           1
          1500.0
        ZeroEnergy[1/cm]            0
        ElectronicLevels[1/cm]      2
          0     3
          1000  2
      End
    GroundEnergy[1/cm]              0
  End
End
";

    #[test]
    fn the_deck_partition_functions_give_the_equilibrium_constant_of_the_thermochemistry() {
        // Both routes describe A + B <=> AB with harmonic vibrations, classical rigid rotors, the same masses and
        // electronic levels; they differ only by the 1 cm-1 cell counting of the deck densities.
        use crate::molecule::{MolType, MoleculeBuilder};
        use crate::thermal::equilibrium::equilibrium_constant_from_thermochemistry;
        let deck = parse_mess_input(ASSOCIATION_DECK).unwrap();
        let m = chemical_activation_model_from_mess(&deck, &MessNetworkSettings::default()).unwrap();
        let species = |name: &str| m.species_partition_functions.iter().find(|s| s.name == name).unwrap();
        let molecule = |name: &str, freq: Vec<f64>, brot: Vec<f64>, mass: f64, levels: Vec<(f64, f64)>, dh0: f64| {
            MoleculeBuilder::new(name.to_string(), MolType::mol)
                .freq(freq)
                .brot(brot)
                .mass(mass)
                .multi(levels[0].1)
                .electronic_levels(levels)
                .dh0(dh0)
                .build()
        };
        assert_eq!(species("AB").kind, DeckSpeciesKind::Well);
        assert_eq!(species("TS").kind, DeckSpeciesKind::Barrier);
        assert_eq!(species("P").kind, DeckSpeciesKind::Bimolecular);
        for t in [500.0, 1500.0] {
            let mut a = molecule("A", vec![1000.0], vec![2.0, 1.0, 0.5], 30.0, vec![(0.0, 2.0)], 0.0);
            let mut b = molecule("B", vec![1500.0], vec![1.5], 32.0, vec![(0.0, 3.0), (1000.0, 2.0)], 0.0);
            let mut ab = molecule("AB", vec![800.0, 1200.0], vec![0.5, 0.2, 0.15], 62.0, vec![(0.0, 2.0), (300.0, 2.0)], -8000.0);
            let thermo =
                equilibrium_constant_from_thermochemistry(&mut [(1.0, &mut a), (1.0, &mut b)], &mut [(1.0, &mut ab)], t, 1.0e5, 0.0)
                    .unwrap();
            // The cells hold the states in ((i-1) dE, i dE] at E = i dE (`microcanonical_builder`): the classical-rotor
            // states lie half a cell low in the Boltzmann factor, so every species with classical rotors has Q low by
            // exp(-dE/2kT), and Q_AB/(Q_A Q_B) keeps one factor exp(+dE/2kT) (0.14% at 500 K with 1 cm-1 cells).
            let half_cell = (0.5 / (KB_CM * t)).exp();
            let k = deck_equilibrium_constant(species("AB"), species("P"), t);
            assert!(
                (k / (thermo.k_c * half_cell) - 1.0).abs() < 1e-4,
                "T = {t}: deck {k:e} vs thermochemistry {:e} cm3 x exp(dE/2kT)",
                thermo.k_c
            );
            // and the reverse ratio
            assert!((deck_equilibrium_constant(species("P"), species("AB"), t) * k - 1.0).abs() < 1e-12);
            // the transition state: Q from the cells = q_el q_rot q_vib of the thermochemistry, x exp(-dE/2kT)
            let mut ts = molecule("TS", vec![700.0], vec![0.45, 0.18, 0.14], 62.0, vec![(0.0, 2.0)], 500.0);
            ts.eval_all_therm_func(t, 1.0e5, 0.0);
            let q_ts = ts.thermo.pfelec * ts.thermo.pfrot * ts.thermo.pfvib;
            assert!((species("TS").partition_function(t) * half_cell / q_ts - 1.0).abs() < 1e-4, "T = {t}");
            assert_eq!(species("TS").ground_energy_cm1, 500.0);
        }
    }


    #[test]
    fn every_bimolecular_species_has_its_high_pressure_association_rate() {
        // Product P of the association deck (no Reactant), formed from AB through the tight TS:
        //   k_inf(A + B -> AB) = (kT/h) Q_TS exp(-(E_TS - E_P)/kT) / Q_P   (Q_P per cm3 with the relative translation),
        // on the cells: sum_E W(E) exp(-E/kT) dE = Q_TS dE/(1 - exp(-dE/kT)) = Q_TS kT (1 + dE/2kT + ...).
        let deck = parse_mess_input(ASSOCIATION_DECK).unwrap();
        let m = chemical_activation_model_from_mess(&deck, &MessNetworkSettings::default()).unwrap();
        assert!(m.entrance_high_pressure_rate.is_none());
        let names: Vec<&str> = m.bimolecular_high_pressure_rates.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(names, vec!["P"]);
        let species = |name: &str| m.species_partition_functions.iter().find(|s| s.name == name).unwrap();
        for t in [500.0, 1500.0] {
            let kt = KB_CM * t;
            let tst = kt / H_PLANCK_CM * species("TS").partition_function(t) * (-(500.0 - 0.0) / kt).exp()
                / species("P").partition_function(t);
            let cells = 1.0 / kt / (1.0 - (-1.0 / kt).exp());
            let k = m.bimolecular_high_pressure_rates[0].1.rate_cm3_s(t);
            assert!((k / (tst * cells) - 1.0).abs() < 1e-9, "T = {t}: {k:e} vs {:e}", tst * cells);
        }
        // With a Reactant, its entry is the entrance rate.
        let deck = parse_mess_input(include_str!("../../examples/c2h3_chemical_activation.inp")).unwrap();
        let m = chemical_activation_model_from_mess(&deck, &MessNetworkSettings::default()).unwrap();
        let (name, rate) = &m.bimolecular_high_pressure_rates[0];
        assert_eq!(name, "P1");
        assert_eq!(rate.rate_cm3_s(1000.0), m.entrance_high_pressure_rate.as_ref().unwrap().rate_cm3_s(1000.0));
    }
}
