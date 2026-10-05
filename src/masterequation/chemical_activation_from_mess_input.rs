//! Chemical-activation network from an input deck in the MESS input format (`mess_input.rs`).
//!
//! Energy grid. All wells share one grain width dE: the value given in the settings, otherwise
//! EnergyStepOverTemperature x k_B x (lowest temperature of the deck), the finest grain of the deck.
//! Every energy of the deck (well ground states, transition states, bimolecular asymptotes) is
//! rounded once to whole grains on the absolute scale of the deck; thresholds are differences of these
//! integers, so that the two directions of an isomerization use the same transition-state grain and
//! rho_a k_ab = rho_b k_ba = W‡/h holds exactly (detailed balance, `chemical_activation_operator.rs`).
//! All grids end at the same absolute top (the settings, otherwise ModelEnergyLimit). Each well grid
//! starts at the lowest grain that contains states (with classical rotors the grain of the ground
//! state itself is empty, G(0) = 0).
//!
//! Microcanonical rate coefficients. Tight transition states: k(E) = W‡(E - E0)/(h rho(E)) (RRKM;
//! PO14 eq. 9; Forst 1973 Sec. 4.5 for the symmetry and degeneracy factors carried by W‡ and rho).
//! Barrierless channels: the inverse Laplace transform of the high-pressure rate coefficient given in
//! the MarXus `InverseLaplaceTransform` block of the barrier (`barrierless::ilt::ilt_barrierless`;
//! Davies, Green, Pilling, Chem. Phys. Lett. 126, 373 (1986)):
//!   association, k_inf in cm3 s-1: W(E - E_th) from the convolved density of the two fragments
//!     (rho_AB(E) = sum_E' rho_A(E') rho_B(E - E') dE), E_th = E(asymptote) + E_inf;
//!   dissociation, k_inf in s-1: W(E - E_th) from the density of the well, E_th = E(well) + E_inf,
//!     which may not lie below the dissociation asymptote.
//! A phase-space-theory core without an ILT block is refused: the barrierless module is not yet
//! connected for chemical activation. Tunneling corrections of k(E) are not implemented yet; a barrier
//! with a `Tunneling` block is refused unless `ignore_tunneling` is set.
//!
//! Collisions: Lennard-Jones parameters combined as sigma = (sigma_1 + sigma_2)/2,
//! eps = sqrt(eps_1 eps_2) (Troe, J. Chem. Phys. 66, 4758 (1977), Sec. III), reduced mass of the two
//! Masses[amu]; <dE_down>(T) = Factor (T/300 K)^Power (Factor is the value at 300 K in this input
//! format); exponential down with the deck's ExponentCutoff. The `Escape` pseudo-first-order rate
//! constant of a well is the bimolecular sink k_c[D] (PO14 eq. 2).

use std::collections::HashMap;

use crate::barrierless::ilt::ilt_barrierless::{ilt_sum_of_states_association, ilt_sum_of_states_dissociation};
use crate::constants::{CM1_TO_KCAL, H_PLANCK_CM, KB_CM};
use crate::utils::atomic_masses::mass_vector_from_symbols_amu;

use super::chemical_activation_network::{
    Channel, ChannelDestination, ChemicalActivationNetwork, CollisionModel, EnergyTransferParameters,
    LennardJonesPair, Well,
};
use super::mess_input::{rotational_constants_from_geometry_cm1, IltDirection, MessBarrierCore, MessDeck, MessSpeciesRrho};
use super::microcanonical_builder::{rrho_density_of_states, rrho_sum_of_states, SpeciesMicroModel};

/// Conversion factor from Epsilons[1/cm] to K (hc/k_B in cm K).
const CM1_TO_KELVIN: f64 = 1.438_776_877;
/// Reference temperature of the energy-transfer Factor in this input format, K.
const ENERGY_TRANSFER_REFERENCE_TEMPERATURE: f64 = 300.0;

/// Choices that the input deck does not fix.
#[derive(Debug, Clone, Default)]
pub struct MessNetworkSettings {
    /// Grain width in cm-1; None: EnergyStepOverTemperature x k_B x lowest temperature.
    pub grain_width_cm1: Option<f64>,
    /// Absolute top of all grids in cm-1 on the energy scale of the deck; None: ModelEnergyLimit.
    pub top_energy_cm1: Option<f64>,
    /// Accept barriers with a `Tunneling` block although tunneling is not included in k(E).
    pub ignore_tunneling: bool,
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

    // Common grain and top of the grids.
    let lowest_temperature = global.temperatures_kelvin.iter().cloned().fold(f64::INFINITY, f64::min);
    let d_e = match settings.grain_width_cm1 {
        Some(width) => width,
        None => {
            global
                .energy_step_over_temperature
                .ok_or("Input deck: no EnergyStepOverTemperature and no grain width in the settings.")?
                * KB_CM
                * lowest_temperature
        }
    };
    if !(d_e > 0.0) || !d_e.is_finite() {
        return Err(format!("Invalid grain width {d_e} cm-1."));
    }
    let top_cm1 = settings
        .top_energy_cm1
        .or(global.model_energy_limit_kcal_mol.map(|e| e / CM1_TO_KCAL))
        .ok_or("Input deck: no ModelEnergyLimit[kcal/mol] and no top energy in the settings.")?;
    let top_grain = (top_cm1 / d_e).floor() as isize;
    let grain_of = |energy_cm1: f64| (energy_cm1 / d_e).round() as isize;
    // Number of grains from absolute grain `from` up to the common top.
    let grains_from = |from: isize, what: &str| -> Result<usize, String> {
        if top_grain - from < 1 {
            return Err(format!("{what} lies at or above the top of the energy grid ({top_cm1} cm-1)."));
        }
        Ok((top_grain - from + 1) as usize)
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

    // Wells: ground-state grain, density from the ground state up to the top, first grain with states.
    struct WellGrid {
        ground: isize,
        rho_from_ground: Vec<f64>,
        first: usize,
    }
    let well_index: HashMap<&str, usize> =
        deck.well_order.iter().enumerate().map(|(i, name)| (name.as_str(), i)).collect();
    let mut grids = Vec::with_capacity(deck.well_order.len());
    for name in &deck.well_order {
        let species = &deck.wells[name];
        let ground = grain_of(species.zero_energy_cm1);
        let n = grains_from(ground, &format!("Well '{name}'"))?;
        let rho_from_ground = rrho_density_of_states(n, d_e, &species_model(species)?)?;
        let first = rho_from_ground
            .iter()
            .position(|r| *r > 0.0)
            .ok_or_else(|| format!("Well '{name}' has no states below the top of the grid."))?;
        grids.push(WellGrid { ground, rho_from_ground, first });
    }

    // k(E) of well w for a channel opening at absolute grain `threshold` with W(E - E_th) = w_sum.
    let rates = |w: usize, threshold: isize, w_sum: &[f64]| -> Vec<f64> {
        let grid = &grids[w];
        (grid.first..grid.rho_from_ground.len())
            .map(|i| {
                let absolute = grid.ground + i as isize;
                if absolute < threshold {
                    0.0
                } else {
                    w_sum[(absolute - threshold) as usize] / (H_PLANCK_CM * grid.rho_from_ground[i])
                }
            })
            .collect()
    };

    let mut channels: Vec<Vec<Channel>> = vec![Vec::new(); grids.len()];
    for barrier in &deck.barriers {
        let name = &barrier.name;
        if barrier.has_tunneling && !settings.ignore_tunneling {
            return Err(format!(
                "Barrier '{name}' has a Tunneling block, but tunneling corrections of k(E) are not implemented \
                 yet; set ignore_tunneling to run without them."
            ));
        }
        let phase_space_core = matches!(barrier.core, MessBarrierCore::PhaseSpaceTheory { .. });
        let tight_sum_of_states = || -> Result<(isize, Vec<f64>), String> {
            if phase_space_core {
                return Err(format!(
                    "Barrier '{name}' is barrierless (phase-space-theory core): give its high-pressure rate \
                     coefficient in an InverseLaplaceTransform block; the barrierless module is not yet \
                     connected for chemical activation."
                ));
            }
            let threshold = grain_of(barrier.rrho.zero_energy_cm1);
            let n = grains_from(threshold, &format!("Barrier '{name}'"))?;
            Ok((threshold, rrho_sum_of_states(n, d_e, &species_model(&barrier.rrho)?)?))
        };

        match (well_index.get(barrier.left.as_str()), well_index.get(barrier.right.as_str())) {
            (Some(&a), Some(&b)) => {
                if barrier.inverse_laplace_transform.is_some() {
                    return Err(format!("Barrier '{name}': an ILT channel between two wells is not supported."));
                }
                let (threshold, w_sum) = tight_sum_of_states()?;
                channels[a].push(Channel {
                    name: name.clone(),
                    destination: ChannelDestination::Well { index: b },
                    rate_constant_s_inv: rates(a, threshold, &w_sum),
                });
                channels[b].push(Channel {
                    name: name.clone(),
                    destination: ChannelDestination::Well { index: a },
                    rate_constant_s_inv: rates(b, threshold, &w_sum),
                });
            }
            (Some(&w), None) | (None, Some(&w)) => {
                let other = if well_index.contains_key(barrier.left.as_str()) { &barrier.right } else { &barrier.left };
                let bimolecular = deck.bimolecular.get(other).ok_or_else(|| {
                    format!("Barrier '{name}' connects to '{other}', which is neither a Well nor a Bimolecular species.")
                })?;
                let (threshold, w_sum) = match &barrier.inverse_laplace_transform {
                    None => tight_sum_of_states()?,
                    Some(ilt) => {
                        let e_inf = ilt.high_pressure_rate.activation_energy_cm1;
                        match ilt.direction {
                            IltDirection::Association => {
                                let threshold = grain_of(bimolecular.ground_energy_cm1 + e_inf);
                                let n = grains_from(threshold, &format!("Barrier '{name}'"))?;
                                let rho_a = rrho_density_of_states(n, d_e, &species_model(&bimolecular.fragment_a)?)?;
                                let rho_b = rrho_density_of_states(n, d_e, &species_model(&bimolecular.fragment_b)?)?;
                                // rho_AB(E) = sum_E' rho_A(E') rho_B(E - E') dE.
                                let rho_ab: Vec<f64> =
                                    (0..n).map(|i| (0..=i).map(|j| rho_a[j] * rho_b[i - j]).sum::<f64>() * d_e).collect();
                                let mass_a: f64 = mass_vector_from_symbols_amu(&bimolecular.fragment_a.geometry_symbols)?.iter().sum();
                                let mass_b: f64 = mass_vector_from_symbols_amu(&bimolecular.fragment_b.geometry_symbols)?.iter().sum();
                                let w_sum = ilt_sum_of_states_association(
                                    &ilt.high_pressure_rate,
                                    &rho_ab,
                                    mass_a * mass_b / (mass_a + mass_b),
                                    d_e,
                                )
                                .map_err(|e| format!("Barrier '{name}': {e}"))?;
                                (threshold, w_sum)
                            }
                            IltDirection::Dissociation => {
                                let grid = &grids[w];
                                let threshold = grain_of(deck.wells[&deck.well_order[w]].zero_energy_cm1 + e_inf);
                                if threshold < grain_of(bimolecular.ground_energy_cm1) {
                                    return Err(format!(
                                        "Barrier '{name}': the ILT threshold E(well) + E_inf lies below the asymptote \
                                         '{other}'; k(E) would be non-zero below the dissociation energy."
                                    ));
                                }
                                let n = grains_from(threshold, &format!("Barrier '{name}'"))?;
                                let w_sum = ilt_sum_of_states_dissociation(&ilt.high_pressure_rate, &grid.rho_from_ground[..n], d_e)
                                    .map_err(|e| format!("Barrier '{name}': {e}"))?;
                                (threshold, w_sum)
                            }
                        }
                    }
                };
                channels[w].push(Channel {
                    name: name.clone(),
                    destination: ChannelDestination::Products { name: other.clone() },
                    rate_constant_s_inv: rates(w, threshold, &w_sum),
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
            bottom_offset_grains: grid.ground + grid.first as isize,
            density_of_states: grid.rho_from_ground[grid.first..].to_vec(),
            channels,
            lennard_jones: lennard_jones.clone(),
            energy_transfer: energy_transfer.clone(),
            bimolecular_sink_s_inv: deck.well_escape_rate_s_inv.get(name).copied().unwrap_or(0.0),
        })
        .collect();
    let network = ChemicalActivationNetwork { grain_width_cm1: d_e, wells };
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

    /// HCO + O2 (R) -> W1 (ILT association, barrierless) <-> W2 (tight) -> OH + CO2 (P, tight);
    /// W2 escapes with 1e5 s-1. Energies in kcal/mol relative to R.
    const DECK: &str = r#"
TemperatureList[K]            300. 500.
PressureList[torr]            10. 760.
EnergyStepOverTemperature     0.2
ModelEnergyLimit[kcal/mol]    40
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

    fn model() -> MessChemicalActivationModel {
        let deck = parse_mess_input(DECK).unwrap();
        chemical_activation_model_from_mess(&deck, &MessNetworkSettings::default()).unwrap()
    }

    fn grain_width() -> f64 {
        0.2 * KB_CM * 300.0
    }

    #[test]
    fn wells_share_one_grid_up_to_the_model_energy_limit() {
        let m = model();
        let d_e = grain_width();
        assert!((m.network.grain_width_cm1 - d_e).abs() < 1e-12);
        let top = (40.0 / CM1_TO_KCAL / d_e).floor() as isize;
        for (w, e_kcal) in [(0, -30.0), (1, -25.0)] {
            let well = &m.network.wells[w];
            assert_eq!(well.bottom_offset_grains + well.grain_count() as isize - 1, top, "well {w}");
            // Grid starts at the ground-state grain or the first grain above it that holds states.
            let ground = (e_kcal / CM1_TO_KCAL / d_e).round() as isize;
            assert!(well.bottom_offset_grains - ground <= 1, "well {w} starts {} grains above its ground state", well.bottom_offset_grains - ground);
            assert!(well.density_of_states.iter().all(|r| *r > 0.0));
        }
        assert_eq!(m.temperatures_kelvin, vec![300.0, 500.0]);
        assert_eq!(m.pressures_torr, vec![10.0, 760.0]);
    }

    #[test]
    fn isomerization_rates_obey_detailed_balance_exactly() {
        let balance = isomerization_detailed_balance(&model().network);
        assert_eq!(balance.len(), 1);
        assert!(balance[0].max_relative_deviation < 1e-12, "{}", balance[0].max_relative_deviation);
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
        // The isomerization opens at the same absolute grain seen from both wells: k = 0 below the
        // transition-state grain; with classical rotors W‡(0) = 0, so the first non-zero k lies in the
        // grain above it.
        let d_e = grain_width();
        let ts = (-5.0 / CM1_TO_KCAL / d_e).round() as isize;
        for (w, c) in [(0, 1), (1, 0)] {
            let well = &m.network.wells[w];
            let first_open = well.channels[c].rate_constant_s_inv.iter().position(|k| *k > 0.0).unwrap() as isize;
            assert_eq!(first_open + well.bottom_offset_grains, ts + 1, "well {w}");
        }
    }

    #[test]
    fn association_ilt_forms_the_entrance_channel_at_the_asymptote() {
        let m = model();
        assert_eq!(m.entrance_channels, vec![(0, 0)]);
        // E_inf = 0: the threshold is the asymptote (absolute grain 0). Both fragments have classical
        // rotors, whose densities vanish in their first grain, so the convolved fragment density and
        // hence W_ILT vanish in the first two grains above the threshold.
        let well = &m.network.wells[0];
        let k = &well.channels[0].rate_constant_s_inv;
        let first_open = k.iter().position(|k| *k > 0.0).unwrap() as isize;
        assert_eq!(first_open + well.bottom_offset_grains, 2);
        assert!(k.iter().all(|k| k.is_finite()));
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
        let deck = parse_mess_input(&deck).unwrap();
        let err = chemical_activation_model_from_mess(&deck, &MessNetworkSettings::default()).unwrap_err();
        assert!(err.contains("B0"), "{err}");
    }

    #[test]
    fn tunneling_must_be_ignored_explicitly() {
        let deck = DECK.replace(
            "      ZeroEnergy[kcal/mol] -5",
            "      Tunneling Eckart\n        ImaginaryFrequency[1/cm] 1500\n        WellDepth[kcal/mol] 25\n        WellDepth[kcal/mol] 20\n      End\n      ZeroEnergy[kcal/mol] -5",
        );
        let deck = parse_mess_input(&deck).unwrap();
        assert!(chemical_activation_model_from_mess(&deck, &MessNetworkSettings::default()).is_err());
        let settings = MessNetworkSettings { ignore_tunneling: true, ..Default::default() };
        assert!(chemical_activation_model_from_mess(&deck, &settings).is_ok());
    }

    #[test]
    fn the_deck_runs_through_the_chemical_activation_driver() {
        let m = model();
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
