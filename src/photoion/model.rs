//! Breakdown curves of a photoionization deck (`deck::PhotoionizationDeck`).
//!
//! - Densities of states of every species on cells covering all ion energies reached, from the species data of the
//!   deck (`chemical_activation_from_mess_input::deck_species_model`: harmonic vibrations, classical rotors, internal
//!   rotors with Kilpatrick-Pitzer reduced moments).
//! - The neutral distribution at the sample temperature (SBB10 eq. 1), cut where the remaining tail is below 1e-15.
//! - At each photon energy the ion distribution (SBB10 eq. 2) and the channels of the ion tree:
//!   - one Fast channel: the ions above the limit dissociate (SBB10 eq. 23);
//!   - statistical channels of the molecular ion: competition within the flight time (SBB10 eqs. 21 and 24);
//!   - the fragment ions formed pass their distribution (SBB10 eq. 5 summed over the dissociating ions) to their own
//!     channels, which must be Fast (SBB10 p. 1235, steps 4-6).
//! - Not available yet, refused with an error: statistical channels of a fragment ion (they need the time at which
//!   the fragment is formed within the flight time), several Fast channels of one ion (their branching needs rates),
//!   and Fast and statistical channels of one ion together.
//!
//! Reference: B. Sztáray, A. Bodi, T. Baer, J. Mass Spectrom. 45, 1233 (2010) (SBB10).

use std::collections::HashMap;

use super::deck::{ChannelRate, PhotoionizationChannel, PhotoionizationDeck};
use super::energy_distributions::{ion_energy_distribution, neutral_thermal_distribution, EV_TO_CM1};
use super::product_energy::{translational_density, ProductPartitioning};
use super::rates::rate_constants;
use crate::constants::KB_CM;
use crate::masterequation::chemical_activation_from_mess_input::deck_species_model;
use crate::masterequation::microcanonical_builder::{rrho_density_of_states, SpeciesMicroModel};
use crate::rrkm::internal_rotor::ReducedMomentModel;

/// Tail of the neutral thermal distribution left out.
const NEUTRAL_TAIL: f64 = 1e-15;

/// Rate constants and product partitionings of the channels, by channel name.
struct Channels<'c> {
    rates: &'c HashMap<String, Vec<f64>>,
    partitionings: &'c HashMap<String, ProductPartitioning>,
}

/// Abundances of the ions at one photon energy, in the order of `PhotoionizationModel::ions`.
#[derive(Debug, Clone, PartialEq)]
pub struct BreakdownRow {
    pub photon_energy_ev: f64,
    pub abundances: Vec<f64>,
}

/// A photoionization deck with its densities of states and the neutral distribution.
#[derive(Debug, Clone)]
pub struct PhotoionizationModel<'a> {
    pub deck: &'a PhotoionizationDeck,
    pub models: HashMap<String, SpeciesMicroModel>,
    /// Densities of states (per cm-1) of every species on `ion_cells` cells.
    pub densities: HashMap<String, Vec<f64>>,
    pub neutral_distribution: Vec<f64>,
    /// Cells covering every ion energy reached at the highest photon energy.
    pub ion_cells: usize,
    /// The molecular ion, then the fragment ions in the order of the channels.
    pub ions: Vec<String>,
}

impl<'a> PhotoionizationModel<'a> {
    pub fn new(deck: &'a PhotoionizationDeck) -> Result<Self, String> {
        let cell = deck.cell_width_cm1;
        let kt = KB_CM * deck.temperature_kelvin;
        // neutral distribution: counted up to 200 kT, cut where the tail is below NEUTRAL_TAIL
        let neutral_species = deck.species(&deck.neutral).ok_or("Photoionization: the neutral is not defined.")?;
        let probe_cells = (200.0 * kt / cell).ceil() as usize + 1;
        let neutral_model = deck_species_model(neutral_species, ReducedMomentModel::default(), probe_cells as f64 * cell)?;
        let full = neutral_thermal_distribution(&rrho_density_of_states(probe_cells, cell, &neutral_model)?, cell, deck.temperature_kelvin)?;
        let mut tail = 0.0;
        let mut neutral_cells = full.len();
        while neutral_cells > 1 && tail + full[neutral_cells - 1] < NEUTRAL_TAIL {
            tail += full[neutral_cells - 1];
            neutral_cells -= 1;
        }
        let kept: f64 = full[..neutral_cells].iter().sum();
        let neutral_distribution: Vec<f64> = full[..neutral_cells].iter().map(|p| p / kept).collect();
        let highest = deck.photon_energies_ev.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let excess = ((highest - deck.ionization_energy_ev) * EV_TO_CM1 / cell).round().max(0.0) as usize;
        let resolution = deck.resolution_fwhm_cm1.map_or(0, |fwhm| (6.0 * fwhm / (2.0 * (2.0 * std::f64::consts::LN_2).sqrt()) / cell).ceil() as usize);
        let ion_cells = neutral_cells + excess + resolution + 1;
        let mut models = HashMap::new();
        let mut densities = HashMap::new();
        for species in &deck.species {
            let model = deck_species_model(species, ReducedMomentModel::default(), ion_cells as f64 * cell)?;
            densities.insert(species.name.clone(), rrho_density_of_states(ion_cells, cell, &model)?);
            models.insert(species.name.clone(), model);
        }
        let mut ions = vec![deck.ion.clone()];
        for c in &deck.channels {
            if !ions.contains(&c.fragment_ion) {
                ions.push(c.fragment_ion.clone());
            }
        }
        Ok(Self { deck, models, densities, neutral_distribution, ion_cells, ions })
    }

    /// Fractional abundances of the ions at every photon energy of the deck.
    pub fn breakdown_curves(&self) -> Result<Vec<BreakdownRow>, String> {
        let deck = self.deck;
        let rates = self.statistical_rates()?;
        let partitionings = self.partitionings()?;
        deck.photon_energies_ev
            .iter()
            .map(|&h_nu| {
                let ion = ion_energy_distribution(
                    &self.neutral_distribution,
                    deck.cell_width_cm1,
                    (h_nu - deck.ionization_energy_ev) * EV_TO_CM1,
                    deck.resolution_fwhm_cm1,
                )?;
                let mut abundances = vec![0.0; self.ions.len()];
                let channels = Channels { rates: &rates, partitionings: &partitionings };
                self.distribute(&deck.ion, deck.ionization_energy_ev, ion, &channels, &mut abundances)?;
                Ok(BreakdownRow { photon_energy_ev: h_nu, abundances })
            })
            .collect()
    }

    /// k(E) on the molecular-ion cells of every statistical channel of the molecular ion, by channel name.
    fn statistical_rates(&self) -> Result<HashMap<String, Vec<f64>>, String> {
        let deck = self.deck;
        let mut rates = HashMap::new();
        for c in &deck.channels {
            let ChannelRate::Statistical(model) = c.rate_model else { continue };
            if c.parent != deck.ion {
                return Err(format!(
                    "Channel '{}': a statistical RateModel for the fragment ion '{}' is not available yet (it needs the \
                     time at which the fragment is formed within the flight time); use Fast.",
                    c.name, c.parent
                ));
            }
            let ts = c.transition_state.as_deref().ok_or_else(|| format!("Channel '{}': no TransitionState.", c.name))?;
            let k = rate_constants(
                model,
                &self.models[&deck.ion],
                &self.models[ts],
                self.onset(c, deck.ionization_energy_ev),
                self.ion_cells,
                deck.cell_width_cm1,
                c.reaction_degeneracy,
            )
            .map_err(|e| format!("Channel '{}': {e}", c.name))?;
            rates.insert(c.name.clone(), k);
        }
        Ok(rates)
    }

    /// Limit of channel c on the energy cells of its parent, formed at `parent_energy_ev` above the neutral.
    fn onset(&self, c: &PhotoionizationChannel, parent_energy_ev: f64) -> usize {
        ((c.appearance_energy_ev - parent_energy_ev) * EV_TO_CM1 / self.deck.cell_width_cm1).round() as usize
    }

    /// The eq. 5 partitioning of every channel, by channel name, for excess energies up to `ion_cells` cells.
    fn partitionings(&self) -> Result<HashMap<String, ProductPartitioning>, String> {
        self.deck
            .channels
            .iter()
            .map(|c| {
                let translation = translational_density(self.ion_cells, self.deck.cell_width_cm1, c.translational_degrees_of_freedom)?;
                let partitioning = ProductPartitioning::new(
                    &self.densities[&c.fragment_ion],
                    &self.densities[&c.neutral_fragment],
                    &translation,
                    self.ion_cells,
                )
                .map_err(|e| format!("Channel '{}': {e}", c.name))?;
                Ok((c.name.clone(), partitioning))
            })
            .collect()
    }

    /// Fragment-ion distribution of channel c from the dissociating weights on the parent's cells (SBB10 eq. 5).
    fn daughter(&self, c: &PhotoionizationChannel, weights: &[f64], onset: usize, channels: &Channels) -> Result<Vec<f64>, String> {
        channels.partitionings[&c.name].daughter(weights, onset).map_err(|e| format!("Channel '{}': {e}", c.name))
    }

    /// Adds the abundances of `name`, formed at `formed_ev` with the weights on its cells, and of its fragments.
    fn distribute(
        &self,
        name: &str,
        formed_ev: f64,
        weights: Vec<f64>,
        channels: &Channels,
        abundances: &mut [f64],
    ) -> Result<(), String> {
        let index = self.ions.iter().position(|i| i == name).expect("ion of the tree");
        let children: Vec<&PhotoionizationChannel> = self.deck.channels.iter().filter(|c| c.parent == name).collect();
        if children.is_empty() {
            abundances[index] += weights.iter().sum::<f64>();
            return Ok(());
        }
        let fast = children.iter().filter(|c| c.rate_model == ChannelRate::Fast).count();
        if fast == children.len() {
            if fast > 1 {
                return Err(format!(
                    "Ion '{name}': {fast} Fast channels; the branching between parallel channels needs a statistical \
                     RateModel."
                ));
            }
            let c = children[0];
            let onset = self.onset(c, formed_ev);
            abundances[index] += weights[..onset.min(weights.len())].iter().sum::<f64>();
            let daughter = self.daughter(c, &weights, onset, channels)?;
            return self.distribute(&c.fragment_ion, c.appearance_energy_ev, daughter, channels, abundances);
        }
        if fast > 0 {
            return Err(format!("Ion '{name}': Fast and statistical channels together are not available yet."));
        }
        // statistical channels of the molecular ion competing within the flight time (SBB10 eqs. 21, 24)
        let tau = self.deck.flight_time_s.ok_or("Photoionization: statistical channels need FlightTime[s].")?;
        let k: Vec<&Vec<f64>> = children.iter().map(|c| &channels.rates[&c.name]).collect();
        let mut channel_weights = vec![vec![0.0; weights.len()]; children.len()];
        for (i, w) in weights.iter().enumerate().filter(|(_, w)| **w != 0.0) {
            let total: f64 = k.iter().map(|k| k[i]).sum();
            let survive = (-total * tau).exp();
            abundances[index] += w * survive;
            if total > 0.0 {
                for (j, kj) in k.iter().enumerate() {
                    channel_weights[j][i] = w * kj[i] / total * (1.0 - survive);
                }
            }
        }
        for (c, w) in children.iter().zip(channel_weights) {
            let daughter = self.daughter(c, &w, self.onset(c, formed_ev), channels)?;
            self.distribute(&c.fragment_ion, c.appearance_energy_ev, daughter, channels, abundances)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::photoion::breakdown::{breakdown_curve, parallel_breakdown_curve, sequential_breakdown, SequentialStep};
    use crate::photoion::deck::parse_photoionization_deck;
    use crate::photoion::deck::tests::DECK;
    use crate::photoion::energy_distributions::{ion_energy_distribution, EV_TO_CM1};
    use crate::photoion::product_energy::translational_density;
    use crate::photoion::rates::{rate_constants, RateModel};

    /// Cell width of the test deck (cm-1).
    const CELL: f64 = 10.0;

    /// The test deck with only the bromine loss, as a fast channel.
    fn fast_bromine_loss() -> String {
        let start = DECK.find("  Channel H2Loss").unwrap();
        let end = DECK[start..].find("  End\n").unwrap() + start + "  End\n".len();
        format!("{}{}", &DECK[..start], &DECK[end..]).replace("    RateModel               RRKM\n", "")
    }

    #[test]
    fn a_single_fast_channel_gives_the_breakdown_curve_of_the_thermal_ion_distribution() {
        let deck = parse_photoionization_deck(&fast_bromine_loss()).unwrap();
        let run = PhotoionizationModel::new(&deck).unwrap();
        let curves = run.breakdown_curves().unwrap();
        let expected = breakdown_curve(&run.neutral_distribution, CELL, 10.307, 11.133, &deck.photon_energies_ev, Some(60.0)).unwrap();
        assert_eq!(run.ions, vec!["EtBr+".to_string(), "C2H5+".to_string()]);
        for (row, point) in curves.iter().zip(&expected) {
            assert_eq!(row.photon_energy_ev, point.photon_energy_ev);
            assert!((row.abundances[0] - point.parent).abs() < 1e-12 && (row.abundances[1] - point.fragment).abs() < 1e-12);
        }
    }

    #[test]
    fn an_rrkm_channel_competes_with_the_flight_time() {
        let text = fast_bromine_loss().replace("    TransitionState         TS1\n", "    RateModel RRKM\n    TransitionState TS1\n");
        let deck = parse_photoionization_deck(&text).unwrap();
        let run = PhotoionizationModel::new(&deck).unwrap();
        let curves = run.breakdown_curves().unwrap();
        let onset = ((11.133 - 10.307) * EV_TO_CM1 / CELL).round() as usize;
        let k = rate_constants(RateModel::Rrkm, &run.models["EtBr+"], &run.models["TS1"], onset, run.ion_cells, CELL, 1.0).unwrap();
        let expected = parallel_breakdown_curve(&run.neutral_distribution, CELL, 10.307, &[&k], 2.0e-5, &deck.photon_energies_ev, Some(60.0)).unwrap();
        for (row, point) in curves.iter().zip(&expected) {
            assert!((row.abundances[0] - point.parent).abs() < 1e-12 && (row.abundances[1] - point.fragments[0]).abs() < 1e-12);
        }
    }

    #[test]
    fn a_fast_sequential_channel_dissociates_the_fragment_ion_further() {
        let text = DECK.replace("    RateModel               RRKM\n", "").replace("PhotonEnergyList[eV]      11.00 11.10 11.20", "PhotonEnergyList[eV] 13.2 13.6 14.0");
        let deck = parse_photoionization_deck(&text).unwrap();
        let run = PhotoionizationModel::new(&deck).unwrap();
        assert_eq!(run.ions, vec!["EtBr+".to_string(), "C2H5+".to_string(), "C2H3+".to_string()]);
        let curves = run.breakdown_curves().unwrap();
        let tr2 = translational_density(run.ion_cells, CELL, 2).unwrap();
        let tr3 = translational_density(run.ion_cells, CELL, 3).unwrap();
        let rho = |name: &str| run.densities[name].as_slice();
        let steps = [
            SequentialStep { onset_cells: ((11.133 - 10.307) * EV_TO_CM1 / CELL).round() as usize, rho_fragment: rho("C2H5+"), rho_neutral: rho("Br"), rho_translation: &tr2 },
            SequentialStep { onset_cells: ((13.40 - 11.133) * EV_TO_CM1 / CELL).round() as usize, rho_fragment: rho("C2H3+"), rho_neutral: rho("H2"), rho_translation: &tr3 },
        ];
        for (row, h_nu) in curves.iter().zip(&deck.photon_energies_ev) {
            let ion = ion_energy_distribution(&run.neutral_distribution, CELL, (h_nu - 10.307) * EV_TO_CM1, Some(60.0)).unwrap();
            let expected = sequential_breakdown(&ion, &steps).unwrap();
            for (a, b) in row.abundances.iter().zip(&expected) {
                assert!((a - b).abs() < 1e-12, "{h_nu} eV: {:?} vs {expected:?}", row.abundances);
            }
            assert!((row.abundances.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        }
        assert!(curves[2].abundances[2] > curves[0].abundances[2]);
    }

    #[test]
    fn combinations_without_a_model_yet_are_refused() {
        let sacm = parse_photoionization_deck(&DECK.replace("RateModel               RRKM", "RateModel SimplifiedSACM")).unwrap();
        assert!(PhotoionizationModel::new(&sacm).unwrap().breakdown_curves().unwrap_err().contains("not available yet"));
        let slow_second = DECK.replace("    TranslationalDegreesOfFreedom 3\n", "    TranslationalDegreesOfFreedom 3\n    RateModel RRKM\n    TransitionState TS1\n");
        let deck = parse_photoionization_deck(&slow_second).unwrap();
        assert!(PhotoionizationModel::new(&deck).unwrap().breakdown_curves().unwrap_err().contains("not available yet"));
        let two_fast = DECK.replace("    RateModel               RRKM\n", "").replace("    Parent                  C2H5+\n", "");
        let deck = parse_photoionization_deck(&two_fast).unwrap();
        assert!(PhotoionizationModel::new(&deck).unwrap().breakdown_curves().unwrap_err().contains("Fast"));
    }
}
