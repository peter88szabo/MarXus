//! Input deck of a photoionization calculation: MESS-format species blocks and a `Photoionization` block.
//!
//! The molecular data of every species (neutral, molecular ion, transition states, fragment ions, neutral fragments)
//! are given as `Species NAME` followed by an RRHO or Atom block with the syntax of the MESS Fragment blocks
//! (`mess_input::parse_species_blocks`: geometry or rotational constants, frequencies, hindered rotors, electronic
//! levels; ZeroEnergy is optional and not used). The energetics are the ionization energy and the 0 K appearance
//! energies of the fragment ions, measured from the ground state of the neutral (Sztáray, Bodi, Baer, J. Mass
//! Spectrom. 45, 1233 (2010), Fig. 2).
//!
//!   Photoionization
//!     Temperature[K]               298          sample temperature
//!     Neutral                      EtBr         species of the neutral precursor
//!     Ion                          EtBr+        species of the molecular ion
//!     IonizationEnergy[eV]         10.307       adiabatic ionization energy
//!     PhotonEnergyList[eV]         11.0 11.1    photon energies, or
//!     PhotonEnergyRange[eV]        11.0 11.2 0.01   first, last, step
//!     EnergyResolution[1/cm]       60           optional: FWHM of the Gaussian photon and electron resolution
//!     FlightTime[s]                2e-5         maximum flight time; required by statistical rate models
//!     CellWidth[1/cm]              1            optional: width of the energy cells (default 1)
//!     Channel NAME
//!       Parent                     EtBr+        optional: the dissociating ion (default: the molecular ion)
//!       FragmentIon                C2H5+
//!       NeutralFragment            Br
//!       AppearanceEnergy[eV]       11.133
//!       RateModel                  Fast         Fast (default), RRKM, PhaseSpaceTheory, SimplifiedSACM
//!       TransitionState            TS1          species of the transition state (RRKM)
//!       ReactionDegeneracy         1            optional, default 1
//!       TranslationalDegreesOfFreedom 2         optional (1, 2 or 3), default 2: partitioning of the excess energy
//!     End
//!   End

use super::rates::RateModel;
use crate::masterequation::mess_input::{first_token, parse_f64, parse_species_blocks, strip_comment, unit_tag, MessSpeciesRrho};
use std::path::Path;

/// Rate treatment of a channel: every ion above the limit dissociates (Fast), or a statistical rate model with the
/// flight time (SBB10 eq. 24).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChannelRate {
    Fast,
    Statistical(RateModel),
}

/// One dissociation channel of an ion.
#[derive(Debug, Clone, PartialEq)]
pub struct PhotoionizationChannel {
    pub name: String,
    pub parent: String,
    pub fragment_ion: String,
    pub neutral_fragment: String,
    pub appearance_energy_ev: f64,
    pub rate_model: ChannelRate,
    pub transition_state: Option<String>,
    pub reaction_degeneracy: f64,
    pub translational_degrees_of_freedom: u32,
}

/// A parsed photoionization deck.
#[derive(Debug, Clone)]
pub struct PhotoionizationDeck {
    pub temperature_kelvin: f64,
    pub neutral: String,
    pub ion: String,
    pub ionization_energy_ev: f64,
    pub photon_energies_ev: Vec<f64>,
    pub resolution_fwhm_cm1: Option<f64>,
    pub flight_time_s: Option<f64>,
    pub cell_width_cm1: f64,
    pub channels: Vec<PhotoionizationChannel>,
    pub species: Vec<MessSpeciesRrho>,
}

impl PhotoionizationDeck {
    pub fn species(&self, name: &str) -> Option<&MessSpeciesRrho> {
        self.species.iter().find(|s| s.name == name)
    }
}

pub fn parse_photoionization_deck_file(path: impl AsRef<Path>) -> Result<PhotoionizationDeck, String> {
    let path = path.as_ref();
    let input = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    parse_photoionization_deck(&input).map_err(|e| format!("{}: {e}", path.display()))
}

pub fn parse_photoionization_deck(input: &str) -> Result<PhotoionizationDeck, String> {
    let species = parse_species_blocks(input)?;
    let lines: Vec<&str> = input.lines().map(strip_comment).collect();
    let start = lines
        .iter()
        .position(|l| first_token(l) == Some("Photoionization"))
        .ok_or("Photoionization deck: no Photoionization block.")?;
    let context = |what: &str| format!("Photoionization block: {what}");
    let (mut temperature, mut neutral, mut ion, mut ionization_energy) = (None, None, None, None);
    let (mut photon_energies, mut resolution, mut flight_time, mut cell_width) = (None, None, None, 1.0);
    let mut channels: Vec<PhotoionizationChannel> = Vec::new();
    let mut i = start + 1;
    loop {
        let line = *lines.get(i).ok_or_else(|| context("missing End"))?;
        i += 1;
        let tokens: Vec<&str> = line.split_whitespace().collect();
        let Some(&key) = tokens.first() else { continue };
        let value = |k: usize| tokens.get(k).copied().ok_or_else(|| context(&format!("malformed line '{line}'")));
        let number = |k: usize| value(k).and_then(|v| parse_f64(v).map_err(|e| context(&e)));
        let unit = unit_tag(line).unwrap_or("");
        let base = key.split('[').next().unwrap_or(key);
        let ev = |what: &str| if unit == "eV" { Ok(()) } else { Err(context(&format!("{what} needs the unit [eV]: '{line}'"))) };
        match base {
            "End" => break,
            "Temperature" => temperature = Some(number(1)?),
            "Neutral" => neutral = Some(value(1)?.to_string()),
            "Ion" => ion = Some(value(1)?.to_string()),
            "IonizationEnergy" => {
                ev("IonizationEnergy")?;
                ionization_energy = Some(number(1)?);
            }
            "PhotonEnergyList" => {
                ev("PhotonEnergyList")?;
                photon_energies = Some((1..tokens.len()).map(number).collect::<Result<Vec<_>, _>>()?);
            }
            "PhotonEnergyRange" => {
                ev("PhotonEnergyRange")?;
                let (first, last, step) = (number(1)?, number(2)?, number(3)?);
                if !(step > 0.0 && last >= first) {
                    return Err(context(&format!("PhotonEnergyRange needs first <= last and a positive step: '{line}'")));
                }
                let n = ((last - first) / step + 1e-9).floor() as usize;
                photon_energies = Some((0..=n).map(|k| first + k as f64 * step).collect());
            }
            "EnergyResolution" => resolution = Some(number(1)?),
            "FlightTime" => flight_time = Some(number(1)?),
            "CellWidth" => cell_width = number(1)?,
            "Channel" => {
                let name = value(1)?.to_string();
                let (channel, next) = parse_channel(&lines, i, &name)?;
                if channels.iter().any(|c| c.name == name) {
                    return Err(context(&format!("Channel '{name}' is defined twice.")));
                }
                channels.push(channel);
                i = next;
            }
            other => {
                return Err(context(&format!(
                    "unknown keyword '{other}' (Temperature[K], Neutral, Ion, IonizationEnergy[eV], PhotonEnergyList[eV], \
                     PhotonEnergyRange[eV], EnergyResolution[1/cm], FlightTime[s], CellWidth[1/cm], Channel)"
                )))
            }
        }
    }
    let missing = |what: &str| context(&format!("missing {what}"));
    let ion = ion.ok_or_else(|| missing("Ion"))?;
    let deck = PhotoionizationDeck {
        temperature_kelvin: temperature.ok_or_else(|| missing("Temperature[K]"))?,
        neutral: neutral.ok_or_else(|| missing("Neutral"))?,
        ion: ion.clone(),
        ionization_energy_ev: ionization_energy.ok_or_else(|| missing("IonizationEnergy[eV]"))?,
        photon_energies_ev: photon_energies.ok_or_else(|| missing("PhotonEnergyList[eV] or PhotonEnergyRange[eV]"))?,
        resolution_fwhm_cm1: resolution,
        flight_time_s: flight_time,
        cell_width_cm1: cell_width,
        channels: channels
            .into_iter()
            .map(|mut c| {
                if c.parent.is_empty() {
                    c.parent = ion.clone();
                }
                c
            })
            .collect(),
        species,
    };
    check_deck(&deck)?;
    Ok(deck)
}

fn parse_channel(lines: &[&str], mut i: usize, name: &str) -> Result<(PhotoionizationChannel, usize), String> {
    let context = |what: &str| format!("Channel '{name}': {what}");
    let mut channel = PhotoionizationChannel {
        name: name.to_string(),
        parent: String::new(),
        fragment_ion: String::new(),
        neutral_fragment: String::new(),
        appearance_energy_ev: f64::NAN,
        rate_model: ChannelRate::Fast,
        transition_state: None,
        reaction_degeneracy: 1.0,
        translational_degrees_of_freedom: 2,
    };
    loop {
        let line = *lines.get(i).ok_or_else(|| context("missing End"))?;
        i += 1;
        let tokens: Vec<&str> = line.split_whitespace().collect();
        let Some(&key) = tokens.first() else { continue };
        let value = || tokens.get(1).copied().ok_or_else(|| context(&format!("malformed line '{line}'")));
        let base = key.split('[').next().unwrap_or(key);
        match base {
            "End" => break,
            "Parent" => channel.parent = value()?.to_string(),
            "FragmentIon" => channel.fragment_ion = value()?.to_string(),
            "NeutralFragment" => channel.neutral_fragment = value()?.to_string(),
            "AppearanceEnergy" => {
                if unit_tag(line) != Some("eV") {
                    return Err(context(&format!("AppearanceEnergy needs the unit [eV]: '{line}'")));
                }
                channel.appearance_energy_ev = parse_f64(value()?).map_err(|e| context(&e))?;
            }
            "RateModel" => {
                channel.rate_model = match value()? {
                    "Fast" => ChannelRate::Fast,
                    "RRKM" => ChannelRate::Statistical(RateModel::Rrkm),
                    "PhaseSpaceTheory" => ChannelRate::Statistical(RateModel::PhaseSpaceTheory),
                    "SimplifiedSACM" => ChannelRate::Statistical(RateModel::SimplifiedStatisticalAdiabaticChannel),
                    other => {
                        return Err(context(&format!(
                            "RateModel '{other}' unknown (Fast, RRKM, PhaseSpaceTheory, SimplifiedSACM)"
                        )))
                    }
                }
            }
            "TransitionState" => channel.transition_state = Some(value()?.to_string()),
            "ReactionDegeneracy" => channel.reaction_degeneracy = parse_f64(value()?).map_err(|e| context(&e))?,
            "TranslationalDegreesOfFreedom" => {
                channel.translational_degrees_of_freedom = match value()?.parse::<u32>() {
                    Ok(d) if (1..=3).contains(&d) => d,
                    _ => return Err(context("TranslationalDegreesOfFreedom must be 1, 2 or 3")),
                }
            }
            other => {
                return Err(context(&format!(
                    "unknown keyword '{other}' (Parent, FragmentIon, NeutralFragment, AppearanceEnergy[eV], RateModel, \
                     TransitionState, ReactionDegeneracy, TranslationalDegreesOfFreedom)"
                )))
            }
        }
    }
    if channel.fragment_ion.is_empty() || channel.neutral_fragment.is_empty() || channel.appearance_energy_ev.is_nan() {
        return Err(context("FragmentIon, NeutralFragment and AppearanceEnergy[eV] are required"));
    }
    Ok((channel, i))
}

/// Species defined, channel parents known, statistical channels with a transition state and the flight time, and
/// appearance energies above the ionization energy and above the appearance energy of the parent ion.
fn check_deck(deck: &PhotoionizationDeck) -> Result<(), String> {
    let defined = |name: &str, role: &str| {
        deck.species(name).map(|_| ()).ok_or_else(|| format!("Photoionization deck: species '{name}' ({role}) is not defined."))
    };
    defined(&deck.neutral, "Neutral")?;
    defined(&deck.ion, "Ion")?;
    for c in &deck.channels {
        let context = |what: &str| format!("Channel '{}': {what}", c.name);
        defined(&c.fragment_ion, "FragmentIon")?;
        defined(&c.neutral_fragment, "NeutralFragment")?;
        let parent_energy = if c.parent == deck.ion {
            deck.ionization_energy_ev
        } else {
            deck.channels
                .iter()
                .find(|p| p.fragment_ion == c.parent)
                .map(|p| p.appearance_energy_ev)
                .ok_or_else(|| context(&format!("Parent '{}' is neither the Ion nor the FragmentIon of a channel", c.parent)))?
        };
        if !(c.appearance_energy_ev > parent_energy) {
            return Err(context(&format!(
                "the appearance energy {} eV lies below the energy {parent_energy} eV at which its parent '{}' is formed",
                c.appearance_energy_ev, c.parent
            )));
        }
        if let ChannelRate::Statistical(_) = c.rate_model {
            let ts = c.transition_state.as_deref().ok_or_else(|| context("a statistical RateModel needs a TransitionState"))?;
            defined(ts, "TransitionState")?;
            if deck.flight_time_s.is_none() {
                return Err(context("a statistical RateModel needs FlightTime[s] in the Photoionization block"));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub(crate) const DECK: &str = r#"
! C2H5Br+ -> C2H5+ + Br, then C2H5+ -> C2H3+ + H2 (test deck; molecular data illustrative)
Photoionization
  Temperature[K]            298
  Neutral                   EtBr
  Ion                       EtBr+
  IonizationEnergy[eV]      10.307
  PhotonEnergyList[eV]      11.00 11.10 11.20
  EnergyResolution[1/cm]    60
  FlightTime[s]             2.0e-5
  CellWidth[1/cm]           10
  Channel BrLoss
    FragmentIon             C2H5+
    NeutralFragment         Br
    AppearanceEnergy[eV]    11.133
    RateModel               RRKM
    TransitionState         TS1
  End
  Channel H2Loss
    Parent                  C2H5+
    FragmentIon             C2H3+
    NeutralFragment         H2
    AppearanceEnergy[eV]    13.40
    TranslationalDegreesOfFreedom 3
  End
End
Species EtBr
  RRHO
    RotationalConstants[1/cm] 3
      0.98 0.13 0.12
    Mass[amu] 108
    Core RigidRotor
      SymmetryFactor 1
    End
    Frequencies[1/cm] 3
      290 960 2950
  End
Species EtBr+
  RRHO
    RotationalConstants[1/cm] 3
      0.95 0.12 0.11
    Mass[amu] 108
    Core RigidRotor
      SymmetryFactor 1
    End
    Frequencies[1/cm] 3
      250 900 2900
    ElectronicLevels[1/cm] 1
      0 2
  End
Species TS1
  RRHO
    RotationalConstants[1/cm] 3
      0.9 0.1 0.09
    Mass[amu] 108
    Core RigidRotor
      SymmetryFactor 1
    End
    Frequencies[1/cm] 2
      200 2900
  End
Species C2H5+
  RRHO
    RotationalConstants[1/cm] 3
      4.2 0.9 0.8
    Mass[amu] 29
    Core RigidRotor
      SymmetryFactor 1
    End
    Frequencies[1/cm] 2
      800 3000
  End
Species C2H3+
  RRHO
    RotationalConstants[1/cm] 3
      5.0 1.0 0.9
    Mass[amu] 27
    Core RigidRotor
      SymmetryFactor 1
    End
    Frequencies[1/cm] 1
      1000
  End
Species Br
  Atom
    Mass[amu] 78.918
    ElectronicLevels[1/cm] 1
      0 4
  End
Species H2
  RRHO
    RotationalConstants[1/cm] 1
      60.8
    Mass[amu] 2.016
    Core RigidRotor
      SymmetryFactor 2
    End
    Frequencies[1/cm] 1
      4401
  End
"#;

    #[test]
    fn the_photoionization_block_gives_the_experiment_and_the_channels() {
        let deck = parse_photoionization_deck(DECK).unwrap();
        assert_eq!((deck.temperature_kelvin, deck.neutral.as_str(), deck.ion.as_str()), (298.0, "EtBr", "EtBr+"));
        assert_eq!(deck.ionization_energy_ev, 10.307);
        assert_eq!(deck.photon_energies_ev, vec![11.0, 11.1, 11.2]);
        assert_eq!((deck.resolution_fwhm_cm1, deck.flight_time_s, deck.cell_width_cm1), (Some(60.0), Some(2.0e-5), 10.0));
        let default_cells = parse_photoionization_deck(&DECK.replace("  CellWidth[1/cm]           10\n", "")).unwrap();
        assert_eq!(default_cells.cell_width_cm1, 1.0);
        assert_eq!(deck.species.len(), 7);
        let [br, h2] = &deck.channels[..] else { panic!("{:?}", deck.channels) };
        assert_eq!((br.name.as_str(), br.parent.as_str(), br.fragment_ion.as_str(), br.neutral_fragment.as_str()), ("BrLoss", "EtBr+", "C2H5+", "Br"));
        assert_eq!((br.appearance_energy_ev, br.rate_model, br.transition_state.as_deref()), (11.133, ChannelRate::Statistical(RateModel::Rrkm), Some("TS1")));
        assert_eq!((br.reaction_degeneracy, br.translational_degrees_of_freedom), (1.0, 2));
        assert_eq!((h2.parent.as_str(), h2.rate_model, h2.translational_degrees_of_freedom), ("C2H5+", ChannelRate::Fast, 3));
    }

    #[test]
    fn a_photon_energy_range_is_expanded_and_the_rate_options_are_read() {
        let deck = DECK.replace("PhotonEnergyList[eV]      11.00 11.10 11.20", "PhotonEnergyRange[eV]     11.0 11.2 0.05");
        let parsed = parse_photoionization_deck(&deck).unwrap();
        let expected = [11.0, 11.05, 11.1, 11.15, 11.2];
        assert_eq!(parsed.photon_energies_ev.len(), 5);
        assert!(parsed.photon_energies_ev.iter().zip(expected).all(|(a, b)| (a - b).abs() < 1e-12));
        for (keyword, model) in [("PhaseSpaceTheory", RateModel::PhaseSpaceTheory), ("SimplifiedSACM", RateModel::SimplifiedStatisticalAdiabaticChannel)] {
            let deck = DECK.replace("RateModel               RRKM", &format!("RateModel {keyword}"));
            assert_eq!(parse_photoionization_deck(&deck).unwrap().channels[0].rate_model, ChannelRate::Statistical(model));
        }
    }

    #[test]
    fn inconsistent_photoionization_decks_are_refused() {
        for (from, to, expected) in [
            ("  FlightTime[s]             2.0e-5\n", "  Flight 2\n", "Flight"),
            ("    TransitionState         TS1\n", "", "TransitionState"),
            ("  FlightTime[s]             2.0e-5\n", "", "FlightTime"),
            ("AppearanceEnergy[eV]    11.133", "AppearanceEnergy[eV]    10.0", "below"),
            ("FragmentIon             C2H3+", "FragmentIon             C2H4+", "C2H4+"),
            ("Parent                  C2H5+", "Parent                  C3H7+", "C3H7+"),
            ("RateModel               RRKM", "RateModel               Troe", "Troe"),
            ("IonizationEnergy[eV]      10.307", "IonizationEnergy[kcal/mol] 237", "[eV]"),
            ("Neutral                   EtBr\n", "", "Neutral"),
        ] {
            assert!(DECK.contains(from), "{from}");
            let err = parse_photoionization_deck(&DECK.replace(from, to)).unwrap_err();
            assert!(err.contains(expected), "expected '{expected}' in: {err}");
        }
    }
}
