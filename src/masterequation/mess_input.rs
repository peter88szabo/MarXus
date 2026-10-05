//! Reader of input decks in the MESS input format (a subset).
//!
//! MarXus reads this format so that existing decks can be used; the chemical-activation network is
//! built from the parsed deck in `chemical_activation_from_mess_input.rs`.
//!
//! Read:
//! - TemperatureList[K], PressureList[torr | atm | bar], EnergyStepOverTemperature,
//!   ModelEnergyLimit[kcal/mol], Reactant
//! - Model -> EnergyRelaxation -> Exponential (Factor[1/cm] at 300 K, Power, ExponentCutoff)
//! - Model -> CollisionFrequency -> LennardJones (Epsilons, Sigmas, Masses of species and bath gas)
//! - Bimolecular (two Fragment blocks, RRHO or Atom with Mass[amu], GroundEnergy), Well (RRHO, Escape pseudo-first-order rate
//!   constant), Barrier (RRHO, Core RigidRotor or PhaseSpaceTheory, presence of a Tunneling block)
//!   - RRHO -> Geometry[angstrom] N (rotational constants via `inertia::get_brot`)
//!   - RRHO -> Core -> SymmetryFactor
//!   - RRHO -> Frequencies[1/cm] N
//!   - RRHO -> ZeroEnergy[kcal/mol | kJ/mol | 1/cm]
//!   - RRHO -> ElectronicLevels[1/cm] N (only the ground-level degeneracy is used)
//! - MarXus extension: `InverseLaplaceTransform ... End` inside the RRHO block of a barrier
//!   (high-pressure rate coefficient of a barrierless channel, see `IltSpecification`).
//!
//!   - RRHO -> Tunneling Eckart (ImaginaryFrequency, two WellDepth values); other models recorded
//! Not read: excited electronic levels, hindered rotors and other model types.

use crate::constants::CM1_TO_KCAL;
use std::collections::HashMap;
use std::path::Path;

use crate::barrierless::ilt::ilt_barrierless::ModifiedArrhenius;
use crate::inertia::inertia::get_brot;
use crate::utils::atomic_masses::mass_vector_from_symbols_amu;

#[derive(Clone, Debug)]
pub struct MessGlobal {
    /// TemperatureList[K].
    pub temperatures_kelvin: Vec<f64>,
    /// PressureList[torr | atm | bar], converted to Torr.
    pub pressures_torr: Vec<f64>,
    /// ExponentCutoff of the exponential-down model: transitions beyond cutoff x <dE_down> are neglected.
    pub exponent_cutoff: Option<f64>,
    /// ExcessEnergyOverTemperature: top of the master-equation grid above the highest barrier, in kT.
    pub excess_energy_over_temperature: Option<f64>,
    pub energy_step_over_temperature: Option<f64>,
    pub model_energy_limit_kcal_mol: Option<f64>,

    pub alpha_factor_cm1: Option<f64>,
    pub alpha_power: Option<f64>,

    pub lj_epsilons_cm1: Option<(f64, f64)>,
    pub lj_sigmas_angstrom: Option<(f64, f64)>,
    pub lj_masses_amu: Option<(f64, f64)>,

    pub reactant_name: Option<String>,
    pub excess_reactant_concentration_cm3: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct MessSpeciesRrho {
    pub name: String,
    pub geometry_symbols: Vec<String>,
    pub geometry_angstrom: Vec<[f64; 3]>,
    pub symmetry_factor: f64,
    pub vibrational_frequencies_cm1: Vec<f64>,
    pub zero_energy_cm1: f64,
    pub electronic_degeneracy_ground: f64,
    /// Mass of an `Atom` fragment (amu); None for RRHO species (mass from the geometry).
    pub atom_mass_amu: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct MessBimolecular {
    pub name: String,
    pub fragment_a: MessSpeciesRrho,
    pub fragment_b: MessSpeciesRrho,
    /// Bimolecular asymptote energy (cm^-1) relative to the same reference used for wells/barriers.
    pub ground_energy_cm1: f64,
}

#[derive(Clone, Debug)]
pub enum MessBarrierCore {
    TightRrho,
    PhaseSpaceTheory {
        fragment_a_geometry_symbols: Vec<String>,
        fragment_a_geometry_angstrom: Vec<[f64; 3]>,
        fragment_b_geometry_symbols: Vec<String>,
        fragment_b_geometry_angstrom: Vec<[f64; 3]>,
        symmetry_operations: f64,
        potential_prefactor_au: f64,
        potential_power_exponent: f64,
    },
}

/// Direction to which the high-pressure rate coefficient of an ILT channel refers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum IltDirection {
    /// Association of the bimolecular side to the well, k_inf in cm3 s-1.
    Association,
    /// Dissociation of the well, k_inf in s-1.
    Dissociation,
}

/// MarXus keyword block `InverseLaplaceTransform ... End` inside the RRHO block of a barrier: k(E) of
/// this barrierless channel follows from the inverse Laplace transform of k_inf(T)
/// (`barrierless::ilt::ilt_barrierless`) instead of the core of the barrier.
///
///   InverseLaplaceTransform
///     Direction                    Association      (or Dissociation)
///     PreExponential[cm^3/s]       6.0e-12          ([1/s] for Dissociation)
///     TemperatureExponent          -0.5
///     ReferenceTemperature[K]      298.0
///     ActivationEnergy[kcal/mol]   0.0              ([1/cm] and [kJ/mol] also accepted)
///   End
#[derive(Clone, Debug, PartialEq)]
pub struct IltSpecification {
    pub direction: IltDirection,
    pub high_pressure_rate: ModifiedArrhenius,
}

#[derive(Clone, Debug)]
pub struct MessBarrier {
    pub name: String,
    pub left: String,
    pub right: String,
    pub rrho: MessSpeciesRrho,
    pub core: MessBarrierCore,
    /// ILT parameters of a barrierless channel (MarXus keyword block), if given.
    pub inverse_laplace_transform: Option<IltSpecification>,
    /// `Tunneling` block of the barrier, if given.
    pub tunneling: Option<TunnelingSpecification>,
}

/// Tunneling model of a barrier.
#[derive(Clone, Debug, PartialEq)]
pub enum TunnelingSpecification {
    /// `Tunneling Eckart`: magnitude of the imaginary frequency and the two well depths (barrier heights
    /// seen from the two sides), the parameters of the Eckart barrier (Miller, J. Am. Chem. Soc. 101,
    /// 6810 (1979), eq. 8).
    Eckart { imaginary_frequency_cm1: f64, well_depths_cm1: [f64; 2] },
    /// Any other tunneling model (not implemented; refused when the network is built).
    Unsupported { model: String },
}

#[derive(Clone, Debug)]
pub struct MessDeck {
    pub global: MessGlobal,
    pub bimolecular: HashMap<String, MessBimolecular>,
    pub wells: HashMap<String, MessSpeciesRrho>,
    pub barriers: Vec<MessBarrier>,
    /// `Escape` pseudo-first-order rate constant of a well (s-1): the bimolecular sink k_c[D].
    pub well_escape_rate_s_inv: HashMap<String, f64>,
    /// Well names in the order of the input deck.
    pub well_order: Vec<String>,
}


fn strip_comment(mut line: &str) -> &str {
    if let Some(idx) = line.find('#') {
        line = &line[..idx];
    }
    if let Some(idx) = line.find('!') {
        line = &line[..idx];
    }
    line.trim()
}

fn first_token(line: &str) -> Option<&str> {
    line.split_whitespace().next()
}

fn parse_f64(raw: &str) -> Result<f64, String> {
    raw.trim()
        .parse::<f64>()
        .map_err(|_| format!("Invalid float '{}'", raw))
}

fn parse_usize(raw: &str) -> Result<usize, String> {
    raw.trim()
        .parse::<usize>()
        .map_err(|_| format!("Invalid integer '{}'", raw))
}

fn parse_key_value_whitespace(line: &str) -> Option<(&str, &str)> {
    let mut it = line.split_whitespace();
    let k = it.next()?;
    let v = it.next()?;
    Some((k, v))
}

fn energy_to_cm1(value: f64, unit_tag: &str) -> Result<f64, String> {
    let u = unit_tag.trim().to_lowercase();
    if u.contains("kcal") {
        Ok(value / CM1_TO_KCAL)
    } else if u.contains("kj") {
        // 1 kcal = 4.184 kJ (thermochemical calorie).
        Ok(value / (CM1_TO_KCAL * 4.184))
    } else if u.contains("1/cm") || u.contains("cm") {
        Ok(value)
    } else {
        Err(format!("Unsupported energy unit tag '{}'", unit_tag))
    }
}

// Keywords that open a block closed by its own `End`. `Barrier`, `Fragment` and `Species` are
// headers without an `End`: they end with the `End` of the model block (RRHO) that follows them.
fn is_block_starter(tok: &str) -> bool {
    matches!(
        tok,
        "Model"
            | "Bimolecular"
            | "Well"
            | "RRHO"
            | "Core"
            | "RigidRotor"
            | "Tunneling"
            | "Escape"
            | "EnergyRelaxation"
            | "CollisionFrequency"
            | "Exponential"
            | "LennardJones"
            | "TimeEvolution"
            | "InverseLaplaceTransform"
            | "Atom"
    )
}

fn collect_block(lines: &[String], start: usize) -> (Vec<String>, usize) {
    let mut depth: i32 = 0;
    let mut out: Vec<String> = Vec::new();
    let mut i = start;

    while i < lines.len() {
        let raw = &lines[i];
        let line = strip_comment(raw);
        if line.is_empty() {
            i += 1;
            continue;
        }

        let tok = first_token(line).unwrap_or("");
        if tok == "End" {
            depth -= 1;
            out.push(line.to_string());
            i += 1;
            if depth <= 0 {
                break;
            }
            continue;
        }

        if is_block_starter(tok) {
            depth += 1;
        }

        out.push(line.to_string());
        i += 1;
    }

    (out, i)
}

fn parse_geometry(
    block: &[String],
    key: &str,
) -> Result<Option<(Vec<String>, Vec<[f64; 3]>)>, String> {
    for (idx, line) in block.iter().enumerate() {
        if line.starts_with(key) {
            let parts: Vec<&str> = line.split_whitespace().collect();
            let n = parts
                .last()
                .ok_or_else(|| format!("Malformed {} line: {}", key, line))?;
            let nat = parse_usize(n)?;
            let mut symbols: Vec<String> = Vec::with_capacity(nat);
            let mut coords: Vec<[f64; 3]> = Vec::with_capacity(nat);
            for j in 0..nat {
                let l = block
                    .get(idx + 1 + j)
                    .ok_or_else(|| format!("Unexpected EOF while reading {}", key))?;
                let fields: Vec<&str> = l.split_whitespace().collect();
                if fields.len() < 4 {
                    return Err(format!("Malformed geometry line: {}", l));
                }
                symbols.push(fields[0].to_string());
                coords.push([
                    parse_f64(fields[1])?,
                    parse_f64(fields[2])?,
                    parse_f64(fields[3])?,
                ]);
            }
            return Ok(Some((symbols, coords)));
        }
    }
    Ok(None)
}

fn parse_frequencies(block: &[String]) -> Result<Vec<f64>, String> {
    for (idx, line) in block.iter().enumerate() {
        if line.starts_with("Frequencies") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            let n = parts
                .last()
                .ok_or_else(|| format!("Malformed Frequencies line: {}", line))?;
            let nfreq = parse_usize(n)?;

            let mut out: Vec<f64> = Vec::with_capacity(nfreq);
            let mut j = idx + 1;
            while j < block.len() && out.len() < nfreq {
                if block[j].starts_with("ZeroEnergy") || block[j].starts_with("ElectronicLevels") {
                    break;
                }
                for tok in block[j].split_whitespace() {
                    if out.len() == nfreq {
                        break;
                    }
                    if let Ok(v) = parse_f64(tok) {
                        out.push(v);
                    }
                }
                j += 1;
            }
            if out.len() != nfreq {
                return Err(format!(
                    "Expected {} frequencies, got {} while parsing",
                    nfreq,
                    out.len()
                ));
            }
            return Ok(out);
        }
    }
    Err("Missing Frequencies[1/cm] block".into())
}

fn parse_symmetry_factor(block: &[String]) -> Result<f64, String> {
    for line in block {
        if line.starts_with("SymmetryFactor") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            let v = parts
                .last()
                .ok_or_else(|| format!("Malformed SymmetryFactor line: {}", line))?;
            return parse_f64(v);
        }
    }
    Ok(1.0)
}

fn parse_zero_energy_cm1(block: &[String]) -> Result<f64, String> {
    for line in block {
        if line.starts_with("ZeroEnergy") {
            // Examples:
            //   ZeroEnergy[kcal/mol]  19.5
            //   ZeroEnergy[1/cm]      0
            let unit_tag = line
                .split('[')
                .nth(1)
                .and_then(|s| s.split(']').next())
                .unwrap_or("cm-1");
            let parts: Vec<&str> = line.split_whitespace().collect();
            let v = parts
                .last()
                .ok_or_else(|| format!("Malformed ZeroEnergy line: {}", line))?;
            return energy_to_cm1(parse_f64(v)?, unit_tag);
        }
    }
    Err("Missing ZeroEnergy[...]".into())
}

fn parse_electronic_degeneracy_ground(block: &[String]) -> Result<f64, String> {
    for (idx, line) in block.iter().enumerate() {
        if line.starts_with("ElectronicLevels") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            let n = parts
                .last()
                .ok_or_else(|| format!("Malformed ElectronicLevels line: {}", line))?;
            let nlev = parse_usize(n)?;
            if nlev == 0 {
                return Ok(1.0);
            }
            let first = block
                .get(idx + 1)
                .ok_or_else(|| "ElectronicLevels count but no following line".to_string())?;
            let fields: Vec<&str> = first.split_whitespace().collect();
            if fields.len() < 2 {
                return Err(format!("Malformed ElectronicLevels entry: {}", first));
            }
            let _e = parse_f64(fields[0])?;
            let g = parse_f64(fields[1])?;
            return Ok(g);
        }
    }
    Ok(1.0)
}

/// A Fragment of a Bimolecular block: an `Atom` (Mass[amu], ElectronicLevels) or an RRHO species.
fn parse_fragment(block: &[String], name: &str) -> Result<MessSpeciesRrho, String> {
    if !block.iter().any(|l| first_token(l) == Some("Atom")) {
        return parse_rrho_species(block, name);
    }
    let mass_line = block
        .iter()
        .find(|l| first_token(l).map_or(false, |t| t.starts_with("Mass")))
        .ok_or_else(|| format!("Atom fragment '{name}' without Mass[amu]"))?;
    let mass = parse_f64(mass_line.split_whitespace().nth(1).ok_or_else(|| format!("Malformed line: {mass_line}"))?)?;
    Ok(MessSpeciesRrho {
        name: name.to_string(),
        geometry_symbols: Vec::new(),
        geometry_angstrom: Vec::new(),
        symmetry_factor: 1.0,
        vibrational_frequencies_cm1: Vec::new(),
        zero_energy_cm1: 0.0,
        electronic_degeneracy_ground: parse_electronic_degeneracy_ground(block)?,
        atom_mass_amu: Some(mass),
    })
}

fn parse_rrho_species(block: &[String], name: &str) -> Result<MessSpeciesRrho, String> {
    parse_rrho_species_impl(block, name, true)
}

// `geometry_required = false` for phase-space-theory barriers: their rotational treatment comes
// from the two fragment geometries of the core, not from a molecular geometry.
fn parse_rrho_species_impl(
    block: &[String],
    name: &str,
    geometry_required: bool,
) -> Result<MessSpeciesRrho, String> {
    let (symbols, coords) = match parse_geometry(block, "Geometry[angstrom]")? {
        Some(geometry) => geometry,
        None if !geometry_required => (Vec::new(), Vec::new()),
        None => return Err(format!("RRHO species '{}' missing Geometry[angstrom]", name)),
    };
    let symmetry_factor = parse_symmetry_factor(block)?;
    let vib = parse_frequencies(block)?;
    let zero_energy_cm1 = parse_zero_energy_cm1(block)?;
    let electronic_degeneracy_ground = parse_electronic_degeneracy_ground(block)?;

    Ok(MessSpeciesRrho {
        atom_mass_amu: None,
        name: name.to_string(),
        geometry_symbols: symbols,
        geometry_angstrom: coords,
        symmetry_factor,
        vibrational_frequencies_cm1: vib,
        zero_energy_cm1,
        electronic_degeneracy_ground,
    })
}

/// Unit tag between square brackets of the first token, e.g. "kcal/mol" for "ZeroEnergy[kcal/mol]".
fn unit_tag(line: &str) -> Option<&str> {
    let first = first_token(line)?;
    first.split('[').nth(1).and_then(|t| t.split(']').next())
}

/// All numbers after the keyword of a list line such as "TemperatureList[K] 300 400 500".
fn parse_list(line: &str) -> Result<Vec<f64>, String> {
    line.split_whitespace().skip(1).map(parse_f64).collect()
}

/// The MarXus `InverseLaplaceTransform ... End` block inside a barrier, if present.
fn parse_inverse_laplace_transform(block: &[String], barrier: &str) -> Result<Option<IltSpecification>, String> {
    let Some(start) = block.iter().position(|l| first_token(l) == Some("InverseLaplaceTransform")) else {
        return Ok(None);
    };
    let context = |what: &str| format!("Barrier '{barrier}', InverseLaplaceTransform: {what}");
    let mut direction = None;
    let mut pre_exponential = None;
    let mut exponent = None;
    let mut reference_temperature = None;
    let mut activation_energy = None;
    for line in &block[start + 1..] {
        let key = first_token(line).unwrap_or("");
        if key == "End" {
            break;
        }
        let value = line.split_whitespace().nth(1).ok_or_else(|| context(&format!("no value in '{line}'")))?;
        if key == "Direction" {
            direction = Some(match value {
                "Association" => IltDirection::Association,
                "Dissociation" => IltDirection::Dissociation,
                other => return Err(context(&format!("unknown Direction '{other}' (Association or Dissociation)"))),
            });
        } else if key.starts_with("PreExponential") {
            pre_exponential = Some((parse_f64(value)?, unit_tag(line).unwrap_or("").to_string()));
        } else if key == "TemperatureExponent" {
            exponent = Some(parse_f64(value)?);
        } else if key.starts_with("ReferenceTemperature") {
            reference_temperature = Some(parse_f64(value)?);
        } else if key.starts_with("ActivationEnergy") {
            let unit = unit_tag(line).ok_or_else(|| context("ActivationEnergy needs a unit, e.g. [kcal/mol]"))?;
            activation_energy = Some(energy_to_cm1(parse_f64(value)?, unit)?);
        } else {
            return Err(context(&format!("unknown keyword '{key}'")));
        }
    }
    let direction = direction.ok_or_else(|| context("missing Direction"))?;
    let (a, unit) = pre_exponential.ok_or_else(|| context("missing PreExponential"))?;
    let expected_unit = match direction {
        IltDirection::Association => "cm^3/s",
        IltDirection::Dissociation => "1/s",
    };
    if unit != expected_unit {
        return Err(context(&format!(
            "PreExponential[{unit}] does not match Direction {direction:?}, which requires PreExponential[{expected_unit}]"
        )));
    }
    Ok(Some(IltSpecification {
        direction,
        high_pressure_rate: ModifiedArrhenius {
            pre_exponential: a,
            temperature_exponent: exponent.ok_or_else(|| context("missing TemperatureExponent"))?,
            reference_temperature_kelvin: reference_temperature.ok_or_else(|| context("missing ReferenceTemperature[K]"))?,
            activation_energy_cm1: activation_energy.ok_or_else(|| context("missing ActivationEnergy"))?,
        },
    }))
}

/// The `Tunneling <model> ... End` block of a barrier, if present.
fn parse_tunneling(block: &[String], barrier: &str) -> Result<Option<TunnelingSpecification>, String> {
    let Some(start) = block.iter().position(|l| first_token(l) == Some("Tunneling")) else {
        return Ok(None);
    };
    let model = block[start].split_whitespace().nth(1).unwrap_or("").to_string();
    if model != "Eckart" {
        return Ok(Some(TunnelingSpecification::Unsupported { model }));
    }
    let mut frequency = None;
    let mut depths = Vec::new();
    for line in &block[start + 1..] {
        let key = first_token(line).unwrap_or("");
        if key == "End" {
            break;
        }
        let value = line
            .split_whitespace()
            .nth(1)
            .ok_or_else(|| format!("Barrier '{barrier}', Tunneling: no value in '{line}'"))?;
        if key.starts_with("ImaginaryFrequency") {
            frequency = Some(parse_f64(value)?.abs());
        } else if key.starts_with("WellDepth") {
            let unit = unit_tag(line).ok_or_else(|| format!("Barrier '{barrier}': WellDepth needs a unit"))?;
            depths.push(energy_to_cm1(parse_f64(value)?, unit)?);
        }
    }
    let imaginary_frequency_cm1 =
        frequency.ok_or_else(|| format!("Barrier '{barrier}', Tunneling Eckart: missing ImaginaryFrequency."))?;
    if depths.len() != 2 || depths.iter().any(|d| !(*d > 0.0)) {
        return Err(format!(
            "Barrier '{barrier}', Tunneling Eckart: two positive WellDepth values are required, got {depths:?}."
        ));
    }
    Ok(Some(TunnelingSpecification::Eckart { imaginary_frequency_cm1, well_depths_cm1: [depths[0], depths[1]] }))
}

fn parse_phasespace_core(block: &[String]) -> Result<Option<MessBarrierCore>, String> {
    // Look for "Core PhaseSpaceTheory" within the barrier RRHO block.
    if !block
        .iter()
        .any(|l| l.starts_with("Core") && l.contains("PhaseSpaceTheory"))
    {
        return Ok(None);
    }

    // In MESS the core contains 2x FragmentGeometry[angstrom] blocks.
    let mut fragments: Vec<(Vec<String>, Vec<[f64; 3]>)> = Vec::new();
    let mut symmetry_operations: Option<f64> = None;
    let mut v0_au: Option<f64> = None;
    let mut n: Option<f64> = None;

    let mut i = 0usize;
    while i < block.len() {
        let line = &block[i];
        if line.starts_with("FragmentGeometry[angstrom]") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            let nat = parse_usize(
                parts
                    .last()
                    .ok_or("Malformed FragmentGeometry".to_string())?,
            )?;
            let mut symbols: Vec<String> = Vec::with_capacity(nat);
            let mut coords: Vec<[f64; 3]> = Vec::with_capacity(nat);
            for j in 0..nat {
                let l = block
                    .get(i + 1 + j)
                    .ok_or_else(|| "Unexpected EOF in FragmentGeometry".to_string())?;
                let fields: Vec<&str> = l.split_whitespace().collect();
                if fields.len() < 4 {
                    return Err(format!("Malformed FragmentGeometry line: {}", l));
                }
                symbols.push(fields[0].to_string());
                coords.push([
                    parse_f64(fields[1])?,
                    parse_f64(fields[2])?,
                    parse_f64(fields[3])?,
                ]);
            }
            fragments.push((symbols, coords));
            i += 1 + nat;
            continue;
        }

        if line.starts_with("SymmetryFactor") {
            symmetry_operations = Some(parse_f64(line.split_whitespace().last().unwrap())?);
        } else if line.starts_with("PotentialPrefactor") {
            v0_au = Some(parse_f64(line.split_whitespace().last().unwrap())?);
        } else if line.starts_with("PotentialPowerExponent") {
            n = Some(parse_f64(line.split_whitespace().last().unwrap())?);
        }
        i += 1;
    }

    if fragments.len() != 2 {
        return Err(
            "PhaseSpaceTheory core must define exactly two FragmentGeometry blocks.".into(),
        );
    }

    Ok(Some(MessBarrierCore::PhaseSpaceTheory {
        fragment_a_geometry_symbols: fragments[0].0.clone(),
        fragment_a_geometry_angstrom: fragments[0].1.clone(),
        fragment_b_geometry_symbols: fragments[1].0.clone(),
        fragment_b_geometry_angstrom: fragments[1].1.clone(),
        symmetry_operations: symmetry_operations.unwrap_or(1.0),
        potential_prefactor_au: v0_au.ok_or("Missing PotentialPrefactor[au]")?,
        potential_power_exponent: n.ok_or("Missing PotentialPowerExponent")?,
    }))
}

pub fn parse_mess_input(input: &str) -> Result<MessDeck, String> {
    let lines: Vec<String> = input.lines().map(|s| s.to_string()).collect();

    let mut global = MessGlobal {
        temperatures_kelvin: Vec::new(),
        pressures_torr: Vec::new(),
        exponent_cutoff: None,
        excess_energy_over_temperature: None,
        energy_step_over_temperature: None,
        model_energy_limit_kcal_mol: None,
        alpha_factor_cm1: None,
        alpha_power: None,
        lj_epsilons_cm1: None,
        lj_sigmas_angstrom: None,
        lj_masses_amu: None,
        reactant_name: None,
        excess_reactant_concentration_cm3: None,
    };

    let mut wells: HashMap<String, MessSpeciesRrho> = HashMap::new();
    let mut bimolecular: HashMap<String, MessBimolecular> = HashMap::new();
    let mut barriers: Vec<MessBarrier> = Vec::new();
    let mut well_escape_rate_s_inv: HashMap<String, f64> = HashMap::new();
    let mut well_order: Vec<String> = Vec::new();

    let mut i = 0usize;
    while i < lines.len() {
        let line = strip_comment(&lines[i]);
        if line.is_empty() {
            i += 1;
            continue;
        }

        // Global scalar parameters
        if let Some((k, v)) = parse_key_value_whitespace(line) {
            match k {
                "TemperatureList[K]" => global.temperatures_kelvin = parse_list(line)?,
                "PressureList[torr]" => global.pressures_torr = parse_list(line)?,
                "PressureList[atm]" => {
                    global.pressures_torr = parse_list(line)?.iter().map(|p| p * 760.0).collect()
                }
                // 1 bar = 1e5 Pa, 1 Torr = 101325/760 Pa.
                "PressureList[bar]" => {
                    global.pressures_torr = parse_list(line)?.iter().map(|p| p * 1.0e5 * 760.0 / 101_325.0).collect()
                }
                "ExponentCutoff" => global.exponent_cutoff = Some(parse_f64(v)?),
                "ExcessEnergyOverTemperature" => global.excess_energy_over_temperature = Some(parse_f64(v)?),
                "EnergyStepOverTemperature" => {
                    global.energy_step_over_temperature = Some(parse_f64(v)?)
                }
                "ModelEnergyLimit[kcal/mol]" => {
                    global.model_energy_limit_kcal_mol = Some(parse_f64(v)?)
                }
                "Reactant" => global.reactant_name = Some(v.to_string()),
                "ExcessReactantConcentration[molecule/cm^3]" => {
                    global.excess_reactant_concentration_cm3 = Some(parse_f64(v)?)
                }
                "Factor[1/cm]" => global.alpha_factor_cm1 = Some(parse_f64(v)?),
                "Power" => global.alpha_power = Some(parse_f64(v)?),
                _ => {}
            }
        }

        // Lennard-Jones model lines (3-list)
        if line.starts_with("Epsilons[1/cm]") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 3 {
                global.lj_epsilons_cm1 = Some((parse_f64(parts[1])?, parse_f64(parts[2])?));
            }
        }
        if line.starts_with("Sigmas[angstrom]") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 3 {
                global.lj_sigmas_angstrom = Some((parse_f64(parts[1])?, parse_f64(parts[2])?));
            }
        }
        if line.starts_with("Masses[amu]") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 3 {
                global.lj_masses_amu = Some((parse_f64(parts[1])?, parse_f64(parts[2])?));
            }
        }

        // Species blocks
        if first_token(line) == Some("Well") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 2 {
                return Err(format!("Malformed Well line: {}", line));
            }
            let name = parts[1].to_string();
            let (block, next) = collect_block(&lines, i);
            let rrho = parse_rrho_species(&block, &name)?;
            for l in &block {
                if first_token(l) == Some("PseudoFirstOrderRateConstant[1/sec]") {
                    let v = l.split_whitespace().nth(1).ok_or_else(|| format!("Malformed escape line: {l}"))?;
                    well_escape_rate_s_inv.insert(name.clone(), parse_f64(v)?);
                }
            }
            if wells.contains_key(&name) {
                return Err(format!("Input deck: well '{name}' defined twice."));
            }
            well_order.push(name.clone());
            wells.insert(name, rrho);
            i = next;
            continue;
        }

        if first_token(line) == Some("Bimolecular") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 2 {
                return Err(format!("Malformed Bimolecular line: {}", line));
            }
            let name = parts[1].to_string();
            let (block, next) = collect_block(&lines, i);

            // Parse the two fragments within the bimolecular block.
            let mut fragment_blocks: Vec<(String, Vec<String>)> = Vec::new();
            let mut j = 0usize;
            while j < block.len() {
                let l = &block[j];
                if first_token(l) == Some("Fragment") {
                    let p: Vec<&str> = l.split_whitespace().collect();
                    if p.len() < 2 {
                        return Err(format!("Malformed Fragment line: {}", l));
                    }
                    let frag_name = p[1].to_string();
                    let (frag_block, j_next) = collect_block(&block, j);
                    fragment_blocks.push((frag_name, frag_block));
                    j = j_next;
                    continue;
                }
                j += 1;
            }
            if fragment_blocks.len() != 2 {
                return Err(format!(
                    "Bimolecular '{}' must contain exactly two Fragment blocks (got {}).",
                    name,
                    fragment_blocks.len()
                ));
            }

            let frag_a = parse_fragment(&fragment_blocks[0].1, &fragment_blocks[0].0)?;
            let frag_b = parse_fragment(&fragment_blocks[1].1, &fragment_blocks[1].0)?;

            let mut ground_energy_cm1: Option<f64> = None;
            for l in &block {
                if first_token(l).map_or(false, |t| t.starts_with("GroundEnergy")) {
                    let unit = unit_tag(l).ok_or_else(|| format!("Bimolecular '{name}': GroundEnergy needs a unit"))?;
                    let value = l.split_whitespace().nth(1).ok_or_else(|| format!("Malformed line: {l}"))?;
                    ground_energy_cm1 = Some(energy_to_cm1(parse_f64(value)?, unit)?);
                }
            }
            let ground_energy_cm1 = ground_energy_cm1
                .ok_or_else(|| format!("Bimolecular '{}' missing GroundEnergy", name))?;

            bimolecular.insert(
                name.clone(),
                MessBimolecular {
                    name,
                    fragment_a: frag_a,
                    fragment_b: frag_b,
                    ground_energy_cm1,
                },
            );

            i = next;
            continue;
        }

        if first_token(line) == Some("Barrier") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 4 {
                return Err(format!("Malformed Barrier line: {}", line));
            }
            let name = parts[1].to_string();
            let left = parts[2].to_string();
            let right = parts[3].to_string();
            let (block, next) = collect_block(&lines, i);

            let core = parse_phasespace_core(&block)?.unwrap_or(MessBarrierCore::TightRrho);
            let geometry_required = matches!(core, MessBarrierCore::TightRrho);
            let rrho = parse_rrho_species_impl(&block, &name, geometry_required)?;
            let inverse_laplace_transform = parse_inverse_laplace_transform(&block, &name)?;
            let tunneling = parse_tunneling(&block, &name)?;

            barriers.push(MessBarrier {
                name,
                left,
                right,
                rrho,
                core,
                inverse_laplace_transform,
                tunneling,
            });

            i = next;
            continue;
        }

        i += 1;
    }

    if global.temperatures_kelvin.is_empty() || global.temperatures_kelvin.iter().any(|t| !(*t > 0.0)) {
        return Err("Input deck: TemperatureList[K] missing or not positive.".into());
    }
    if global.pressures_torr.is_empty() || global.pressures_torr.iter().any(|p| !(*p > 0.0)) {
        return Err("Input deck: PressureList[torr | atm | bar] missing or not positive.".into());
    }

    Ok(MessDeck {
        global,
        bimolecular,
        wells,
        barriers,
        well_escape_rate_s_inv,
        well_order,
    })
}

/// Convenience wrapper for parsing a MESS input from a file on disk.
pub fn parse_mess_input_file(path: impl AsRef<Path>) -> Result<MessDeck, String> {
    let text = std::fs::read_to_string(path.as_ref())
        .map_err(|e| format!("Failed to read MESS input file: {e}"))?;
    parse_mess_input(&text)
}

#[cfg(test)]
mod tests {
    use super::{parse_mess_input, MessBarrierCore};

    #[test]
    fn bimolecular_fragment_headers_do_not_require_end_blocks() {
        // In MESS, `Fragment <name>` is a header for the following RRHO block and does not
        // necessarily have its own `End`. Our block collector must therefore not treat
        // `Fragment` as a depth-increasing block starter.
        let deck = r#"
TemperatureList[K] 300.
PressureList[torr] 760.
EnergyStepOverTemperature 0.2
Model
  CollisionFrequency
    LennardJones
      Epsilons[1/cm] 417.0 33.4
      Sigmas[angstrom] 6.5 3.9
      Masses[amu] 149 28
    End
  End

  Bimolecular R
    Fragment A
      RRHO
        Geometry[angstrom] 1
        H 0 0 0
        Core RigidRotor
          SymmetryFactor 1.0
        End
        Frequencies[1/cm] 1
        100.0
        ZeroEnergy[1/cm] 0
        ElectronicLevels[1/cm] 1
          0 1
      End
    Fragment B
      RRHO
        Geometry[angstrom] 1
        H 0 0 0
        Core RigidRotor
          SymmetryFactor 1.0
        End
        Frequencies[1/cm] 1
        200.0
        ZeroEnergy[1/cm] 0
        ElectronicLevels[1/cm] 1
          0 1
      End
    GroundEnergy[kcal/mol] 0.0
  End
End
"#;

        let parsed = parse_mess_input(deck).expect("should parse");
        let r = parsed
            .bimolecular
            .get("R")
            .expect("should have bimolecular R");
        assert_eq!(r.fragment_a.name, "A");
        assert_eq!(r.fragment_b.name, "B");
        assert_eq!(r.fragment_a.vibrational_frequencies_cm1.len(), 1);
        assert_eq!(r.fragment_b.vibrational_frequencies_cm1.len(), 1);
    }

    #[test]
    fn well_with_escape_constant_block_parses() {
        // Many real MESS decks include `Escape Constant ... End` inside a Well block.
        // The block collector must not terminate the Well early when it sees that nested `End`.
        let deck = r#"
TemperatureList[K] 300.
PressureList[torr] 760.
EnergyStepOverTemperature 0.2
Model
  CollisionFrequency
    LennardJones
      Epsilons[1/cm] 417.0 33.4
      Sigmas[angstrom] 6.5 3.9
      Masses[amu] 149 28
    End
  End

Well W1
  Escape Constant
    PseudoFirstOrderRateConstant[1/sec]  2.5E7
  End
  Species
    RRHO
      Geometry[angstrom] 1
      H 0 0 0
      Core RigidRotor
        SymmetryFactor 1.0
      End
      Frequencies[1/cm] 1
      100.0
      ZeroEnergy[1/cm] 0
      ElectronicLevels[1/cm] 1
        0 1
    End
  End
End
"#;

        let parsed = parse_mess_input(deck).expect("should parse");
        let w1 = parsed.wells.get("W1").expect("should have W1");
        assert_eq!(w1.geometry_symbols.len(), 1);
        assert_eq!(w1.vibrational_frequencies_cm1.len(), 1);
    }

    #[test]
    fn consecutive_barriers_are_all_read() {
        // In this input format a `Barrier <name> <left> <right>` header has no `End` of its own:
        // the barrier ends with the `End` of the model block (RRHO) that follows it.
        let deck = r#"
TemperatureList[K] 300.
PressureList[torr] 760.
Model
Well W1
  Species
    RRHO
      Geometry[angstrom] 1
      H 0 0 0
      Core RigidRotor
        SymmetryFactor 1.0
      End
      Frequencies[1/cm] 1
      100.0
      ZeroEnergy[1/cm] 0
      ElectronicLevels[1/cm] 1
        0 1
    End
  End
Barrier B1 W1 P1
  RRHO
    Geometry[angstrom] 1
    H 0 0 0
    Core RigidRotor
      SymmetryFactor 1.0
    End
    Frequencies[1/cm] 1
    300.0
    ZeroEnergy[1/cm] 1000
    ElectronicLevels[1/cm] 1
      0 1
  End
Barrier B2 W1 P2
  RRHO
    Geometry[angstrom] 1
    H 0 0 0
    Core RigidRotor
      SymmetryFactor 1.0
    End
    Frequencies[1/cm] 1
    400.0
    ZeroEnergy[1/cm] 2000
    ElectronicLevels[1/cm] 1
      0 1
  End
End
"#;

        let parsed = parse_mess_input(deck).expect("should parse");
        let names: Vec<&str> = parsed.barriers.iter().map(|b| b.name.as_str()).collect();
        assert_eq!(names, vec!["B1", "B2"]);
        assert_eq!(parsed.barriers[1].right, "P2");
        assert_eq!(parsed.barriers[1].rrho.vibrational_frequencies_cm1, vec![400.0]);
    }

    #[test]
    fn phase_space_barrier_without_molecular_geometry_parses() {
        // A phase-space-theory barrier is described by the two fragment geometries of its core;
        // it has no molecular Geometry[angstrom] block.
        let deck = r#"
TemperatureList[K] 300.
PressureList[torr] 760.
Model
Barrier B0 R W1
  RRHO
    Stoichiometry C1O2
    Core PhaseSpaceTheory
      FragmentGeometry[angstrom] 1
      C 0.0 0.0 0.0
      FragmentGeometry[angstrom] 2
      O 0.0 0.0 0.0
      O 0.0 0.0 1.2
      SymmetryFactor 2.0
      PotentialPrefactor[au] 2.4
      PotentialPowerExponent 6.
    End
    Frequencies[1/cm] 1
    1585.0
    ZeroEnergy[kcal/mol] 0.0
    ElectronicLevels[1/cm] 1
      0 3
  End
End
"#;

        let parsed = parse_mess_input(deck).expect("should parse");
        assert_eq!(parsed.barriers.len(), 1);
        assert!(parsed.barriers[0].rrho.geometry_symbols.is_empty());
        assert!(matches!(parsed.barriers[0].core, MessBarrierCore::PhaseSpaceTheory { .. }));
    }

    #[test]
    fn temperature_and_pressure_lists_are_read_in_torr() {
        let deck = "TemperatureList[K] 300. 400. 500.\nPressureList[atm] 1 10\nModel\nEnd\n";
        let parsed = parse_mess_input(deck).expect("should parse");
        assert_eq!(parsed.global.temperatures_kelvin, vec![300.0, 400.0, 500.0]);
        assert_eq!(parsed.global.pressures_torr, vec![760.0, 7600.0]);
        let deck = "TemperatureList[K] 300.\nPressureList[bar] 1\nModel\nEnd\n";
        let parsed = parse_mess_input(deck).expect("should parse");
        assert!((parsed.global.pressures_torr[0] - 750.061_682_704).abs() < 1e-6);
    }

    #[test]
    fn wells_keep_the_order_of_the_input_deck() {
        let second = ONE_WELL_DECK.replace("Well W1", "Well A0");
        let deck = format!("{}{}", ONE_WELL_DECK.trim_end().trim_end_matches("End"), &second[second.find("Well A0").unwrap()..]);
        let parsed = parse_mess_input(&deck).expect("should parse");
        assert_eq!(parsed.well_order, vec!["W1".to_string(), "A0".to_string()]);
    }

    #[test]
    fn atom_fragments_and_ground_energy_in_wavenumbers_are_read() {
        let deck = r#"
TemperatureList[K] 1000.
PressureList[atm] 1.
Model
  Bimolecular P1
    Fragment C2H2
      RRHO
        Geometry[angstrom] 2
        C 0 0 -0.6
        C 0 0 0.6
        Core RigidRotor
          SymmetryFactor 2
        End
        Frequencies[1/cm] 1
        2000.0
        ZeroEnergy[1/cm] 0
        ElectronicLevels[1/cm] 1
          0 1
      End
    Fragment H
      Atom
        Mass[amu]    1
        ElectronicLevels[1/cm]          1
                0       2
      End
    GroundEnergy[1/cm]                  350.0
  End
End
"#;
        let parsed = parse_mess_input(deck).expect("should parse");
        let p1 = &parsed.bimolecular["P1"];
        assert_eq!(p1.ground_energy_cm1, 350.0);
        let h = &p1.fragment_b;
        assert_eq!(h.atom_mass_amu, Some(1.0));
        assert_eq!(h.electronic_degeneracy_ground, 2.0);
        assert!(h.vibrational_frequencies_cm1.is_empty() && h.geometry_symbols.is_empty());
        assert_eq!(p1.fragment_a.atom_mass_amu, None);
    }

    #[test]
    fn keywords_that_begin_like_block_names_are_not_blocks() {
        // `WellCutoff` is a global keyword, not a `Well` block.
        let deck = ONE_WELL_DECK.replace("PressureList[torr] 760.\n", "PressureList[torr] 760.\nWellCutoff 10\n");
        let parsed = parse_mess_input(&deck).expect("should parse");
        assert_eq!(parsed.well_order, vec!["W1".to_string()]);
    }

    #[test]
    fn excess_energy_over_temperature_is_read() {
        let deck = ONE_WELL_DECK.replace("PressureList[torr] 760.\n", "PressureList[torr] 760.\nExcessEnergyOverTemperature 40\n");
        let parsed = parse_mess_input(&deck).expect("should parse");
        assert_eq!(parsed.global.excess_energy_over_temperature, Some(40.0));
    }

    #[test]
    fn exponent_cutoff_and_well_escape_rate_are_read() {
        let deck = ONE_WELL_DECK.replace("Well W1\n", "Well W1\n  Escape Constant\n    PseudoFirstOrderRateConstant[1/sec]  2.5E7\n  End\n");
        let parsed = parse_mess_input(&deck).expect("should parse");
        assert_eq!(parsed.global.exponent_cutoff, Some(15.0));
        assert_eq!(parsed.well_escape_rate_s_inv.get("W1"), Some(&2.5e7));
    }

    const ILT_BARRIER_DECK: &str = r#"
TemperatureList[K] 300.
PressureList[torr] 760.
Model
Barrier B0 R W1
  RRHO
    Stoichiometry C1O2
    Core PhaseSpaceTheory
      FragmentGeometry[angstrom] 1
      C 0.0 0.0 0.0
      FragmentGeometry[angstrom] 2
      O 0.0 0.0 0.0
      O 0.0 0.0 1.2
      SymmetryFactor 2.0
      PotentialPrefactor[au] 2.4
      PotentialPowerExponent 6.
    End
    InverseLaplaceTransform
      Direction                    Association
      PreExponential[cm^3/s]       6.0e-12
      TemperatureExponent          -0.5
      ReferenceTemperature[K]      298.0
      ActivationEnergy[kcal/mol]   0.1
    End
    Frequencies[1/cm] 1
    1585.0
    ZeroEnergy[kcal/mol] 0.0
    ElectronicLevels[1/cm] 1
      0 3
  End
Barrier B1 W1 P1
  RRHO
    Geometry[angstrom] 1
    H 0 0 0
    Core RigidRotor
      SymmetryFactor 1.0
    End
    Tunneling Eckart
      ImaginaryFrequency[1/cm] 1500
      WellDepth[kcal/mol] 20
      WellDepth[kcal/mol] 25
    End
    Frequencies[1/cm] 1
    300.0
    ZeroEnergy[1/cm] 1000
    ElectronicLevels[1/cm] 1
      0 1
  End
End
"#;

    #[test]
    fn inverse_laplace_transform_block_is_read_inside_a_barrier() {
        use crate::barrierless::ilt::ilt_barrierless::ModifiedArrhenius;
        let parsed = parse_mess_input(ILT_BARRIER_DECK).expect("should parse");
        assert_eq!(parsed.barriers.len(), 2);
        let ilt = parsed.barriers[0].inverse_laplace_transform.as_ref().expect("ILT block");
        assert_eq!(ilt.direction, super::IltDirection::Association);
        let expected = ModifiedArrhenius {
            pre_exponential: 6.0e-12,
            temperature_exponent: -0.5,
            reference_temperature_kelvin: 298.0,
            activation_energy_cm1: 0.1 / crate::constants::CM1_TO_KCAL,
        };
        assert_eq!(ilt.high_pressure_rate, expected);
        assert_eq!(parsed.barriers[0].rrho.vibrational_frequencies_cm1, vec![1585.0]);
        assert!(parsed.barriers[1].inverse_laplace_transform.is_none());
    }

    #[test]
    fn eckart_tunneling_parameters_are_read() {
        let parsed = parse_mess_input(ILT_BARRIER_DECK).expect("should parse");
        assert_eq!(parsed.barriers[0].tunneling, None);
        let kcal = crate::constants::CM1_TO_KCAL;
        assert_eq!(
            parsed.barriers[1].tunneling,
            Some(super::TunnelingSpecification::Eckart {
                imaginary_frequency_cm1: 1500.0,
                well_depths_cm1: [20.0 / kcal, 25.0 / kcal],
            })
        );
        assert_eq!(parsed.barriers[1].rrho.vibrational_frequencies_cm1, vec![300.0]);
    }

    #[test]
    fn other_tunneling_models_are_recorded_as_unsupported() {
        let deck = ILT_BARRIER_DECK.replace("Tunneling Eckart", "Tunneling Read");
        let parsed = parse_mess_input(&deck).expect("should parse");
        assert_eq!(parsed.barriers[1].tunneling, Some(super::TunnelingSpecification::Unsupported { model: "Read".into() }));
        let deck = ILT_BARRIER_DECK.replace("      WellDepth[kcal/mol] 25\n", "");
        assert!(parse_mess_input(&deck).is_err(), "an Eckart block needs two WellDepth values");
    }

    #[test]
    fn inverse_laplace_transform_units_must_match_the_direction() {
        let deck = ILT_BARRIER_DECK.replace("PreExponential[cm^3/s]", "PreExponential[1/s]");
        assert!(parse_mess_input(&deck).is_err());
        let deck = ILT_BARRIER_DECK.replace("Direction                    Association", "Direction Dissociation");
        assert!(parse_mess_input(&deck).is_err());
        let deck = ILT_BARRIER_DECK
            .replace("Direction                    Association", "Direction Dissociation")
            .replace("PreExponential[cm^3/s]", "PreExponential[1/s]");
        let parsed = parse_mess_input(&deck).expect("dissociation with 1/s parses");
        assert_eq!(parsed.barriers[0].inverse_laplace_transform.as_ref().unwrap().direction, super::IltDirection::Dissociation);
    }

    const ONE_WELL_DECK: &str = r#"
TemperatureList[K] 300.
PressureList[torr] 760.
Model
  EnergyRelaxation
    Exponential
      Factor[1/cm] 200.0
      Power 0.85
      ExponentCutoff 15
  End
  CollisionFrequency
    LennardJones
      Epsilons[1/cm] 417.0 33.4
      Sigmas[angstrom] 6.5 3.9
      Masses[amu] 149 28
  End
Well W1
  Species
    RRHO
      Geometry[angstrom] 1
      H 0 0 0
      Core RigidRotor
        SymmetryFactor 1.0
      End
      Frequencies[1/cm] 1
      100.0
      ZeroEnergy[1/cm] 0
      ElectronicLevels[1/cm] 1
        0 1
    End
  End
End
"#;

}

pub(crate) fn rotational_constants_from_geometry_cm1(
    symbols: &[String],
    coords_angstrom: &[[f64; 3]],
) -> Result<Vec<f64>, String> {
    if symbols.is_empty() || coords_angstrom.is_empty() || symbols.len() != coords_angstrom.len() {
        return Err("Geometry must have matching non-empty symbols and coordinates.".into());
    }

    if symbols.len() == 1 {
        return Ok(Vec::new());
    }

    let masses_amu = mass_vector_from_symbols_amu(symbols)?;
    let coords: Vec<[f64; 3]> = coords_angstrom.to_vec();
    let brot = get_brot(&coords, &masses_amu);

    let mut finite: Vec<f64> = brot
        .into_iter()
        .filter(|b| b.is_finite() && *b > 0.0 && *b < 1.0e6)
        .collect();

    // Normalize common representations:
    // - if inertia returns two equal components for linear rotor, prefer [B]
    if finite.len() == 2 {
        let b = 0.5 * (finite[0] + finite[1]);
        return Ok(vec![b]);
    }

    if finite.len() >= 3 {
        finite.truncate(3);
        return Ok(finite);
    }

    Err("Failed to compute usable rotational constants from geometry.".into())
}
