//! Minimal MESS input parser (subset) and conversion to MarXus multiwell inputs.
//!
//! Scope (intentionally limited)
//! ----------------------------
//! This parser is designed to cover the common subset needed to:
//! - read wells and barriers from a MESS deck
//! - compute RRHO microcanonical inputs (ρ(E), W‡(E), k(E)) using existing MarXus routines
//! - build a MarXus multiwell network topology (`Vec<WellDefinition>`) and per-channel models
//!
//! Supported (for the provided decks):
//! - TemperatureList[K], PressureList[torr], EnergyStepOverTemperature, ModelEnergyLimit[kcal/mol]
//! - Model -> EnergyRelaxation -> Exponential (Factor/Power)
//! - Model -> CollisionFrequency -> LennardJones (Epsilons/Sigmas/Masses; uses the first/second as
//!   "species/bath", and constructs a single reduced-mass + sigma/epsilon set)
//! - Species blocks: Bimolecular, Well, Barrier
//!   - RRHO -> Geometry[angstrom] N (computes rotational constants via `inertia::get_brot`)
//!   - RRHO -> Core -> RigidRotor -> SymmetryFactor
//!   - RRHO -> Frequencies[1/cm] N
//!   - RRHO -> ZeroEnergy[kcal/mol] or ZeroEnergy[1/cm]
//!   - RRHO -> ElectronicLevels[1/cm] N (uses the *ground-level degeneracy* as a simple factor)
//!   - Barrier -> Core PhaseSpaceTheory (fragment geometries + V0 + n + symmetry)
//!
//! Not yet supported:
//! - multiple temperatures/pressures (we read the first value)
//! - full electronic level handling (excited electronic energies as shifted state manifolds)
//! - tunneling models (Eckart/Wigner); currently ignored in microcanonical building
//! - advanced MESS options like hindered rotors, adiabatic channels, etc.
//!
//! The goal is to let you *start* from a MESS deck and obtain a consistent multiwell micro model.

use std::collections::HashMap;
use std::path::Path;

use crate::barrierless::phasespace::phase_space_theory::PhaseSpaceTheoryModel;
use crate::barrierless::phasespace::types::{
    CaptureFragment, CaptureFragmentRotorModel, PhaseSpaceTheoryInput,
};
use crate::inertia::inertia::get_brot;
use crate::masterequation::graph_utils::DisjointSetUnion;
use crate::masterequation::microcanonical_builder::{
    ChannelMicroModel, SpeciesMicroModel, TransitionStateModel,
};
use crate::masterequation::reaction_network::{
    CollisionKernelImplementation, CollisionModelParams, MasterEquationSettings, ReactionChannel,
    WellDefinition,
};
use crate::utils::atomic_masses::mass_vector_from_symbols_amu;

const KB_CM1_PER_K: f64 = 0.695_034_76;
const KCAL_PER_MOL_PER_CM1: f64 = 2.859_144e-3;

#[derive(Clone, Debug)]
pub struct MessBuildOptions {
    /// Override ΔE (cm^-1). If None, computed from EnergyStepOverTemperature: ΔE = step_over_T * kB * T.
    pub energy_grain_width_cm1: Option<f64>,
    /// Maximum energy (kcal/mol) included above each well minimum (local grain axis).
    pub max_energy_kcal_mol: f64,
    /// Collision band half width (grains).
    pub collision_band_half_width: usize,
    /// Nonreactive grains at bottom (k(E)=0 assumed below).
    pub nonreactive_grain_count: usize,
}

impl Default for MessBuildOptions {
    fn default() -> Self {
        Self {
            energy_grain_width_cm1: Some(20.0),
            max_energy_kcal_mol: 50.0,
            collision_band_half_width: 20,
            nonreactive_grain_count: 10,
        }
    }
}

#[derive(Clone, Debug)]
pub struct MessGlobal {
    pub temperature_kelvin: f64,
    pub pressure_torr: f64,
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

#[derive(Clone, Debug)]
pub struct MessBarrier {
    pub name: String,
    pub left: String,
    pub right: String,
    pub rrho: MessSpeciesRrho,
    pub core: MessBarrierCore,
}

#[derive(Clone, Debug)]
pub struct MessDeck {
    pub global: MessGlobal,
    pub bimolecular: HashMap<String, MessBimolecular>,
    pub wells: HashMap<String, MessSpeciesRrho>,
    pub barriers: Vec<MessBarrier>,
}

pub struct MessBuiltNetwork {
    pub settings: MasterEquationSettings,
    pub wells: Vec<WellDefinition>,
    pub well_models_by_name: HashMap<String, SpeciesMicroModel>,
    pub channel_models_by_name: HashMap<String, ChannelMicroModel>,
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
        Ok(value / KCAL_PER_MOL_PER_CM1)
    } else if u.contains("1/cm") || u.contains("cm") {
        Ok(value)
    } else {
        Err(format!("Unsupported energy unit tag '{}'", unit_tag))
    }
}

fn is_block_starter(tok: &str) -> bool {
    matches!(
        tok,
        "Model"
            | "Bimolecular"
            | "Well"
            | "Barrier"
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

fn parse_rrho_species(block: &[String], name: &str) -> Result<MessSpeciesRrho, String> {
    let (symbols, coords) = parse_geometry(block, "Geometry[angstrom]")?
        .ok_or_else(|| format!("RRHO species '{}' missing Geometry[angstrom]", name))?;
    let symmetry_factor = parse_symmetry_factor(block)?;
    let vib = parse_frequencies(block)?;
    let zero_energy_cm1 = parse_zero_energy_cm1(block)?;
    let electronic_degeneracy_ground = parse_electronic_degeneracy_ground(block)?;

    Ok(MessSpeciesRrho {
        name: name.to_string(),
        geometry_symbols: symbols,
        geometry_angstrom: coords,
        symmetry_factor,
        vibrational_frequencies_cm1: vib,
        zero_energy_cm1,
        electronic_degeneracy_ground,
    })
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
        temperature_kelvin: 0.0,
        pressure_torr: 0.0,
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
                "TemperatureList[K]" => global.temperature_kelvin = parse_f64(v)?,
                "PressureList[torr]" => global.pressure_torr = parse_f64(v)?,
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
        if line.starts_with("Well") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 2 {
                return Err(format!("Malformed Well line: {}", line));
            }
            let name = parts[1].to_string();
            let (block, next) = collect_block(&lines, i);
            let rrho = parse_rrho_species(&block, &name)?;
            wells.insert(name, rrho);
            i = next;
            continue;
        }

        if line.starts_with("Bimolecular") {
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
                if l.starts_with("Fragment") {
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

            let frag_a = parse_rrho_species(&fragment_blocks[0].1, &fragment_blocks[0].0)?;
            let frag_b = parse_rrho_species(&fragment_blocks[1].1, &fragment_blocks[1].0)?;

            let mut ground_energy_cm1: Option<f64> = None;
            for l in &block {
                if l.starts_with("GroundEnergy[kcal/mol]") {
                    let parts: Vec<&str> = l.split_whitespace().collect();
                    ground_energy_cm1 = Some(energy_to_cm1(
                        parse_f64(parts.last().unwrap())?,
                        "kcal/mol",
                    )?);
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

        if line.starts_with("Barrier") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 4 {
                return Err(format!("Malformed Barrier line: {}", line));
            }
            let name = parts[1].to_string();
            let left = parts[2].to_string();
            let right = parts[3].to_string();
            let (block, next) = collect_block(&lines, i);

            let rrho = parse_rrho_species(&block, &name)?;
            let core = parse_phasespace_core(&block)?.unwrap_or(MessBarrierCore::TightRrho);

            barriers.push(MessBarrier {
                name,
                left,
                right,
                rrho,
                core,
            });

            i = next;
            continue;
        }

        i += 1;
    }

    if global.temperature_kelvin <= 0.0 {
        return Err("MESS input missing TemperatureList[K]".into());
    }
    if global.pressure_torr <= 0.0 {
        return Err("MESS input missing PressureList[torr]".into());
    }

    Ok(MessDeck {
        global,
        bimolecular,
        wells,
        barriers,
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
    use super::parse_mess_input;

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
}

impl MessDeck {
    /// Convert to a MarXus multiwell network + per-well and per-channel micro models.
    ///
    /// This builds:
    /// - `WellDefinition` list containing ONLY MESS `Well` species
    /// - Channels for:
    ///   - well<->well barriers: creates both directions
    ///   - well->sink barriers: creates a sink channel from the well
    ///
    /// Bimolecular species (including the MESS "Reactant") are treated as sinks (not wells).
    pub fn build_multiwell_network(
        &self,
        options: MessBuildOptions,
    ) -> Result<MessBuiltNetwork, String> {
        let temperature = self.global.temperature_kelvin;
        let pressure_torr = self.global.pressure_torr;

        let energy_grain_width_cm1 = if let Some(v) = options.energy_grain_width_cm1 {
            v
        } else {
            let step_over_t = self.global.energy_step_over_temperature.ok_or(
                "Missing EnergyStepOverTemperature (or set energy_grain_width_cm1 override).",
            )?;
            step_over_t * KB_CM1_PER_K * temperature
        };
        if energy_grain_width_cm1 <= 0.0 || !energy_grain_width_cm1.is_finite() {
            return Err("Computed energy grain width is invalid.".into());
        }

        // Global collision model
        let (m1, m2) = self.global.lj_masses_amu.unwrap_or((0.0, 0.0));
        let reduced_mass_amu = if m1 > 0.0 && m2 > 0.0 {
            (m1 * m2) / (m1 + m2)
        } else {
            0.0
        };
        let (eps1_cm1, _eps2_cm1) = self.global.lj_epsilons_cm1.unwrap_or((0.0, 0.0));
        let epsilon_kelvin = eps1_cm1 * 1.438_776_877;
        let (sigma1_a, _sigma2_a) = self.global.lj_sigmas_angstrom.unwrap_or((0.0, 0.0));

        let default_collision = CollisionModelParams {
            lennard_jones_sigma_angstrom: sigma1_a,
            lennard_jones_epsilon_kelvin: epsilon_kelvin,
            reduced_mass_amu,
            alpha_at_1000K_cm1: self.global.alpha_factor_cm1.unwrap_or(200.0),
            alpha_temperature_exponent: self.global.alpha_power.unwrap_or(0.85),
        };

        let settings = MasterEquationSettings {
            temperature_kelvin: temperature,
            pressure_torr,
            boltzmann_constant_wavenumber_per_kelvin: KB_CM1_PER_K,
            collision_band_half_width: options.collision_band_half_width,
            collision_kernel_implementation: CollisionKernelImplementation::Mess,
            outgoing_rate_threshold: 0.0,
            internal_rate_threshold: 0.0,
            bathgas_number_density_prefactor: 3.262e16,
            mean_speed_prefactor: 1.0e4,
            enforce_interwell_detailed_balance: false,
            linear_solver: crate::masterequation::reaction_network::MultiwellLinearSolver::Gmres,
            krylov_tolerance: 1e-10,
            krylov_max_iter: 8000,
            gmres_restart: 50,
        };

        // Order wells deterministically (sorted by name).
        let mut well_names: Vec<String> = self.wells.keys().cloned().collect();
        well_names.sort();

        // Compute alignment offsets from the MESS absolute energies (relative reference):
        // offset_w = round( E_well_min / ΔE ), so that absolute_energy_cm1(w, i) = (i + offset_w)*ΔE
        let mut wells_out: Vec<WellDefinition> = Vec::with_capacity(well_names.len());
        for name in &well_names {
            let sp = &self.wells[name];
            let offset = (sp.zero_energy_cm1 / energy_grain_width_cm1).round() as isize;

            let max_energy_cm1 = options.max_energy_kcal_mol / KCAL_PER_MOL_PER_CM1;
            let n_bins = (max_energy_cm1 / energy_grain_width_cm1).ceil().max(1.0) as usize + 1;

            wells_out.push(WellDefinition {
                well_name: name.clone(),
                energy_grain_width_cm1: energy_grain_width_cm1,
                lowest_included_grain_index: 0,
                one_past_highest_included_grain_index: n_bins,
                alignment_offset_in_grains: offset,
                nonreactive_grain_count: options.nonreactive_grain_count,
                collision_params: default_collision.clone(),
                channels: Vec::new(),
            });
        }

        let well_index_by_name: HashMap<String, usize> = wells_out
            .iter()
            .enumerate()
            .map(|(i, w)| (w.well_name.clone(), i))
            .collect();

        // Build micro models for wells.
        let mut well_models_by_name: HashMap<String, SpeciesMicroModel> = HashMap::new();
        for name in &well_names {
            let sp = &self.wells[name];
            let brot = rotational_constants_from_geometry_cm1(
                &sp.geometry_symbols,
                &sp.geometry_angstrom,
            )?;
            well_models_by_name.insert(
                name.clone(),
                SpeciesMicroModel {
                    name: name.clone(),
                    vibrational_frequencies_cm1: sp.vibrational_frequencies_cm1.clone(),
                    rotational_constants_cm1: brot,
                    symmetry_number: sp.symmetry_factor,
                    chirality_number: 1.0,
                    electronic_degeneracy: sp.electronic_degeneracy_ground,
                },
            );
        }

        // Validate barrier connectivity and build per-species adjacency for efficient construction.
        let mut barrier_indices_by_species: HashMap<String, Vec<usize>> = HashMap::new();
        let mut first_barrier_by_unordered_pair: HashMap<(String, String), String> = HashMap::new();
        for (barrier_index, barrier) in self.barriers.iter().enumerate() {
            if barrier.left == barrier.right {
                return Err(format!(
                    "Barrier '{}' connects '{}' to itself; this is not supported.",
                    barrier.name, barrier.left
                ));
            }

            // MESS disallows multiple barriers connecting the same unordered pair of species.
            let (a, b) = if barrier.left <= barrier.right {
                (barrier.left.clone(), barrier.right.clone())
            } else {
                (barrier.right.clone(), barrier.left.clone())
            };
            let key = (a, b);
            if let Some(existing) = first_barrier_by_unordered_pair.get(&key) {
                return Err(format!(
                    "Multiple barriers connect the same species pair ({} <-> {}): '{}' and '{}'.",
                    key.0, key.1, existing, barrier.name
                ));
            }
            first_barrier_by_unordered_pair.insert(key, barrier.name.clone());

            barrier_indices_by_species
                .entry(barrier.left.clone())
                .or_default()
                .push(barrier_index);
            barrier_indices_by_species
                .entry(barrier.right.clone())
                .or_default()
                .push(barrier_index);
        }

        // Ensure the well graph (using only well-well barriers) is connected, like MESS does.
        if wells_out.len() > 1 {
            let mut dsu = DisjointSetUnion::new(wells_out.len());
            for barrier in &self.barriers {
                let left_well = well_index_by_name.get(&barrier.left).copied();
                let right_well = well_index_by_name.get(&barrier.right).copied();
                if let (Some(i), Some(j)) = (left_well, right_well) {
                    dsu.union(i, j);
                }
            }
            let roots = dsu.component_roots();
            if roots.len() > 1 {
                return Err(format!(
                    "Wells are not all connected by well–well barriers (found {} disconnected components).",
                    roots.len()
                ));
            }
        }

        // Build channels and per-channel micro models keyed by channel name.
        let mut channel_models_by_name: HashMap<String, ChannelMicroModel> = HashMap::new();

        for from_idx in 0..wells_out.len() {
            let from_name = wells_out[from_idx].well_name.clone();

            // Find all barriers touching this well (either well<->well or well->sink).
            let touching = barrier_indices_by_species
                .get(&from_name)
                .cloned()
                .unwrap_or_default();
            for barrier_index in touching {
                let b = &self.barriers[barrier_index];
                let other = if b.left == from_name {
                    b.right.clone()
                } else {
                    b.left.clone()
                };

                let connected_well_index = well_index_by_name.get(&other).copied();

                // Add channel in topology.
                let ch_name = format!("{}_to_{}", from_name, other);
                wells_out[from_idx].channels.push(ReactionChannel {
                    name: ch_name.clone(),
                    connected_well_index,
                });

                // Threshold energy relative to the FROM well minimum:
                // E0 = E_TS - E_from_well
                let e_from = self.wells.get(&from_name).ok_or_else(|| {
                    format!(
                        "Internal error: well '{}' not found while building barrier thresholds",
                        from_name
                    )
                })?;
                let e0_cm1 = (b.rrho.zero_energy_cm1 - e_from.zero_energy_cm1).max(0.0);

                // Build TS model:
                let ts_model = match &b.core {
                    MessBarrierCore::TightRrho => {
                        let brot_ts = rotational_constants_from_geometry_cm1(
                            &b.rrho.geometry_symbols,
                            &b.rrho.geometry_angstrom,
                        )?;
                        TransitionStateModel::TightRRHO {
                            species: SpeciesMicroModel {
                                name: b.name.clone(),
                                vibrational_frequencies_cm1: b
                                    .rrho
                                    .vibrational_frequencies_cm1
                                    .clone(),
                                rotational_constants_cm1: brot_ts,
                                symmetry_number: b.rrho.symmetry_factor,
                                chirality_number: 1.0,
                                electronic_degeneracy: b.rrho.electronic_degeneracy_ground,
                            },
                        }
                    }
                    MessBarrierCore::PhaseSpaceTheory {
                        fragment_a_geometry_symbols,
                        fragment_a_geometry_angstrom,
                        fragment_b_geometry_symbols,
                        fragment_b_geometry_angstrom,
                        symmetry_operations,
                        potential_prefactor_au,
                        potential_power_exponent,
                    } => {
                        let pst = PhaseSpaceTheoryModel::new(PhaseSpaceTheoryInput {
                            fragment_a: CaptureFragment {
                                mass_amu: None,
                                rotor: CaptureFragmentRotorModel::GeometryAngstrom {
                                    symbols: fragment_a_geometry_symbols.clone(),
                                    coordinates_angstrom: fragment_a_geometry_angstrom.clone(),
                                },
                            },
                            fragment_b: CaptureFragment {
                                mass_amu: None,
                                rotor: CaptureFragmentRotorModel::GeometryAngstrom {
                                    symbols: fragment_b_geometry_symbols.clone(),
                                    coordinates_angstrom: fragment_b_geometry_angstrom.clone(),
                                },
                            },
                            symmetry_operations: *symmetry_operations,
                            potential_prefactor_au: *potential_prefactor_au,
                            potential_power_exponent: *potential_power_exponent,
                        })?;

                        TransitionStateModel::PhaseSpaceTheoryRRHO {
                            pst_core: pst,
                            vibrational_frequencies_cm1: b.rrho.vibrational_frequencies_cm1.clone(),
                            electronic_degeneracy: b.rrho.electronic_degeneracy_ground,
                        }
                    }
                };

                channel_models_by_name.insert(
                    ch_name,
                    ChannelMicroModel {
                        threshold_energy_cm1: e0_cm1,
                        transition_state: ts_model,
                    },
                );
            }
        }

        Ok(MessBuiltNetwork {
            settings,
            wells: wells_out,
            well_models_by_name,
            channel_models_by_name,
        })
    }
}

fn rotational_constants_from_geometry_cm1(
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
