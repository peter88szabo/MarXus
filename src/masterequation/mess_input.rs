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
//! - MarXus extension: `MarXus ... End` in the header, the solution method and its settings
//!   (`parse_marxus_header`, `solution_method.rs`).
//!
//!   - RRHO -> Tunneling Eckart (ImaginaryFrequency, two WellDepth values); other models recorded
//! Not read: excited electronic levels and other model types. Internal rotors (`Rotor` blocks) are refused with
//! an error.

use crate::constants::CM1_TO_KCAL;
use crate::rrkm::internal_rotor::{
    potential_from_equidistant_points, potential_from_fourier_expansion, TorsionalPotential, DEFAULT_BASIS_SIZE, MAX_BASIS_SIZE,
};
use std::collections::HashMap;
use std::path::Path;

use crate::barrierless::ilt::ilt_barrierless::ModifiedArrhenius;
use crate::barrierless::phasespace::types::PstTstLevel;
use crate::inertia::inertia::get_brot;
use crate::utils::atomic_masses::mass_vector_from_symbols_amu;

use super::solution_method::{
    chemical_subspace_criterion_from_keyword, collision_integral_from_keyword, rotor_reduced_moment_from_keyword,
    eigen_solver_from_keyword, initial_state_from_keyword, integrator_from_keyword, SolutionMethod,
    SolutionSettings, STEADY_STATE_KEYWORD_REPLACED,
};

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
    /// ChemicalEigenvalueMax: CSE species merging, chemical eigenvalues <= this x the lowest relaxation
    /// eigenvalue (only MESS's absolute mode, 0 < value < 1).
    pub chemical_eigenvalue_max: Option<f64>,
    /// WellProjectionThreshold: primary wells of the CSE species partition.
    pub well_projection_threshold: Option<f64>,

    pub alpha_factor_cm1: Option<f64>,
    pub alpha_power: Option<f64>,

    pub lj_epsilons_cm1: Option<(f64, f64)>,
    pub lj_sigmas_angstrom: Option<(f64, f64)>,
    pub lj_masses_amu: Option<(f64, f64)>,

    pub reactant_name: Option<String>,
    pub excess_reactant_concentration_cm3: Option<f64>,

    /// Solution method and its settings from the MarXus header block (`parse_marxus_header`).
    pub solution: SolutionSettings,
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
    /// All electronic levels (energy above the lowest level in cm-1, degeneracy), from `ElectronicLevels[unit] N`; the
    /// first is (0, electronic_degeneracy_ground). Without the block: [(0, 1)].
    pub electronic_levels: Vec<(f64, f64)>,
    /// Mass of an `Atom` fragment (amu); None for RRHO species (mass from the geometry).
    pub atom_mass_amu: Option<f64>,
    /// Rotational constants (cm-1) given in place of a geometry (`RotationalConstants[1/cm] N`, N = 3, or 1 for a
    /// linear rotor; MarXus extension of the deck format for species data given that way).
    pub rotational_constants_cm1: Option<Vec<f64>>,
    /// Mass (amu) of an RRHO species given by its rotational constants (`Mass[amu]`).
    pub mass_amu: Option<f64>,
    /// One-dimensional internal rotors (`Rotor Hindered`, `Rotor Free`); their torsions are not among the
    /// Frequencies.
    pub internal_rotors: Vec<MessInternalRotor>,
}

/// Kind of an internal rotor of the deck.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MessRotorKind {
    Hindered,
    Free,
}

/// One-dimensional internal rotor of an RRHO species, MESS syntax (one block per rotor, beside the Core):
///
///   Rotor Hindered                  (or Free: no potential)
///     HamiltonSizeMin  999          (optional: smallest Fourier basis, odd; default 999)
///     HamiltonSizeMax  1999         (optional: largest Fourier basis, odd; default 1999)
///     Group            5 6 7        (rotating atoms, 1-based in the geometry)
///     Axis             1 2          (two atoms on the rotation axis)
///     Symmetry         3            (rotor symmetry number; default 1)
///     Potential[kcal/mol] N         (N equidistant points on [0, 360/Symmetry), the first a minimum;
///       V_1 ... V_N                  or FourierExpansion[kcal/mol] n followed by n lines "index value",
///   End                              the coefficients c_0, a_1, b_1, a_2, ... in that order)
///
/// MarXus extension: `RotationalConstant[1/cm] B` gives B of the rotor, in place of the reduced moment from the
/// geometry; it is required for a species given by its rotational constants. `GridSize` and `ThermalPowerMax` are
/// accepted and have no effect.
#[derive(Clone, Debug, PartialEq)]
pub struct MessInternalRotor {
    pub kind: MessRotorKind,
    /// Rotating atoms, 0-based.
    pub group: Vec<usize>,
    /// The two atoms on the axis, 0-based.
    pub axis: (usize, usize),
    pub symmetry: u32,
    /// Torsional potential (cm-1) in x = Symmetry phi; zero for a free rotor.
    pub potential: TorsionalPotential,
    pub hamilton_size_min: usize,
    pub hamilton_size_max: usize,
    /// Rotational constant (cm-1) given in the deck.
    pub rotational_constant_cm1: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct MessBimolecular {
    pub name: String,
    pub fragment_a: MessSpeciesRrho,
    pub fragment_b: MessSpeciesRrho,
    /// Bimolecular asymptote energy (cm^-1) relative to the same reference used for wells/barriers.
    pub ground_energy_cm1: f64,
    /// MarXus `FragmentWell <name>`: the fragment (named as the well) that is a well of the network; the other fragment
    /// is a partner in excess (reports/nonthermal_sources_design.md, Section 15.1).
    pub fragment_well: Option<String>,
    /// MarXus `PartnerConcentration[molecule/cm^3]`: the concentration of that partner (or of the excess fragment of a
    /// lumped state).
    pub partner_concentration_cm3: Option<f64>,
    /// MarXus `LumpedState`: the pair is one thermal state of the master equation, pseudo-first-order in its excess
    /// fragment (reports/nonthermal_sources_design.md, Section 15.2).
    pub lumped_state: bool,
    /// MarXus `ExcessFragment <name>`: the fragment in excess of a lumped state.
    pub excess_fragment: Option<String>,
}

/// MarXus block `FragmentEnergy <kind> ... End` in a barrier to a bimolecular species with a fragment well: the energy
/// partitioning P(e | X) of the fragment (`fragment_partition.rs`).
///
///   FragmentEnergy TwoPieceGaussian
///     Mu[1/cm]          -2685.46  0.50158       (intercept [unit], gradient in X)
///     SigmaLeft[1/cm]   -132.008  0.053485
///     SigmaRight[1/cm]  -1122.11  0.13401
///   End
///   FragmentEnergy ModifiedPrior   (Order, TemperatureExponent, ReferenceTemperature[K])
///   FragmentEnergy Prior
#[derive(Clone, Debug, PartialEq)]
pub enum FragmentEnergySpecification {
    Prior,
    ModifiedPrior { order: f64, temperature_exponent: f64, reference_temperature_kelvin: f64 },
    /// [intercept in cm-1, gradient] of mu, sigma_L and sigma_R.
    TwoPieceGaussian { mu: [f64; 2], sigma_left: [f64; 2], sigma_right: [f64; 2] },
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
        /// `TSTLevel` (T, E, EJ, J=0); default EJ, as in MESS.
        tst_level: PstTstLevel,
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
    /// Energy partitioning of the fragment well (MarXus `FragmentEnergy` block), if given.
    pub fragment_energy: Option<FragmentEnergySpecification>,
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
    /// MarXus block in a Well: its own collision parameters in place of the global model.
    pub well_collision: HashMap<String, WellCollisionOverride>,
    /// Well names in the order of the input deck.
    pub well_order: Vec<String>,
    /// `Dummy` bimolecular species (as in MESS): products without molecular data, in deck order.
    pub dummy_bimolecular: Vec<String>,
    /// The `Preparation ... End` block (a MarXus extension: initial population, sources, bath history), comments
    /// stripped; interpreted by `preparation_input::parse_preparation`.
    pub preparation_block: Option<Vec<String>>,
}


/// MarXus block `MarXus ... End` in a Well: collision parameters of that well in place of the global ones of the Model
/// (MESS has one model for all wells). Each given pair replaces the global pair (bath, complex):
///
///   MarXus
///     Factor[1/cm]              98.3      (<dE_down> at the reference temperature)
///     Power                     1.0       (<dE_down> ~ (T/T_ref)^Power)
///     ReferenceTemperature[K]   295       (MarXus; the global model refers to 300 K)
///     Epsilons[1/cm]            57.0 150.2
///     Sigmas[angstrom]          3.74 4.6
///     Masses[amu]               28.0 57.0
///   End
#[derive(Clone, Debug, Default, PartialEq)]
pub struct WellCollisionOverride {
    pub factor_cm1: Option<f64>,
    pub power: Option<f64>,
    pub reference_temperature_kelvin: Option<f64>,
    pub epsilons_cm1: Option<(f64, f64)>,
    pub sigmas_angstrom: Option<(f64, f64)>,
    pub masses_amu: Option<(f64, f64)>,
}

/// The MarXus block of a Well (its lines removed from the block) and the parameters in it.
fn split_well_collision_block(block: &[String], well: &str) -> Result<(Vec<String>, Option<WellCollisionOverride>), String> {
    let Some(start) = block.iter().position(|l| first_token(l) == Some("MarXus")) else {
        return Ok((block.to_vec(), None));
    };
    let end = block[start..].iter().position(|l| first_token(l) == Some("End")).map(|k| start + k).ok_or_else(|| format!("Well '{well}': MarXus block without End."))?;
    let context = |what: &str| format!("Well '{well}', MarXus block: {what}");
    let mut o = WellCollisionOverride::default();
    for line in &block[start + 1..end] {
        let key = first_token(line).unwrap_or("");
        let values = line.split_whitespace().skip(1).map(parse_f64).collect::<Result<Vec<f64>, String>>()?;
        let one = || values.first().copied().ok_or_else(|| context(&format!("{key} needs a value")));
        let two = || if values.len() == 2 { Ok((values[0], values[1])) } else { Err(context(&format!("{key} needs two values (bath, complex)"))) };
        match key {
            "Factor[1/cm]" => o.factor_cm1 = Some(one()?),
            "Power" => o.power = Some(one()?),
            "ReferenceTemperature[K]" => o.reference_temperature_kelvin = Some(one()?),
            "Epsilons[1/cm]" => o.epsilons_cm1 = Some(two()?),
            "Sigmas[angstrom]" => o.sigmas_angstrom = Some(two()?),
            "Masses[amu]" => o.masses_amu = Some(two()?),
            other => {
                return Err(context(&format!(
                    "unknown keyword '{other}' (Factor[1/cm], Power, ReferenceTemperature[K], Epsilons[1/cm], Sigmas[angstrom], Masses[amu])"
                )))
            }
        }
    }
    let rest = block[..start].iter().chain(&block[end + 1..]).cloned().collect();
    Ok((rest, Some(o)))
}

pub(crate) fn strip_comment(mut line: &str) -> &str {
    if let Some(idx) = line.find('#') {
        line = &line[..idx];
    }
    if let Some(idx) = line.find('!') {
        line = &line[..idx];
    }
    line.trim()
}

pub(crate) fn first_token(line: &str) -> Option<&str> {
    line.split_whitespace().next()
}

pub(crate) fn parse_f64(raw: &str) -> Result<f64, String> {
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

pub(crate) fn energy_to_cm1(value: f64, unit_tag: &str) -> Result<f64, String> {
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
            | "FragmentEnergy"
            | "MarXus"
            | "Atom"
            | "Rotor"
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

/// `ElectronicLevels[unit] N` followed by N lines "energy degeneracy" (unit 1/cm by default, also kcal/mol, kJ/mol):
/// the levels sorted by energy; the lowest must lie at 0 (the zero of the species). Without the block: [(0, 1)].
fn parse_electronic_levels(block: &[String], name: &str) -> Result<Vec<(f64, f64)>, String> {
    let Some(idx) = block.iter().position(|l| l.starts_with("ElectronicLevels")) else {
        return Ok(vec![(0.0, 1.0)]);
    };
    let line = &block[idx];
    let context = |what: &str| format!("Species '{name}', ElectronicLevels: {what}");
    let unit = unit_tag(line).unwrap_or("1/cm");
    let n = parse_usize(line.split_whitespace().last().ok_or_else(|| context("malformed line"))?)?;
    if n == 0 {
        return Ok(vec![(0.0, 1.0)]);
    }
    let mut levels = Vec::with_capacity(n);
    for k in 0..n {
        let entry = block.get(idx + 1 + k).ok_or_else(|| context(&format!("{n} levels expected, {k} found")))?;
        let fields: Vec<&str> = entry.split_whitespace().collect();
        let (Some(e), Some(g)) = (fields.first().and_then(|v| parse_f64(v).ok()), fields.get(1).and_then(|v| parse_f64(v).ok())) else {
            return Err(context(&format!("{n} levels expected; line '{entry}' is not 'energy degeneracy'")));
        };
        if !(g > 0.0) {
            return Err(context(&format!("degeneracy {g} must be positive")));
        }
        levels.push((energy_to_cm1(e, unit)?, g));
    }
    levels.sort_by(|a, b| a.0.total_cmp(&b.0));
    if levels[0].0 != 0.0 {
        return Err(context(&format!("the lowest level lies at {} 1/cm; it must lie at 0, the zero of the species", levels[0].0)));
    }
    Ok(levels)
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
        electronic_levels: parse_electronic_levels(block, name)?,
        atom_mass_amu: Some(mass),
        rotational_constants_cm1: None,
        mass_amu: None,
        internal_rotors: Vec::new(),
    })
}

fn parse_rrho_species(block: &[String], name: &str) -> Result<MessSpeciesRrho, String> {
    parse_rrho_species_impl(block, name, true, false)
}

// `geometry_required = false` for phase-space-theory barriers: their rotational treatment comes
// from the two fragment geometries of the core, not from a molecular geometry.
// `frequencies_optional = true` for barriers given by an inverse Laplace transform: their k(E) comes from k_inf(T)
// and the reactant states, not from the barrier's own RRHO data.
fn parse_rrho_species_impl(
    block: &[String],
    name: &str,
    geometry_required: bool,
    frequencies_optional: bool,
) -> Result<MessSpeciesRrho, String> {
    let (block, rotor_blocks) = split_rotor_blocks(block, name)?;
    let block = &block[..];
    let rotational_constants_cm1 = parse_rotational_constants(block, name)?;
    let mass_amu = parse_rrho_mass(block, name)?;
    let (symbols, coords) = match (parse_geometry(block, "Geometry[angstrom]")?, &rotational_constants_cm1) {
        (Some(_), Some(_)) => {
            return Err(format!(
                "RRHO species '{name}': give either Geometry[angstrom] or RotationalConstants[1/cm], not both."
            ))
        }
        (Some(geometry), None) => geometry,
        (None, Some(_)) => {
            if mass_amu.is_none() {
                return Err(format!("RRHO species '{name}': RotationalConstants[1/cm] needs Mass[amu]."));
            }
            (Vec::new(), Vec::new())
        }
        (None, None) if !geometry_required => (Vec::new(), Vec::new()),
        (None, None) => {
            return Err(format!("RRHO species '{}' missing Geometry[angstrom] (or RotationalConstants[1/cm])", name))
        }
    };
    let symmetry_factor = parse_symmetry_factor(block)?;
    let vib = if frequencies_optional && !block.iter().any(|l| l.starts_with("Frequencies")) {
        Vec::new()
    } else {
        parse_frequencies(block)?
    };
    let zero_energy_cm1 = parse_zero_energy_cm1(block)?;
    let electronic_degeneracy_ground = parse_electronic_degeneracy_ground(block)?;
    let electronic_levels = parse_electronic_levels(block, name)?;
    let atoms = (!symbols.is_empty()).then_some(symbols.len());
    let internal_rotors = rotor_blocks
        .iter()
        .enumerate()
        .map(|(r, rotor)| parse_rotor(rotor, &format!("RRHO species '{name}', Rotor {}", r + 1), atoms))
        .collect::<Result<Vec<_>, _>>()?;

    Ok(MessSpeciesRrho {
        atom_mass_amu: None,
        name: name.to_string(),
        geometry_symbols: symbols,
        geometry_angstrom: coords,
        symmetry_factor,
        vibrational_frequencies_cm1: vib,
        zero_energy_cm1,
        electronic_degeneracy_ground,
        electronic_levels,
        rotational_constants_cm1,
        mass_amu,
        internal_rotors,
    })
}

/// The `Rotor ... End` blocks of an RRHO block, and the RRHO block without them.
fn split_rotor_blocks(block: &[String], name: &str) -> Result<(Vec<String>, Vec<Vec<String>>), String> {
    let (mut rest, mut rotors) = (Vec::new(), Vec::new());
    let mut lines = block.iter();
    while let Some(line) = lines.next() {
        if first_token(line) != Some("Rotor") {
            rest.push(line.clone());
            continue;
        }
        let mut rotor = vec![line.clone()];
        loop {
            let next = lines.next().ok_or_else(|| format!("RRHO species '{name}': '{line}' without End."))?;
            rotor.push(next.clone());
            if first_token(next) == Some("End") {
                break;
            }
        }
        rotors.push(rotor);
    }
    Ok((rest, rotors))
}

/// The data lines after `block[idx]`: the following lines that start with a number.
fn data_lines(block: &[String], idx: usize) -> &[String] {
    let count = block[idx + 1..].iter().take_while(|l| first_token(l).map_or(false, |t| parse_f64(t).is_ok())).count();
    &block[idx + 1..idx + 1 + count]
}

/// One `Rotor Hindered` or `Rotor Free` block (`MessInternalRotor`); `atoms` is the size of the species
/// geometry, None without a geometry.
fn parse_rotor(block: &[String], context: &str, atoms: Option<usize>) -> Result<MessInternalRotor, String> {
    let kind = match block[0].split_whitespace().nth(1) {
        Some("Hindered") => MessRotorKind::Hindered,
        Some("Free") => MessRotorKind::Free,
        other => {
            return Err(format!(
                "{context}: rotor type '{}' is not supported (Hindered or Free).",
                other.unwrap_or("")
            ))
        }
    };
    let index = |token: &str| -> Result<usize, String> {
        match parse_usize(token) {
            Ok(i) if i >= 1 => Ok(i - 1),
            _ => Err(format!("{context}: atom index '{token}' must be a positive integer (1-based).")),
        }
    };
    let (mut group, mut axis, mut symmetry) = (Vec::new(), None, 1_u32);
    let (mut points, mut fourier) = (None, None);
    let (mut size_min, mut size_max, mut rotational_constant_cm1) = (DEFAULT_BASIS_SIZE, MAX_BASIS_SIZE, None);
    let mut idx = 1;
    while idx < block.len() {
        let line = &block[idx];
        let tokens: Vec<&str> = line.split_whitespace().collect();
        let value = |k: usize| tokens.get(k).copied().ok_or_else(|| format!("{context}: malformed line '{line}'."));
        let mut data = 0;
        match tokens[0] {
            "End" => break,
            "Group" => {
                group = tokens[1..].iter().map(|t| index(t)).collect::<Result<Vec<_>, _>>()?;
                if (1..group.len()).any(|i| group[..i].contains(&group[i])) {
                    return Err(format!("{context}: an atom appears twice in the Group."));
                }
            }
            "Axis" => {
                let (a, b) = (index(value(1)?)?, index(value(2)?)?);
                if a == b {
                    return Err(format!("{context}: the two Axis atoms must differ."));
                }
                axis = Some((a, b));
            }
            "Symmetry" => {
                symmetry = match parse_usize(value(1)?) {
                    Ok(s) if s >= 1 => s as u32,
                    _ => return Err(format!("{context}: Symmetry must be a positive integer.")),
                };
            }
            "HamiltonSizeMin" => size_min = parse_usize(value(1)?)?,
            "HamiltonSizeMax" => size_max = parse_usize(value(1)?)?,
            "GridSize" | "ThermalPowerMax" => {}
            "RotationalConstant[1/cm]" => {
                let b = parse_f64(value(1)?)?;
                if !(b > 0.0) {
                    return Err(format!("{context}: RotationalConstant[1/cm] must be positive."));
                }
                rotational_constant_cm1 = Some(b);
            }
            key if key.starts_with("Potential[") || key.starts_with("FourierExpansion[") => {
                if kind == MessRotorKind::Free {
                    return Err(format!("{context}: a Free rotor has no potential ('{key}')."));
                }
                if points.is_some() || fourier.is_some() {
                    return Err(format!("{context}: the potential is given twice."));
                }
                let unit = unit_tag(line).unwrap_or("");
                let n = parse_usize(value(1)?)?;
                if n == 0 {
                    return Err(format!("{context}: {key} needs at least one value."));
                }
                let lines = data_lines(block, idx);
                let values = if key.starts_with("Potential[") {
                    // N values over one or more lines; the rest of the line with the N-th value is not read (as in
                    // the MESS reader)
                    let mut values = Vec::with_capacity(n);
                    for l in lines {
                        if values.len() == n {
                            break;
                        }
                        data += 1;
                        for token in l.split_whitespace().take(n - values.len()) {
                            values.push(parse_f64(token)?);
                        }
                    }
                    values
                } else {
                    // n lines "index value": the index is not used, the values are taken in order
                    data = lines.len().min(n);
                    lines[..data]
                        .iter()
                        .map(|l| match l.split_whitespace().collect::<Vec<_>>()[..] {
                            [_, v] => parse_f64(v),
                            _ => Err(format!("{context}: {key} line '{l}' needs 'index value'.")),
                        })
                        .collect::<Result<Vec<_>, _>>()?
                };
                if values.len() != n {
                    return Err(format!("{context}: {key} expects {n} values, got {}.", values.len()));
                }
                let values = values.iter().map(|v| energy_to_cm1(*v, unit)).collect::<Result<Vec<_>, _>>()?;
                if key.starts_with("Potential[") {
                    points = Some(values);
                } else {
                    fourier = Some(values);
                }
            }
            key => {
                return Err(format!(
                    "{context}: keyword '{key}' is not supported in a Rotor block (Group, Axis, Symmetry, \
                     Potential[unit], FourierExpansion[kcal/mol], HamiltonSizeMin, HamiltonSizeMax, \
                     RotationalConstant[1/cm], GridSize, ThermalPowerMax)."
                ))
            }
        }
        idx += 1 + data;
    }
    if group.is_empty() {
        return Err(format!("{context}: Group (the rotating atoms) is missing."));
    }
    let axis = axis.ok_or_else(|| format!("{context}: Axis (two atoms) is missing."))?;
    if group.contains(&axis.0) || group.contains(&axis.1) {
        return Err(format!("{context}: the Group must not contain an axis atom."));
    }
    if let Some(n) = atoms {
        if group.iter().chain([&axis.0, &axis.1]).any(|&a| a >= n) {
            return Err(format!("{context}: Group or Axis atom beyond the {n} atoms of the geometry."));
        }
    } else if rotational_constant_cm1.is_none() {
        return Err(format!(
            "{context}: the species has no geometry for the reduced moment; give RotationalConstant[1/cm] in the Rotor block."
        ));
    }
    if size_min % 2 == 0 || size_max % 2 == 0 || size_min > size_max {
        return Err(format!(
            "{context}: HamiltonSizeMin ({size_min}) and HamiltonSizeMax ({size_max}) must be odd, the minimum not above the maximum."
        ));
    }
    let potential = match (kind, points, fourier) {
        (MessRotorKind::Free, _, _) => TorsionalPotential { constant: 0.0, cosine: Vec::new(), sine: Vec::new() },
        (MessRotorKind::Hindered, Some(points), None) => {
            let n = points.len();
            if n > 1 && points[1] + points[n - 1] - 2.0 * points[0] < 0.0 {
                return Err(format!("{context}: the first point of the Potential must be a minimum (it is at angle 0)."));
            }
            potential_from_equidistant_points(&points)?
        }
        (MessRotorKind::Hindered, None, Some(coefficients)) => potential_from_fourier_expansion(&coefficients)?,
        _ => {
            return Err(format!(
                "{context}: a Hindered rotor needs its potential, Potential[unit] or FourierExpansion[kcal/mol]."
            ))
        }
    };
    Ok(MessInternalRotor {
        kind,
        group,
        axis,
        symmetry,
        potential,
        hamilton_size_min: size_min,
        hamilton_size_max: size_max,
        rotational_constant_cm1,
    })
}

/// `RotationalConstants[1/cm] N` followed by N values (3, or 1 for a linear rotor), in place of a geometry.
fn parse_rotational_constants(block: &[String], name: &str) -> Result<Option<Vec<f64>>, String> {
    let Some(idx) = block.iter().position(|l| first_token(l).map_or(false, |t| t.starts_with("RotationalConstants")))
    else {
        return Ok(None);
    };
    let line = &block[idx];
    if unit_tag(line) != Some("1/cm") {
        return Err(format!("RRHO species '{name}': RotationalConstants needs the unit [1/cm]: {line}"));
    }
    let n = parse_usize(line.split_whitespace().nth(1).ok_or_else(|| format!("Malformed line: {line}"))?)?;
    if n != 1 && n != 3 {
        return Err(format!("RRHO species '{name}': RotationalConstants[1/cm] needs 3 values, or 1 for a linear rotor (got {n})."));
    }
    let mut values = Vec::with_capacity(n);
    for l in &block[idx + 1..] {
        let tokens: Vec<f64> = l.split_whitespace().map_while(|t| parse_f64(t).ok()).collect();
        if tokens.is_empty() {
            break;
        }
        values.extend(tokens);
        if values.len() >= n {
            break;
        }
    }
    if values.len() != n || values.iter().any(|b| !(*b > 0.0)) {
        return Err(format!("RRHO species '{name}': RotationalConstants[1/cm] expects {n} positive values, got {values:?}."));
    }
    Ok(Some(values))
}

/// `Mass[amu] m` of an RRHO species given by its rotational constants.
fn parse_rrho_mass(block: &[String], name: &str) -> Result<Option<f64>, String> {
    let Some(line) = block.iter().find(|l| first_token(l).map_or(false, |t| t.starts_with("Mass["))) else {
        return Ok(None);
    };
    if unit_tag(line) != Some("amu") {
        return Err(format!("RRHO species '{name}': Mass needs the unit [amu]: {line}"));
    }
    let m = parse_f64(line.split_whitespace().nth(1).ok_or_else(|| format!("Malformed line: {line}"))?)?;
    if !(m > 0.0) {
        return Err(format!("RRHO species '{name}': Mass[amu] must be positive."));
    }
    Ok(Some(m))
}

/// Unit tag between square brackets of the first token, e.g. "kcal/mol" for "ZeroEnergy[kcal/mol]".
pub(crate) fn unit_tag(line: &str) -> Option<&str> {
    let first = first_token(line)?;
    first.split('[').nth(1).and_then(|t| t.split(']').next())
}

/// The MarXus header block, a MarXus extension of the deck header: the solution method and its settings
/// (`solution_method.rs`). Every keyword is optional; the command line overrides them.
///
///   MarXus
///     Method                              SteadyStateOlzmann (SteadyStateAbsorbingBarrier, CSE or TimeIntegration; required)
///     AbsorbingBarrierBelowThreshold[kT]  10                 (intermediate steady state)
///     EigenSolver                         InverseIteration   (FullDecomposition or Lapack)
///     SumRuleTolerance                    1.5e-2             (thermal eigenpair of the final steady state)
///     Integrator                          Rodas4             (time integration: Rodas4, Rodas3, Ros4, Ros3, Ros2)
///     InitialState                        Pulse              (time integration: Pulse or Continuous)
///     TimeRange[s]                        1e-12  1e2         (time integration: first and last output time)
///     TimesPerDecade                      4                  (time integration: output times per decade)
///     IntegrationTolerance                1e-6               (time integration: relative tolerance)
///     CompareWithCse                      1e-2               (time integration of a Preparation in a constant bath:
///                                                             compare with the CSE description in time; tolerance)
///     NCores                              8                  (cores of the run: the (T, p) conditions are
///                                                             computed in batches of up to NCores at a time)
///     CollisionIntegral                   Neufeld            (Lennard-Jones Omega(2,2)*: Neufeld or Troe)
///     RotorReducedMoment                  Pitzer             (internal rotors from a geometry: Pitzer or BondAxis)
///     ChemicalSubspaceCriterion           RelaxationProjection  (CSE merging: RelaxationProjection, as MESS's direct
///                                                             method, or EigenvalueRatio)
///   End
fn parse_marxus_header(block: &[String]) -> Result<SolutionSettings, String> {
    let context = |what: &str| format!("MarXus header block: {what}");
    let mut settings = SolutionSettings::default();
    let mut seen: Vec<&str> = Vec::new();
    for line in block.iter().skip(1) {
        let key = first_token(line).unwrap_or("");
        if key == "End" {
            break;
        }
        if seen.contains(&key) {
            return Err(context(&format!("{key} given twice")));
        }
        seen.push(key);
        let value = line.split_whitespace().nth(1).ok_or_else(|| context(&format!("no value in '{line}'")))?;
        match key {
            "Method" => settings.method = Some(SolutionMethod::from_keyword(value).map_err(|e| context(&e))?),
            "SteadyState" => return Err(context(&format!("the keyword SteadyState no longer exists: {STEADY_STATE_KEYWORD_REPLACED}."))),
            "AbsorbingBarrierBelowThreshold[kT]" => settings.absorbing_barrier_kt = Some(parse_f64(value)?),
            "EigenSolver" => settings.eigen_solver = Some(eigen_solver_from_keyword(value).map_err(|e| context(&e))?),
            "SumRuleTolerance" => settings.sum_rule_tolerance = Some(parse_f64(value)?),
            "Integrator" => settings.integrator = Some(integrator_from_keyword(value).map_err(|e| context(&e))?),
            "InitialState" => settings.initial_state = Some(initial_state_from_keyword(value).map_err(|e| context(&e))?),
            "TimeRange[s]" => {
                let last = line
                    .split_whitespace()
                    .nth(2)
                    .ok_or_else(|| context("TimeRange[s] needs two values, the first and the last output time"))?;
                settings.time_range_s = Some((parse_f64(value)?, parse_f64(last)?));
            }
            "TimesPerDecade" => settings.times_per_decade = Some(parse_usize(value)?),
            "IntegrationTolerance" => settings.integration_tolerance = Some(parse_f64(value)?),
            "CompareWithCse" => settings.cse_comparison_tolerance = Some(parse_f64(value)?),
            "NCores" => {
                let cores = parse_usize(value).map_err(|e| context(&format!("NCores: {e}")))?;
                if cores == 0 {
                    return Err(context("NCores: the number of cores must be at least 1"));
                }
                settings.cores = Some(cores);
            }
            "CollisionIntegral" => settings.collision_integral = Some(collision_integral_from_keyword(value).map_err(|e| context(&e))?),
            "RotorReducedMoment" => {
                settings.rotor_reduced_moment = Some(rotor_reduced_moment_from_keyword(value).map_err(|e| context(&e))?)
            }
            "ChemicalSubspaceCriterion" => {
                settings.chemical_subspace_criterion =
                    Some(chemical_subspace_criterion_from_keyword(value).map_err(|e| context(&e))?)
            }
            _ => {
                return Err(context(&format!(
                    "unknown keyword '{key}' (Method, AbsorbingBarrierBelowThreshold[kT], EigenSolver, \
                     SumRuleTolerance, Integrator, InitialState, TimeRange[s], TimesPerDecade, IntegrationTolerance, \
                     CompareWithCse, NCores, CollisionIntegral, RotorReducedMoment, ChemicalSubspaceCriterion)"
                )))
            }
        }
    }
    Ok(settings)
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

/// The MarXus `FragmentEnergy <kind> ... End` block inside a barrier, if present.
fn parse_fragment_energy(block: &[String], barrier: &str) -> Result<Option<FragmentEnergySpecification>, String> {
    let Some(start) = block.iter().position(|l| first_token(l) == Some("FragmentEnergy")) else {
        return Ok(None);
    };
    let context = |what: &str| format!("Barrier '{barrier}', FragmentEnergy: {what}");
    let kind = block[start].split_whitespace().nth(1).unwrap_or("").to_string();
    let mut keys: Vec<(String, Vec<f64>, Option<String>)> = Vec::new();
    for line in &block[start + 1..] {
        let key = first_token(line).unwrap_or("");
        if key == "End" {
            break;
        }
        let values = line.split_whitespace().skip(1).map(parse_f64).collect::<Result<Vec<f64>, String>>()?;
        keys.push((key.split('[').next().unwrap_or("").to_string(), values, unit_tag(line).map(|u| u.to_string())));
    }
    let find = |name: &str| keys.iter().find(|(k, _, _)| k == name).ok_or_else(|| context(&format!("missing {name}")));
    let allowed = |names: &[&str]| -> Result<(), String> {
        match keys.iter().find(|(k, _, _)| !names.contains(&k.as_str())) {
            Some((k, _, _)) => Err(context(&format!("unknown keyword '{k}' for {kind} ({})", names.join(", ")))),
            None => Ok(()),
        }
    };
    let scalar = |name: &str| -> Result<f64, String> {
        let (_, v, _) = find(name)?;
        v.first().copied().ok_or_else(|| context(&format!("{name} needs a value")))
    };
    // Intercept in its unit (cm-1 after conversion), gradient without unit.
    let line = |name: &str| -> Result<[f64; 2], String> {
        let (_, v, unit) = find(name)?;
        if v.len() != 2 {
            return Err(context(&format!("{name} needs an intercept and a gradient")));
        }
        let unit = unit.as_deref().ok_or_else(|| context(&format!("{name} needs an energy unit, e.g. {name}[1/cm]")))?;
        Ok([energy_to_cm1(v[0], unit)?, v[1]])
    };
    let spec = match kind.as_str() {
        "Prior" => {
            allowed(&[])?;
            FragmentEnergySpecification::Prior
        }
        "ModifiedPrior" => {
            allowed(&["Order", "TemperatureExponent", "ReferenceTemperature"])?;
            FragmentEnergySpecification::ModifiedPrior {
                order: scalar("Order")?,
                temperature_exponent: scalar("TemperatureExponent")?,
                reference_temperature_kelvin: scalar("ReferenceTemperature")?,
            }
        }
        "TwoPieceGaussian" => {
            allowed(&["Mu", "SigmaLeft", "SigmaRight"])?;
            FragmentEnergySpecification::TwoPieceGaussian { mu: line("Mu")?, sigma_left: line("SigmaLeft")?, sigma_right: line("SigmaRight")? }
        }
        other => return Err(context(&format!("unknown kind '{other}' (Prior, ModifiedPrior, TwoPieceGaussian)"))),
    };
    Ok(Some(spec))
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
    let mut tst_level = PstTstLevel::default();

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
        } else if first_token(line) == Some("TSTLevel") {
            tst_level = match line.split_whitespace().nth(1) {
                Some("T") => PstTstLevel::T,
                Some("E") => PstTstLevel::E,
                Some("EJ") => PstTstLevel::EJ,
                Some("J=0") => PstTstLevel::J0,
                other => return Err(format!("PhaseSpaceTheory TSTLevel: unknown level {other:?} (T, E, EJ, J=0).")),
            };
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
        tst_level,
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
        chemical_eigenvalue_max: None,
        well_projection_threshold: None,
        alpha_factor_cm1: None,
        alpha_power: None,
        lj_epsilons_cm1: None,
        lj_sigmas_angstrom: None,
        lj_masses_amu: None,
        reactant_name: None,
        excess_reactant_concentration_cm3: None,
        solution: SolutionSettings::default(),
    };

    let mut wells: HashMap<String, MessSpeciesRrho> = HashMap::new();
    let mut bimolecular: HashMap<String, MessBimolecular> = HashMap::new();
    let mut barriers: Vec<MessBarrier> = Vec::new();
    let mut well_escape_rate_s_inv: HashMap<String, f64> = HashMap::new();
    let mut well_collision: HashMap<String, WellCollisionOverride> = HashMap::new();
    let mut well_order: Vec<String> = Vec::new();
    let mut dummy_bimolecular: Vec<String> = Vec::new();
    let mut marxus_header_read = false;
    let mut preparation_block: Option<Vec<String>> = None;

    let mut i = 0usize;
    while i < lines.len() {
        let line = strip_comment(&lines[i]);
        if line.is_empty() {
            i += 1;
            continue;
        }

        // MarXus header block: solution method and its settings.
        if first_token(line) == Some("MarXus") {
            if marxus_header_read {
                return Err("Input deck: the MarXus header block is given twice.".into());
            }
            let (block, next) = collect_block(&lines, i);
            global.solution = parse_marxus_header(&block)?;
            marxus_header_read = true;
            i = next;
            continue;
        }

        // Preparation block: initial population, sources and bath history (`preparation_input.rs`).
        if first_token(line) == Some("Preparation") {
            if preparation_block.is_some() {
                return Err("Input deck: the Preparation block is given twice.".into());
            }
            let (block, next) = super::preparation_input::collect_preparation_block(&lines, i)?;
            preparation_block = Some(block);
            i = next;
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
                "ChemicalEigenvalueMax" => {
                    let x = parse_f64(v)?;
                    if !(x > 0.0 && x < 1.0) {
                        return Err(format!(
                            "ChemicalEigenvalueMax {x}: MarXus accepts 0 < value < 1, read as the relaxational projection \
                             1 - F_ne <= value (MESS direct method; the default), or with `ChemicalSubspaceCriterion \
                             EigenvalueRatio` in the MarXus block as Lambda <= value x Lambda_(N+1). MESS reads a value > 1 \
                             as Lambda_(N+1)/Lambda >= value: give EigenvalueRatio with 1/value. Values below 0 are not \
                             implemented."
                        ));
                    }
                    global.chemical_eigenvalue_max = Some(x);
                }
                "WellProjectionThreshold" => global.well_projection_threshold = Some(parse_f64(v)?),
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
            let (block, collision) = split_well_collision_block(&block, &name)?;
            if let Some(collision) = collision {
                well_collision.insert(name.clone(), collision);
            }
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
            // `Dummy` (as in MESS): a product without molecular data; the block has no End.
            let following = (i + 1..lines.len()).find(|&k| !strip_comment(&lines[k]).is_empty());
            if let Some(k) = following.filter(|&k| first_token(strip_comment(&lines[k])) == Some("Dummy")) {
                dummy_bimolecular.push(name);
                i = k + 1;
                continue;
            }
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
            let mut fragment_well = None;
            let mut partner_concentration_cm3 = None;
            let mut lumped_state = false;
            let mut excess_fragment = None;
            for l in &block {
                let key = first_token(l).unwrap_or("");
                let value = || l.split_whitespace().nth(1).ok_or_else(|| format!("Bimolecular '{name}': no value in '{l}'"));
                if key == "LumpedState" {
                    lumped_state = true;
                } else if key == "ExcessFragment" {
                    excess_fragment = Some(value()?.to_string());
                } else if key == "FragmentWell" {
                    fragment_well = Some(value()?.to_string());
                } else if key.starts_with("PartnerConcentration") {
                    if unit_tag(l) != Some("molecule/cm^3") {
                        return Err(format!("Bimolecular '{name}': PartnerConcentration needs the unit [molecule/cm^3]."));
                    }
                    partner_concentration_cm3 = Some(parse_f64(value()?)?);
                }
            }
            if lumped_state {
                if fragment_well.is_some() {
                    return Err(format!("Bimolecular '{name}': LumpedState and FragmentWell exclude each other."));
                }
                if excess_fragment.is_none() || partner_concentration_cm3.is_none() {
                    return Err(format!(
                        "Bimolecular '{name}': LumpedState needs ExcessFragment <name> and PartnerConcentration[molecule/cm^3]."
                    ));
                }
            } else if excess_fragment.is_some() {
                return Err(format!("Bimolecular '{name}': ExcessFragment belongs to a LumpedState."));
            } else if fragment_well.is_some() != partner_concentration_cm3.is_some() {
                return Err(format!("Bimolecular '{name}': FragmentWell and PartnerConcentration[molecule/cm^3] go together."));
            }

            bimolecular.insert(
                name.clone(),
                MessBimolecular {
                    name,
                    fragment_a: frag_a,
                    fragment_b: frag_b,
                    ground_energy_cm1,
                    fragment_well,
                    partner_concentration_cm3,
                    lumped_state,
                    excess_fragment,
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
            // A barrier given by an inverse Laplace transform takes its k(E) from k_inf(T) and the fragments.
            let has_ilt = block.iter().any(|l| first_token(l) == Some("InverseLaplaceTransform"));
            let geometry_required = matches!(core, MessBarrierCore::TightRrho) && !has_ilt;
            let rrho = parse_rrho_species_impl(&block, &name, geometry_required, has_ilt)?;
            let inverse_laplace_transform = parse_inverse_laplace_transform(&block, &name)?;
            let fragment_energy = parse_fragment_energy(&block, &name)?;
            let tunneling = parse_tunneling(&block, &name)?;

            barriers.push(MessBarrier {
                name,
                left,
                right,
                rrho,
                core,
                inverse_laplace_transform,
                fragment_energy,
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
        well_collision,
        well_order,
        dummy_bimolecular,
        preparation_block,
    })
}

/// Convenience wrapper for parsing a MESS input from a file on disk.
/// `Species NAME` lines at the top level of a deck, each followed by an RRHO or Atom block with the syntax of the
/// Fragment blocks (MarXus extension, used by the photoionization decks of `photoion::deck`): the molecular data of
/// species outside a network. ZeroEnergy is optional (0 when absent); names must be unique. The `Species` line inside a
/// MESS Well block has no name and is not read here.
pub fn parse_species_blocks(input: &str) -> Result<Vec<MessSpeciesRrho>, String> {
    let lines: Vec<String> = input.lines().map(|s| s.to_string()).collect();
    let mut species: Vec<MessSpeciesRrho> = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        let tokens: Vec<&str> = strip_comment(&lines[i]).split_whitespace().collect();
        let ["Species", name] = tokens[..] else {
            i += 1;
            continue;
        };
        let mut j = i + 1;
        while j < lines.len() && strip_comment(&lines[j]).is_empty() {
            j += 1;
        }
        let opener = lines.get(j).and_then(|l| first_token(strip_comment(l))).unwrap_or("");
        if opener != "RRHO" && opener != "Atom" {
            return Err(format!("Species '{name}': an RRHO or Atom block must follow the Species line."));
        }
        if species.iter().any(|s| s.name == name) {
            return Err(format!("Species '{name}' is defined twice."));
        }
        let (mut block, next) = collect_block(&lines, j);
        let parsed = if opener == "Atom" {
            parse_fragment(&block, name)?
        } else {
            if !block.iter().any(|l| l.starts_with("ZeroEnergy")) {
                block.push("ZeroEnergy[1/cm] 0".to_string());
            }
            parse_rrho_species(&block, name)?
        };
        species.push(parsed);
        i = next;
    }
    Ok(species)
}

pub fn parse_mess_input_file(path: impl AsRef<Path>) -> Result<MessDeck, String> {
    let text = std::fs::read_to_string(path.as_ref())
        .map_err(|e| format!("Failed to read MESS input file: {e}"))?;
    parse_mess_input(&text)
}

#[cfg(test)]
mod tests {
    use super::{parse_mess_input, parse_species_blocks, MessBarrierCore, MessRotorKind, PstTstLevel, SolutionSettings, CM1_TO_KCAL};

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
        // Without TSTLevel the MESS default EJ applies.
        assert!(matches!(
            parsed.barriers[0].core,
            MessBarrierCore::PhaseSpaceTheory { tst_level: PstTstLevel::EJ, .. }
        ));
        // TSTLevel is read: E, EJ, T, J=0.
        for (keyword, level) in [("E", PstTstLevel::E), ("EJ", PstTstLevel::EJ), ("T", PstTstLevel::T), ("J=0", PstTstLevel::J0)] {
            let with_level = deck.replace("PotentialPowerExponent 6.", &format!("PotentialPowerExponent 6.\n      TSTLevel {keyword}"));
            let parsed = parse_mess_input(&with_level).expect("should parse");
            match &parsed.barriers[0].core {
                MessBarrierCore::PhaseSpaceTheory { tst_level, .. } => assert_eq!(*tst_level, level, "{keyword}"),
                _ => panic!("not a phase-space core"),
            }
        }
        let unknown = deck.replace("PotentialPowerExponent 6.", "PotentialPowerExponent 6.\n      TSTLevel X");
        assert!(parse_mess_input(&unknown).is_err());
    }

    const MARXUS_HEADER_DECK: &str = "TemperatureList[K] 300.\nPressureList[torr] 760\n\
MarXus\n  Method SteadyStateOlzmann   ! the thermal eigenpair is part of it\n  EigenSolver Lapack\n  \
SumRuleTolerance 2e-2\n  AbsorbingBarrierBelowThreshold[kT] 5\nEnd\nModel\nEnd\n";

    #[test]
    fn the_marxus_header_block_gives_the_solution_settings() {
        use super::super::chemical_activation_eigen::EigenSolver;
        use super::super::solution_method::SolutionMethod;
        let parsed = parse_mess_input(MARXUS_HEADER_DECK).expect("should parse");
        assert_eq!(
            parsed.global.solution,
            SolutionSettings {
                method: Some(SolutionMethod::SteadyStateOlzmann),
                absorbing_barrier_kt: Some(5.0),
                eigen_solver: Some(EigenSolver::FullDecompositionLapack),
                sum_rule_tolerance: Some(0.02),
                ..Default::default()
            }
        );
        // The header keywords around the block are still read.
        assert_eq!(parsed.global.pressures_torr, vec![760.0]);
        let cse = parse_mess_input(
            &MARXUS_HEADER_DECK.replace("Method SteadyStateOlzmann", "Method CSE"),
        )
        .unwrap();
        assert_eq!(
            cse.global.solution.method,
            Some(SolutionMethod::ChemicallySignificantEigenvalues)
        );
    }

    #[test]
    fn the_cse_merging_thresholds_are_read_from_the_header() {
        // MESS keywords: ChemicalEigenvalueMax (chemical eigenvalues <= this x the lowest relaxation eigenvalue)
        // and WellProjectionThreshold (primary wells of the partition), both used by the CSE species merging.
        let deck = MARXUS_HEADER_DECK.replace(
            "PressureList[torr] 760\n",
            "PressureList[torr] 760\nChemicalEigenvalueMax 0.2\nWellProjectionThreshold 0.3\n",
        );
        let g = parse_mess_input(&deck).unwrap().global;
        assert_eq!(g.chemical_eigenvalue_max, Some(0.2));
        assert_eq!(g.well_projection_threshold, Some(0.3));
        let absent = parse_mess_input(MARXUS_HEADER_DECK).unwrap().global;
        assert_eq!(absent.chemical_eigenvalue_max, None);
        // Only MESS's absolute threshold, 0 < value < 1, is implemented; the other MESS modes are refused.
        for bad in ["ChemicalEigenvalueMax 2", "ChemicalEigenvalueMax -0.5", "ChemicalEigenvalueMax 0"] {
            let deck = MARXUS_HEADER_DECK.replace("PressureList[torr] 760\n", &format!("PressureList[torr] 760\n{bad}\n"));
            let err = parse_mess_input(&deck).unwrap_err();
            assert!(err.contains("ChemicalEigenvalueMax"), "{bad}: {err}");
        }
    }

    #[test]
    fn the_marxus_header_block_gives_the_number_of_cores() {
        let deck = MARXUS_HEADER_DECK.replace("EigenSolver Lapack", "EigenSolver Lapack\n  NCores 8");
        assert_eq!(parse_mess_input(&deck).unwrap().global.solution.cores, Some(8));
        assert_eq!(parse_mess_input(MARXUS_HEADER_DECK).unwrap().global.solution.cores, None);
        for bad in ["NCores 0", "NCores 2.5", "NCores -1"] {
            let deck = MARXUS_HEADER_DECK.replace("EigenSolver Lapack", &format!("EigenSolver Lapack\n  {bad}"));
            let err = parse_mess_input(&deck).unwrap_err();
            assert!(err.contains("NCores"), "{bad}: {err}");
        }
    }

    #[test]
    fn the_marxus_header_block_gives_the_chemical_subspace_criterion() {
        use crate::masterequation::chemically_significant_eigenvalues::ChemicalSubspaceCriterion;
        let deck = MARXUS_HEADER_DECK.replace("EigenSolver Lapack", "EigenSolver Lapack\n  ChemicalSubspaceCriterion EigenvalueRatio");
        assert_eq!(parse_mess_input(&deck).unwrap().global.solution.chemical_subspace_criterion, Some(ChemicalSubspaceCriterion::EigenvalueRatio));
        assert_eq!(parse_mess_input(MARXUS_HEADER_DECK).unwrap().global.solution.chemical_subspace_criterion, None);
        let bad = MARXUS_HEADER_DECK.replace("EigenSolver Lapack", "EigenSolver Lapack\n  ChemicalSubspaceCriterion Gap");
        assert!(parse_mess_input(&bad).unwrap_err().contains("ChemicalSubspaceCriterion"));
    }

    #[test]
    fn the_marxus_header_block_gives_the_rotor_reduced_moment() {
        use crate::rrkm::internal_rotor::ReducedMomentModel;
        let deck = MARXUS_HEADER_DECK.replace("EigenSolver Lapack", "EigenSolver Lapack\n  RotorReducedMoment BondAxis");
        assert_eq!(parse_mess_input(&deck).unwrap().global.solution.rotor_reduced_moment, Some(ReducedMomentModel::BondAxis));
        assert_eq!(parse_mess_input(MARXUS_HEADER_DECK).unwrap().global.solution.rotor_reduced_moment, None);
        let bad = MARXUS_HEADER_DECK.replace("EigenSolver Lapack", "EigenSolver Lapack\n  RotorReducedMoment Smith");
        assert!(parse_mess_input(&bad).unwrap_err().contains("RotorReducedMoment"));
    }

    #[test]
    fn the_marxus_header_block_gives_the_collision_integral() {
        use super::super::collisional_relaxation::CollisionIntegral;
        let deck = MARXUS_HEADER_DECK.replace("EigenSolver Lapack", "EigenSolver Lapack\n  CollisionIntegral Neufeld");
        assert_eq!(parse_mess_input(&deck).unwrap().global.solution.collision_integral, Some(CollisionIntegral::Neufeld1972));
        assert_eq!(parse_mess_input(MARXUS_HEADER_DECK).unwrap().global.solution.collision_integral, None);
        let bad = MARXUS_HEADER_DECK.replace("EigenSolver Lapack", "EigenSolver Lapack\n  CollisionIntegral Smith");
        assert!(parse_mess_input(&bad).unwrap_err().contains("CollisionIntegral"));
    }

    #[test]
    fn the_marxus_header_block_gives_the_time_integration_settings() {
        use super::super::direct_time_integration::InitialState;
        use super::super::solution_method::SolutionMethod;
        use crate::numeric::integrators::rosenbrock_methods::RosenbrockMethod;
        let deck = "TemperatureList[K] 300.\nPressureList[torr] 760\nMarXus\n  Method TimeIntegration\n  Integrator Ros4\n  \
InitialState Continuous\n  TimeRange[s] 1e-10 1e1\n  TimesPerDecade 3\n  IntegrationTolerance 1e-7\nEnd\nModel\nEnd\n";
        let s = parse_mess_input(deck).unwrap().global.solution;
        assert_eq!(s.method, Some(SolutionMethod::TimeIntegration));
        assert_eq!(s.integrator, Some(RosenbrockMethod::Ros4));
        assert_eq!(s.initial_state, Some(InitialState::ContinuousFormation));
        assert_eq!(s.time_range_s, Some((1e-10, 1e1)));
        assert_eq!(s.times_per_decade, Some(3));
        assert_eq!(s.integration_tolerance, Some(1e-7));
        let one_value = deck.replace("TimeRange[s] 1e-10 1e1", "TimeRange[s] 1e-10");
        assert!(parse_mess_input(&one_value)
            .unwrap_err()
            .contains("TimeRange"));
        assert_eq!(s.cse_comparison_tolerance, None);
        let compared = deck.replace("TimesPerDecade 3", "TimesPerDecade 3\n  CompareWithCse 1e-3");
        assert_eq!(parse_mess_input(&compared).unwrap().global.solution.cse_comparison_tolerance, Some(1e-3));
    }

    #[test]
    fn a_deck_without_a_marxus_header_block_leaves_the_solution_settings_unset() {
        let parsed = parse_mess_input("TemperatureList[K] 300.\nPressureList[torr] 760\nModel\nEnd\n").unwrap();
        assert_eq!(parsed.global.solution, SolutionSettings::default());
    }

    #[test]
    fn the_marxus_header_block_refuses_unknown_keywords_values_and_repetitions() {
        let unknown = MARXUS_HEADER_DECK.replace("EigenSolver Lapack", "Solver Lapack");
        assert!(parse_mess_input(&unknown).unwrap_err().contains("MarXus"));
        // The removed keyword SteadyState names the two steady-state methods that replace it.
        let old_keyword = MARXUS_HEADER_DECK.replace(
            "EigenSolver Lapack",
            "EigenSolver Lapack\n  SteadyState Final",
        );
        assert!(parse_mess_input(&old_keyword)
            .unwrap_err()
            .contains("SteadyStateAbsorbingBarrier"));
        let eigenvalue =
            MARXUS_HEADER_DECK.replace("Method SteadyStateOlzmann", "Method Eigenvalue");
        assert!(parse_mess_input(&eigenvalue)
            .unwrap_err()
            .contains("final steady state"));
        let repeated = MARXUS_HEADER_DECK.replace(
            "EigenSolver Lapack",
            "EigenSolver Lapack\n  EigenSolver Full",
        );
        assert!(parse_mess_input(&repeated)
            .unwrap_err()
            .contains("EigenSolver"));
        let two_blocks =
            MARXUS_HEADER_DECK.replace("Model\n", "MarXus\n  Method CSE\nEnd\nModel\n");
        assert!(parse_mess_input(&two_blocks)
            .unwrap_err()
            .contains("MarXus"));
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
    fn rotational_constants_and_mass_replace_the_geometry_of_an_rrho_species() {
        let deck = "TemperatureList[K] 300.\nPressureList[torr] 760\nModel\n  Well W1\n    Species\n      RRHO\n        \
RotationalConstants[1/cm] 3\n          0.1081 0.161\n          0.3105\n        Mass[amu] 75\n        Core RigidRotor\n          \
SymmetryFactor 1\n        End\n        Frequencies[1/cm] 1\n          500\n        ZeroEnergy[kcal/mol] -30\n      End\n  End\nEnd\n";
        let w = &parse_mess_input(deck).unwrap().wells["W1"];
        assert_eq!(w.rotational_constants_cm1, Some(vec![0.1081, 0.161, 0.3105]));
        assert_eq!(w.mass_amu, Some(75.0));
        assert!(w.geometry_symbols.is_empty());
        // One constant: a linear rotor.
        let linear = deck.replace("RotationalConstants[1/cm] 3\n          0.1081 0.161\n          0.3105", "RotationalConstants[1/cm] 1\n 1.449");
        assert_eq!(parse_mess_input(&linear).unwrap().wells["W1"].rotational_constants_cm1, Some(vec![1.449]));
        // Two constants, missing values, no mass, or both a geometry and constants: errors.
        for bad in [
            deck.replace("RotationalConstants[1/cm] 3\n          0.1081 0.161\n          0.3105", "RotationalConstants[1/cm] 2\n 1 2"),
            deck.replace("          0.3105\n", ""),
            deck.replace("        Mass[amu] 75\n", ""),
            deck.replace("        Mass[amu] 75\n", "        Mass[amu] 75\n        Geometry[angstrom] 1\n        H 0 0 0\n"),
        ] {
            assert!(parse_mess_input(&bad).is_err(), "{bad}");
        }
    }

    /// Deck with one well W1 whose RRHO block holds `species` (geometry or rotational constants, rotors).
    fn rotor_deck(species: &str) -> String {
        format!(
            "TemperatureList[K] 300.\nPressureList[torr] 760\nModel\n  Well W1\n    Species\n      RRHO\n{species}\n        \
             Frequencies[1/cm] 2\n          500 1200\n        ZeroEnergy[kcal/mol] -30\n      End\n  End\nEnd\n"
        )
    }

    const ETHANE_CORE: &str = "        Geometry[angstrom] 8\n          C 0 0 -0.765\n          C 0 0 0.765\n          \
H 1.02 0 -1.16\n          H -0.51 0.883 -1.16\n          H -0.51 -0.883 -1.16\n          H -1.02 0 1.16\n          \
H 0.51 -0.883 1.16\n          H 0.51 0.883 1.16\n        Core RigidRotor\n          SymmetryFactor 6\n        End";

    #[test]
    fn species_blocks_give_the_molecular_data_of_a_deck_without_a_network() {
        // MarXus extension (photoionization decks): `Species NAME` followed by an RRHO or Atom block, the syntax of
        // the Fragment blocks; ZeroEnergy is optional (0 when absent).
        let deck = "Photoionization\n  Temperature[K] 298\nEnd\nSpecies EtBr   ! neutral\n  RRHO\n    Geometry[angstrom] 3\n      \
C 0 0 0\n      C 0 0 1.5\n      Br 1.9 0 0\n    Core RigidRotor\n      SymmetryFactor 1\n    End\n    Frequencies[1/cm] 2\n      \
290 960\n    ElectronicLevels[1/cm] 1\n      0 1\n  End\nSpecies Br\n  Atom\n    Mass[amu] 78.918\n    \
ElectronicLevels[1/cm] 1\n      0 4\n  End\nSpecies TS\n  RRHO\n    RotationalConstants[1/cm] 3\n      0.9 0.1 0.09\n    \
Mass[amu] 108.0\n    Core RigidRotor\n      SymmetryFactor 1\n    End\n    Frequencies[1/cm] 1\n      700\n    \
ZeroEnergy[kcal/mol] 30\n  End\n";
        let species = parse_species_blocks(deck).unwrap();
        let names: Vec<&str> = species.iter().map(|s| s.name.as_str()).collect();
        assert_eq!(names, ["EtBr", "Br", "TS"]);
        assert_eq!(species[0].vibrational_frequencies_cm1, vec![290.0, 960.0]);
        assert_eq!((species[0].zero_energy_cm1, species[0].geometry_symbols.len()), (0.0, 3));
        assert_eq!((species[1].atom_mass_amu, species[1].electronic_degeneracy_ground), (Some(78.918), 4.0));
        assert!((species[2].zero_energy_cm1 - 30.0 / CM1_TO_KCAL).abs() < 1e-9);
        let twice = format!("{deck}Species Br\n  Atom\n    Mass[amu] 79\n  End\n");
        assert!(parse_species_blocks(&twice).unwrap_err().contains("Br"));
        assert!(parse_species_blocks("Species X\nPhotoionization\nEnd\n").unwrap_err().contains("X"));
    }

    #[test]
    fn every_electronic_level_is_read() {
        let deck = |levels: &str| rotor_deck(&format!("{ETHANE_CORE}\n        ElectronicLevels[1/cm] {levels}"));
        let w1 = &parse_mess_input(&deck("3\n          0 2\n          139.7 2\n          1000 4")).unwrap().wells["W1"];
        assert_eq!(w1.electronic_levels, vec![(0.0, 2.0), (139.7, 2.0), (1000.0, 4.0)]);
        assert_eq!(w1.electronic_degeneracy_ground, 2.0);
        // the lowest level is the zero of the species; a missing level line is an error
        assert!(parse_mess_input(&deck("2\n          50 2\n          139.7 2")).unwrap_err().contains("ElectronicLevels"));
        assert!(parse_mess_input(&deck("2\n          0 2")).unwrap_err().contains("ElectronicLevels"));
        // an Atom fragment, O(3P): 3P2, 3P1, 3P0
        let o = parse_species_blocks("Species O\n  Atom\n    Mass[amu] 15.995\n    ElectronicLevels[1/cm] 3\n      0 5\n      158.3 3\n      227 1\n  End\n").unwrap();
        assert_eq!(o[0].electronic_levels, vec![(0.0, 5.0), (158.3, 3.0), (227.0, 1.0)]);
        // without the block: one level with degeneracy 1
        let plain = &parse_mess_input(&rotor_deck(ETHANE_CORE)).unwrap().wells["W1"];
        assert_eq!(plain.electronic_levels, vec![(0.0, 1.0)]);
    }

    #[test]
    fn hindered_and_free_rotor_blocks_are_read() {
        use crate::rrkm::internal_rotor::{potential_from_equidistant_points, TorsionalPotential};
        let deck = rotor_deck(&format!(
            "{ETHANE_CORE}\n        Rotor Hindered   ! CH3\n          Group 3 4 5\n          Axis 1 2\n          Symmetry 3\n          \
             Potential[kcal/mol] 2\n          0. 2.45\n        End\n        Rotor Free\n          HamiltonSizeMin 11\n          \
             HamiltonSizeMax 21\n          Group 6\n          Axis 2 1\n        End"
        ));
        let w1 = &parse_mess_input(&deck).unwrap().wells["W1"];
        assert_eq!(w1.vibrational_frequencies_cm1, vec![500.0, 1200.0]);
        assert_eq!(w1.symmetry_factor, 6.0);
        assert_eq!(w1.geometry_symbols.len(), 8);
        let [hindered, free] = &w1.internal_rotors[..] else { panic!("{:?}", w1.internal_rotors) };
        assert_eq!(hindered.kind, MessRotorKind::Hindered);
        assert_eq!((hindered.group.clone(), hindered.axis, hindered.symmetry), (vec![2, 3, 4], (0, 1), 3));
        assert_eq!(hindered.potential, potential_from_equidistant_points(&[0.0, 2.45 / CM1_TO_KCAL]).unwrap());
        assert_eq!((hindered.hamilton_size_min, hindered.hamilton_size_max), (999, 1999));
        assert_eq!(hindered.rotational_constant_cm1, None);
        assert_eq!(free.kind, MessRotorKind::Free);
        assert_eq!((free.group.clone(), free.axis, free.symmetry), (vec![5], (1, 0), 1));
        assert_eq!(free.potential, TorsionalPotential { constant: 0.0, cosine: Vec::new(), sine: Vec::new() });
        assert_eq!((free.hamilton_size_min, free.hamilton_size_max), (11, 21));
    }

    #[test]
    fn potential_values_are_read_over_lines_and_the_rest_of_the_last_line_is_ignored() {
        // As in the MESS reader: N values over one or more lines; the rest of the line with the N-th value is not read.
        use crate::rrkm::internal_rotor::potential_from_equidistant_points;
        let deck = |values: &str| {
            rotor_deck(&format!(
                "{ETHANE_CORE}\n        Rotor Hindered\n          Group 3 4 5\n          Axis 1 2\n          Symmetry 3\n          \
                 Potential[kcal/mol] 4\n{values}\n        End"
            ))
        };
        let expected = potential_from_equidistant_points(&[0.0, 1.0 / CM1_TO_KCAL, 2.0 / CM1_TO_KCAL, 1.0 / CM1_TO_KCAL]).unwrap();
        for values in ["          0. 1.\n          2. 1.", "          0. 1. 2. 1. 0.5 0.7"] {
            assert_eq!(parse_mess_input(&deck(values)).unwrap().wells["W1"].internal_rotors[0].potential, expected);
        }
        let short = parse_mess_input(&deck("          0. 1. 2.")).unwrap_err();
        assert!(short.contains("W1") && short.contains("4"), "{short}");
        let extra_line = parse_mess_input(&deck("          0. 1. 2. 1.\n          0.5")).unwrap_err();
        assert!(extra_line.contains("W1"), "{extra_line}");
    }

    #[test]
    fn fourier_expansion_lines_give_the_coefficients_in_order() {
        use crate::rrkm::internal_rotor::potential_from_fourier_expansion;
        let deck = rotor_deck(&format!(
            "{ETHANE_CORE}\n        Rotor Hindered\n          Group 3 4 5\n          Axis 1 2\n          Symmetry 3\n          \
             FourierExpansion[kcal/mol] 3\n          0 1.2\n          1 -1.2\n          2 0.3\n        End"
        ));
        let w1 = &parse_mess_input(&deck).unwrap().wells["W1"];
        let expected = potential_from_fourier_expansion(&[1.2 / CM1_TO_KCAL, -1.2 / CM1_TO_KCAL, 0.3 / CM1_TO_KCAL]).unwrap();
        assert_eq!(w1.internal_rotors[0].potential, expected);
    }

    #[test]
    fn a_rotor_of_a_species_given_by_rotational_constants_needs_its_rotational_constant() {
        // MarXus extension: RotationalConstant[1/cm] of the rotor, for species without a geometry.
        let constants = "        RotationalConstants[1/cm] 3\n          0.1081 0.161 0.3105\n        Mass[amu] 75\n        \
Core RigidRotor\n          SymmetryFactor 1\n        End";
        let rotor = "\n        Rotor Hindered\n          Group 5 6 7\n          Axis 1 2\n          Symmetry 3\n          \
Potential[kcal/mol] 2\n          0. 2.45";
        let given = rotor_deck(&format!("{constants}{rotor}\n          RotationalConstant[1/cm] 5.6\n        End"));
        assert_eq!(parse_mess_input(&given).unwrap().wells["W1"].internal_rotors[0].rotational_constant_cm1, Some(5.6));
        let err = parse_mess_input(&rotor_deck(&format!("{constants}{rotor}\n        End"))).unwrap_err();
        assert!(err.contains("W1") && err.contains("RotationalConstant[1/cm]"), "{err}");
    }

    #[test]
    fn rotor_blocks_with_errors_are_refused() {
        let hindered = |body: &str| rotor_deck(&format!("{ETHANE_CORE}\n        Rotor Hindered\n{body}\n        End"));
        let potential = "          Potential[kcal/mol] 2\n          0. 2.45";
        for (deck, expected) in [
            (rotor_deck(&format!("{ETHANE_CORE}\n        Rotor Umbrella\n          Group 3\n          Axis 1 2\n        End")), "Umbrella"),
            (hindered(&format!("          Group 3 4 5\n          Axis 1 2\n          Symmetry 3\n          Potential[kcal/mol] 3\n          1.0 0.0 0.5")), "minimum"),
            (hindered(&format!("          Axis 1 2\n{potential}")), "Group"),
            (hindered(&format!("          Group 3 4 5\n{potential}")), "Axis"),
            (hindered(&format!("          Group 2 3\n          Axis 1 2\n{potential}")), "axis"),
            (hindered(&format!("          Group 3 4 9\n          Axis 1 2\n{potential}")), "8 atoms"),
            (hindered("          Group 3 4 5\n          Axis 1 2"), "potential"),
            (hindered(&format!("          Group 3 4 5\n          Axis 1 2\n          LevelEnergyMax[kcal/mol] 20\n{potential}")), "LevelEnergyMax"),
            (hindered(&format!("          HamiltonSizeMin 100\n          Group 3 4 5\n          Axis 1 2\n{potential}")), "odd"),
            (rotor_deck(&format!("{ETHANE_CORE}\n        Rotor Free\n          Group 3\n          Axis 1 2\n{potential}\n        End")), "Free"),
        ] {
            let err = parse_mess_input(&deck).unwrap_err();
            assert!(err.contains("W1") && err.contains(expected), "expected '{expected}' in: {err}");
        }
    }

    #[test]
    fn a_barrier_given_by_an_inverse_laplace_transform_needs_no_geometry() {
        // The adapter takes the association k(E) of such a barrier from k_inf(T) and the fragments only.
        let start = ILT_BARRIER_DECK.find("    Core PhaseSpaceTheory").unwrap();
        let end = ILT_BARRIER_DECK.find("    InverseLaplaceTransform").unwrap();
        let deck = format!("{}{}", &ILT_BARRIER_DECK[..start], &ILT_BARRIER_DECK[end..]);
        let parsed = parse_mess_input(&deck).expect("an ILT barrier without a geometry");
        assert!(parsed.barriers[0].inverse_laplace_transform.is_some());
        // Nor frequencies.
        let deck = deck.replacen("    Frequencies[1/cm] 1\n    1585.0\n", "", 1);
        let parsed = parse_mess_input(&deck).expect("an ILT barrier without frequencies");
        assert!(parsed.barriers[0].rrho.vibrational_frequencies_cm1.is_empty());
    }

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

    #[test]
    fn fragment_wells_and_energy_partitionings_are_read() {
        // MarXus keywords: `FragmentWell`, `PartnerConcentration[molecule/cm^3]` in a Bimolecular block, and the
        // `FragmentEnergy <kind> ... End` block in a barrier (reports/nonthermal_sources_design.md, Section 15.1).
        let deck = "TemperatureList[K] 300.\nPressureList[torr] 760\nModel\n  Bimolecular P\n    Fragment A\n      RRHO\n        \
Geometry[angstrom] 1\n        H 0 0 0\n        Core RigidRotor\n          SymmetryFactor 1\n        End\n        Frequencies[1/cm] 0\n        \
ZeroEnergy[1/cm] 0\n        ElectronicLevels[1/cm] 1\n          0 1\n      End\n    Fragment B\n      RRHO\n        Geometry[angstrom] 1\n        \
H 0 0 0\n        Core RigidRotor\n          SymmetryFactor 1\n        End\n        Frequencies[1/cm] 0\n        ZeroEnergy[1/cm] 0\n        \
ElectronicLevels[1/cm] 1\n          0 1\n      End\n    GroundEnergy[kcal/mol] 0.0\n    FragmentWell B\n    \
PartnerConcentration[molecule/cm^3] 2.5e14\n  End\n  Barrier X W P\n    RRHO\n      Geometry[angstrom] 1\n      H 0 0 0\n      Core RigidRotor\n        \
SymmetryFactor 1\n      End\n      FragmentEnergy TwoPieceGaussian\n        Mu[kJ/mol]          -32.0  0.5\n        SigmaLeft[1/cm]     -132.0  0.05\n        \
SigmaRight[1/cm]    -1122.0  0.13\n      End\n      Frequencies[1/cm] 0\n      ZeroEnergy[kcal/mol] 5\n      ElectronicLevels[1/cm] 1\n        0 1\n    End\nEnd\n";
        let parsed = parse_mess_input(deck).unwrap();
        let p = &parsed.bimolecular["P"];
        assert_eq!((p.fragment_well.as_deref(), p.partner_concentration_cm3), (Some("B"), Some(2.5e14)));
        match parsed.barriers[0].fragment_energy.as_ref().expect("FragmentEnergy") {
            super::FragmentEnergySpecification::TwoPieceGaussian { mu, sigma_left, sigma_right } => {
                assert!((mu[0] - super::energy_to_cm1(-32.0, "kJ/mol").unwrap()).abs() < 1e-9 && mu[1] == 0.5);
                assert_eq!((sigma_left, sigma_right), (&[-132.0, 0.05], &[-1122.0, 0.13]));
            }
            other => panic!("{other:?}"),
        }
        let modified = deck.replace(
            "FragmentEnergy TwoPieceGaussian\n        Mu[kJ/mol]          -32.0  0.5\n        SigmaLeft[1/cm]     -132.0  0.05\n        SigmaRight[1/cm]    -1122.0  0.13\n",
            "FragmentEnergy ModifiedPrior\n        Order 0.27\n        TemperatureExponent 0\n        ReferenceTemperature[K] 298\n",
        );
        assert_eq!(
            parse_mess_input(&modified).unwrap().barriers[0].fragment_energy,
            Some(super::FragmentEnergySpecification::ModifiedPrior { order: 0.27, temperature_exponent: 0.0, reference_temperature_kelvin: 298.0 })
        );
        let unknown = deck.replace("FragmentEnergy TwoPieceGaussian", "FragmentEnergy Uniform");
        assert!(parse_mess_input(&unknown).unwrap_err().contains("FragmentEnergy"));
    }

    #[test]
    fn a_well_can_have_its_own_collision_parameters() {
        let deck = ONE_WELL_DECK.to_string();
        let with = deck.replace(
            "\nWell W1\n",
            "\nWell W1\n  MarXus\n    Factor[1/cm] 98.3\n    Power 1.0\n    ReferenceTemperature[K] 295\n    Epsilons[1/cm] 57.0 150.2\n    Sigmas[angstrom] 3.74 4.6\n    Masses[amu] 28.0 57.0\n  End\n",
        );
        assert_ne!(with, deck);
        let parsed = parse_mess_input(&with).unwrap();
        let o = &parsed.well_collision["W1"];
        assert_eq!((o.factor_cm1, o.power, o.reference_temperature_kelvin), (Some(98.3), Some(1.0), Some(295.0)));
        assert_eq!((o.epsilons_cm1, o.sigmas_angstrom, o.masses_amu), (Some((57.0, 150.2)), Some((3.74, 4.6)), Some((28.0, 57.0))));
        // The species of the well is read as before (the Masses line is not taken for the species).
        assert_eq!(parsed.wells["W1"].zero_energy_cm1, parse_mess_input(&deck).unwrap().wells["W1"].zero_energy_cm1);
        assert!(parse_mess_input(&deck).unwrap().well_collision.is_empty());
        let bad = with.replace("    Power 1.0\n", "    Exponent 1.0\n");
        assert_ne!(bad, with);
        assert!(parse_mess_input(&bad).unwrap_err().contains("Exponent"));
    }

    #[test]
    fn a_lumped_reactant_state_is_read_with_its_excess_fragment() {
        // MarXus keywords in a Bimolecular block: `LumpedState`, `ExcessFragment <name>`, `PartnerConcentration`
        // (reports/nonthermal_sources_design.md, Section 15.2).
        let deck = "TemperatureList[K] 300.\nPressureList[torr] 760\nModel\n  Bimolecular R\n    Fragment A\n      RRHO\n        \
Geometry[angstrom] 1\n        H 0 0 0\n        Core RigidRotor\n          SymmetryFactor 1\n        End\n        Frequencies[1/cm] 0\n        \
ZeroEnergy[1/cm] 0\n        ElectronicLevels[1/cm] 1\n          0 1\n      End\n    Fragment X\n      RRHO\n        Geometry[angstrom] 1\n        \
H 0 0 0\n        Core RigidRotor\n          SymmetryFactor 1\n        End\n        Frequencies[1/cm] 0\n        ZeroEnergy[1/cm] 0\n        \
ElectronicLevels[1/cm] 1\n          0 1\n      End\n    GroundEnergy[kcal/mol] 0.0\n    LumpedState\n    ExcessFragment X\n    \
PartnerConcentration[molecule/cm^3] 3e15\n  End\nEnd\n";
        let r = &parse_mess_input(deck).unwrap().bimolecular["R"];
        assert_eq!((r.lumped_state, r.excess_fragment.as_deref(), r.partner_concentration_cm3), (true, Some("X"), Some(3e15)));
        let without = deck.replace("    ExcessFragment X\n", "");
        assert!(parse_mess_input(&without).unwrap_err().contains("ExcessFragment"));
    }
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
