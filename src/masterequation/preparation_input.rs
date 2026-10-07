//! Deck syntax of a prepared experiment: the top-level `Preparation ... End` block, a MarXus extension of the deck
//! format (design: reports/nonthermal_sources_design.md, Section 8; design note N,
//! papers/Reactant_flux_initiation/MarXus_Nonthermal_Sources.tex, App. A).
//!
//!   Preparation
//!     InitialPopulation
//!       Amount                      1.0
//!       Distribution Thermal
//!         Well                      W1
//!         PreparationTemperature[K] 300
//!       End
//!     End
//!     Source laser
//!       Profile Rectangular
//!         Amount                    1.0
//!         Start[s]                  0
//!         Stop[s]                   1e-8
//!       End
//!       Distribution Gaussian
//!         Well                      W1
//!         Centre[kcal/mol]          25.0
//!         Width[1/cm]               200
//!       End
//!     End
//!     Bath
//!       Segment Start[s] 0     Temperature[K] 300   Pressure[torr] 760
//!       Segment Start[s] 1e-4  Temperature[K] 1500  Pressure[torr] 760
//!     End
//!   End
//!
//! Distributions (`Distribution <kind>`): `Thermal` (Well, PreparationTemperature[K]); `Gaussian` (Well, Centre[unit],
//! Width[unit] = the standard deviation, Representation Density|PerStateWeight); `Tabulated` (Well, File, Representation BinMass|Density|
//! PerStateWeight, EnergyUnit 1/cm|kcal/mol|kJ/mol, Support Complete|Truncate; file lines "lower upper value");
//! `SingleEnergy` (Well, Energy[unit]); `ThermalEntrance` (the entrance flux of the deck's Reactant at the bath
//! temperature); `Mixture` (sub-blocks `Component <weight>` with one distribution each). Every distribution may have
//! `EnergyReference AboveWellGround|Absolute` (default AboveWellGround: above the ZeroEnergy of the well; Absolute: the
//! energy scale of the deck) and `Shift[unit]` (a shift of the whole distribution, e.g. by a photon energy).
//! Profiles (`Profile <kind>`): `Impulse` (Time[s], Amount), `Rectangular` (Start[s], Stop[s], Amount), `Gaussian`
//! (Centre[s], Width[s] = the standard deviation, Amount), `Feed` (Start[s], optional Stop[s], Rate[1/s]), `PrecursorDecay` (Start[s],
//! FormationRate[1/s], TotalLossRate[1/s], PrecursorAmount), `Tabulated` (File; lines "time rate"), `Train`
//! (sub-blocks `Profile <kind>`). Units of energies: [1/cm], [kcal/mol], [kJ/mol]; pressures: [torr], [atm], [bar].

use std::path::Path;

use super::mess_input::{energy_to_cm1, first_token, parse_f64, strip_comment, unit_tag};
use super::chemical_activation_from_mess_input::MessChemicalActivationModel;
use super::chemical_activation_network::Conditions;
use super::chemical_activation_sources::thermal_entrance_source;
use super::mess_input::MessDeck;
use super::prepared_distributions::{
    gaussian, mixture, shifted, tabulated, thermal, EnergyReference, GrainDistribution, Representation, Support, Tabulated,
};
use super::prepared_time_integration::{BathSegment, InitialPopulation, Preparation, SourceChannel};
use super::source_profiles::TimeProfile;
use crate::constants::KB_CM;

/// Keywords that open a sub-block of the Preparation block (closed by `End`).
const BLOCK_KEYWORDS: [&str; 7] = ["Preparation", "InitialPopulation", "Source", "Profile", "Distribution", "Bath", "Component"];

/// A distribution as given in the deck (energies in cm-1); built on the grid by `preparation_from_deck`.
#[derive(Debug, Clone, PartialEq)]
pub enum DistributionSpec {
    Thermal { well: String, temperature_kelvin: f64, shift_cm1: f64 },
    Gaussian { well: String, centre_cm1: f64, width_cm1: f64, representation: Representation, reference: EnergyReference, shift_cm1: f64 },
    Tabulated {
        well: String,
        bins_cm1: Vec<(f64, f64)>,
        values: Vec<f64>,
        representation: Representation,
        reference: EnergyReference,
        support: Support,
        shift_cm1: f64,
    },
    SingleEnergy { well: String, energy_cm1: f64, reference: EnergyReference, shift_cm1: f64 },
    /// The thermal entrance flux of the deck's Reactant at the bath temperature (PO14 eqs. 7, 9).
    ThermalEntrance,
    Mixture(Vec<(f64, DistributionSpec)>),
}

#[derive(Debug, Clone, PartialEq)]
pub struct InitialSpec {
    pub amount: f64,
    pub distribution: DistributionSpec,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SourceSpec {
    pub name: String,
    pub profile: TimeProfile,
    pub distribution: DistributionSpec,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BathSegmentSpec {
    pub start_s: f64,
    pub temperature_kelvin: f64,
    pub pressure_torr: f64,
}

/// The Preparation block of a deck.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PreparationSpec {
    pub initial: Option<InitialSpec>,
    pub sources: Vec<SourceSpec>,
    /// None: the bath of every condition of the deck (TemperatureList x PressureList).
    pub bath: Option<Vec<BathSegmentSpec>>,
}

/// The lines of the `Preparation` block starting at line `start` (comments stripped, empty lines dropped, the opening
/// `Preparation` and its closing `End` included) and the index of the line after it.
pub fn collect_preparation_block(lines: &[String], start: usize) -> Result<(Vec<String>, usize), String> {
    let mut depth = 0usize;
    let mut out = Vec::new();
    for (i, raw) in lines.iter().enumerate().skip(start) {
        let line = strip_comment(raw).trim();
        if line.is_empty() {
            continue;
        }
        let tok = first_token(line).unwrap_or("");
        if tok == "End" {
            out.push(line.to_string());
            depth = depth.checked_sub(1).ok_or("Preparation block: an `End` without a block.")?;
            if depth == 0 {
                return Ok((out, i + 1));
            }
            continue;
        }
        if BLOCK_KEYWORDS.contains(&tok) {
            depth += 1;
        }
        out.push(line.to_string());
    }
    Err("Preparation block: a block is not closed with `End`.".into())
}

/// A sub-block: its keyword with arguments, its key lines (tokens) and its sub-blocks.
#[derive(Debug, Clone, Default)]
struct Block {
    keyword: String,
    args: Vec<String>,
    keys: Vec<Vec<String>>,
    children: Vec<Block>,
}

fn tree(lines: &[String]) -> Result<Block, String> {
    let mut stack: Vec<Block> = Vec::new();
    for line in lines {
        let tokens: Vec<String> = line.split_whitespace().map(str::to_string).collect();
        let tok = tokens[0].as_str();
        if tok == "End" {
            let done = stack.pop().ok_or("Preparation block: unmatched `End`.")?;
            match stack.last_mut() {
                Some(parent) => parent.children.push(done),
                None => return Ok(done),
            }
        } else if BLOCK_KEYWORDS.contains(&tok) {
            stack.push(Block { keyword: tok.to_string(), args: tokens[1..].to_vec(), ..Block::default() });
        } else {
            stack.last_mut().ok_or("Preparation block: a key outside the block.")?.keys.push(tokens);
        }
    }
    Err("Preparation block: a block is not closed with `End`.".into())
}

impl Block {
    fn context(&self) -> String {
        format!("`{} {}`", self.keyword, self.args.join(" "))
    }

    /// The key lines whose first token, without its unit tag, is `name`; refused if given twice.
    fn key(&self, name: &str) -> Result<Option<&Vec<String>>, String> {
        let found: Vec<&Vec<String>> = self.keys.iter().filter(|k| k[0].split('[').next() == Some(name)).collect();
        match found.len() {
            0 => Ok(None),
            1 => Ok(Some(found[0])),
            _ => Err(format!("{}: `{name}` is given twice.", self.context())),
        }
    }

    fn value(&self, name: &str) -> Result<Option<String>, String> {
        match self.key(name)? {
            None => Ok(None),
            Some(tokens) => tokens.get(1).cloned().map(Some).ok_or_else(|| format!("{}: `{name}` has no value.", self.context())),
        }
    }

    fn required(&self, name: &str) -> Result<String, String> {
        self.value(name)?.ok_or_else(|| format!("{}: `{name}` is required.", self.context()))
    }

    fn number(&self, name: &str) -> Result<Option<f64>, String> {
        self.value(name)?.map(|v| parse_f64(&v).map_err(|e| format!("{}: `{name}`: {e}", self.context()))).transpose()
    }

    fn required_number(&self, name: &str) -> Result<f64, String> {
        self.number(name)?.ok_or_else(|| format!("{}: `{name}` is required.", self.context()))
    }

    /// An energy `name[unit]` in cm-1.
    fn energy(&self, name: &str) -> Result<Option<f64>, String> {
        match self.key(name)? {
            None => Ok(None),
            Some(tokens) => {
                let unit = unit_tag(&tokens[0]).ok_or_else(|| format!("{}: `{name}` needs a unit, e.g. {name}[1/cm].", self.context()))?;
                let value = parse_f64(tokens.get(1).ok_or_else(|| format!("{}: `{name}` has no value.", self.context()))?)?;
                energy_to_cm1(value, unit).map(Some)
            }
        }
    }

    fn required_energy(&self, name: &str) -> Result<f64, String> {
        self.energy(name)?.ok_or_else(|| format!("{}: `{name}[unit]` is required.", self.context()))
    }

    /// Refuses keys and sub-blocks that are not in `keys` / `blocks`.
    fn only(&self, keys: &[&str], blocks: &[&str]) -> Result<(), String> {
        for k in &self.keys {
            let name = k[0].split('[').next().unwrap_or("");
            if !keys.contains(&name) {
                return Err(format!("{}: unknown keyword `{}` (allowed: {}).", self.context(), k[0], keys.join(", ")));
            }
        }
        for c in &self.children {
            if !blocks.contains(&c.keyword.as_str()) {
                return Err(format!("{}: unexpected block `{}`.", self.context(), c.keyword));
            }
        }
        Ok(())
    }

    fn children_named(&self, keyword: &str) -> Vec<&Block> {
        self.children.iter().filter(|c| c.keyword == keyword).collect()
    }

    fn single_child(&self, keyword: &str) -> Result<&Block, String> {
        let found = self.children_named(keyword);
        match found.len() {
            1 => Ok(found[0]),
            0 => Err(format!("{}: a `{keyword} ...` block is required.", self.context())),
            _ => Err(format!("{}: more than one `{keyword}` block.", self.context())),
        }
    }

    fn kind(&self) -> Result<&str, String> {
        self.args.first().map(String::as_str).ok_or_else(|| format!("`{}` needs a kind, e.g. `{} Thermal`.", self.keyword, self.keyword))
    }
}

fn keyword_of<T: Copy>(value: &str, table: &[(&str, T)], what: &str) -> Result<T, String> {
    let wanted = value.to_lowercase();
    table
        .iter()
        .find(|(k, _)| k.to_lowercase() == wanted)
        .map(|(_, v)| *v)
        .ok_or_else(|| format!("{what} `{value}`: allowed {}.", table.iter().map(|(k, _)| *k).collect::<Vec<_>>().join(", ")))
}

fn representation(b: &Block, default: Representation) -> Result<Representation, String> {
    match b.value("Representation")? {
        None => Ok(default),
        Some(v) => keyword_of(
            &v,
            &[("BinMass", Representation::BinMass), ("Density", Representation::Density), ("PerStateWeight", Representation::PerStateWeight)],
            "Representation",
        ),
    }
}

fn reference(b: &Block) -> Result<EnergyReference, String> {
    match b.value("EnergyReference")? {
        None => Ok(EnergyReference::AboveWellGround),
        Some(v) => keyword_of(&v, &[("AboveWellGround", EnergyReference::AboveWellGround), ("Absolute", EnergyReference::Absolute)], "EnergyReference"),
    }
}

fn distribution(b: &Block, base: Option<&Path>) -> Result<DistributionSpec, String> {
    let common = ["Well", "EnergyReference", "Shift"];
    let shift = b.energy("Shift")?.unwrap_or(0.0);
    let with = |extra: &[&'static str]| -> Vec<&'static str> { common.iter().chain(extra).copied().collect() };
    match b.kind()? {
        "Thermal" => {
            b.only(&with(&["PreparationTemperature"]), &[])?;
            Ok(DistributionSpec::Thermal { well: b.required("Well")?, temperature_kelvin: b.required_number("PreparationTemperature")?, shift_cm1: shift })
        }
        "Gaussian" => {
            b.only(&with(&["Centre", "Width", "Representation"]), &[])?;
            Ok(DistributionSpec::Gaussian {
                well: b.required("Well")?,
                centre_cm1: b.required_energy("Centre")?,
                width_cm1: b.required_energy("Width")?,
                representation: representation(b, Representation::Density)?,
                reference: reference(b)?,
                shift_cm1: shift,
            })
        }
        "Tabulated" => {
            b.only(&with(&["File", "Representation", "EnergyUnit", "Support"]), &[])?;
            let unit = b.value("EnergyUnit")?.unwrap_or_else(|| "1/cm".into());
            let factor = energy_to_cm1(1.0, &unit)?;
            let file = b.required("File")?;
            let path = base.map_or_else(|| Path::new(&file).to_path_buf(), |d| d.join(&file));
            let (bins_cm1, values) = read_distribution_file(&path, factor)?;
            let support = match b.value("Support")? {
                None => Support::Complete,
                Some(v) => keyword_of(&v, &[("Complete", Support::Complete), ("Truncate", Support::Truncate)], "Support")?,
            };
            Ok(DistributionSpec::Tabulated {
                well: b.required("Well")?,
                bins_cm1,
                values,
                representation: representation(b, Representation::BinMass)?,
                reference: reference(b)?,
                support,
                shift_cm1: shift,
            })
        }
        "SingleEnergy" => {
            b.only(&with(&["Energy"]), &[])?;
            Ok(DistributionSpec::SingleEnergy { well: b.required("Well")?, energy_cm1: b.required_energy("Energy")?, reference: reference(b)?, shift_cm1: shift })
        }
        "ThermalEntrance" => {
            b.only(&[], &[])?;
            Ok(DistributionSpec::ThermalEntrance)
        }
        "Mixture" => {
            b.only(&[], &["Component"])?;
            let mut components = Vec::new();
            for c in b.children_named("Component") {
                c.only(&[], &["Distribution"])?;
                let weight = parse_f64(c.args.first().ok_or_else(|| "`Component` needs its weight, e.g. `Component 0.3`.".to_string())?)?;
                components.push((weight, distribution(c.single_child("Distribution")?, base)?));
            }
            if components.is_empty() {
                return Err("`Distribution Mixture` needs `Component <weight>` blocks.".into());
            }
            Ok(DistributionSpec::Mixture(components))
        }
        other => Err(format!(
            "`Distribution {other}`: unknown kind (Thermal, Gaussian, Tabulated, SingleEnergy, ThermalEntrance, Mixture)."
        )),
    }
}

fn profile(b: &Block, base: Option<&Path>) -> Result<TimeProfile, String> {
    let p = match b.kind()? {
        "Impulse" => {
            b.only(&["Time", "Amount"], &[])?;
            TimeProfile::Impulse { time_s: b.required_number("Time")?, amount: b.required_number("Amount")? }
        }
        "Rectangular" => {
            b.only(&["Start", "Stop", "Amount"], &[])?;
            TimeProfile::Rectangular { start_s: b.required_number("Start")?, end_s: b.required_number("Stop")?, amount: b.required_number("Amount")? }
        }
        "Gaussian" => {
            b.only(&["Centre", "Width", "Amount"], &[])?;
            TimeProfile::Gaussian { centre_s: b.required_number("Centre")?, sigma_s: b.required_number("Width")?, amount: b.required_number("Amount")? }
        }
        "Feed" => {
            b.only(&["Start", "Stop", "Rate"], &[])?;
            TimeProfile::Feed { start_s: b.number("Start")?.unwrap_or(0.0), end_s: b.number("Stop")?, rate: b.required_number("Rate")? }
        }
        "PrecursorDecay" => {
            b.only(&["Start", "FormationRate", "TotalLossRate", "PrecursorAmount"], &[])?;
            TimeProfile::PrecursorDecay {
                start_s: b.number("Start")?.unwrap_or(0.0),
                formation_rate_s_inv: b.required_number("FormationRate")?,
                total_loss_rate_s_inv: b.required_number("TotalLossRate")?,
                precursor_amount: b.required_number("PrecursorAmount")?,
            }
        }
        "Tabulated" => {
            b.only(&["File"], &[])?;
            let file = b.required("File")?;
            read_profile_file(&base.map_or_else(|| Path::new(&file).to_path_buf(), |d| d.join(&file)))?
        }
        "Train" => {
            b.only(&[], &["Profile"])?;
            TimeProfile::Train(b.children_named("Profile").into_iter().map(|c| profile(c, base)).collect::<Result<_, _>>()?)
        }
        other => {
            return Err(format!("`Profile {other}`: unknown kind (Impulse, Rectangular, Gaussian, Feed, PrecursorDecay, Tabulated, Train)."))
        }
    };
    p.validate().map_err(|e| format!("{}: {e}", b.context()))?;
    Ok(p)
}

fn pressure_torr(tokens: &[String], k: usize) -> Result<f64, String> {
    let unit = unit_tag(&tokens[k]).unwrap_or("").to_lowercase();
    let value = parse_f64(tokens.get(k + 1).ok_or("Bath segment: a pressure needs a value.")?)?;
    match unit.as_str() {
        "torr" => Ok(value),
        "atm" => Ok(value * 760.0),
        "bar" => Ok(value * 1.0e5 * 760.0 / 101_325.0),
        _ => Err(format!("Bath segment: pressure unit `{unit}` (torr, atm or bar).")),
    }
}

/// Interprets a Preparation block (as returned by `collect_preparation_block`); relative file names are read from the
/// current directory (`parse_preparation_in`: from a given directory).
pub fn parse_preparation(lines: &[String]) -> Result<PreparationSpec, String> {
    parse_preparation_in(lines, None)
}

/// `parse_preparation` with files relative to `base`.
pub fn parse_preparation_in(lines: &[String], base: Option<&Path>) -> Result<PreparationSpec, String> {
    let root = tree(lines)?;
    if root.keyword != "Preparation" {
        return Err("A Preparation block must start with `Preparation`.".into());
    }
    root.only(&[], &["InitialPopulation", "Source", "Bath"])?;
    let mut spec = PreparationSpec::default();
    let initial = root.children_named("InitialPopulation");
    if initial.len() > 1 {
        return Err("Preparation: more than one InitialPopulation block.".into());
    }
    if let Some(b) = initial.first() {
        b.only(&["Amount"], &["Distribution"])?;
        spec.initial = Some(InitialSpec { amount: b.number("Amount")?.unwrap_or(1.0), distribution: distribution(b.single_child("Distribution")?, base)? });
    }
    for b in root.children_named("Source") {
        b.only(&[], &["Profile", "Distribution"])?;
        let name = b.args.first().cloned().ok_or("`Source` needs a name, e.g. `Source laser`.")?;
        if spec.sources.iter().any(|s| s.name == name) {
            return Err(format!("Preparation: two sources named `{name}`."));
        }
        spec.sources.push(SourceSpec { name, profile: profile(b.single_child("Profile")?, base)?, distribution: distribution(b.single_child("Distribution")?, base)? });
    }
    let baths = root.children_named("Bath");
    if baths.len() > 1 {
        return Err("Preparation: more than one Bath block.".into());
    }
    if let Some(b) = baths.first() {
        b.only(&["Segment"], &[])?;
        let mut segments = Vec::new();
        for tokens in &b.keys {
            let find = |name: &str| tokens.iter().position(|t| t.split('[').next() == Some(name));
            let (s, t, p) = (find("Start"), find("Temperature"), find("Pressure"));
            let (Some(s), Some(t), Some(p)) = (s, t, p) else {
                return Err(format!("Bath segment `{}`: needs Start[s], Temperature[K] and Pressure[unit].", tokens.join(" ")));
            };
            let number = |k: usize| parse_f64(tokens.get(k + 1).map(String::as_str).unwrap_or(""));
            segments.push(BathSegmentSpec { start_s: number(s)?, temperature_kelvin: number(t)?, pressure_torr: pressure_torr(tokens, p)? });
        }
        if segments.is_empty() || segments[0].start_s != 0.0 || segments.windows(2).any(|w| w[1].start_s <= w[0].start_s) {
            return Err("Bath: the segments must start at 0 s and follow in increasing order.".into());
        }
        spec.bath = Some(segments);
    }
    if spec.initial.is_none() && spec.sources.is_empty() {
        return Err("Preparation: neither an InitialPopulation nor a Source is given.".into());
    }
    Ok(spec)
}

/// The library preparation (`prepared_time_integration::Preparation`) of `spec` on the network of `model` (built
/// from `deck`). Energies above the ground of a well are placed with its ZeroEnergy on the absolute scale of the
/// deck. Without a Bath block, the bath is `conditions` for all times. `ThermalEntrance` is the entrance flux of the
/// deck's Reactant at the temperature of the first bath segment.
pub fn preparation_from_deck(
    spec: &PreparationSpec,
    deck: &MessDeck,
    model: &MessChemicalActivationModel,
    conditions: &Conditions,
) -> Result<Preparation, String> {
    let bath: Vec<BathSegment> = match &spec.bath {
        Some(segments) => segments
            .iter()
            .map(|s| BathSegment { start_s: s.start_s, conditions: Conditions { temperature_kelvin: s.temperature_kelvin, pressure_torr: s.pressure_torr } })
            .collect(),
        None => vec![BathSegment { start_s: 0.0, conditions: conditions.clone() }],
    };
    let kt = KB_CM * bath[0].conditions.temperature_kelvin;
    let initial = match &spec.initial {
        Some(i) => Some(InitialPopulation { amount: i.amount, distribution: build_distribution(&i.distribution, deck, model, kt)? }),
        None => None,
    };
    let channels = spec
        .sources
        .iter()
        .map(|s| {
            Ok(SourceChannel {
                name: s.name.clone(),
                distribution: build_distribution(&s.distribution, deck, model, kt).map_err(|e| format!("Source '{}': {e}", s.name))?,
                profile: s.profile.clone(),
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    Ok(Preparation { initial, channels, bath })
}

/// One distribution of the deck on the grid of the network.
pub fn build_distribution(spec: &DistributionSpec, deck: &MessDeck, model: &MessChemicalActivationModel, kt_cm1: f64) -> Result<GrainDistribution, String> {
    let network = &model.network;
    let well_index = |name: &str| -> Result<usize, String> {
        network.wells.iter().position(|w| w.name == name).ok_or_else(|| format!("Distribution: `{name}` is not a well of the deck."))
    };
    // Energy on the absolute scale of the deck.
    let absolute = |well: &str, energy: f64, reference: EnergyReference| -> Result<f64, String> {
        Ok(match reference {
            EnergyReference::Absolute => energy,
            EnergyReference::AboveWellGround => deck.wells.get(well).ok_or_else(|| format!("Distribution: `{well}` is not a well of the deck."))?.zero_energy_cm1 + energy,
        })
    };
    let shift = |d: GrainDistribution, shift_cm1: f64| -> Result<GrainDistribution, String> {
        if shift_cm1 == 0.0 {
            Ok(d)
        } else {
            shifted(network, &d, shift_cm1, Support::Complete)
        }
    };
    match spec {
        DistributionSpec::Thermal { well, temperature_kelvin, shift_cm1 } => shift(thermal(network, well_index(well)?, *temperature_kelvin)?, *shift_cm1),
        DistributionSpec::Gaussian { well, centre_cm1, width_cm1, representation, reference, shift_cm1 } => {
            let w = well_index(well)?;
            shift(gaussian(network, w, absolute(well, *centre_cm1, *reference)?, *width_cm1, *representation, EnergyReference::Absolute)?, *shift_cm1)
        }
        DistributionSpec::Tabulated { well, bins_cm1, values, representation, reference, support, shift_cm1 } => {
            let w = well_index(well)?;
            let offset = absolute(well, 0.0, *reference)?;
            let table = Tabulated {
                bins: bins_cm1.iter().map(|&(lo, hi)| (lo + offset, hi + offset)).collect(),
                values: values.clone(),
                representation: *representation,
                reference: EnergyReference::Absolute,
            };
            shift(tabulated(network, w, &table, *support)?, *shift_cm1)
        }
        DistributionSpec::SingleEnergy { well, energy_cm1, reference, shift_cm1 } => {
            let w = well_index(well)?;
            let e = absolute(well, *energy_cm1, *reference)?;
            let table = Tabulated { bins: vec![(e, e)], values: vec![1.0], representation: Representation::BinMass, reference: EnergyReference::Absolute };
            shift(tabulated(network, w, &table, Support::Complete)?, *shift_cm1)
        }
        DistributionSpec::ThermalEntrance => {
            if model.entrance_channels.is_empty() {
                return Err("Distribution ThermalEntrance: the deck has no Reactant with channels into the wells.".into());
            }
            let mass = thermal_entrance_source(network, &model.entrance_channels, kt_cm1)?;
            Ok(GrainDistribution { mass, lost_fraction: 0.0, description: "thermal entrance flux of the Reactant".into() })
        }
        DistributionSpec::Mixture(components) => {
            let built = components.iter().map(|(a, c)| Ok((*a, build_distribution(c, deck, model, kt_cm1)?))).collect::<Result<Vec<_>, String>>()?;
            mixture(&built)
        }
    }
}

fn numbers_of(line: &str) -> Result<Vec<f64>, String> {
    line.split(|c: char| c.is_whitespace() || c == ',').filter(|t| !t.is_empty()).map(parse_f64).collect()
}

/// A distribution file: lines "lower upper value" (whitespace or commas; `#` comments), energies multiplied by
/// `energy_factor` (to cm-1).
pub fn read_distribution_file(path: &Path, energy_factor: f64) -> Result<(Vec<(f64, f64)>, Vec<f64>), String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("Distribution file {}: {e}", path.display()))?;
    let (mut bins, mut values) = (Vec::new(), Vec::new());
    for line in text.lines().map(|l| l.split('#').next().unwrap_or("").trim()).filter(|l| !l.is_empty()) {
        let x = numbers_of(line).map_err(|e| format!("Distribution file {}: {e}", path.display()))?;
        if x.len() != 3 {
            return Err(format!("Distribution file {}: line `{line}` needs lower, upper and value.", path.display()));
        }
        bins.push((x[0] * energy_factor, x[1] * energy_factor));
        values.push(x[2]);
    }
    Ok((bins, values))
}

/// A rate file: lines "time rate" (s, 1/s).
pub fn read_profile_file(path: &Path) -> Result<TimeProfile, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("Profile file {}: {e}", path.display()))?;
    let (mut times_s, mut rates) = (Vec::new(), Vec::new());
    for line in text.lines().map(|l| l.split('#').next().unwrap_or("").trim()).filter(|l| !l.is_empty()) {
        let x = numbers_of(line).map_err(|e| format!("Profile file {}: {e}", path.display()))?;
        if x.len() != 2 {
            return Err(format!("Profile file {}: line `{line}` needs time and rate.", path.display()));
        }
        times_s.push(x[0]);
        rates.push(x[1]);
    }
    let p = TimeProfile::Tabulated { times_s, rates };
    p.validate()?;
    Ok(p)
}

#[cfg(test)]
mod tests {
    use super::*;

    const BLOCK: &str = "
Preparation
  InitialPopulation
    Amount                      0.5
    Distribution Thermal
      Well                      W1
      PreparationTemperature[K] 300
    End
  End
  Source laser                  ! photolysis
    Profile Rectangular
      Amount                    1.0
      Start[s]                  0
      Stop[s]                   1e-8
    End
    Distribution Gaussian
      Well                      W1
      Centre[kcal/mol]          25.0
      Width[1/cm]               200
      Representation            PerStateWeight
    End
  End
  Source feed
    Profile Train
      Profile Impulse
        Time[s]                 1e-6
        Amount                  0.25
      End
      Profile Feed
        Start[s]                1e-6
        Rate[1/s]               10
      End
    End
    Distribution Mixture
      Component 0.4
        Distribution SingleEnergy
          Well                  W1
          Energy[1/cm]          12000
          EnergyReference       Absolute
        End
      End
      Component 0.6
        Distribution ThermalEntrance
        End
      End
    End
  End
  Bath
    Segment Start[s] 0     Temperature[K] 300   Pressure[torr] 760
    Segment Start[s] 1e-4  Temperature[K] 1500  Pressure[atm]  2
  End
End
";

    fn block_lines(text: &str) -> Vec<String> {
        text.lines().map(|l| l.to_string()).collect()
    }

    #[test]
    fn a_preparation_block_is_read_into_its_parts() {
        let lines = block_lines(BLOCK);
        let start = lines.iter().position(|l| l.trim() == "Preparation").unwrap();
        let (block, next) = collect_preparation_block(&lines, start).unwrap();
        assert_eq!(next, lines.len());
        let spec = parse_preparation(&block).unwrap();
        let initial = spec.initial.unwrap();
        assert_eq!(initial.amount, 0.5);
        assert_eq!(initial.distribution, DistributionSpec::Thermal { well: "W1".into(), temperature_kelvin: 300.0, shift_cm1: 0.0 });
        assert_eq!(spec.sources.len(), 2);
        let laser = &spec.sources[0];
        assert_eq!(laser.name, "laser");
        assert_eq!(laser.profile, TimeProfile::Rectangular { start_s: 0.0, end_s: 1e-8, amount: 1.0 });
        match &laser.distribution {
            DistributionSpec::Gaussian { well, centre_cm1, width_cm1, representation, reference, .. } => {
                assert_eq!(well, "W1");
                assert!((centre_cm1 - 25.0 / crate::constants::CM1_TO_KCAL).abs() < 1e-9);
                assert_eq!((*width_cm1, *representation, *reference), (200.0, Representation::PerStateWeight, EnergyReference::AboveWellGround));
            }
            other => panic!("{other:?}"),
        }
        let feed = &spec.sources[1];
        assert_eq!(
            feed.profile,
            TimeProfile::Train(vec![
                TimeProfile::Impulse { time_s: 1e-6, amount: 0.25 },
                TimeProfile::Feed { start_s: 1e-6, end_s: None, rate: 10.0 },
            ])
        );
        match &feed.distribution {
            DistributionSpec::Mixture(components) => {
                assert_eq!(components.len(), 2);
                assert_eq!(components[0].0, 0.4);
                assert_eq!(components[0].1, DistributionSpec::SingleEnergy { well: "W1".into(), energy_cm1: 12000.0, reference: EnergyReference::Absolute, shift_cm1: 0.0 });
                assert_eq!(components[1].1, DistributionSpec::ThermalEntrance);
            }
            other => panic!("{other:?}"),
        }
        let bath = spec.bath.unwrap();
        assert_eq!(bath.len(), 2);
        assert_eq!((bath[1].start_s, bath[1].temperature_kelvin), (1e-4, 1500.0));
        assert!((bath[1].pressure_torr - 1520.0).abs() < 1e-9);
    }

    #[test]
    fn errors_in_a_preparation_block_are_reported() {
        let cases = [
            ("Preparation\n  Source s\n    Profile Impulse\n      Time[s] 0\n      Amount 1\n    End\n  End\nEnd", "Distribution"),
            ("Preparation\n  Source s\n    Profile Impulse\n      Time[s] 0\n      Amount 1\n    End\n    Distribution Thermal\n      Well W\n    End\n  End\nEnd", "PreparationTemperature"),
            ("Preparation\n  Source s\n    Profile Sawtooth\n    End\n  End\nEnd", "Sawtooth"),
            ("Preparation\n  Colour blue\nEnd", "Colour"),
            ("Preparation\n  InitialPopulation\n    Amount 1\n", "End"),
            ("Preparation\n  Source s\n    Profile Impulse\n      Time[s] 0\n      Time[s] 1\n      Amount 1\n    End\n  End\nEnd", "twice"),
        ];
        for (text, expected) in cases {
            let lines = block_lines(text);
            let result = collect_preparation_block(&lines, 0).and_then(|(block, _)| parse_preparation(&block));
            let e = result.expect_err(text);
            assert!(e.contains(expected), "{text}\n-> {e}");
        }
    }

    #[test]
    fn tabulated_files_are_read_with_comments_and_units() {
        let dir = std::env::temp_dir().join(format!("marxus_preparation_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("dist.dat"), "# lower upper value\n0 1 0.5\n1, 2, 1.5\n").unwrap();
        std::fs::write(dir.join("rate.dat"), "# t R\n0 0\n1e-6 2\n2e-6 0\n").unwrap();
        let (bins, values) = read_distribution_file(&dir.join("dist.dat"), 1.0 / crate::constants::CM1_TO_KCAL).unwrap();
        assert_eq!(values, vec![0.5, 1.5]);
        assert!((bins[1].1 - 2.0 / crate::constants::CM1_TO_KCAL).abs() < 1e-9);
        assert_eq!(read_profile_file(&dir.join("rate.dat")).unwrap(), TimeProfile::Tabulated { times_s: vec![0.0, 1e-6, 2e-6], rates: vec![0.0, 2.0, 0.0] });
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
