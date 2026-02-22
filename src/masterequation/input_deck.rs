use std::collections::HashMap;

use super::reaction_network::{
    CollisionKernelImplementation, CollisionModelParams, MasterEquationSettings, ReactionChannel,
    MultiwellLinearSolver, WellDefinition,
};

macro_rules! bail {
    ($($arg:tt)*) => {
        return Err(format!($($arg)*))
    };
}

/// Energy unit choices for input values.
#[derive(Clone, Copy, Debug)]
pub enum EnergyUnit {
    Cm1,
    KcalMol,
    KjMol,
}

impl EnergyUnit {
    fn from_str(raw: &str) -> Result<Self, String> {
        let s = raw.trim().to_lowercase();
        match s.as_str() {
            "cm-1" | "cm^-1" | "cm1" => Ok(EnergyUnit::Cm1),
            "kcal/mol" | "kcal" => Ok(EnergyUnit::KcalMol),
            "kj/mol" | "kj" => Ok(EnergyUnit::KjMol),
            _ => bail!("Unknown energy unit: {}", raw),
        }
    }

    /// Convert an energy value in this unit to cm^-1.
    pub fn to_cm1(self, value: f64) -> f64 {
        match self {
            EnergyUnit::Cm1 => value,
            EnergyUnit::KcalMol => value / 2.859_144, // 1 cm^-1 = 2.859144e-3 kcal/mol
            EnergyUnit::KjMol => value / 0.011_962_66, // 1 cm^-1 = 0.01196266 kJ/mol
        }
    }
}

/// Species energy definition (exactly one of Eelec or dH0).
#[derive(Clone, Debug)]
pub enum SpeciesEnergy {
    Eelec(f64),
    DH0(f64),
}

/// Species kind: a well (single species) or a bimolecular parent that only names fragments.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SpeciesKind {
    Well,
    BimolecularParent,
}

/// Parsed per-species data block.
#[derive(Clone, Debug)]
pub struct SpeciesDefinition {
    pub name: String,
    pub kind: SpeciesKind,

    /// For bimolecular parent: names of the two fragments.
    pub fragments: Option<(String, String)>,

    /// Geometry (Å) if provided.
    pub geometry_angstrom: Option<Vec<f64>>,

    /// Rotational constants (cm^-1) if provided.
    pub brot_cm1: Option<Vec<f64>>,

    /// Linear or nonlinear flag.
    pub is_linear: Option<bool>,

    /// Vibrational frequencies (cm^-1).
    pub vibrational_frequencies_cm1: Vec<f64>,

    /// Species energy (exactly one field).
    pub energy: Option<SpeciesEnergy>,

    /// Electronic degeneracy.
    pub electronic_degeneracy: Option<f64>,

    /// Spin multiplicity.
    pub spin_multiplicity: Option<f64>,

    /// Symmetry number.
    pub symmetry_number: Option<f64>,

    /// Chirality number.
    pub chirality_number: Option<f64>,

    /// Optional collision-parameter overrides (per well).
    pub collision_params_override: Option<CollisionModelParams>,

    /// Optional grain/offset overrides (per well).
    pub energy_grain_width_cm1: Option<f64>,
    pub lowest_included_grain_index: Option<usize>,
    pub one_past_highest_included_grain_index: Option<usize>,
    pub alignment_offset_in_grains: Option<isize>,
    pub nonreactive_grain_count: Option<usize>,
}

/// Global settings parsed from the input file.
#[derive(Clone, Debug)]
pub struct MasterEquationInputSettings {
    pub temperature_kelvin: f64,
    pub pressure_torr: f64,
    pub boltzmann_constant_wavenumber_per_kelvin: f64,
    pub collision_band_half_width: usize,
    /// Collision-kernel implementation for collisional energy transfer: "spd" or "mess".
    pub collision_kernel_implementation: CollisionKernelImplementation,
    pub outgoing_rate_threshold: f64,
    pub internal_rate_threshold: f64,
    pub bathgas_number_density_prefactor: f64,
    pub mean_speed_prefactor: f64,
    pub enforce_interwell_detailed_balance: bool,
    pub linear_solver: MultiwellLinearSolver,
    pub krylov_tolerance: f64,
    pub krylov_max_iter: usize,
    pub gmres_restart: usize,

    /// Default collision parameters (used unless species overrides).
    pub default_collision_params: CollisionModelParams,

    /// Default grain/offset parameters (used unless species overrides).
    pub default_energy_grain_width_cm1: f64,
    pub default_lowest_included_grain_index: usize,
    pub default_one_past_highest_included_grain_index: usize,
    pub default_alignment_offset_in_grains: isize,
    pub default_nonreactive_grain_count: usize,

    /// Energy unit for Eelec/dH0 fields.
    pub energy_unit: EnergyUnit,
}

/// A reversible reaction link parsed from the reactions section.
#[derive(Clone, Debug)]
pub struct ReactionLink {
    pub left: String,
    pub right: String,
}

/// Full parsed input.
#[derive(Clone, Debug)]
pub struct MasterEquationInput {
    pub settings: MasterEquationInputSettings,
    pub reactions: Vec<ReactionLink>,
    pub species: Vec<SpeciesDefinition>,
}

impl MasterEquationInput {
    /// Build ME settings and wells (with channels) from parsed input.
    pub fn build_settings_and_wells(
        &self,
    ) -> Result<(MasterEquationSettings, Vec<WellDefinition>), String> {
        let settings = MasterEquationSettings {
            temperature_kelvin: self.settings.temperature_kelvin,
            pressure_torr: self.settings.pressure_torr,
            boltzmann_constant_wavenumber_per_kelvin: self
                .settings
                .boltzmann_constant_wavenumber_per_kelvin,
            collision_band_half_width: self.settings.collision_band_half_width,
            collision_kernel_implementation: self.settings.collision_kernel_implementation,
            outgoing_rate_threshold: self.settings.outgoing_rate_threshold,
            internal_rate_threshold: self.settings.internal_rate_threshold,
            bathgas_number_density_prefactor: self.settings.bathgas_number_density_prefactor,
            mean_speed_prefactor: self.settings.mean_speed_prefactor,
            enforce_interwell_detailed_balance: self.settings.enforce_interwell_detailed_balance,
            linear_solver: self.settings.linear_solver,
            krylov_tolerance: self.settings.krylov_tolerance,
            krylov_max_iter: self.settings.krylov_max_iter,
            gmres_restart: self.settings.gmres_restart,
        };

        let mut wells: Vec<WellDefinition> = Vec::new();
        let mut well_index_by_name: HashMap<String, usize> = HashMap::new();

        for sp in &self.species {
            if sp.kind != SpeciesKind::Well {
                continue;
            }

            let collision_params = sp
                .collision_params_override
                .clone()
                .unwrap_or_else(|| self.settings.default_collision_params.clone());

            let well = WellDefinition {
                well_name: sp.name.clone(),
                energy_grain_width_cm1: sp
                    .energy_grain_width_cm1
                    .unwrap_or(self.settings.default_energy_grain_width_cm1),
                lowest_included_grain_index: sp
                    .lowest_included_grain_index
                    .unwrap_or(self.settings.default_lowest_included_grain_index),
                one_past_highest_included_grain_index: sp
                    .one_past_highest_included_grain_index
                    .unwrap_or(self.settings.default_one_past_highest_included_grain_index),
                alignment_offset_in_grains: sp
                    .alignment_offset_in_grains
                    .unwrap_or(self.settings.default_alignment_offset_in_grains),
                nonreactive_grain_count: sp
                    .nonreactive_grain_count
                    .unwrap_or(self.settings.default_nonreactive_grain_count),
                collision_params,
                channels: Vec::new(),
            };

            let idx = wells.len();
            wells.push(well);
            well_index_by_name.insert(sp.name.clone(), idx);
        }

        // Build channels from reactions (reversible).
        for link in &self.reactions {
            let left_is_well = well_index_by_name.contains_key(&link.left);
            let right_is_well = well_index_by_name.contains_key(&link.right);

            if !left_is_well && !right_is_well {
                bail!(
                    "Reaction links non-well species on both sides: {} <--> {}",
                    link.left,
                    link.right
                );
            }

            if left_is_well {
                let from_idx = well_index_by_name[&link.left];
                let connected = if right_is_well {
                    Some(well_index_by_name[&link.right])
                } else {
                    None
                };
                wells[from_idx].channels.push(ReactionChannel {
                    name: format!("{}_to_{}", link.left, link.right),
                    connected_well_index: connected,
                });
            }

            if right_is_well {
                let from_idx = well_index_by_name[&link.right];
                let connected = if left_is_well {
                    Some(well_index_by_name[&link.left])
                } else {
                    None
                };
                wells[from_idx].channels.push(ReactionChannel {
                    name: format!("{}_to_{}", link.right, link.left),
                    connected_well_index: connected,
                });
            }
        }

        Ok((settings, wells))
    }
}

/// Parse a full master-equation input string using the custom text format.
pub fn parse_master_equation_input(input: &str) -> Result<MasterEquationInput, String> {
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum Section {
        None,
        Parameters,
        Reactions,
        Species,
    }

    let mut section = Section::None;

    let mut raw_parameters: HashMap<String, String> = HashMap::new();
    let mut reactions: Vec<ReactionLink> = Vec::new();
    let mut species: Vec<SpeciesDefinition> = Vec::new();

    let mut current_species: Option<SpeciesDefinition> = None;

    for (line_no, raw_line) in input.lines().enumerate() {
        let line_number = line_no + 1;
        let mut line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(idx) = line.find('#') {
            line = line[..idx].trim();
        }
        if line.is_empty() {
            continue;
        }

        let lower = line.to_lowercase();
        if lower == "[parameters]" {
            if let Some(sp) = current_species.take() {
                species.push(sp);
            }
            section = Section::Parameters;
            continue;
        }
        if lower == "[reactions]" {
            if let Some(sp) = current_species.take() {
                species.push(sp);
            }
            section = Section::Reactions;
            continue;
        }
        if lower == "[species]" {
            if let Some(sp) = current_species.take() {
                species.push(sp);
            }
            section = Section::Species;
            continue;
        }

        match section {
            Section::Parameters => {
                let (key, value) = split_key_value(line, line_number)?;
                raw_parameters.insert(key, value);
            }
            Section::Reactions => {
                let reaction = parse_reaction_line(line, line_number)?;
                reactions.push(reaction);
            }
            Section::Species => {
                if line.ends_with(':') {
                    if let Some(sp) = current_species.take() {
                        species.push(sp);
                    }
                    let name = line.trim_end_matches(':').trim().to_string();
                    if name.is_empty() {
                        bail!("Empty species name at line {}", line_number);
                    }
                    current_species = Some(SpeciesDefinition {
                        name,
                        kind: SpeciesKind::Well,
                        fragments: None,
                        geometry_angstrom: None,
                        brot_cm1: None,
                        is_linear: None,
                        vibrational_frequencies_cm1: Vec::new(),
                        energy: None,
                        electronic_degeneracy: None,
                        spin_multiplicity: None,
                        symmetry_number: None,
                        chirality_number: None,
                        collision_params_override: None,
                        energy_grain_width_cm1: None,
                        lowest_included_grain_index: None,
                        one_past_highest_included_grain_index: None,
                        alignment_offset_in_grains: None,
                        nonreactive_grain_count: None,
                    });
                    continue;
                }

                let sp = current_species
                    .as_mut()
                    .ok_or_else(|| format!("Species key outside block at line {}", line_number))?;
                parse_species_key_value(sp, line, line_number)?;
            }
            Section::None => {
                bail!("Line outside any section at {}: {}", line_number, line);
            }
        }
    }

    if let Some(sp) = current_species.take() {
        species.push(sp);
    }

    let settings = build_settings_from_parameters(&raw_parameters)?;
    validate_species_blocks(&species)?;

    Ok(MasterEquationInput {
        settings,
        reactions,
        species,
    })
}

fn split_key_value(line: &str, line_number: usize) -> Result<(String, String), String> {
    let parts: Vec<&str> = line.splitn(2, '=').collect();
    if parts.len() != 2 {
        bail!("Expected key=value at line {}", line_number);
    }
    Ok((parts[0].trim().to_string(), parts[1].trim().to_string()))
}

fn parse_reaction_line(line: &str, line_number: usize) -> Result<ReactionLink, String> {
    let cleaned = line.trim().trim_end_matches(';').trim().to_string();
    if cleaned.is_empty() {
        bail!("Empty reaction line at {}", line_number);
    }
    if !cleaned.contains("<-->") {
        bail!("Reactions must use '<-->' at line {}", line_number);
    }
    let parts: Vec<&str> = cleaned.split("<-->").collect();
    if parts.len() != 2 {
        bail!("Invalid reaction format at line {}", line_number);
    }
    let left = parts[0].trim();
    let right = parts[1].trim();
    if left.is_empty() || right.is_empty() {
        bail!("Invalid reaction endpoints at line {}", line_number);
    }
    Ok(ReactionLink {
        left: left.to_string(),
        right: right.to_string(),
    })
}

fn parse_species_key_value(
    sp: &mut SpeciesDefinition,
    line: &str,
    line_number: usize,
) -> Result<(), String> {
    let (key, value) = split_key_value(line, line_number)?;
    let key_lower = key.to_lowercase();

    match key_lower.as_str() {
        "well" => {
            let flag = parse_bool(&value, line_number)?;
            if flag {
                sp.kind = SpeciesKind::Well;
            }
            Ok(())
        }
        "bimol" => {
            let flag = parse_bool(&value, line_number)?;
            if flag {
                sp.kind = SpeciesKind::BimolecularParent;
            }
            Ok(())
        }
        "fragments" => {
            let parts: Vec<&str> = value.split('+').collect();
            if parts.len() != 2 {
                bail!("Fragments must be 'A + B' at line {}", line_number);
            }
            let a = parts[0].trim().to_string();
            let b = parts[1].trim().to_string();
            if a.is_empty() || b.is_empty() {
                bail!("Invalid fragment names at line {}", line_number);
            }
            sp.fragments = Some((a, b));
            Ok(())
        }
        "geometry" => {
            sp.geometry_angstrom = Some(parse_list_f64(&value, line_number)?);
            Ok(())
        }
        "brot" => {
            sp.brot_cm1 = Some(parse_list_f64(&value, line_number)?);
            Ok(())
        }
        "linear" => {
            sp.is_linear = Some(parse_bool(&value, line_number)?);
            Ok(())
        }
        "vib" | "vibrations" | "frequencies" => {
            sp.vibrational_frequencies_cm1 = parse_list_f64(&value, line_number)?;
            Ok(())
        }
        "eelec" => {
            if sp.energy.is_some() {
                bail!("Both Eelec and dH0 provided at line {}", line_number);
            }
            sp.energy = Some(SpeciesEnergy::Eelec(parse_f64(&value, line_number)?));
            Ok(())
        }
        "dh0" => {
            if sp.energy.is_some() {
                bail!("Both Eelec and dH0 provided at line {}", line_number);
            }
            sp.energy = Some(SpeciesEnergy::DH0(parse_f64(&value, line_number)?));
            Ok(())
        }
        "electronic_degeneracy" => {
            sp.electronic_degeneracy = Some(parse_f64(&value, line_number)?);
            Ok(())
        }
        "spin_multiplicity" => {
            sp.spin_multiplicity = Some(parse_f64(&value, line_number)?);
            Ok(())
        }
        "symmetry_number" => {
            sp.symmetry_number = Some(parse_f64(&value, line_number)?);
            Ok(())
        }
        "chirality_number" => {
            sp.chirality_number = Some(parse_f64(&value, line_number)?);
            Ok(())
        }
        "collision_sigma" => {
            let sigma = parse_f64(&value, line_number)?;
            let mut params = sp
                .collision_params_override
                .clone()
                .unwrap_or(CollisionModelParams {
                    lennard_jones_sigma_angstrom: sigma,
                    lennard_jones_epsilon_kelvin: 0.0,
                    reduced_mass_amu: 0.0,
                    alpha_at_1000K_cm1: 0.0,
                    alpha_temperature_exponent: 0.0,
                });
            params.lennard_jones_sigma_angstrom = sigma;
            sp.collision_params_override = Some(params);
            Ok(())
        }
        "collision_epsilon" => {
            let eps = parse_f64(&value, line_number)?;
            let mut params = sp
                .collision_params_override
                .clone()
                .unwrap_or(CollisionModelParams {
                    lennard_jones_sigma_angstrom: 0.0,
                    lennard_jones_epsilon_kelvin: eps,
                    reduced_mass_amu: 0.0,
                    alpha_at_1000K_cm1: 0.0,
                    alpha_temperature_exponent: 0.0,
                });
            params.lennard_jones_epsilon_kelvin = eps;
            sp.collision_params_override = Some(params);
            Ok(())
        }
        "collision_reduced_mass" => {
            let mu = parse_f64(&value, line_number)?;
            let mut params = sp
                .collision_params_override
                .clone()
                .unwrap_or(CollisionModelParams {
                    lennard_jones_sigma_angstrom: 0.0,
                    lennard_jones_epsilon_kelvin: 0.0,
                    reduced_mass_amu: mu,
                    alpha_at_1000K_cm1: 0.0,
                    alpha_temperature_exponent: 0.0,
                });
            params.reduced_mass_amu = mu;
            sp.collision_params_override = Some(params);
            Ok(())
        }
        "collision_alpha_1000k" => {
            let alpha = parse_f64(&value, line_number)?;
            let mut params = sp
                .collision_params_override
                .clone()
                .unwrap_or(CollisionModelParams {
                    lennard_jones_sigma_angstrom: 0.0,
                    lennard_jones_epsilon_kelvin: 0.0,
                    reduced_mass_amu: 0.0,
                    alpha_at_1000K_cm1: alpha,
                    alpha_temperature_exponent: 0.0,
                });
            params.alpha_at_1000K_cm1 = alpha;
            sp.collision_params_override = Some(params);
            Ok(())
        }
        "collision_alpha_exponent" => {
            let alpha = parse_f64(&value, line_number)?;
            let mut params = sp
                .collision_params_override
                .clone()
                .unwrap_or(CollisionModelParams {
                    lennard_jones_sigma_angstrom: 0.0,
                    lennard_jones_epsilon_kelvin: 0.0,
                    reduced_mass_amu: 0.0,
                    alpha_at_1000K_cm1: 0.0,
                    alpha_temperature_exponent: alpha,
                });
            params.alpha_temperature_exponent = alpha;
            sp.collision_params_override = Some(params);
            Ok(())
        }
        "energy_grain_width" => {
            sp.energy_grain_width_cm1 = Some(parse_f64(&value, line_number)?);
            Ok(())
        }
        "lowest_grain" => {
            sp.lowest_included_grain_index = Some(parse_usize(&value, line_number)?);
            Ok(())
        }
        "highest_grain_exclusive" => {
            sp.one_past_highest_included_grain_index = Some(parse_usize(&value, line_number)?);
            Ok(())
        }
        "alignment_offset" => {
            sp.alignment_offset_in_grains = Some(parse_isize(&value, line_number)?);
            Ok(())
        }
        "nonreactive_grains" => {
            sp.nonreactive_grain_count = Some(parse_usize(&value, line_number)?);
            Ok(())
        }
        _ => bail!("Unknown species key '{}' at line {}", key, line_number),
    }
}

fn build_settings_from_parameters(
    raw: &HashMap<String, String>,
) -> Result<MasterEquationInputSettings, String> {
    let get = |k: &str| {
        raw.get(k)
            .cloned()
            .ok_or_else(|| format!("Missing parameter: {}", k))
    };

    let temperature_kelvin = parse_f64(&get("temperature_kelvin")?, 0)?;
    let pressure_torr = parse_f64(&get("pressure_torr")?, 0)?;
    let boltzmann_constant_wavenumber_per_kelvin =
        parse_f64(&get("boltzmann_constant_wavenumber_per_kelvin")?, 0)?;
    let collision_band_half_width = parse_usize(&get("collision_band_half_width")?, 0)?;
    let collision_kernel_implementation = raw
        .get("collision_kernel_implementation")
        .map(|s| s.trim().to_lowercase())
        .unwrap_or_else(|| "spd".to_string());
    let collision_kernel_implementation = match collision_kernel_implementation.as_str() {
        "spd" => CollisionKernelImplementation::Spd,
        "mess" => CollisionKernelImplementation::Mess,
        other => bail!(
            "Unknown collision_kernel_implementation='{}'. Use 'spd' or 'mess'.",
            other
        ),
    };
    let outgoing_rate_threshold = parse_f64(&get("outgoing_rate_threshold")?, 0)?;
    let internal_rate_threshold = parse_f64(&get("internal_rate_threshold")?, 0)?;
    let bathgas_number_density_prefactor = parse_f64(&get("bathgas_number_density_prefactor")?, 0)?;
    let mean_speed_prefactor = parse_f64(&get("mean_speed_prefactor")?, 0)?;
    let enforce_interwell_detailed_balance = raw
        .get("enforce_interwell_detailed_balance")
        .map(|s| s.trim().to_lowercase())
        .map(|s| matches!(s.as_str(), "true" | "1" | "yes" | "y" | "on"))
        .unwrap_or(false);

    let linear_solver = raw
        .get("linear_solver")
        .map(|s| s.trim().to_lowercase())
        .unwrap_or_else(|| "direct".to_string());
    let linear_solver = match linear_solver.as_str() {
        "direct" => MultiwellLinearSolver::Direct,
        "gmres" => MultiwellLinearSolver::Gmres,
        "bicgstab" | "bicg" => MultiwellLinearSolver::BiCgStab,
        other => bail!(
            "Unknown linear_solver='{}'. Use direct | gmres | bicgstab.",
            other
        ),
    };

    let krylov_tolerance = raw
        .get("krylov_tolerance")
        .map(|s| parse_f64(s, 0))
        .transpose()?
        .unwrap_or(1e-10);
    let krylov_max_iter = raw
        .get("krylov_max_iter")
        .map(|s| parse_usize(s, 0))
        .transpose()?
        .unwrap_or(8000);
    let gmres_restart = raw
        .get("gmres_restart")
        .map(|s| parse_usize(s, 0))
        .transpose()?
        .unwrap_or(50);

    let default_collision_params = CollisionModelParams {
        lennard_jones_sigma_angstrom: parse_f64(&get("default_collision_sigma_angstrom")?, 0)?,
        lennard_jones_epsilon_kelvin: parse_f64(&get("default_collision_epsilon_kelvin")?, 0)?,
        reduced_mass_amu: parse_f64(&get("default_collision_reduced_mass_amu")?, 0)?,
        alpha_at_1000K_cm1: parse_f64(&get("default_collision_alpha_1000k_cm1")?, 0)?,
        alpha_temperature_exponent: parse_f64(&get("default_collision_alpha_exponent")?, 0)?,
    };

    let default_energy_grain_width_cm1 = parse_f64(&get("default_energy_grain_width_cm1")?, 0)?;
    let default_lowest_included_grain_index =
        parse_usize(&get("default_lowest_included_grain_index")?, 0)?;
    let default_one_past_highest_included_grain_index =
        parse_usize(&get("default_one_past_highest_included_grain_index")?, 0)?;
    let default_alignment_offset_in_grains =
        parse_isize(&get("default_alignment_offset_in_grains")?, 0)?;
    let default_nonreactive_grain_count = parse_usize(&get("default_nonreactive_grain_count")?, 0)?;

    let energy_unit = EnergyUnit::from_str(&get("energy_unit")?)?;

    Ok(MasterEquationInputSettings {
        temperature_kelvin,
        pressure_torr,
        boltzmann_constant_wavenumber_per_kelvin,
        collision_band_half_width,
        collision_kernel_implementation,
        outgoing_rate_threshold,
        internal_rate_threshold,
        bathgas_number_density_prefactor,
        mean_speed_prefactor,
        enforce_interwell_detailed_balance,
        linear_solver,
        krylov_tolerance,
        krylov_max_iter,
        gmres_restart,
        default_collision_params,
        default_energy_grain_width_cm1,
        default_lowest_included_grain_index,
        default_one_past_highest_included_grain_index,
        default_alignment_offset_in_grains,
        default_nonreactive_grain_count,
        energy_unit,
    })
}

fn validate_species_blocks(species: &[SpeciesDefinition]) -> Result<(), String> {
    for sp in species {
        match sp.kind {
            SpeciesKind::Well => {
                if sp.fragments.is_some() {
                    bail!("Fragments defined for well '{}'", sp.name);
                }
                if sp.geometry_angstrom.is_none() && sp.brot_cm1.is_none() {
                    bail!("Species '{}' must define geometry or brot", sp.name);
                }
                if sp.geometry_angstrom.is_some() && sp.brot_cm1.is_some() {
                    bail!("Species '{}' cannot define both geometry and brot", sp.name);
                }
                if sp.is_linear.is_none() {
                    bail!("Species '{}' must define linear=true/false", sp.name);
                }
                if sp.vibrational_frequencies_cm1.is_empty() {
                    bail!("Species '{}' must define vibrational frequencies", sp.name);
                }
                if sp.energy.is_none() {
                    bail!("Species '{}' must define Eelec or dH0", sp.name);
                }
                if sp.electronic_degeneracy.is_none()
                    || sp.spin_multiplicity.is_none()
                    || sp.symmetry_number.is_none()
                    || sp.chirality_number.is_none()
                {
                    bail!("Species '{}' must define electronic_degeneracy, spin_multiplicity, symmetry_number, chirality_number", sp.name);
                }
                if let Some(params) = &sp.collision_params_override {
                    if params.lennard_jones_sigma_angstrom <= 0.0
                        || params.lennard_jones_epsilon_kelvin <= 0.0
                        || params.reduced_mass_amu <= 0.0
                        || params.alpha_at_1000K_cm1 <= 0.0
                    {
                        bail!(
                            "Species '{}' collision overrides must define all collision_* fields",
                            sp.name
                        );
                    }
                }
            }
            SpeciesKind::BimolecularParent => {
                if sp.fragments.is_none() {
                    bail!("Bimolecular species '{}' must define fragments", sp.name);
                }
                if sp.geometry_angstrom.is_some()
                    || sp.brot_cm1.is_some()
                    || sp.is_linear.is_some()
                    || !sp.vibrational_frequencies_cm1.is_empty()
                    || sp.energy.is_some()
                    || sp.electronic_degeneracy.is_some()
                    || sp.spin_multiplicity.is_some()
                    || sp.symmetry_number.is_some()
                    || sp.chirality_number.is_some()
                {
                    bail!(
                        "Bimolecular species '{}' must only define fragments",
                        sp.name
                    );
                }
            }
        }
    }
    Ok(())
}

fn parse_bool(raw: &str, line_number: usize) -> Result<bool, String> {
    match raw.trim().to_lowercase().as_str() {
        "true" | "yes" | "1" => Ok(true),
        "false" | "no" | "0" => Ok(false),
        _ => bail!("Invalid boolean '{}' at line {}", raw, line_number),
    }
}

fn parse_f64(raw: &str, line_number: usize) -> Result<f64, String> {
    raw.trim()
        .parse::<f64>()
        .map_err(|_| format!("Invalid number '{}' at line {}", raw, line_number))
}

fn parse_usize(raw: &str, line_number: usize) -> Result<usize, String> {
    raw.trim()
        .parse::<usize>()
        .map_err(|_| format!("Invalid integer '{}' at line {}", raw, line_number))
}

fn parse_isize(raw: &str, line_number: usize) -> Result<isize, String> {
    raw.trim()
        .parse::<isize>()
        .map_err(|_| format!("Invalid integer '{}' at line {}", raw, line_number))
}

fn parse_list_f64(raw: &str, line_number: usize) -> Result<Vec<f64>, String> {
    let mut s = raw.trim().to_string();
    if s.starts_with('[') && s.ends_with(']') {
        s = s[1..s.len() - 1].to_string();
    }
    if s.trim().is_empty() {
        return Ok(Vec::new());
    }
    let mut out = Vec::new();
    for part in s.split(',') {
        let v = parse_f64(part.trim(), line_number)?;
        out.push(v);
    }
    Ok(out)
}
