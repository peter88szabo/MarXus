use std::collections::HashMap;

use super::input_deck::{
    MasterEquationInputSettings, ReactionLink, SpeciesDefinition, SpeciesKind,
};
use super::reaction_network::{MasterEquationSettings, ReactionChannel, WellDefinition};

pub fn build_master_equation_settings(
    settings: &MasterEquationInputSettings,
) -> MasterEquationSettings {
    MasterEquationSettings {
        temperature_kelvin: settings.temperature_kelvin,
        pressure_torr: settings.pressure_torr,
        boltzmann_constant_wavenumber_per_kelvin: settings.boltzmann_constant_wavenumber_per_kelvin,
        collision_band_half_width: settings.collision_band_half_width,
        collision_kernel_implementation: settings.collision_kernel_implementation,
        outgoing_rate_threshold: settings.outgoing_rate_threshold,
        internal_rate_threshold: settings.internal_rate_threshold,
        bathgas_number_density_prefactor: settings.bathgas_number_density_prefactor,
        mean_speed_prefactor: settings.mean_speed_prefactor,
        enforce_interwell_detailed_balance: settings.enforce_interwell_detailed_balance,
        linear_solver: settings.linear_solver,
        krylov_tolerance: settings.krylov_tolerance,
        krylov_max_iter: settings.krylov_max_iter,
        gmres_restart: settings.gmres_restart,
    }
}

pub fn build_wells_and_channels(
    species: &[SpeciesDefinition],
    reactions: &[ReactionLink],
    settings: &MasterEquationInputSettings,
) -> Result<Vec<WellDefinition>, String> {
    let mut wells: Vec<WellDefinition> = Vec::new();
    let mut well_index_by_name: HashMap<String, usize> = HashMap::new();

    for sp in species {
        if sp.kind != SpeciesKind::Well {
            continue;
        }

        let collision_params = sp
            .collision_params_override
            .clone()
            .unwrap_or_else(|| settings.default_collision_params.clone());

        let well = WellDefinition {
            well_name: sp.name.clone(),
            energy_grain_width_cm1: sp
                .energy_grain_width_cm1
                .unwrap_or(settings.default_energy_grain_width_cm1),
            lowest_included_grain_index: sp
                .lowest_included_grain_index
                .unwrap_or(settings.default_lowest_included_grain_index),
            one_past_highest_included_grain_index: sp
                .one_past_highest_included_grain_index
                .unwrap_or(settings.default_one_past_highest_included_grain_index),
            alignment_offset_in_grains: sp
                .alignment_offset_in_grains
                .unwrap_or(settings.default_alignment_offset_in_grains),
            nonreactive_grain_count: sp
                .nonreactive_grain_count
                .unwrap_or(settings.default_nonreactive_grain_count),
            collision_params,
            channels: Vec::new(),
        };

        let idx = wells.len();
        wells.push(well);
        well_index_by_name.insert(sp.name.clone(), idx);
    }

    for link in reactions {
        let left_is_well = well_index_by_name.contains_key(&link.left);
        let right_is_well = well_index_by_name.contains_key(&link.right);

        if !left_is_well && !right_is_well {
            return Err(format!(
                "Reaction links non-well species on both sides: {} <--> {}",
                link.left, link.right
            ));
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

    Ok(wells)
}
