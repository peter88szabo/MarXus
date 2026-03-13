//! Network-level microcanonical builder glue.
//!
//! The low-level microcanonical builder (`microcanonical_builder`) works on *ordered arrays*:
//!   - `wells: &[WellDefinition]`
//!   - `well_models: &[SpeciesMicroModel]` (same order as `wells`)
//!   - `channel_micro_models: &[Vec<ChannelMicroModel>]` (same order as `wells`, and per-well
//!     ordering matches `WellDefinition.channels`)
//!
//! In practice, user inputs (including MESS-like decks) are usually keyed by *names*:
//!   - wells by name (e.g., "G2", "G3", ...)
//!   - channels by label (e.g., "G2_to_G3" or "R_to_G2")
//!
//! This module provides a thin adapter that:
//!   1) walks the already-built `WellDefinition` network topology
//!   2) looks up the corresponding micro models by name
//!   3) produces the ordered arrays required by the microcanonical builder
//!   4) returns the computed `MicrocanonicalNetworkData`
//!
//! It intentionally does not implement any physics model itself.

use std::collections::HashMap;

use crate::masterequation::microcanonical_builder::{
    build_microcanonical_network_data, ChannelMicroModel, MicrocanonicalNetworkData,
    SpeciesMicroModel,
};
use crate::masterequation::reaction_network::WellDefinition;

/// Build microcanonical arrays for a named network.
///
/// Inputs:
/// - `wells`: the network topology with ordered wells + ordered channels per well
/// - `well_models_by_name`: micro models keyed by `WellDefinition.well_name`
/// - `channel_models_by_name`: channel micro models keyed by `ReactionChannel.name`
///
/// Errors if any well or channel is missing a micro model.
pub fn build_microcanonical_network_data_from_named_models(
    wells: &[WellDefinition],
    well_models_by_name: &HashMap<String, SpeciesMicroModel>,
    channel_models_by_name: &HashMap<String, ChannelMicroModel>,
) -> Result<MicrocanonicalNetworkData, String> {
    let mut well_models: Vec<SpeciesMicroModel> = Vec::with_capacity(wells.len());
    let mut channel_micro_models: Vec<Vec<ChannelMicroModel>> = Vec::with_capacity(wells.len());

    for well in wells {
        let well_model = well_models_by_name.get(&well.well_name).ok_or_else(|| {
            format!(
                "Missing micro model for well '{}'. Provide a `SpeciesMicroModel` keyed by this name.",
                well.well_name
            )
        })?;
        well_models.push(well_model.clone());

        let mut per_well_channels: Vec<ChannelMicroModel> = Vec::with_capacity(well.channels.len());
        for ch in &well.channels {
            let ch_model = channel_models_by_name.get(&ch.name).ok_or_else(|| {
                format!(
                    "Missing micro model for channel '{}' (from well '{}'). Provide a `ChannelMicroModel` keyed by this channel name.",
                    ch.name, well.well_name
                )
            })?;
            per_well_channels.push(ch_model.clone());
        }
        channel_micro_models.push(per_well_channels);
    }

    build_microcanonical_network_data(wells, &well_models, &channel_micro_models)
}
