use super::reaction_network::MicrocanonicalProvider;

/// Minimal placeholder microcanonical provider for multi-well examples.
///
/// This is not a physical model; it only exists to make examples runnable.
pub struct PlaceholderMicroData;

impl MicrocanonicalProvider for PlaceholderMicroData {
    fn density_of_states(&self, well_index: usize, local_grain_index: usize) -> f64 {
        let w = (well_index as f64) + 1.0;
        let i = local_grain_index as f64;
        // Positive, increasing ρ(E); keep variation moderate across wells.
        w * (1.0 + 0.01 * i + 1e-6 * i * i)
    }

    fn microcanonical_rate(
        &self,
        well_index: usize,
        channel_index: usize,
        local_grain_index: usize,
    ) -> f64 {
        let w = (well_index as f64) + 1.0;
        let c = (channel_index as f64) + 1.0;
        let i = local_grain_index as f64;
        // Simple increasing k(E) with well/channel scaling.
        (w * c) * (1.0e-4 * i).max(0.0)
    }
}

