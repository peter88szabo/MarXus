use super::energy_grained_me::MicrocanonicalGridData;

/// Minimal placeholder microcanonical provider for energy-grained examples.
///
/// This is not a physical model; it only exists to make examples runnable.
pub struct PlaceholderMicrocanonicalData {
    pub channel_count: usize,
}

impl MicrocanonicalGridData for PlaceholderMicrocanonicalData {
    fn density_of_states(&self, energy_bin_index: usize) -> f64 {
        let i = energy_bin_index as f64;
        // Positive, mildly increasing ρ(E) to avoid pathological ratios.
        1.0 + 0.01 * i + 1e-6 * i * i
    }

    fn microcanonical_rate_for_channel(
        &self,
        channel_index: usize,
        energy_bin_index: usize,
    ) -> f64 {
        if channel_index >= self.channel_count {
            return 0.0;
        }
        let i = energy_bin_index as f64;
        // Simple increasing loss rate with energy; channel-dependent scaling.
        let scale = 1.0 + (channel_index as f64);
        scale * (1.0e-4 * i).max(0.0)
    }

    fn number_of_unimolecular_channels(&self) -> usize {
        self.channel_count
    }
}

