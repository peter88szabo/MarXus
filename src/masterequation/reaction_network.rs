/// A single reaction channel from a given well.
/// If `connected_well_index` is Some(j), this channel is treated as an internal isomerization
/// into well j at aligned energy grains; otherwise it is an outgoing sink.
#[derive(Clone, Debug)]
pub struct ReactionChannel {
    pub name: String,
    pub connected_well_index: Option<usize>,
}

/// Parameters needed to compute collision frequency and energy transfer per well.
#[derive(Clone, Debug)]
pub struct CollisionModelParams {
    pub lennard_jones_sigma_angstrom: f64,
    pub lennard_jones_epsilon_kelvin: f64,
    pub reduced_mass_amu: f64,

    /// Energy-transfer parameter: alpha(T) = alpha_at_1000K * (T/1000)^alpha_temperature_exponent
    pub alpha_at_1000K_cm1: f64,
    pub alpha_temperature_exponent: f64,
}

/// A well with its energy grid, microcanonical data sources, and reaction channels.
#[derive(Clone, Debug)]
pub struct WellDefinition {
    pub well_name: String,

    /// Energy grain width ΔE in cm^-1
    pub energy_grain_width_cm1: f64,

    /// Lowest included local grain index (truncation)
    pub lowest_included_grain_index: usize,

    /// One past the highest included local grain index
    pub one_past_highest_included_grain_index: usize,

    /// Alignment offset in grain units (used to align absolute energy between wells)
    pub alignment_offset_in_grains: isize,

    /// Number of non-reactive grains at the bottom (k(E)=0 assumed below this)
    pub nonreactive_grain_count: usize,

    pub collision_params: CollisionModelParams,

    pub channels: Vec<ReactionChannel>,
}

/// Collision-kernel implementation choice for collisional energy transfer.
///
/// This is intended to match the single-well `CollisionKernelModel` options so that
/// users can run consistent falloff calculations in both single- and multi-well modes.
#[derive(Clone, Debug, Copy, PartialEq, Eq)]
pub enum CollisionKernelImplementation {
    /// MarXus' symmetric/SPD-style exponential-down kernel with a reactive-region
    /// probability cap (implemented as an effective collision-frequency scaling).
    Spd,

    /// MESS-style exponential-down kernel normalization (per-source normalization
    /// with an explicit ΔE=0 self-weight).
    Mess,
}

/// Global settings (temperature/pressure, banding, thresholds, etc.).
#[derive(Clone, Debug)]
pub struct MasterEquationSettings {
    /// Temperature in K
    pub temperature_kelvin: f64,
    /// Pressure in Torr
    pub pressure_torr: f64,

    /// Boltzmann constant in (cm^-1)/K
    pub boltzmann_constant_wavenumber_per_kelvin: f64,

    /// Band half-width (in grains) for collision transitions; keep small for speed.
    pub collision_band_half_width: usize,

    /// Collision-kernel implementation used for collisional energy transfer.
    pub collision_kernel_implementation: CollisionKernelImplementation,

    /// Zero outgoing microcanonical rates smaller than this threshold (s^-1)
    pub outgoing_rate_threshold: f64,

    /// If internal forward/backward rate is smaller than this, zero both directions (s^-1)
    pub internal_rate_threshold: f64,

    /// Bath-gas constant used for number density conversion.
    /// For SSUMES-like units you likely want RTCMMLC ~ (k_B in appropriate units),
    /// but here we keep it as a configurable constant.
    pub bathgas_number_density_prefactor: f64,

    /// Mean relative speed prefactor (model constant), gives v̄ = prefactor * sqrt(T / μ)
    pub mean_speed_prefactor: f64,
}

/// Trait that supplies microcanonical data ρ(E) and k_c(E).
/// Replace the placeholder implementation with your real data provider.
pub trait MicrocanonicalProvider {
    /// Density of states ρ at (well_index, local_grain_index)
    fn density_of_states(&self, well_index: usize, local_grain_index: usize) -> f64;

    /// Microcanonical rate constant k_c(E) at (well_index, channel_index, local_grain_index), in s^-1
    fn microcanonical_rate(
        &self,
        well_index: usize,
        channel_index: usize,
        local_grain_index: usize,
    ) -> f64;
}

/// Chemical activation definition: inject population into one well via a chosen channel.
#[derive(Clone, Debug)]
pub struct ChemicalActivationDefinition {
    pub activated_well_index: usize,
    pub recombination_channel_index: usize,
}

/// Results container
pub struct SolutionResults {
    pub steady_state_population: Vec<f64>,
    pub per_well_per_channel_rates: Vec<Vec<f64>>,
    pub total_outgoing_rate_constant: f64,
}
