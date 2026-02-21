/// Trait that supplies microcanonical data on the discretized energy grid.
pub trait MicrocanonicalGridData {
    /// Density of states ρ(E_i) at energy-bin index i.
    fn density_of_states(&self, energy_bin_index: usize) -> f64;

    /// Microcanonical unimolecular rate constant for a given channel at bin i, in s^-1.
    fn microcanonical_rate_for_channel(&self, channel_index: usize, energy_bin_index: usize)
        -> f64;

    /// Number of unimolecular channels included in k_loss(E) = sum over channels.
    fn number_of_unimolecular_channels(&self) -> usize;
}

/// Energy grid definition: uniform bins.
#[derive(Clone, Debug)]
pub struct EnergyGrid {
    /// Number of energy bins.
    pub number_of_bins: usize,

    /// Bin width ΔE in cm^-1.
    pub bin_width_wavenumber: f64,

    /// Optional bin-zero energy offset in cm^-1.
    pub energy_origin_wavenumber: f64,
}

impl EnergyGrid {
    /// Energy at bin value E_i (cm^-1), with i=0..n-1:
    /// E_i = E0 + i * ΔE
    pub fn energy_at_bin(&self, energy_bin_index: usize) -> f64 {
        self.energy_origin_wavenumber + (energy_bin_index as f64) * self.bin_width_wavenumber
    }

    /// Convert a wavenumber energy shift (cm^-1) into an integer bin shift.
    pub fn energy_shift_to_bin_shift(&self, energy_shift_wavenumber: f64) -> Result<isize, String> {
        if self.bin_width_wavenumber <= 0.0 {
            return Err("Energy grid bin width must be positive.".into());
        }

        let exact = energy_shift_wavenumber / self.bin_width_wavenumber;
        let rounded = exact.round();
        let mismatch = (exact - rounded).abs() * self.bin_width_wavenumber;

        // Many input energies (barriers, threshold energies, reaction exothermicities)
        // are not exactly commensurate with the chosen ΔE. In that case we map to the
        // nearest bin shift (round) as long as the mismatch is <= ΔE/2 (nearest-bin rule).
        //
        // This is important for matching MESS-style decks, where E0 and ΔE are specified
        // independently (and are rarely exact multiples).
        let tol = 1e-8 * self.bin_width_wavenumber;
        let half_bin = 0.5 * self.bin_width_wavenumber + tol;
        if mismatch > half_bin {
            return Err(format!(
                "Energy shift {} cm^-1 is too far from the nearest bin multiple of ΔE={} cm^-1 (mismatch {} cm^-1).",
                energy_shift_wavenumber, self.bin_width_wavenumber, mismatch
            ));
        }

        Ok(rounded as isize)
    }
}

/// Physical + numerical settings for the Olzmann-style ME.
#[derive(Clone, Debug)]
pub struct OlzmannMasterEquationSettings {
    /// Temperature (K)
    pub temperature_kelvin: f64,

    /// Boltzmann constant in (cm^-1)/K.
    pub boltzmann_constant_wavenumber_per_kelvin: f64,

    /// Collision frequency ω (s^-1).
    pub collision_frequency_per_second: f64,

    /// Pseudo-first-order capture loss rate k_capture (s^-1), e.g. k_c [D].
    pub pseudo_first_order_capture_loss_per_second: f64,

    /// Base downward stepladder probability per collision.
    pub stepladder_base_downward_probability: f64,

    /// Collision-kernel model used to assemble the collisional part of J.
    pub collision_kernel_model: CollisionKernelModel,
}

#[derive(Clone, Debug)]
pub enum CollisionKernelModel {
    /// Original nearest-neighbor stepladder operator.
    Stepladder,

    /// Banded exponential-down kernel assembled in the symmetric/SPD style used by MarXus.
    ///
    /// This kernel enforces detailed balance in the band weights and then uses a global
    /// collision-frequency scaling safeguard in the reactive region so that implied
    /// per-collision leave probabilities do not exceed 1.
    ExponentialBanded {
        /// Mean downward energy transfer <DeltaE_down> in cm^-1.
        mean_downstep_wavenumber: f64,

        /// Exponent cutoff from MESS-like parameters, used to limit bandwidth.
        exponent_cutoff: f64,
    },

    /// Banded exponential-down kernel assembled in a MESS-like way:
    ///
    /// - Construct unnormalized band weights using an exponential-down shape and
    ///   an equilibrium correction for upward transitions.
    /// - Row-normalize per source grain (including a self-weight at ΔE=0) to form
    ///   a per-collision probability kernel.
    /// - Use ω as the physical collision frequency with no additional global cap.
    ///
    /// This is provided as an optional alternative to `ExponentialBanded` so users
    /// can select which normalization convention they want from an input deck.
    ExponentialBandedMess {
        /// Mean downward energy transfer <DeltaE_down> in cm^-1.
        mean_downstep_wavenumber: f64,

        /// Exponent cutoff from MESS-like parameters, used to limit bandwidth.
        exponent_cutoff: f64,
    },
}

/// Chemical-activation source construction variants.
#[derive(Clone, Debug)]
pub enum SourceConstructionChoice {
    /// Thermal reactants using TS sum-of-states and threshold E0.
    ThermalTransitionStateSumOfStates {
        threshold_wavenumber: f64,
        transition_state_sum_of_states: Vec<f64>,
    },

    /// Nonthermal reactants via shifted convolution.
    NonthermalShiftedConvolution {
        threshold_wavenumber: f64,
        reactant_a_distribution: Vec<f64>,
        reactant_b_distribution: Vec<f64>,
    },

    /// Shift approximation for consecutive activation.
    ShiftApproximation {
        reactant_distribution: Vec<f64>,
        reaction_energy_wavenumber: f64,
        partner_mean_energy_wavenumber: f64,
    },
}

/// Outputs from a steady-state CA solve.
pub struct SteadyStateChemicalActivationResult {
    /// Steady-state bin populations N_ss (unnormalized, scales with formation flux).
    pub steady_state_bin_populations: Vec<f64>,

    /// Normalized steady-state distribution over bins.
    pub normalized_steady_state_distribution: Vec<f64>,

    /// Channel-resolved chemically activated rate constants.
    pub chemically_activated_rate_constants: Vec<f64>,

    /// Total unimolecular loss rate averaged over the steady-state distribution.
    pub total_unimolecular_loss_rate_constant: f64,
}

/// Outputs from thermal (no-source) unimolecular dissociation in a single well.
pub struct ThermalDissociationResult {
    /// Canonical Boltzmann-like population over bins:
    /// n_i ∝ ρ(E_i) exp(-E_i / k_B T), normalized to sum 1.
    pub normalized_thermal_distribution: Vec<f64>,

    /// Channel-resolved thermal dissociation rate constants:
    /// k_th(channel) = sum_i n_i * k_channel(E_i).
    pub thermal_rate_constants: Vec<f64>,

    /// Total unimolecular thermal loss rate:
    /// k_loss_total = sum_i n_i * sum_channels k_channel(E_i).
    pub total_unimolecular_loss_rate_constant: f64,
}
