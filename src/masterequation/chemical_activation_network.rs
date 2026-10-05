//! Network definition of the steady-state chemical-activation (CA) master equation.
//!
//! The formulation follows Olzmann and co-workers:
//!   dN/dt = R F - J N,   J = omega (I - P) + K + k_c[D] I                     (PO14 eq. 2; O02 eqs. 5-6)
//! with N the grain populations of the energized intermediates, F the normalized nascent (source)
//! distribution, omega the collision frequency, P the collisional transition probabilities, K the
//! diagonal matrix of microcanonical rate coefficients and k_c[D] a pseudo-first-order bimolecular
//! loss of the intermediate. Several wells (isomers) are coupled by microcanonical isomerization rates
//! at equal absolute energies, as in multiwell steady-state solvers.
//!
//! References used throughout the chemical-activation files:
//!   O91  Olzmann, Gebhardt, Scherzer, Int. J. Chem. Kinet. 23, 825 (1991)
//!   O02  Olzmann, Phys. Chem. Chem. Phys. 4, 3614 (2002)
//!   GO10 Gonzalez-Garcia, Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010)
//!   PO14 Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014)
//!   T77  Troe, J. Chem. Phys. 66, 4758 (1977)
//!   R19  Robertson, Comprehensive Chemical Kinetics 43 (2019)
//!   PR03 Pilling, Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003)
//!   CD07 Carstensen, Dean, Comprehensive Chemical Kinetics 42 (2007)
//!
//! Energy grid: all wells share one grain width dE. Grain i of a well lies at E_i = i dE above that
//! well's bottom; on the common (absolute) energy scale it lies at (i + offset) dE with an integer
//! offset per well, so that grains of different wells at the same absolute energy coincide exactly.

/// Temperature dependence of the mean energy transferred in deactivating collisions,
/// <dE_down>(T) = <dE_down>(T_ref) (T/T_ref)^n.
#[derive(Debug, Clone)]
pub struct EnergyTransferParameters {
    /// <dE_down> at the reference temperature, cm-1.
    pub mean_down_at_reference_cm1: f64,
    /// Reference temperature T_ref, K (stated explicitly: input formats differ, e.g. 300 K or 1000 K).
    pub reference_temperature_kelvin: f64,
    /// Temperature exponent n.
    pub temperature_exponent: f64,
}

impl EnergyTransferParameters {
    /// <dE_down>(T) in cm-1.
    pub fn mean_down_cm1(&self, temperature_kelvin: f64) -> f64 {
        self.mean_down_at_reference_cm1
            * (temperature_kelvin / self.reference_temperature_kelvin).powf(self.temperature_exponent)
    }
}

/// Lennard-Jones parameters of the intermediate-bath gas pair (already combined:
/// sigma_AM = (sigma_A + sigma_M)/2, eps_AM = sqrt(eps_A eps_M), T77 Sec. III).
#[derive(Debug, Clone)]
pub struct LennardJonesPair {
    pub sigma_angstrom: f64,
    pub epsilon_kelvin: f64,
    pub reduced_mass_amu: f64,
}

/// Where the flux of a unimolecular channel goes.
#[derive(Debug, Clone, PartialEq)]
pub enum ChannelDestination {
    /// Bimolecular (or other) products: an irreversible exit of the network.
    Products { name: String },
    /// Another well of the network (isomerization), entered at the same absolute energy.
    Well { index: usize },
}

/// A unimolecular channel of a well with its microcanonical rate coefficients.
#[derive(Debug, Clone)]
pub struct Channel {
    pub name: String,
    pub destination: ChannelDestination,
    /// Classical reaction threshold (grain of the well). With tunneling, k(E) > 0 already below it; the
    /// absorbing barrier of the intermediate steady state refers to the classical threshold. None: the
    /// threshold is the first grain with k > 0.
    pub threshold_grain: Option<usize>,
    /// k(E_i) in s-1 for every grain i of the well (zero below the threshold).
    pub rate_constant_s_inv: Vec<f64>,
}

/// One well (isomer) of the network.
#[derive(Debug, Clone)]
pub struct Well {
    pub name: String,
    /// Absolute energy of grain 0 (the well bottom) in whole grains.
    pub bottom_offset_grains: isize,
    /// Density of states rho(E_i) in states per cm-1 for every grain i of the well.
    pub density_of_states: Vec<f64>,
    pub channels: Vec<Channel>,
    pub lennard_jones: LennardJonesPair,
    pub energy_transfer: EnergyTransferParameters,
    /// Pseudo-first-order bimolecular loss k_c[D] of the intermediate in s-1, energy independent
    /// (O02 eq. 5, GO10 eq. 6, PO14 eqs. 1-2). Zero if absent.
    pub bimolecular_sink_s_inv: f64,
}

impl Well {
    pub fn grain_count(&self) -> usize {
        self.density_of_states.len()
    }

    /// Lowest reaction threshold over all channels: the explicit classical threshold of a channel, or its
    /// first grain with k > 0.
    pub fn lowest_threshold_grain(&self) -> Option<usize> {
        self.channels
            .iter()
            .filter_map(|c| c.threshold_grain.or_else(|| c.rate_constant_s_inv.iter().position(|k| *k > 0.0)))
            .min()
    }
}

/// The chemical-activation network: wells on a common grain grid.
#[derive(Debug, Clone)]
pub struct ChemicalActivationNetwork {
    /// Grain width dE in cm-1, common to all wells.
    pub grain_width_cm1: f64,
    pub wells: Vec<Well>,
}

impl ChemicalActivationNetwork {
    /// Absolute energy (common scale) of grain `grain` of well `well`, cm-1.
    pub fn absolute_energy_cm1(&self, well: usize, grain: usize) -> f64 {
        (grain as isize + self.wells[well].bottom_offset_grains) as f64 * self.grain_width_cm1
    }

    /// Grain of well `to` at the absolute energy of grain `grain` of well `from` (may be negative or
    /// beyond the grid of `to`; the caller decides what that means).
    pub fn aligned_grain(&self, from: usize, grain: usize, to: usize) -> isize {
        grain as isize + self.wells[from].bottom_offset_grains - self.wells[to].bottom_offset_grains
    }

    /// Consistency checks of the input data.
    pub fn validate(&self) -> Result<(), String> {
        if !(self.grain_width_cm1 > 0.0) {
            return Err("The grain width must be positive.".into());
        }
        if self.wells.is_empty() {
            return Err("The network has no wells.".into());
        }
        for (w, well) in self.wells.iter().enumerate() {
            let n = well.grain_count();
            if n == 0 {
                return Err(format!("Well '{}' has an empty energy grid.", well.name));
            }
            if well.density_of_states.iter().any(|r| !(*r > 0.0) || !r.is_finite()) {
                return Err(format!(
                    "Well '{}': the density of states must be positive and finite in every grain.",
                    well.name
                ));
            }
            if !(well.bimolecular_sink_s_inv >= 0.0) {
                return Err(format!("Well '{}': the bimolecular sink must be >= 0.", well.name));
            }
            let lj = &well.lennard_jones;
            if !(lj.sigma_angstrom > 0.0 && lj.epsilon_kelvin > 0.0 && lj.reduced_mass_amu > 0.0) {
                return Err(format!("Well '{}': Lennard-Jones parameters must be positive.", well.name));
            }
            let et = &well.energy_transfer;
            if !(et.mean_down_at_reference_cm1 > 0.0 && et.reference_temperature_kelvin > 0.0) {
                return Err(format!(
                    "Well '{}': <dE_down> and its reference temperature must be positive.",
                    well.name
                ));
            }
            for channel in &well.channels {
                if channel.rate_constant_s_inv.len() != n {
                    return Err(format!(
                        "Well '{}', channel '{}': {} rate coefficients for {} grains.",
                        well.name,
                        channel.name,
                        channel.rate_constant_s_inv.len(),
                        n
                    ));
                }
                if channel.threshold_grain.map_or(false, |t| t >= n) {
                    return Err(format!(
                        "Well '{}', channel '{}': threshold grain {:?} beyond the grid ({n} grains).",
                        well.name, channel.name, channel.threshold_grain
                    ));
                }
                if channel.rate_constant_s_inv.iter().any(|k| !(*k >= 0.0) || !k.is_finite()) {
                    return Err(format!(
                        "Well '{}', channel '{}': rate coefficients must be finite and >= 0.",
                        well.name, channel.name
                    ));
                }
                if let ChannelDestination::Well { index } = channel.destination {
                    if index >= self.wells.len() || index == w {
                        return Err(format!(
                            "Well '{}', channel '{}': invalid destination well {}.",
                            well.name, channel.name, index
                        ));
                    }
                }
            }
        }
        Ok(())
    }
}

/// Collisional energy-transfer model (implemented in `collision_kernels.rs`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CollisionModel {
    /// Exponential down with exact normalization (R19 eq. 4.16); transitions up to
    /// `cutoff_in_mean_down` x <dE_down> are kept (beyond, exp(-cutoff) is neglected).
    ExponentialDown { cutoff_in_mean_down: f64 },
    /// Olzmann stepladder (O91 eqs. 13-18) with step size dE_SL = <dE_down>(T): the step
    /// "represents the average amount of energy transferred in down collisions" (GO10, text before eq. 16;
    /// eq. 16 relates it to the average over up and down collisions, <dE> = dE_SL tanh(dE_SL/(2 F_E kT))).
    Stepladder,
}

/// Position of the absorbing barrier of the intermediate steady state.
#[derive(Debug, Clone, PartialEq)]
pub enum AbsorbingBarrier {
    /// `kt_multiple` k_BT below the lowest reaction threshold of each well; the literature places the
    /// absorbing boundary about 10 k_BT below the reaction threshold (PR03; CD07 p. 125), the default.
    /// For wells that are shallow compared with 10 k_BT plus their thermal width the user may choose a
    /// smaller distance; the stabilization then depends on that choice.
    BelowLowestThreshold { kt_multiple: f64 },
    /// Explicit barrier grain for every well (grains below it are absorbing).
    AtGrains(Vec<usize>),
}

impl Default for AbsorbingBarrier {
    fn default() -> Self {
        AbsorbingBarrier::BelowLowestThreshold { kt_multiple: 10.0 }
    }
}

/// Which steady state of the chemical-activation problem is solved.
#[derive(Debug, Clone, PartialEq)]
pub enum SteadyState {
    /// Final (asymptotic) steady state: no absorbing barrier, the stabilized intermediates remain in
    /// the population and react thermally; "there is no more net stabilization" (GO10 p. 12295;
    /// PO14 p. 238).
    Final,
    /// Intermediate steady state: grains below an absorbing barrier are removed and the flux into
    /// them is the stabilization; "implemented by introducing a lower absorbing barrier into the
    /// master equation" (GO10 p. 12295; O02 p. 3616).
    Intermediate { barrier: AbsorbingBarrier },
    /// Rate coefficients from the eigenvalues and eigenvectors of J instead of a steady state (for wells
    /// that are shallow compared with the absorbing-barrier distance plus their thermal width). Not
    /// available yet: selecting it is reported as an error, never replaced by another method.
    EigenvalueAnalysis,
}

/// Temperature and bath-gas pressure.
#[derive(Debug, Clone, Copy)]
pub struct Conditions {
    pub temperature_kelvin: f64,
    pub pressure_torr: f64,
}

/// Model options of a chemical-activation calculation.
#[derive(Debug, Clone)]
pub struct ChemicalActivationOptions {
    pub collision_model: CollisionModel,
    pub steady_state: SteadyState,
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Well with a smooth density of states and one product channel opening at `threshold`.
    pub(crate) fn test_well(name: &str, grains: usize, offset: isize, threshold: usize) -> Well {
        Well {
            name: name.into(),
            bottom_offset_grains: offset,
            density_of_states: (0..grains).map(|i| (1.0 + 0.05 * i as f64).powi(8)).collect(),
            channels: vec![Channel {
                name: format!("{name}-products"),
                destination: ChannelDestination::Products { name: "P".into() },
                threshold_grain: None,
                rate_constant_s_inv: (0..grains)
                    .map(|i| if i >= threshold { 1.0e6 * ((i - threshold) as f64 + 1.0) } else { 0.0 })
                    .collect(),
            }],
            lennard_jones: LennardJonesPair { sigma_angstrom: 4.5, epsilon_kelvin: 300.0, reduced_mass_amu: 20.0 },
            energy_transfer: EnergyTransferParameters {
                mean_down_at_reference_cm1: 200.0,
                reference_temperature_kelvin: 300.0,
                temperature_exponent: 0.85,
            },
            bimolecular_sink_s_inv: 0.0,
        }
    }

    #[test]
    fn aligned_grain_maps_equal_absolute_energies() {
        let network = ChemicalActivationNetwork {
            grain_width_cm1: 10.0,
            wells: vec![test_well("A", 300, 0, 200), test_well("B", 350, -50, 250)],
        };
        // Grain 100 of A lies at 1000 cm-1; well B's bottom is 500 cm-1 lower, so that is grain 150 of B.
        assert_eq!(network.aligned_grain(0, 100, 1), 150);
        assert_eq!(network.aligned_grain(1, 150, 0), 100);
        assert_eq!(network.absolute_energy_cm1(1, 150), 1000.0);
    }

    #[test]
    fn lowest_threshold_is_the_first_grain_with_an_open_channel() {
        let mut well = test_well("A", 300, 0, 200);
        well.channels.push(Channel {
            name: "A-low".into(),
            destination: ChannelDestination::Products { name: "Q".into() },
            threshold_grain: None,
            rate_constant_s_inv: (0..300).map(|i| if i >= 180 { 1.0 } else { 0.0 }).collect(),
        });
        assert_eq!(well.lowest_threshold_grain(), Some(180));
    }

    #[test]
    fn an_explicit_classical_threshold_replaces_the_first_open_grain() {
        // Tunneling opens the channel at grain 200, its classical threshold is grain 250.
        let mut well = test_well("A", 300, 0, 200);
        well.channels[0].threshold_grain = Some(250);
        assert_eq!(well.lowest_threshold_grain(), Some(250));
        // A second channel without explicit threshold, open from grain 230, is lower.
        well.channels.push(Channel {
            name: "A-second".into(),
            destination: ChannelDestination::Products { name: "Q".into() },
            threshold_grain: None,
            rate_constant_s_inv: (0..300).map(|i| if i >= 230 { 1.0 } else { 0.0 }).collect(),
        });
        assert_eq!(well.lowest_threshold_grain(), Some(230));
        let mut network = ChemicalActivationNetwork { grain_width_cm1: 10.0, wells: vec![well] };
        network.wells[0].channels[0].threshold_grain = Some(300);
        assert!(network.validate().is_err(), "threshold beyond the grid");
    }

    #[test]
    fn mean_down_follows_the_power_law_from_its_reference_temperature() {
        let et = EnergyTransferParameters {
            mean_down_at_reference_cm1: 200.0,
            reference_temperature_kelvin: 300.0,
            temperature_exponent: 0.85,
        };
        assert!((et.mean_down_cm1(300.0) - 200.0).abs() < 1e-12);
        assert!((et.mean_down_cm1(600.0) - 200.0 * 2.0_f64.powf(0.85)).abs() < 1e-12);
    }

    #[test]
    fn validation_rejects_self_isomerization_and_rate_arrays_of_wrong_length() {
        let mut network = ChemicalActivationNetwork { grain_width_cm1: 10.0, wells: vec![test_well("A", 300, 0, 200)] };
        assert!(network.validate().is_ok());
        network.wells[0].channels[0].destination = ChannelDestination::Well { index: 0 };
        assert!(network.validate().is_err());
        network.wells[0].channels[0].destination = ChannelDestination::Products { name: "P".into() };
        network.wells[0].channels[0].rate_constant_s_inv.pop();
        assert!(network.validate().is_err());
    }
}
