#[derive(Clone, Debug)]
pub struct ExponentialDownParameters {
    pub factor_cm1: f64,
    pub power: f64,
    pub exponent_cutoff: f64,
}

#[derive(Clone, Debug)]
pub struct LennardJonesParameters {
    pub epsilon_1_cm1: f64,
    pub epsilon_2_cm1: f64,
    pub sigma_1_angstrom: f64,
    pub sigma_2_angstrom: f64,
    pub mass_1_amu: f64,
    pub mass_2_amu: f64,
}

#[derive(Clone, Debug)]
pub struct DerivedSteadyStateCollisionParameters {
    pub collision_frequency_per_second: f64,
    pub stepladder_base_downward_probability: f64,
    pub average_downstep_cm1: f64,
}

pub fn average_downstep_cm1(
    temperature_kelvin: f64,
    exp_down: &ExponentialDownParameters,
) -> Result<f64, String> {
    if temperature_kelvin <= 0.0 {
        return Err("temperature_kelvin must be > 0.".into());
    }
    if exp_down.factor_cm1 <= 0.0 {
        return Err("EnergyRelaxation.Exponential Factor must be > 0.".into());
    }
    Ok(exp_down.factor_cm1 * (temperature_kelvin / 300.0).powf(exp_down.power))
}

pub fn derive_stepladder_downward_probability(
    temperature_kelvin: f64,
    energy_bin_width_cm1: f64,
    exp_down: &ExponentialDownParameters,
) -> Result<f64, String> {
    if temperature_kelvin <= 0.0 {
        return Err("temperature_kelvin must be > 0.".into());
    }
    if energy_bin_width_cm1 <= 0.0 {
        return Err("energy_bin_width_cm1 must be > 0.".into());
    }
    if exp_down.factor_cm1 <= 0.0 {
        return Err("EnergyRelaxation.Exponential Factor must be > 0.".into());
    }

    let average_downstep_cm1 = average_downstep_cm1(temperature_kelvin, exp_down)?;
    let exponent = (energy_bin_width_cm1 / average_downstep_cm1)
        .max(0.0)
        .min(exp_down.exponent_cutoff.max(0.0));

    Ok((1.0 - (-exponent).exp()).clamp(1e-12, 1.0))
}

pub fn derive_steady_state_collision_parameters(
    temperature_kelvin: f64,
    pressure_atm: f64,
    energy_bin_width_cm1: f64,
    exp_down: &ExponentialDownParameters,
    lj: &LennardJonesParameters,
) -> Result<DerivedSteadyStateCollisionParameters, String> {
    if temperature_kelvin <= 0.0 {
        return Err("temperature_kelvin must be > 0.".into());
    }
    if pressure_atm <= 0.0 {
        return Err("pressure_atm must be > 0.".into());
    }
    if energy_bin_width_cm1 <= 0.0 {
        return Err("energy_bin_width_cm1 must be > 0.".into());
    }
    if exp_down.factor_cm1 <= 0.0 {
        return Err("EnergyRelaxation.Exponential Factor must be > 0.".into());
    }

    let average_downstep_cm1 = average_downstep_cm1(temperature_kelvin, exp_down)?;
    let stepladder_base_downward_probability =
        derive_stepladder_downward_probability(temperature_kelvin, energy_bin_width_cm1, exp_down)?;

    let collision_frequency_per_second =
        estimate_lj_collision_frequency_per_second(temperature_kelvin, pressure_atm, lj)?;

    Ok(DerivedSteadyStateCollisionParameters {
        collision_frequency_per_second,
        stepladder_base_downward_probability,
        average_downstep_cm1,
    })
}

pub fn estimate_lj_collision_frequency_per_second(
    temperature_kelvin: f64,
    pressure_atm: f64,
    lj: &LennardJonesParameters,
) -> Result<f64, String> {
    if lj.sigma_1_angstrom <= 0.0
        || lj.sigma_2_angstrom <= 0.0
        || lj.mass_1_amu <= 0.0
        || lj.mass_2_amu <= 0.0
    {
        return Err("Lennard-Jones sigmas and masses must be > 0.".into());
    }

    const ATM_TO_PA: f64 = 101_325.0;
    const KB_J_PER_K: f64 = 1.380_649e-23;
    const AMU_TO_KG: f64 = 1.660_539_066_60e-27;
    const CM1_TO_K: f64 = 1.438_776_877;

    let sigma_m = 0.5 * (lj.sigma_1_angstrom + lj.sigma_2_angstrom) * 1.0e-10;
    let epsilon_cm1 = (lj.epsilon_1_cm1 * lj.epsilon_2_cm1).sqrt().max(1e-12);
    let epsilon_over_k_kelvin = epsilon_cm1 * CM1_TO_K;
    let t_star = (temperature_kelvin / epsilon_over_k_kelvin).max(1e-6);

    let omega_22 = 1.16145 * t_star.powf(-0.14874)
        + 0.52487 * (-0.77320 * t_star).exp()
        + 2.16178 * (-2.437887 * t_star).exp();

    let m1 = lj.mass_1_amu * AMU_TO_KG;
    let m2 = lj.mass_2_amu * AMU_TO_KG;
    let reduced_mass = (m1 * m2) / (m1 + m2);

    let relative_speed =
        (8.0 * KB_J_PER_K * temperature_kelvin / (std::f64::consts::PI * reduced_mass)).sqrt();
    let collision_rate_constant_m3_s =
        std::f64::consts::PI * sigma_m * sigma_m * relative_speed * omega_22;

    let pressure_pa = pressure_atm * ATM_TO_PA;
    let number_density_m3 = pressure_pa / (KB_J_PER_K * temperature_kelvin);

    Ok(collision_rate_constant_m3_s * number_density_m3)
}
