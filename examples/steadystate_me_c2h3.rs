mod cases;

use MarXus::masterequation::api::{
    run_master_equation, EnergyGrainedSteadyStateInput, EnergyGrainedSteadyStateMode,
    EnergyGrainedSteadyStateOutput,
};
use MarXus::masterequation::collisional_energy_transfer::{
    derive_steady_state_collision_parameters, ExponentialDownParameters, LennardJonesParameters,
};
use MarXus::masterequation::energy_grained_me::{
    CollisionKernelModel, EnergyGrid, OlzmannMasterEquationSettings, SourceConstructionChoice,
};
use MarXus::masterequation::high_pressure_limit::{
    eckart_tunneling_kappa, high_pressure_tst_with_thermo, EckartTunnelingInput,
    ReactionMolecularity,
};
use MarXus::masterequation::report::{
    print_high_pressure_limit, print_reaction_header, print_source_threshold,
    print_species_species_table, print_steady_state_summary, ReactionEnergetics,
};
use MarXus::masterequation::text_input::{
    get_f64, get_f64_or, get_string, get_string_or, get_usize, get_vec_f64,
    parse_sectioned_kv,
};

fn main() -> Result<(), String> {
    let cfg = parse_sectioned_kv(include_str!("steadystate_me_c2h3_input.txt"))?;

    let temperature_kelvin = get_f64(&cfg, "conditions", "temperature_kelvin")?;
    let pressure_atm = get_f64(&cfg, "conditions", "pressure_atm")?;

    let grid = EnergyGrid {
        number_of_bins: get_usize(&cfg, "parameters", "number_of_bins")?,
        bin_width_wavenumber: get_f64(&cfg, "parameters", "bin_width_wavenumber")?,
        energy_origin_wavenumber: get_f64(&cfg, "parameters", "energy_origin_wavenumber")?,
    };

    let exp = ExponentialDownParameters {
        factor_cm1: get_f64(&cfg, "collision", "energy_relax_factor_cm1")?,
        power: get_f64(&cfg, "collision", "energy_relax_power")?,
        exponent_cutoff: get_f64(&cfg, "collision", "energy_relax_exponent_cutoff")?,
    };
    let eps = get_vec_f64(&cfg, "collision", "lennard_jones_epsilons_cm1")?;
    let sig = get_vec_f64(&cfg, "collision", "lennard_jones_sigmas_angstrom")?;
    let masses = get_vec_f64(&cfg, "collision", "lennard_jones_masses_amu")?;
    if eps.len() != 2 || sig.len() != 2 || masses.len() != 2 {
        return Err("collision LJ epsilon/sigma/mass lists must each have exactly 2 values.".into());
    }
    let lj = LennardJonesParameters {
        epsilon_1_cm1: eps[0],
        epsilon_2_cm1: eps[1],
        sigma_1_angstrom: sig[0],
        sigma_2_angstrom: sig[1],
        mass_1_amu: masses[0],
        mass_2_amu: masses[1],
    };
    let derived = derive_steady_state_collision_parameters(
        temperature_kelvin,
        pressure_atm,
        grid.bin_width_wavenumber,
        &exp,
        &lj,
    )?;
    println!(
        "Collision model exponential_down_lj_omega (banded): omega={:.6e} s^-1, p_down={:.6e}, <DeltaE_down>={:.3} cm^-1",
        derived.collision_frequency_per_second,
        derived.stepladder_base_downward_probability,
        derived.average_downstep_cm1
    );

    let base_settings = OlzmannMasterEquationSettings {
        temperature_kelvin: get_f64(&cfg, "parameters", "temperature_kelvin")?,
        boltzmann_constant_wavenumber_per_kelvin: get_f64(
            &cfg,
            "parameters",
            "boltzmann_constant_wavenumber_per_kelvin",
        )?,
        collision_frequency_per_second: derived.collision_frequency_per_second,
        pseudo_first_order_capture_loss_per_second: get_f64(
            &cfg,
            "parameters",
            "pseudo_first_order_capture_loss_per_second",
        )?,
        stepladder_base_downward_probability: derived.stepladder_base_downward_probability,
        // NOTE: this example intentionally runs *both* kernel implementations (SPD and MESS-like)
        // back-to-back for comparison, so the actual kernel choice is set per run below.
        collision_kernel_model: CollisionKernelModel::Stepladder,
    };

    let e = ReactionEnergetics {
        well_name: get_string(&cfg, "species", "well_name")?,
        product_name: get_string(&cfg, "species", "product_name")?,
        barrier_name: get_string(&cfg, "species", "barrier_name")?,
        barrier_label: get_string(&cfg, "species", "barrier_label")?,
        well_dh0_kcal_mol: get_f64(&cfg, "species", "well_dh0_kcal_mol")?,
        product_dh0_kcal_mol: get_f64_or(&cfg, "species", "product_dh0_kcal_mol", 0.0)?,
        barrier_dh0_kcal_mol: get_f64(&cfg, "species", "barrier_dh0_kcal_mol")?,
    };

    print_reaction_header(temperature_kelvin, pressure_atm, &e);

    let tunneling_model = get_string_or(&cfg, "tunneling", "model", "none");
    let tunneling_kappa = if tunneling_model == "eckart" {
        eckart_tunneling_kappa(EckartTunnelingInput {
            temperature_kelvin,
            imaginary_frequency_cm1: get_f64(&cfg, "tunneling", "imaginary_frequency_cm1")?,
            forward_barrier_kcal_mol: get_f64_or(
                &cfg,
                "tunneling",
                "forward_barrier_kcal_mol",
                e.ea_forward_kcal_mol(),
            )?,
            reverse_barrier_kcal_mol: get_f64_or(
                &cfg,
                "tunneling",
                "reverse_barrier_kcal_mol",
                e.ea_reverse_kcal_mol(),
            )?,
            integration_step_kcal_mol: get_f64_or(&cfg, "tunneling", "integration_step_kcal_mol", 0.05)?,
            integration_max_kcal_mol: get_f64_or(&cfg, "tunneling", "integration_max_kcal_mol", 120.0)?,
        })?
    } else {
        get_f64_or(&cfg, "high_pressure", "tunneling_kappa", 1.0)?
    };

    let freq_cutoff_cm1 = get_f64_or(&cfg, "high_pressure", "freq_cutoff_cm1", 100.0)?;
    let (mut w1_mol, mut ts_mol) = cases::c2h3_molecules::build_c2h3_w1_and_ts_molecules(
        e.well_dh0_kcal_mol,
        e.barrier_dh0_kcal_mol,
    );
    let hp_forward = high_pressure_tst_with_thermo(
        &mut w1_mol,
        None,
        &mut ts_mol,
        temperature_kelvin,
        pressure_atm,
        freq_cutoff_cm1,
        tunneling_kappa,
        ReactionMolecularity::Unimolecular,
    )?;
    let (mut c2h2_mol, mut h_mol) =
        cases::c2h3_molecules::build_p1_reactant_molecules(e.product_dh0_kcal_mol);
    let hp_reverse = high_pressure_tst_with_thermo(
        &mut c2h2_mol,
        Some(&mut h_mol),
        &mut ts_mol,
        temperature_kelvin,
        pressure_atm,
        freq_cutoff_cm1,
        tunneling_kappa,
        ReactionMolecularity::Bimolecular,
    )?;

    print_high_pressure_limit(
        &e,
        &hp_forward,
        &hp_reverse,
        &tunneling_model,
        tunneling_kappa,
    );

    let provider = get_string(&cfg, "microcanonical", "provider")?;
    if provider != "rrkm_c2h3" {
        return Err("This example expects microcanonical.provider=rrkm_c2h3".into());
    }
    let micro =
        cases::c2h3_rrkm::build_c2h3_rrkm_microcanonical_data(&grid, e.ea_forward_kcal_mol())?;

    let formation_flux = get_f64(&cfg, "parameters", "formation_flux")?;
    let threshold_wavenumber = get_f64(&cfg, "source_thermal_ts", "threshold_wavenumber")?;
    let source = SourceConstructionChoice::ThermalTransitionStateSumOfStates {
        threshold_wavenumber,
        transition_state_sum_of_states: cases::c2h3_rrkm::build_c2h3_rrkm_ts_sum_of_states(&grid),
    };

    print_source_threshold(threshold_wavenumber);

    // -------------------------------------------------------------------------
    // 1) Run with MarXus SPD/symmetric kernel (with reactive-region probability cap)
    // -------------------------------------------------------------------------
    let mut settings_spd = base_settings.clone();
    settings_spd.collision_kernel_model = CollisionKernelModel::ExponentialBanded {
        mean_downstep_wavenumber: derived.average_downstep_cm1,
        exponent_cutoff: exp.exponent_cutoff,
    };

    let out_spd = run_master_equation(
        EnergyGrainedSteadyStateInput {
            energy_grid: grid.clone(),
            settings: settings_spd,
            mode: EnergyGrainedSteadyStateMode::ChemicalActivation {
                formation_flux,
                source_choice: source.clone(),
            },
        },
        &micro,
    )?;
    let res_spd = match out_spd {
        EnergyGrainedSteadyStateOutput::ChemicalActivation(r) => r,
        _ => return Err("Unexpected run result variant (spd).".into()),
    };

    // -------------------------------------------------------------------------
    // 2) Run with MESS-like normalized banded kernel (row-normalized per energy)
    // -------------------------------------------------------------------------
    let mut settings_mess = base_settings.clone();
    settings_mess.collision_kernel_model = CollisionKernelModel::ExponentialBandedMess {
        mean_downstep_wavenumber: derived.average_downstep_cm1,
        exponent_cutoff: exp.exponent_cutoff,
    };

    let out_mess = run_master_equation(
        EnergyGrainedSteadyStateInput {
            energy_grid: grid,
            settings: settings_mess,
            mode: EnergyGrainedSteadyStateMode::ChemicalActivation {
                formation_flux,
                source_choice: source,
            },
        },
        &micro,
    )?;
    let res_mess = match out_mess {
        EnergyGrainedSteadyStateOutput::ChemicalActivation(r) => r,
        _ => return Err("Unexpected run result variant (mess).".into()),
    };

    // -------------------------------------------------------------------------
    // Printing: show both results next to each other, AND print the full
    // "Species-Species Rate Tables" block for BOTH kernel runs (SPD + MESS-like).
    // -------------------------------------------------------------------------
    println!();
    println!("Master-equation steady-state results (two collision-kernel implementations):");

    let k_spd = res_spd
        .chemically_activated_rate_constants
        .first()
        .copied()
        .unwrap_or(0.0);
    let k_mess = res_mess
        .chemically_activated_rate_constants
        .first()
        .copied()
        .unwrap_or(0.0);

    let red_spd = if hp_forward.rate_constant > 0.0 {
        k_spd / hp_forward.rate_constant
    } else {
        0.0
    };
    let red_mess = if hp_forward.rate_constant > 0.0 {
        k_mess / hp_forward.rate_constant
    } else {
        0.0
    };

    let krev_spd = hp_reverse.rate_constant * red_spd;
    let krev_mess = hp_reverse.rate_constant * red_mess;

    println!(
        "{:<20} {:>18} {:>18}",
        "Kernel",
        "SPD (MarXus)",
        "MESS-like"
    );
    println!(
        "{:<20} {:>18.6e} {:>18.6e}",
        "Total loss (s^-1)",
        res_spd.total_unimolecular_loss_rate_constant,
        res_mess.total_unimolecular_loss_rate_constant
    );
    println!(
        "{:<20} {:>18.6e} {:>18.6e}",
        "W1->P1 @1 atm",
        k_spd,
        k_mess
    );
    println!(
        "{:<20} {:>18.6e} {:>18.6e}",
        "P1->W1 @1 atm",
        krev_spd,
        krev_mess
    );

    // -------------------------------------------------------------------------
    // Comparison against MESS (NO TUNNELING) reference output.
    //
    // IMPORTANT: these reference numbers are only meaningful when tunneling is disabled
    // (tunneling.model = none). They correspond to the user-provided MESS printout:
    //
    // Temperature = 1e+03 K, Pressure = 1 atm
    // High Pressure:
    //   W1->P1 = 2.41e+05 s^-1
    //   P1->W1 = 3.88e-11 cm^3 molecule^-1 s^-1
    // 1 atm:
    //   W1->P1 = 1.46e+04 s^-1
    //   P1->W1 = 2.34e-12 cm^3 molecule^-1 s^-1
    // -------------------------------------------------------------------------
    if tunneling_model != "none" {
        println!();
        println!(
            "NOTE: MESS %error comparison is skipped because tunneling.model='{}'. Set tunneling.model='none' to compare.",
            tunneling_model
        );
    } else {
        let mess_hp_w1_p1 = 2.41e5_f64;
        let mess_hp_p1_w1 = 3.88e-11_f64;
        let mess_1atm_w1_p1 = 1.46e4_f64;
        let mess_1atm_p1_w1 = 2.34e-12_f64;

        let perr = |our: f64, reference: f64| -> f64 {
            if reference == 0.0 {
                0.0
            } else {
                100.0 * (our - reference) / reference
            }
        };

        println!();
        println!("Comparison vs MESS reference (no tunneling):");
        println!(
            "{:<28} {:>14} {:>14} {:>10} {:>14} {:>10}",
            "Quantity", "MESS", "SPD", "%err", "MESS-like", "%err"
        );
        println!(
            "{:<28} {:>14.6e} {:>14.6e} {:>10.3} {:>14.6e} {:>10.3}",
            "k_inf  W1->P1 (s^-1)",
            mess_hp_w1_p1,
            hp_forward.rate_constant,
            perr(hp_forward.rate_constant, mess_hp_w1_p1),
            hp_forward.rate_constant,
            perr(hp_forward.rate_constant, mess_hp_w1_p1),
        );
        println!(
            "{:<28} {:>14.6e} {:>14.6e} {:>10.3} {:>14.6e} {:>10.3}",
            "k_inf  P1->W1 (cm^3)",
            mess_hp_p1_w1,
            hp_reverse.rate_constant,
            perr(hp_reverse.rate_constant, mess_hp_p1_w1),
            hp_reverse.rate_constant,
            perr(hp_reverse.rate_constant, mess_hp_p1_w1),
        );
        println!(
            "{:<28} {:>14.6e} {:>14.6e} {:>10.3} {:>14.6e} {:>10.3}",
            "k(1 atm) W1->P1 (s^-1)",
            mess_1atm_w1_p1,
            k_spd,
            perr(k_spd, mess_1atm_w1_p1),
            k_mess,
            perr(k_mess, mess_1atm_w1_p1),
        );
        println!(
            "{:<28} {:>14.6e} {:>14.6e} {:>10.3} {:>14.6e} {:>10.3}",
            "k(1 atm) P1->W1 (cm^3)",
            mess_1atm_p1_w1,
            krev_spd,
            perr(krev_spd, mess_1atm_p1_w1),
            krev_mess,
            perr(krev_mess, mess_1atm_p1_w1),
        );
    }

    // -------------------------------------------------------------------------
    // Full tables: SPD kernel run
    // -------------------------------------------------------------------------
    println!();
    println!("================================================================================");
    println!("FULL OUTPUT: SPD (MarXus) collision kernel");
    println!("================================================================================");
    print_steady_state_summary(
        res_spd.total_unimolecular_loss_rate_constant,
        &res_spd.chemically_activated_rate_constants,
    );
    print_species_species_table(
        temperature_kelvin,
        pressure_atm,
        &e,
        hp_forward.rate_constant,
        hp_reverse.rate_constant,
        k_spd,
        krev_spd,
    );

    // -------------------------------------------------------------------------
    // Full tables: MESS-like kernel run
    // -------------------------------------------------------------------------
    println!();
    println!("================================================================================");
    println!("FULL OUTPUT: MESS-like normalized collision kernel");
    println!("================================================================================");
    print_steady_state_summary(
        res_mess.total_unimolecular_loss_rate_constant,
        &res_mess.chemically_activated_rate_constants,
    );
    print_species_species_table(
        temperature_kelvin,
        pressure_atm,
        &e,
        hp_forward.rate_constant,
        hp_reverse.rate_constant,
        k_mess,
        krev_mess,
    );

    Ok(())
}
