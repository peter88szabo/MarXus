use MarXus::masterequation::api::{
    run_energy_grained_steady_state, EnergyGrainedSteadyStateInput, EnergyGrainedSteadyStateMode,
    EnergyGrainedSteadyStateOutput,
};
use MarXus::masterequation::energy_grained_me::{
    CollisionKernelModel, EnergyGrid, OlzmannMasterEquationSettings, SourceConstructionChoice,
};
use MarXus::masterequation::placeholder_microcanonical::PlaceholderMicrocanonicalData;

fn main() -> Result<(), String> {
    let grid = EnergyGrid {
        number_of_bins: 400,
        bin_width_wavenumber: 10.0,
        energy_origin_wavenumber: 0.0,
    };

    let settings = OlzmannMasterEquationSettings {
        temperature_kelvin: 298.0,
        boltzmann_constant_wavenumber_per_kelvin: 0.69503476,
        collision_frequency_per_second: 1.0e9,
        pseudo_first_order_capture_loss_per_second: 0.0,
        stepladder_base_downward_probability: 0.30,
        collision_kernel_model: CollisionKernelModel::Stepladder,
    };

    let micro = PlaceholderMicrocanonicalData { channel_count: 1 };

    let input = EnergyGrainedSteadyStateInput {
        energy_grid: grid,
        settings,
        mode: EnergyGrainedSteadyStateMode::ChemicalActivation {
            formation_flux: 1.0,
            source_choice: SourceConstructionChoice::ThermalTransitionStateSumOfStates {
                threshold_wavenumber: 500.0,
                transition_state_sum_of_states: {
                    let mut w_ts = vec![0.0; 200];
                    for (k, value) in w_ts.iter_mut().enumerate() {
                        *value = ((k + 1) as f64).powi(2);
                    }
                    w_ts
                },
            },
        },
    };

    let out = run_energy_grained_steady_state(input, &micro)?;
    match out {
        EnergyGrainedSteadyStateOutput::ChemicalActivation(r) => {
            println!("Total loss = {:.6e} s^-1", r.total_unimolecular_loss_rate_constant);
            println!("Channel k  = {:?}", r.chemically_activated_rate_constants);
        }
        _ => return Err("Unexpected run result variant.".into()),
    }

    Ok(())
}
