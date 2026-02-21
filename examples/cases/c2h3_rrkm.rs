use MarXus::inertia::inertia::get_brot;
use MarXus::masterequation::energy_grained_me::{EnergyGrid, MicrocanonicalGridData};
use MarXus::rrkm::rrkm_rate::get_kE;
use MarXus::rrkm::sum_and_density::get_rovib_WE_or_rhoE;

pub struct C2h3RrkmMicrocanonicalData {
    rho_e_w1: Vec<f64>,
    k_e_w1_to_p1: Vec<f64>,
}

impl MicrocanonicalGridData for C2h3RrkmMicrocanonicalData {
    fn density_of_states(&self, energy_bin_index: usize) -> f64 {
        self.rho_e_w1[energy_bin_index]
    }

    fn microcanonical_rate_for_channel(
        &self,
        channel_index: usize,
        energy_bin_index: usize,
    ) -> f64 {
        if channel_index == 0 {
            self.k_e_w1_to_p1[energy_bin_index]
        } else {
            0.0
        }
    }

    fn number_of_unimolecular_channels(&self) -> usize {
        1
    }
}

pub fn build_c2h3_rrkm_microcanonical_data(
    grid: &EnergyGrid,
    ea_forward_kcal_mol: f64,
) -> Result<C2h3RrkmMicrocanonicalData, String> {
    if grid.number_of_bins < 2 {
        return Err("number_of_bins must be >= 2 for RRKM provider.".into());
    }

    let d_e = grid.bin_width_wavenumber;
    let nebin = grid.number_of_bins - 1;
    const KCAL_PER_CM1: f64 = 2.85914e-3;
    const MASS_C: f64 = 12.0;
    const MASS_H: f64 = 1.0;

    let w1_xyz = vec![
        [0.0, 0.0229607853, -0.6194779279],
        [0.0, -0.0846208019, 0.6897727967],
        [0.0, 0.9977802761, -1.1068370147],
        [0.0, -0.8465433753, -1.2673164024],
        [0.0, 0.5835275291, 1.5364927734],
    ];
    let w1_massvec = vec![MASS_C, MASS_C, MASS_H, MASS_H, MASS_H];
    let w1_brot = get_brot(&w1_xyz, &w1_massvec).to_vec();
    let w1_freq = vec![
        722.40, 810.05, 914.06, 1067.96, 1392.07, 1616.66, 3073.17, 3179.01, 3249.24,
    ];

    let ts_xyz = vec![
        [0.0, 0.1183768384, -0.5562510626],
        [0.0, -0.0058980147, 0.6558697182],
        [0.0, 0.5386263624, -1.5351015600],
        [0.0, -0.2821730201, 1.6836028161],
        [0.0, -1.5967941880, -1.3355954010],
    ];
    let ts_massvec = vec![MASS_C, MASS_C, MASS_H, MASS_H, MASS_H];
    let ts_brot = get_brot(&ts_xyz, &ts_massvec).to_vec();
    let ts_freq = vec![
        426.55, 601.17, 623.87, 745.02, 848.44, 1924.52, 3392.98, 3472.32,
    ];

    let d_h0_cm1 = ea_forward_kcal_mol / KCAL_PER_CM1;
    let k_e = get_kE(
        nebin,
        d_e,
        ts_freq.len(),
        w1_freq.len(),
        &ts_freq,
        &w1_freq,
        3,
        3,
        &ts_brot,
        &w1_brot,
        1.0,
        1.0,
        d_h0_cm1,
    );

    let freq_bin_w1: Vec<usize> = w1_freq.iter().map(|w| (w / d_e + 0.5) as usize).collect();
    let rho_e = get_rovib_WE_or_rhoE(
        "den".to_string(),
        w1_freq.len(),
        nebin,
        d_e,
        3,
        &freq_bin_w1,
        &w1_brot,
    );

    Ok(C2h3RrkmMicrocanonicalData {
        rho_e_w1: rho_e,
        k_e_w1_to_p1: k_e,
    })
}

pub fn build_c2h3_rrkm_ts_sum_of_states(grid: &EnergyGrid) -> Vec<f64> {
    let d_e = grid.bin_width_wavenumber;
    let nebin = grid.number_of_bins - 1;
    const MASS_C: f64 = 12.0;
    const MASS_H: f64 = 1.0;

    let ts_xyz = vec![
        [0.0, 0.1183768384, -0.5562510626],
        [0.0, -0.0058980147, 0.6558697182],
        [0.0, 0.5386263624, -1.5351015600],
        [0.0, -0.2821730201, 1.6836028161],
        [0.0, -1.5967941880, -1.3355954010],
    ];
    let ts_massvec = vec![MASS_C, MASS_C, MASS_H, MASS_H, MASS_H];
    let ts_brot = get_brot(&ts_xyz, &ts_massvec).to_vec();
    let ts_freq = vec![
        426.55, 601.17, 623.87, 745.02, 848.44, 1924.52, 3392.98, 3472.32,
    ];
    let freq_bin_ts: Vec<usize> = ts_freq.iter().map(|w| (w / d_e + 0.5) as usize).collect();

    get_rovib_WE_or_rhoE(
        "sum".to_string(),
        ts_freq.len(),
        nebin,
        d_e,
        3,
        &freq_bin_ts,
        &ts_brot,
    )
}
