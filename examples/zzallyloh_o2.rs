use MarXus::barrierless::sacm::{
    channel_eigval_interpol_from_input, compute_thermal_rates_debug, format_correlation_table,
    ProductRotorCase, ReactantRotorCase, SacmInterpolationInput, SacmThermalInput,
};

fn main() {
    let use_interpolation = true;

    // ZZAllylOH + O2 collisions; data mapped from legacy SACM Fortran input.
    let temperatures = (280..=330).step_by(5).map(|t| t as f64).collect();

    let vib_reactant = vec![
        3537.49, 3238.21, 3120.89, 3091.88, 3043.10, 3036.24, 2996.47, 2984.04, 2926.42, 1750.07,
        1647.54, 1519.50, 1470.94, 1444.56, 1440.59, 1440.19, 1422.76, 1385.61, 1367.18, 1332.61,
        1269.97, 1225.52, 1181.68, 1101.72, 1077.70, 1056.16, 1047.02, 995.46, 980.77, 957.87,
        869.24, 818.66, 734.81, 634.41, 565.17, 534.03, 461.15, 386.02, 337.04, 291.10, 252.85,
        224.83, 201.26, 191.78, 169.51, 110.81, 108.21, 81.71, 155.31, 42.21, 34.51,
    ];

    let vib_frag1 = vec![
        3716.44, 3115.66, 3093.44, 3076.52, 3049.20, 3014.44, 2994.15, 2932.40, 1774.17, 1684.67,
        1449.41, 1444.40, 1435.69, 1405.87, 1376.42, 1374.62, 1351.36, 1306.48, 1267.74, 1183.07,
        1113.29, 1083.96, 1059.20, 1019.78, 1016.47, 973.55, 921.22, 860.88, 827.60, 635.45,
        497.52, 447.05, 429.84, 379.45, 271.25, 229.82, 221.49, 167.96, 149.51, 129.45, 100.19,
        24.67,
    ];

    let vib_frag2 = vec![3590.82, 1416.91, 1218.12];

    // QSTAR input from the Fortran channel list (eps,x pairs).
    let qstar_eps_base = vec![
        3716.44, 3590.82, 3115.66, 3093.44, 3076.52, 3049.20, 3014.44, 2994.15, 2932.40, 1774.17,
        1684.67, 1449.41, 1444.40, 1435.69, 1416.91, 1405.87, 1376.42, 1374.62, 1351.36, 1306.48,
        1267.74, 1218.12, 1183.07, 1113.29, 1083.96, 1059.20, 1019.78, 1016.47, 973.55, 921.22,
        860.88, 827.60, 635.45, 497.52, 447.05, 429.84, 379.45, 271.25, 229.82, 221.49, 167.96,
        149.51, 129.45, 100.19, 24.67, 20.94481, 1.15778, 1.09713, 0.105636, 0.040089, 0.033272,
    ];
    let qstar_x_base = vec![
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
    ];

    // EPSE/XE list from the paired input (reactant eigenvalues for SPOL/QSTAR).
    let reactant_eps = vec![
        3537.49, 3238.21, 3120.89, 3091.88, 3043.10, 3036.24, 2996.47, 2984.04, 2926.42, 1750.07,
        1647.54, 1519.50, 1470.94, 1444.56, 1440.59, 1440.19, 1422.76, 1385.61, 1367.18, 1332.61,
        1269.97, 1225.52, 1181.68, 1101.72, 1077.70, 1056.16, 1047.02, 995.46, 980.77, 957.87,
        869.24, 818.66, 734.81, 634.41, 565.17, 534.03, 461.15, 386.02, 337.04, 291.10, 252.85,
        224.83, 201.26, 191.78, 169.51, 110.81, 108.21, 81.71, 42.21, 34.51, 0.068433,
    ];
    let reactant_x = vec![1.0; 50]
        .into_iter()
        .chain([0.5].into_iter())
        .collect::<Vec<f64>>();

    let product_eps = vec![
        3716.44, 3590.82, 3115.66, 3093.44, 3076.52, 3049.20, 3014.44, 2994.15, 2932.40, 1774.17,
        1684.67, 1449.41, 1444.40, 1435.69, 1416.91, 1405.87, 1376.42, 1374.62, 1351.36, 1306.48,
        1267.74, 1218.12, 1183.07, 1113.29, 1083.96, 1059.20, 1019.78, 1016.47, 973.55, 921.22,
        860.88, 827.60, 635.45, 497.52, 447.05, 429.84, 379.45, 271.25, 229.82, 221.49, 167.96,
        149.51, 129.45, 100.19, 24.67, 20.94481, 1.15778, 1.09713, 0.105636, 0.040089, 0.033272,
    ];
    let product_x = vec![
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
    ];

    let reacoord_freq = 155.31;
    let beta_eff = 1.524;
    let alpha_over_beta = 0.35;

    let reactant_zpe = 0.5 * vib_reactant.iter().sum::<f64>();
    let product_zpe = 0.5 * (vib_frag1.iter().sum::<f64>() + vib_frag2.iter().sum::<f64>());
    let h0k = 3217.75;
    let dmorse = h0k + reactant_zpe - product_zpe;
    let zpe_shift = reactant_zpe - product_zpe - reacoord_freq / 2.0;

    let red_mass = 116.0 * 33.0 / (116.0 + 33.0);
    let beta = 0.1218 * (red_mass / dmorse).sqrt() * reacoord_freq;
    let alpha_over_beta_eff = alpha_over_beta * beta / beta_eff;
    let dmoref = (beta / beta_eff).powi(2) * dmorse;

    let nfg = reactant_eps
        .len()
        .min(product_eps.len())
        .min(reactant_x.len())
        .min(product_x.len());
    let mut sum = 0.0;
    let mut count = 0;
    for i in 0..nfg.saturating_sub(1) {
        let x_r = reactant_x[i];
        let x_p = product_x[i];
        if x_p > 0.51 {
            continue;
        }
        if x_r < 0.51 {
            continue;
        }
        sum += reactant_eps[i];
        count += 1;
    }
    let c3 = if count == 0 {
        1.0
    } else {
        let avg = sum / count as f64;
        4.0 - 1.29 * (avg / dmoref).ln()
    };

    let (qstar_eps, qstar_x, channel) = if use_interpolation {
        let interpol = SacmInterpolationInput {
            alfa_beta: alpha_over_beta_eff,
            dmorse: dmoref,
            freq_react: reactant_eps.clone(),
            freq_prod: product_eps.clone(),
            x_react: reactant_x.clone(),
            x_prod: product_x.clone(),
        };
        let channel = channel_eigval_interpol_from_input(&interpol);
        (channel.freq_ts.clone(), channel.x_ts.clone(), Some(channel))
    } else {
        (qstar_eps_base, qstar_x_base, None)
    };

    let qe = 4.058;
    let centrifugal_a1 = 2.0 / (beta * qe);
    let centrifugal_a2 = 1.0 / (beta * qe * beta * qe);

    let input = SacmThermalInput {
        temperatures,
        max_temperature: 330.0,
        rotor_case_reactant: ReactantRotorCase::ProlateSymmetricTop,
        rotor_case_products: ProductRotorCase::SphericalSpherical,
        alpha_over_beta,
        anisotropy_scale: c3 * beta / beta_eff,
        transition_symmetry: 1.0,
        reactant_symmetry: 1.0,
        frag1_symmetry: 1.0,
        frag2_symmetry: 1.0,
        dissociation_energy: dmorse,
        product_zpe,
        zpe_shift,
        enthalpy_0k: h0k,
        base_rot_const: (0.025008 + 0.019026) / 2.0,
        centrifugal_a1,
        centrifugal_a2,
        external_rot_a: 0.068433,
        external_rot_b: 0.025008,
        external_rot_c: 0.019026,
        frag_rot_b1: 0.052,
        frag_rot_b2: 2.985,
        red_mass,
        vib_reactant,
        vib_frag1,
        vib_frag2,
        internal_rot_reactant: vec![],
        internal_rot_frag1: vec![],
        internal_rot_frag2: vec![],
        hindered_barrier_reactant: 0.0,
        hindered_barrier_frag1: 0.0,
        hindered_barrier_frag2: 0.0,
        hindered_count_reactant: 0.0,
        hindered_count_frag1: 0.0,
        hindered_count_frag2: 0.0,
        electronic_q_reactant: 2.0,
        electronic_q_frag1: 1.0,
        electronic_q_frag2: 2.0,
        qstar_eps,
        qstar_x,
    };

    let rates = compute_thermal_rates_debug(&input);
    println!("ZZAllylOH + O2 thermal capture rates (Fortran-style input)");
    println!(
        "SPOL interpolation: {}",
        if use_interpolation {
            "enabled"
        } else {
            "disabled"
        }
    );
    if let Some(channel) = channel {
        println!(
            "{}",
            format_correlation_table(
                &reactant_eps,
                &reactant_x,
                &channel.freq_ts,
                &channel.x_ts,
                &product_eps,
                &product_x,
            )
        );
    }
    println!("T/K     k_diss     k_recomb      K_eq");
    for rate in &rates {
        println!(
            "{:<6.0} {:>10.3e} {:>10.3e} {:>10.3e}",
            rate.temperature,
            rate.dissociation_rate,
            rate.recombination_rate,
            rate.equilibrium_constant
        );
    }

    if let Some(rate) = rates.first() {
        println!("\nDebug terms at T = {:.0} K:", rate.temperature);
        println!("KT      = {:.6}", rate.kt);
        println!("QCENT   = {:.6e}", rate.qcent);
        println!("ENULL(0)= {:.6}", rate.enull_0);
        println!("ENULL ladder (J, E0):");
        for (j, e0) in rate.enull.iter().enumerate() {
            println!("  {:>3}  {:.6}", j, e0);
        }
        println!(
            "QROT (EXTERNAL): {:>12.4e} {:>12.4e} {:>12.4e}",
            rate.qera, rate.qerb, rate.qerc
        );
        println!(
            "QROT (INTERNAL): {:>12.4e} {:>12.4e} {:>12.4e}",
            rate.qira, rate.qirb, rate.qirc
        );
        println!(
            "QVIB:          {:>12.4e} {:>12.4e} {:>12.4e}",
            rate.qviba, rate.qvibb, rate.qvibc
        );
        println!(
            "QEL:           {:>12.4e} {:>12.4e} {:>12.4e}",
            rate.qela, rate.qelb, rate.qelc
        );
        println!("QVRA    = {:.6e}", rate.qvra);
        println!("QVRB    = {:.6e}", rate.qvrb);
        println!("QVRC    = {:.6e}", rate.qvrc);
        println!("QSTAR   = {:.6e}", rate.qstarp);
        println!("QPROD   = {:.6e}", rate.qprod);
    }
}
