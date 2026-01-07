#![allow(non_snake_case)]
mod inertia;
mod molecule;
mod numeric;
mod rrkm;
mod thermal;
mod utils;

use crate::inertia::inertia::get_brot;
use crate::molecule::MolType;
use crate::molecule::MoleculeBuilder;
use crate::rrkm::rrkm_rate::get_kE;
use crate::utils::print_rates::print_rrkm_rates;
use crate::utils::time::format_duration;
use std::time::Instant;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fn main() {
    let logo = r#"
================================================================

    M     M     AAAAA     RRRRR     X   X     U   U     SSSSS
    MM   MM    A     A    R    R     X X      U   U    S
    M M M M    AAAAAAA    RRRRR       X       U   U     SSSS
    M  M  M    A     A    R   R      X X      U   U         S
    M     M    A     A    R    R    X   X     UUUUU    SSSSS

 Molecular Statistical Physics for Kinetics and Thermochemistry

                       Version: 0.1
                      11. Sept. 2024

 Author: Peter Szabo
 Email: peter88szabo@gmail.com

 Features:
    - Thermochemistry:
        + Partition functions
        + Thermodynamical energy functions (U, H, F, G, S, Cv, Cp)
        + Equilbrium constants

    - Kinetics:
        + Canonical-TST high-pressure reaction rate cofficients
        + RRKM: E-dependent rate constants (μ-canonical TST)
        + Sum and Density of states (E- or E,J-resolved)
        + Barrierless reactions:
            * Phase Space Theory (PST) with 1D arbitrary potential
            * Statistical Adiabatic Channel Method (SACM)
            * Inverse Laplace Transform
            * Microcanonical Variational TST
            * Canonical Variational TST
        + Sumberged barrier with a pre-reaction vdW complex:
          μ-canonical J-resolved 2-TST treatment

================================================================
"#;
    println!("{}", logo);

    let start_time = Instant::now();

    //==================================================
    // Input for complex and TS structures to run RRKM
    //==================================================
    // Energy grain
    //--------------------------------------------------
    let dE = 10.0; //cm-1
    let Emax = 2000_000.0; // cm-1
    let nebin = (Emax / dE + 0.5) as usize;

    println!();
    println!("ΔE (energy grain):              {:>12.1} cm-1", dE);
    println!("Emax:                           {:>12.1} cm-1 \n", Emax);
    //--------------------------------------------------
    // Vibrations
    //--------------------------------------------------
    // R4 reaction:
    //let omega_ts = [ 184.88,  319.96,  340.12,  395.92,  452.98,
    //                 492.44,  541.35,  769.56,  829.18,  895.49,
    //                 929.80,  980.90, 1009.13, 1014.58, 1045.27,
    //                1052.21, 1072.72, 1131.33, 1260.49, 1277.77,
    //                1300.86, 1362.49, 1433.14, 1478.45, 1535.70,
    //                1550.28, 1641.89, 3046.95, 3174.54, 3177.29,
    //                3179.81, 3201.11, 3257.75, 3268.34, 3276.94];

    //R4 reaction:
    //let omega_cpx = [  66.49,   96.36,  134.98,  291.92,  357.88,
    //                  431.61,  532.36,  680.46,  720.05,  887.43,
    //                  903.05,  937.96,  987.74,  996.72, 1003.24,
    //                 1049.34, 1134.90, 1170.06, 1220.56, 1313.31,
    //                 1323.55, 1355.86, 1392.40, 1437.39, 1465.22,
    //                 1504.04, 1742.44, 1745.17, 3053.49, 3128.03,
    //                 3166.41, 3175.20, 3187.22, 3203.11, 3257.68, 3299.33];

    // R2 cis reaction:
    //let omega_cpx = [   25.16,  181.86,  193.60,  202.13,  404.57,
    //                   451.61,  622.42,  647.48,  779.21,  815.43,
    //                   889.19,  961.60,  986.09,  989.69, 1028.92,
    //                  1038.35, 1046.69, 1084.58, 1107.00, 1267.75,
    //                  1329.42, 1334.31, 1423.25, 1465.76, 1471.87,
    //                  1717.11, 1728.79, 1733.48, 3149.81, 3153.54,
    //                  3158.61, 3175.43, 3187.75, 3195.44, 3274.53, 3289.20];

    // R2 trans reaction:
    //let omega_cpx = [
    //    109.05, 154.18, 175.97, 339.01, 355.77, 400.01, 612.88, 687.80, 746.10, 857.79, 903.69,
    //    965.21, 968.72, 970.14, 1005.69, 1038.83, 1051.69, 1106.72, 1216.26, 1286.91, 1310.84,
    //    1357.69, 1393.49, 1443.74, 1491.42, 1662.13, 1716.88, 1738.19, 3171.49, 3173.58, 3174.80,
    //    3192.26, 3197.05, 3212.98, 3266.05, 3266.53,
    //];

   // delta_HAPLD at r2SCAN-3c:
    let omega_cpx = [
           46.2, 105.1,  121.2, 153.4, 178.0, 209.99, 250.56, 299.19, 353.39, 372.40, 427.02, 592.9, 625.18, 806.45, 835.75, 898.66, 979.98,
           993.58, 1014.72, 1054.21, 1065.77, 1096.93, 1214.43, 1258.53, 1350.98, 1360.89, 1384.07, 1404.62, 1426.19, 1461.15, 1476.33,
           1524.63, 1682.97, 1742.44, 2921.70, 2994.68, 3022.35, 3056.12, 3094.06, 3124.87, 3148.63, 3754.08];

    // R2 reacion:
    //let omega_ts = [
    //    63.06, 357.13, 364.21, 521.00, 578.21, 688.13, 697.88, 840.98, 850.93, 898.00, 949.24,
    //    973.29, 1008.88, 1024.17, 1034.57, 1043.84, 1091.47, 1131.18, 1215.54, 1274.62, 1288.67,
    //    1420.26, 1465.22, 1513.71, 1530.23, 1581.74, 1640.52, 3156.38, 3158.14, 3173.90, 3178.81,
    //    3185.86, 3202.64, 3278.43, 3330.59,
    //];

    //HPALD to PerHemiac:
    let omega_ts = [
           71.48, 171.12, 203.86, 243.88, 301.37, 309.48, 410.68, 464.99, 540.68, 584.94, 657.23, 697.48,
           816.21, 829.65, 883.69, 932.80, 988.12, 1038.11, 1060.47, 1067.66, 1109.99, 1177.52, 1212.71,
           1283.31, 1311.58, 1351.25, 1369.83, 1393.17, 1452.92, 1455.01, 1464.96, 1485.77, 1718.65,
           1948.12, 2920.77, 2976.47, 3007.09, 3071.55, 3078.75, 3117.02, 3142.40 ];

    let nvib_cpx = omega_cpx.len();
    let nvib_ts = omega_ts.len();

    //use std::iter::Sum;
    let ZPE_cpx = 0.5 * omega_cpx.iter().sum::<f64>(); //0.5*(omega_cpx.iter().sum());
    let ZPE_ts = 0.5 * omega_ts.iter().sum::<f64>(); //0.5*(omega_cpx.iter().sum());

    println!(
        "ZPE of complex:                 {:>12.1} cm-1 {:>12.2} kcal/mol",
        ZPE_cpx,
        ZPE_cpx * 2.85914e-3
    );
    println!(
        "ZPE of TS:                      {:>12.1} cm-1 {:>12.2} kcal/mol\n",
        ZPE_ts,
        ZPE_ts * 2.85914e-3
    );

    //--------------------------------------------------
    // Rotations
    //--------------------------------------------------
    //R2 Complex Rotational constants in cm-1:     0.183752     0.096944     0.067233
    //let Brot_cpx = [0.183752, 0.096944, 0.067233];
    //HPALD:
    let Brot_cpx = [0.1245141, 0.03617214, 0.02993249];

    //R2 TS Rotational constants in cm-1:     0.169426     0.140868     0.081368
    //let Brot_ts = [0.169426, 0.140868, 0.081368];
      
    //HPALD to PerHemiac
    let Brot_ts = [0.09791913, 0.06980448, 0.04425807];

    let nrot_cpx = Brot_cpx.len();
    let nrot_ts = Brot_ts.len();
    //--------------------------------------------------
    // E0 = reaction energy (only electronic no ZPE)
    // ΔH0 = 0K heat of formation (including ZPEs)
    //
    // ΔH0 = E0 + ZPE_ts - ZPE_cpx
    //--------------------------------------------------
    // R4 reaction let Ezero = 10952.58; // in cm-1
    // let Ezero = 7546.117; // R2 cis reacton in cm-1
    //let Ezero = 10671.434; // R2 trans reacton in cm-1
    let Ezero = 241965.0; // HPALD to Per-Hemiacetal reacton in cm-1
    let dH0 = Ezero + ZPE_ts - ZPE_cpx;
    println!(
        "E0  (reaction energy at 0K):    {:>12.1} cm-1 {:>12.2} kcal/mol   (no ZPE only pure electronic)",
        Ezero, Ezero * 2.85914e-3
    );
    println!(
        "ΔH0 (reaction enthalpy at 0K):  {:>12.1} cm-1 {:>12.2} kcal/mol   (ΔH0 = E0 + ZPE_ts - ZPE_cpx)\n",
        dH0, dH0 * 2.85914e-3
    );

    if Ezero >= Emax {
        println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        println! {"Code stoped beacuse Emax is too small comapred to reaction energy (E0)"};
        println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        panic! {"!!!! Increase Emax !!!!"};
    }
    //--------------------------------------------------
    // Symmetry numbers
    //--------------------------------------------------
    let sigma_cpx = 1.0;
    let sigma_ts = 1.0;
    println!("symm_num_cpx:                   {:>10.0}", sigma_cpx);
    println!("symm_num_ts:                    {:>10.0} \n", sigma_ts);
    //==================================================
    // End of Input
    //==================================================

    //Minimum energy (including ZPE) of the reaction as integer energy bin
    let nbin_dH0 = (dH0 / dE + 0.5) as usize;

    // Compute RRKM rate constant
    let kRRKM = get_kE(
        nebin, dE, nvib_ts, nvib_cpx, &omega_ts, &omega_cpx, nrot_ts, nrot_cpx, &Brot_ts,
        &Brot_cpx, sigma_ts, sigma_cpx, dH0,
    );

    let fname = "microcanonical_rates_new.dat";

    print_rrkm_rates(fname, nbin_dH0, nebin, dE, dH0, &kRRKM);

    let end_time = Instant::now();

    let duration = end_time - start_time;
    let duration_seconds = duration.as_secs_f64();

    println!("--------------------------------------------------------");
    println!("RRKM calculation finished \n");
    //println!("Data written into file: microcanonical_rate.dat \n");

    println!("Computational time: {:>15.4} sec", duration_seconds);
    println!(
        "Computational time:          {} hh:mm:ss\n",
        format_duration(duration)
    );

    //============================================================================

    println!("");
    let name = "Water".to_string();
    let moltype = MolType::mol;

    let mut water = MoleculeBuilder::new(name, moltype)
        .freq(vec![440.0, 1600.0, 3600.0])
        .brot(vec![10.0, 10.0, 20.0])
        .dh0(199.9)
        .multi(3.0)
        .chiral(22.0)
        .symnum(6.0)
        .mass(16.0)
        .build();

    println!("ene is provided");
    println!("symnum: {:?}", water.symnum);
    println!("multi: {:?}", water.multi);
    println!("chiral: {:?}", water.chiral);
    println!("zpe: {:?}", water.zpe);
    println!("ene: {:?}", water.ene); // Output: 199.9
    println!("dh0: {:?}", water.dh0); // Output: ene + zpe
    println!("freq: {:?}", water.freq);
    println!("brot: {:?}", water.brot);

    println!("\nThermo");
    println!("gtot: {:?}", water.thermo.gtot);
    println!("srot: {:?}", water.thermo.srot);
    println!("hvib: {:?}", water.thermo.hvib);
    println!("pfvib: {:?}", water.thermo.pfvib);

    //water.calculate_entropy();
    //water.calculate_internal_energy();
}
