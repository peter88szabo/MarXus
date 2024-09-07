#![allow(non_snake_case)]
mod utils;
mod numeric;
mod rrkm;

use crate::rrkm::rrkm_rate::get_kE;
use crate::utils::time::format_duration;
use crate::utils::print_rates::print_rrkm_rates;
use std::time::{Instant};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fn main() {
    let start_time = Instant::now();

    //==================================================
    // Input for complex and TS structures to run RRKM
    //==================================================
    // Energy grain
    //--------------------------------------------------
    let dE = 10.0; //cm-1
    let Emax = 200_000.0; // cm-1
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
    let omega_cpx = [  109.05,  154.18,  175.97,  339.01,  355.77,
                       400.01,  612.88,  687.80,  746.10,  857.79,
                       903.69,  965.21,  968.72,  970.14, 1005.69,
                      1038.83, 1051.69, 1106.72, 1216.26, 1286.91,
                      1310.84, 1357.69, 1393.49, 1443.74, 1491.42,
                      1662.13, 1716.88, 1738.19, 3171.49, 3173.58,
                      3174.80, 3192.26, 3197.05, 3212.98, 3266.05, 3266.53]; 

    // R2 reacion:     
    let omega_ts = [   63.06,  357.13,  364.21,  521.00,  578.21,
                      688.13,  697.88,  840.98,  850.93,  898.00,
                      949.24,  973.29, 1008.88, 1024.17, 1034.57,
                     1043.84, 1091.47, 1131.18, 1215.54, 1274.62,
                     1288.67, 1420.26, 1465.22, 1513.71, 1530.23,
                     1581.74, 1640.52, 3156.38, 3158.14, 3173.90,
                     3178.81, 3185.86, 3202.64, 3278.43, 3330.59];

    let nvib_cpx = omega_cpx.len();
    let nvib_ts = omega_ts.len();

    //use std::iter::Sum;
    let ZPE_cpx = 0.5 * omega_cpx.iter().sum::<f64>(); //0.5*(omega_cpx.iter().sum());
    let ZPE_ts = 0.5 * omega_ts.iter().sum::<f64>(); //0.5*(omega_cpx.iter().sum());

    println!("ZPE of complex:                 {:>12.1} cm-1 {:>12.2} kcal/mol", ZPE_cpx, ZPE_cpx * 2.85914e-3);
    println!("ZPE of TS:                      {:>12.1} cm-1 {:>12.2} kcal/mol\n", ZPE_ts, ZPE_ts * 2.85914e-3);
    //--------------------------------------------------
    // Rotations
    //--------------------------------------------------
   //R2 Complex Rotational constants in cm-1:     0.183752     0.096944     0.067233
    let Brot_cpx = [0.183752, 0.096944, 0.067233];

   //R2 TS Rotational constants in cm-1:     0.169426     0.140868     0.081368 
    let Brot_ts = [0.169426,     0.140868,     0.081368];
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
    let Ezero = 10671.434; // R2 trans reacton in cm-1
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
    println!("Data written into file: microcanonical_rate.dat \n");

    println!("Computational time: {:>15.4} sec", duration_seconds);
    println!(
        "Computational time:          {} hh:mm:ss\n",
        format_duration(duration)
    );
}
