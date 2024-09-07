#![allow(non_snake_case)]
use std::fs;
use std::io::Write;

pub fn print_rrkm_rates(file_name: &str, nbin_dH0: usize, nebin: usize, dE: f64, dH0: f64, kRRKM: &[f64]) {
    let mut file = fs::File::create(file_name).unwrap();

    // Write the header line
    writeln!(
        &mut file,
        "{:>15} {:>15} {:>15} {:>15}",
        "E[cm-1]", "E-dH0[cm-1]", "rate(s-1)", "lifetime(s)"
    ).unwrap();

    // Write the data
    for i in nbin_dH0..=nebin {
        let ene = i as f64 * dE;
        writeln!(
            &mut file,
            "{:>15.1} {:>15.1} {:>15.6e} {:>15.6e}",
            ene, ene - dH0, kRRKM[i], 1.0 / kRRKM[i]
        ).unwrap();
    }
}

