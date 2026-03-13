use super::locate_centrifugal_barrier::centrifugal_barrier_bisect;
use super::morsepot::morse_value_and_derivative;

#[test]
fn test_locatingcentrifug() {
    let de = 0.02;
    let beta = 1.2;
    let re = 3.0;

    let mu = 1000.0;

    let r_lo = 0.0005;
    let r_hi = 20.0;

    println!("\nJ        r_max(bohr)        V_max(Ha)");

    for j in 0..=60 {
        match centrifugal_barrier_bisect(
            j,
            mu,
            |r| morse_value_and_derivative(r, de, beta, re),
            de,
            r_lo,
            r_hi,
        ) {
            Some((rmax, vmax)) => {
                println!("{:3}      {:12.6}      {:14.10}", j, rmax, vmax);
            }
            None => {
                println!("{:3}      {:12}      {}", j, "-", "NO_BARRIER");
            }
        }
    }
}
