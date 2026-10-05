//! Canonical Eckart tunneling correction factors kappa(T) of every barrier of an input deck in the MESS input
//! format that has a `Tunneling Eckart` block.
//!
//!   cargo run --release --example eckart_kappa_from_deck -- deck.inp [T1 T2 ...]
//!
//! kappa(T) = beta exp(beta V_f) integral P(E - V_f) exp(-beta E) dE with the exact Eckart transmission
//! probability P (Miller, J. Am. Chem. Soc. 101, 6810 (1979), eq. 8; `tunneling::eckart`), V_f and V_b the
//! two well depths of the deck and the magnitude of its imaginary frequency. This is the factor by which
//! tunneling raises the high-pressure rate coefficient of the barrier.
//!
//! The same factor of the MESS semiclassical Eckart model (`tunneling::mess_eckart_tunneling`) is given in
//! the column `kappa_mess`.
//!
//! Default temperatures: 100, 200, ..., 2000 K. Output: CSV `barrier,T[K],kappa,kappa_mess`.

use MarXus::masterequation::mess_input::{parse_mess_input_file, TunnelingSpecification};
use MarXus::tunneling::mess_eckart_tunneling::MessEckartTunneling;
use MarXus::tunneling::tunneling::eckart;

/// Boltzmann constant in cm-1/K (CODATA 2018: 0.695 034 800 cm-1/K).
const KB_CM1_PER_K: f64 = 0.695_034_800;

fn main() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let path = args.next().ok_or("usage: eckart_kappa_from_deck deck.inp [T1 T2 ...]")?;
    let mut temperatures: Vec<f64> =
        args.map(|a| a.parse::<f64>().map_err(|_| format!("invalid temperature '{a}'"))).collect::<Result<_, _>>()?;
    if temperatures.is_empty() {
        temperatures = (1..=20).map(|i| 100.0 * i as f64).collect();
    }
    let deck = parse_mess_input_file(&path)?;
    println!("barrier,T[K],kappa,kappa_mess");
    for barrier in &deck.barriers {
        if let Some(TunnelingSpecification::Eckart { imaginary_frequency_cm1, well_depths_cm1 }) = &barrier.tunneling {
            let mess = MessEckartTunneling::new(*imaginary_frequency_cm1, *well_depths_cm1)?;
            for &t in &temperatures {
                let kt = KB_CM1_PER_K * t;
                // Integration from the forward asymptote to 40 kT above the higher of the two asymptotes, in
                // steps of 0.1 cm-1 (Simpson rule).
                let e_max = well_depths_cm1[0].max(well_depths_cm1[1]) + 40.0 * kt;
                let kappa = eckart(1.0 / kt, *imaginary_frequency_cm1, well_depths_cm1[0], well_depths_cm1[1], 0.1, e_max);
                println!("{},{},{:.6e},{:.6e}", barrier.name, t, kappa, mess.canonical_factor(kt));
            }
        }
    }
    Ok(())
}
