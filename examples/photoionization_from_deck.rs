//! Breakdown curves of dissociative photoionization for an input deck (`photoion::deck`): statistical,
//! microcanonical modelling of threshold photoelectron photoion coincidence (PEPICO) experiments without collisions.
//!
//!   cargo run --release --example photoionization_from_deck -- deck.inp [--csv FILE]
//!
//! The model follows B. Sztáray, A. Bodi, T. Baer, J. Mass Spectrom. 45, 1233 (2010) [SBB10]:
//! - thermal energy distribution of the neutral (eq. 1), transposed onto the ion at each photon energy (eq. 2) and
//!   convolved with the energy resolution;
//! - fast dissociation (eq. 23), or RRKM rate constants (eq. 6) competing within the flight time (eqs. 21, 24);
//! - statistical partitioning of the excess energy between fragment ion, neutral fragment and translation (eq. 5),
//!   passed on to sequential fast dissociations.
//! Rate models PhaseSpaceTheory and SimplifiedSACM are accepted in the deck and refused at run time (not available
//! yet).
//!
//! Options:
//! --csv FILE   breakdown curves as CSV: photon energy (eV) and the abundance of every ion
//! Unknown options are refused.

use MarXus::photoion::deck::{parse_photoionization_deck_file, ChannelRate, PhotoionizationDeck};
use MarXus::photoion::energy_distributions::EV_TO_CM1;
use MarXus::photoion::model::{BreakdownRow, PhotoionizationModel};
use std::io::Write;

fn main() -> Result<(), String> {
    let mut deck_path: Option<String> = None;
    let mut csv_path: Option<String> = None;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--csv" => csv_path = Some(args.next().ok_or("--csv needs a file name")?),
            option if option.starts_with("--") => return Err(format!("unknown option '{option}' (--csv)")),
            path if deck_path.is_none() => deck_path = Some(path.to_string()),
            extra => return Err(format!("unexpected argument '{extra}'")),
        }
    }
    let deck_path = deck_path.ok_or("usage: photoionization_from_deck deck.inp [--csv FILE]")?;
    let deck = parse_photoionization_deck_file(&deck_path)?;
    let model = PhotoionizationModel::new(&deck)?;
    let curves = model.breakdown_curves()?;
    print_report(&deck_path, &deck, &model, &curves);
    if let Some(path) = csv_path {
        write_csv(&path, &model, &curves).map_err(|e| format!("{path}: {e}"))?;
    }
    Ok(())
}

fn print_report(path: &str, deck: &PhotoionizationDeck, model: &PhotoionizationModel, curves: &[BreakdownRow]) {
    let rule = "=".repeat(100);
    println!("{rule}\n MarXus dissociative photoionization (statistical, collision-free)\n{rule}");
    println!(" Deck:               {path}");
    println!(" Model:              Sztáray, Bodi, Baer, J. Mass Spectrom. 45, 1233 (2010), eqs. 1, 2, 5, 6, 21, 23, 24");
    println!(" Neutral, ion:       {} -> {}  (IE = {} eV), T = {} K", deck.neutral, deck.ion, deck.ionization_energy_ev, deck.temperature_kelvin);
    println!(
        " Energy grid:        cells of {} cm-1; {} neutral cells (thermal tail below 1e-15 left out), {} ion cells",
        deck.cell_width_cm1,
        model.neutral_distribution.len(),
        model.ion_cells
    );
    match deck.resolution_fwhm_cm1 {
        Some(fwhm) => println!(" Resolution:         Gaussian, FWHM {fwhm} cm-1 (photon and electron energy resolution)"),
        None => println!(" Resolution:         none"),
    }
    if let Some(tau) = deck.flight_time_s {
        println!(" Flight time:        {tau:e} s");
    }
    println!("\n Channels (appearance energies from the ground state of the neutral):");
    println!("   {:<14} {:>10} {:>10} {:>10} {:>12} {:>14}  {}", "name", "parent", "fragment", "neutral", "E0 [eV]", "E0 - IE [cm-1]", "rate");
    for c in &deck.channels {
        let rate = match c.rate_model {
            ChannelRate::Fast => "Fast".to_string(),
            ChannelRate::Statistical(m) => format!("{m:?}, TS {}", c.transition_state.as_deref().unwrap_or("-")),
        };
        println!(
            "   {:<14} {:>10} {:>10} {:>10} {:>12.4} {:>14.1}  {rate}; translational d.o.f. {}",
            c.name,
            c.parent,
            c.fragment_ion,
            c.neutral_fragment,
            c.appearance_energy_ev,
            (c.appearance_energy_ev - deck.ionization_energy_ev) * EV_TO_CM1,
            c.translational_degrees_of_freedom
        );
    }
    println!("\n Breakdown curves (fractional ion abundances):");
    print!("   {:>10}", "h nu [eV]");
    for ion in &model.ions {
        print!(" {:>12}", ion);
    }
    println!();
    for row in curves {
        print!("   {:>10.4}", row.photon_energy_ev);
        for a in &row.abundances {
            print!(" {:>12.6}", a);
        }
        println!();
    }
}

fn write_csv(path: &str, model: &PhotoionizationModel, curves: &[BreakdownRow]) -> std::io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "photon_energy_eV,{}", model.ions.join(","))?;
    for row in curves {
        let values: Vec<String> = row.abundances.iter().map(|a| format!("{a:.10e}")).collect();
        writeln!(f, "{},{}", row.photon_energy_ev, values.join(","))?;
    }
    Ok(())
}
