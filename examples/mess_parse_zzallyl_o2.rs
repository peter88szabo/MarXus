use MarXus::masterequation::mess_input::{parse_mess_input, MessBarrierCore};

fn main() -> Result<(), String> {
    // Parse the embedded input deck (MESS input format) and list what was read.
    let mess = include_str!("mess_zzallyl_o2_case1_excerpt.inp");
    let deck = parse_mess_input(mess)?;

    println!("Temperatures [K]: {:?}", deck.global.temperatures_kelvin);
    println!("Pressures [Torr]: {:?}", deck.global.pressures_torr);
    println!("Reactant: {:?}", deck.global.reactant_name);

    println!("\nWells:");
    for name in &deck.well_order {
        let w = &deck.wells[name];
        println!(
            "  {:<6} E0 = {:>10.1} cm-1  frequencies = {:>3}  escape = {:?} s-1",
            name,
            w.zero_energy_cm1,
            w.vibrational_frequencies_cm1.len(),
            deck.well_escape_rate_s_inv.get(name)
        );
    }

    println!("\nBarriers:");
    for b in &deck.barriers {
        let core = match b.core {
            MessBarrierCore::TightRrho => "tight",
            MessBarrierCore::PhaseSpaceTheory { .. } => "phase-space theory",
        };
        println!(
            "  {:<6} {:>4} <-> {:<4} E0 = {:>10.1} cm-1  core: {core}  ILT: {}  tunneling: {}",
            b.name,
            b.left,
            b.right,
            b.rrho.zero_energy_cm1,
            b.inverse_laplace_transform.is_some(),
            b.tunneling.is_some()
        );
    }
    Ok(())
}
