use MarXus::masterequation::mess_input::{parse_mess_input, MessBuildOptions};

fn main() -> Result<(), String> {
    // Minimal smoke test: parse the embedded MESS deck and build the MarXus multiwell network.
    //
    // This example does not solve the master equation yet; it only demonstrates the new
    // "read MESS input" pathway and the geometry->Brot computation.
    let mess = include_str!("mess_zzallyl_o2_case1_excerpt.inp");
    let deck = parse_mess_input(mess)?;

    let built = deck.build_multiwell_network(MessBuildOptions {
        energy_grain_width_cm1: Some(20.0),
        max_energy_kcal_mol: 50.0,
        collision_band_half_width: 20,
        nonreactive_grain_count: 10,
    })?;

    println!("Parsed MESS deck:");
    println!(
        "T = {:.1} K, P = {:.1} torr, wells = {}, barriers = {}",
        deck.global.temperature_kelvin,
        deck.global.pressure_torr,
        built.wells.len(),
        deck.barriers.len()
    );

    println!();
    println!("Wells + channels (MarXus network):");
    for (i, w) in built.wells.iter().enumerate() {
        println!(
            "  well[{}] {:<4}  ΔE={:>6.1}  bins={}  offset={}  channels={}",
            i,
            w.well_name,
            w.energy_grain_width_cm1,
            w.one_past_highest_included_grain_index,
            w.alignment_offset_in_grains,
            w.channels.len()
        );
        for ch in &w.channels {
            println!("    - {}", ch.name);
        }
    }

    println!();
    println!("Micro models:");
    println!("  wells: {}", built.well_models_by_name.len());
    println!("  channels: {}", built.channel_models_by_name.len());

    Ok(())
}
