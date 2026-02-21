use MarXus::inertia::inertia::get_brot;
use MarXus::molecule::{MolType, MoleculeBuilder, MoleculeStruct};

pub fn build_c2h3_w1_and_ts_molecules(
    well_dh0_kcal_mol: f64,
    barrier_dh0_kcal_mol: f64,
) -> (MoleculeStruct, MoleculeStruct) {
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
    let w1_mass_total: f64 = w1_massvec.iter().sum();
    let w1_brot = get_brot(&w1_xyz, &w1_massvec).to_vec();
    let w1_freq = vec![
        722.40, 810.05, 914.06, 1067.96, 1392.07, 1616.66, 3073.17, 3179.01, 3249.24,
    ];
    let w1_dh0_cm1 = well_dh0_kcal_mol / KCAL_PER_CM1;

    let ts_xyz = vec![
        [0.0, 0.1183768384, -0.5562510626],
        [0.0, -0.0058980147, 0.6558697182],
        [0.0, 0.5386263624, -1.5351015600],
        [0.0, -0.2821730201, 1.6836028161],
        [0.0, -1.5967941880, -1.3355954010],
    ];
    let ts_massvec = vec![MASS_C, MASS_C, MASS_H, MASS_H, MASS_H];
    let ts_mass_total: f64 = ts_massvec.iter().sum();
    let ts_brot = get_brot(&ts_xyz, &ts_massvec).to_vec();
    let ts_freq = vec![
        426.55, 601.17, 623.87, 745.02, 848.44, 1924.52, 3392.98, 3472.32,
    ];
    let ts_dh0_cm1 = barrier_dh0_kcal_mol / KCAL_PER_CM1;

    let w1 = MoleculeBuilder::new("W1".to_string(), MolType::reactant)
        .nlin(false)
        .mass(w1_mass_total)
        .brot(w1_brot)
        .freq(w1_freq)
        .dh0(w1_dh0_cm1)
        .multi(2.0)
        .symnum(1.0)
        .chiral(1.0)
        .build();

    let ts = MoleculeBuilder::new("B1".to_string(), MolType::ts)
        .nlin(false)
        .mass(ts_mass_total)
        .brot(ts_brot)
        .freq(ts_freq)
        .dh0(ts_dh0_cm1)
        .multi(2.0)
        .symnum(1.0)
        .chiral(1.0)
        .build();

    (w1, ts)
}

pub fn build_p1_reactant_molecules(product_dh0_kcal_mol: f64) -> (MoleculeStruct, MoleculeStruct) {
    const KCAL_PER_CM1: f64 = 2.85914e-3;
    const MASS_C: f64 = 12.0;
    const MASS_H: f64 = 1.0;

    let c2h2_massvec = vec![MASS_C, MASS_C, MASS_H, MASS_H];
    // Use stable linear-rotor constants directly to avoid singular inertia-axis artifacts.
    let c2h2_brot = vec![1.176, 1.176, 0.0];
    let c2h2_freq = vec![595.45, 595.45, 746.55, 746.55, 2004.77, 3409.28, 3501.09];
    let c2h2_mass_total: f64 = c2h2_massvec.iter().sum();
    let c2h2_dh0_cm1 = product_dh0_kcal_mol / KCAL_PER_CM1;

    let c2h2 = MoleculeBuilder::new("C2H2".to_string(), MolType::reactant)
        .nlin(true)
        .mass(c2h2_mass_total)
        .brot(c2h2_brot)
        .freq(c2h2_freq)
        .dh0(c2h2_dh0_cm1)
        .multi(1.0)
        .symnum(2.0)
        .chiral(1.0)
        .build();

    let mut h = MoleculeStruct::default();
    h.name = "H".to_string();
    h.moltype = MolType::reactant;
    h.mass = 1.0;
    h.multi = 2.0;
    h.symnum = 1.0;
    h.chiral = 1.0;
    h.dh0 = 0.0;
    h.ene = 0.0;

    (c2h2, h)
}

