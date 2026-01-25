use MarXus::barrierless::sacm::{
    compute_thermal_rates, ProductRotorCase, ReactantRotorCase, SacmThermalInput,
};

fn main() {
    let temperatures = vec![300.0, 500.0, 1000.0, 1500.0];

    // Literature values provided by the user.
    let oh_freq_cm1 = 3737.76;
    let oh_req_ang = 0.969;
    let co_freq_cm1 = 2169.814;
    let co_req_ang = 1.128;

    let oh_brot_cm1 = diatomic_brot_cm1(15.999, 1.008, oh_req_ang);
    let co_brot_cm1 = diatomic_brot_cm1(12.0, 15.999, co_req_ang);

    // Morse well depth and disappearing-mode frequency (from user).
    let de_kj_mol = 150.0;
    let dissociation_energy_cm1 = de_kj_mol * 83.593;
    let morse_mode_cm1 = 250.0;

    // PST-like: no anisotropy reduction (baseline for comparison).
    let pst_input = SacmThermalInput {
        temperatures: temperatures.clone(),
        max_temperature: 1500.0,
        rotor_case_reactant: ReactantRotorCase::Linear,
        rotor_case_products: ProductRotorCase::LinearLinear,
        alpha_over_beta: 4.0,
        anisotropy_scale: 0.0,
        transition_symmetry: 1.0,
        reactant_symmetry: 1.0,
        frag1_symmetry: 1.0,
        frag2_symmetry: 1.0,
        dissociation_energy: dissociation_energy_cm1,
        product_zpe: 0.0,
        zpe_shift: 0.0,
        enthalpy_0k: 0.0,
        base_rot_const: co_brot_cm1,
        centrifugal_a1: 0.0,
        centrifugal_a2: 0.0,
        external_rot_a: 0.0,
        external_rot_b: co_brot_cm1,
        external_rot_c: 0.0,
        frag_rot_b1: co_brot_cm1,
        frag_rot_b2: oh_brot_cm1,
        red_mass: reduced_mass_amu(28.0, 17.0),
        vib_reactant: vec![],
        vib_frag1: vec![co_freq_cm1],
        vib_frag2: vec![oh_freq_cm1],
        internal_rot_reactant: vec![],
        internal_rot_frag1: vec![],
        internal_rot_frag2: vec![],
        hindered_barrier_reactant: 0.0,
        hindered_barrier_frag1: 0.0,
        hindered_barrier_frag2: 0.0,
        hindered_count_reactant: 0.0,
        hindered_count_frag1: 0.0,
        hindered_count_frag2: 0.0,
        electronic_q_reactant: 1.0,
        electronic_q_frag1: 1.0,
        electronic_q_frag2: 1.0,
        qstar_eps: vec![morse_mode_cm1],
        qstar_x: vec![1.0],
    };

    // SACM-like: include anisotropy reduction.
    let mut sacm_input = pst_input.clone();
    sacm_input.anisotropy_scale = 1.0;

    let pst_rates = compute_thermal_rates(&pst_input);
    let sacm_rates = compute_thermal_rates(&sacm_input);

    println!("CO + OH thermal capture rates (illustrative inputs)");
    println!("SPOL interpolation: disabled");
    println!("T/K     k_diss(PST)    k_recomb(PST)    k_diss(SACM)   k_recomb(SACM)");

    for (pst, sacm) in pst_rates.iter().zip(sacm_rates.iter()) {
        println!(
            "{:<6.0} {:>12.3e} {:>12.3e} {:>12.3e} {:>12.3e}",
            pst.temperature,
            pst.dissociation_rate,
            pst.recombination_rate,
            sacm.dissociation_rate,
            sacm.recombination_rate,
        );
    }
}

fn reduced_mass_amu(m1: f64, m2: f64) -> f64 {
    m1 * m2 / (m1 + m2)
}

fn diatomic_brot_cm1(m1_amu: f64, m2_amu: f64, bond_angstrom: f64) -> f64 {
    let mu = reduced_mass_amu(m1_amu, m2_amu);
    // B (cm^-1) = h / (8 pi^2 c I), with I in amu*Å^2 and the constant folded below.
    let inertia = mu * bond_angstrom * bond_angstrom;
    16.857629 / inertia
}
