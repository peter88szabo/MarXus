// thermofuncs.rs
use crate::constants::{
    AMU_TO_ELECTRON_MASS, AU_TO_KCAL, AU_TO_KJ, BOLTZMANN_SI, CLIGHT_SI, CM1_TO_HARTREE, CM1_TO_K,
    CM1_TO_KCAL, HPLANCK_AU, PASCAL_TO_AU, PI, PI_SQ, PLANCK_SI, RGAS_AU, RGAS_SI, TWO_PI,
};
use crate::molecule::MoleculeStruct;
const HPLANCK_AU_SQ: f64 = HPLANCK_AU * HPLANCK_AU;

// -----------------------------
// SI constants needed ONLY for Grimme free-rotor entropy (dimensionless inside ln)
// -----------------------------
// 1 Hartree/mol = AU_TO_KJ kJ/mol = AU_TO_KJ*1000 J/mol
const J_PER_HARTREE_PER_MOL: f64 = AU_TO_KJ * 1000.0;

// Grimme default average moment of inertia (kg*m^2)
const GRIMME_BAV_SI: f64 = 1.0e-44;

// Damping exponent in Grimme qRRHO
const GRIMME_ALPHA: f64 = 4.0;

/// Electronic contributions per mole (Hartree, Hartree/K), from levels (energy above the ground level in cm-1,
/// degeneracy), with x_j = eps_j / (k_B T):
///   q = sum_j g_j e^(-x_j),  U = H = RT <x>,  F = G = -RT ln q,  S = (U - F)/T,  Cv = Cp = R (<x^2> - <x>^2),
/// where <f> = sum_j g_j f(x_j) e^(-x_j) / q (canonical averages over the levels).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElectronicContributions {
    pub partition_function: f64,
    pub internal_energy: f64,
    pub free_energy: f64,
    pub entropy: f64,
    pub heat_capacity: f64,
}

pub fn electronic_contributions(levels: &[(f64, f64)], temp: f64) -> ElectronicContributions {
    let rt = RGAS_AU * temp;
    let (mut q, mut sum_x, mut sum_x2) = (0.0, 0.0, 0.0);
    for &(eps, g) in levels {
        let x = CM1_TO_K * eps / temp;
        let w = g * (-x).exp();
        q += w;
        sum_x += w * x;
        sum_x2 += w * x * x;
    }
    let (mean, mean_sq) = (sum_x / q, sum_x2 / q);
    let internal_energy = rt * mean;
    let free_energy = -rt * q.ln();
    ElectronicContributions {
        partition_function: q,
        internal_energy,
        free_energy,
        entropy: (internal_energy - free_energy) / temp,
        heat_capacity: RGAS_AU * (mean_sq - mean * mean),
    }
}

#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
#[allow(unused_variables)]
#[allow(dead_code)]
impl MoleculeStruct {
    pub fn eval_all_therm_func(&mut self, temp: f64, pressure: f64, freq_cutoff: f64) {
        self.all_electronic(temp);
        self.all_translation(pressure, temp);
        self.all_rotations(temp);
        self.all_vibrations(temp, freq_cutoff);

        self.thermo.utherm =
            self.thermo.uelec + self.thermo.utrans + self.thermo.urot + self.thermo.uvib;

        self.thermo.htherm =
            self.thermo.helec + self.thermo.htrans + self.thermo.hrot + self.thermo.hvib;

        self.thermo.stherm =
            self.thermo.selec + self.thermo.strans + self.thermo.srot + self.thermo.svib;

        self.thermo.ftherm =
            self.thermo.felec + self.thermo.ftrans + self.thermo.frot + self.thermo.fvib;

        self.thermo.gtherm =
            self.thermo.gelec + self.thermo.gtrans + self.thermo.grot + self.thermo.gvib;

        self.thermo.cvtherm =
            self.thermo.cvelec + self.thermo.cvtrans + self.thermo.cvrot + self.thermo.cvvib;

        self.thermo.cptherm =
            self.thermo.cpelec + self.thermo.cptrans + self.thermo.cprot + self.thermo.cpvib;

        self.thermo.utot = self.thermo.utherm + self.dh0 * CM1_TO_HARTREE;
        self.thermo.htot = self.thermo.htherm + self.dh0 * CM1_TO_HARTREE;
        self.thermo.ftot = self.thermo.ftherm + self.dh0 * CM1_TO_HARTREE;
        self.thermo.gtot = self.thermo.gtherm + self.dh0 * CM1_TO_HARTREE;

        self.thermo.stot = self.thermo.stherm;
        self.thermo.cvtot = self.thermo.cvtherm;
        self.thermo.cptot = self.thermo.cptherm;

        self.thermo.pftot =
            self.thermo.pfelec * self.thermo.pftrans * self.thermo.pfrot * self.thermo.pfvib;
    }

    // -----------------------------------------------------------------------------------------
    // Electronic: from the electronic levels (`electronic_contributions`); without levels the ground level alone
    // with degeneracy `multi` (q = multi, U = H = Cv = 0).
    fn all_electronic(&mut self, temp: f64) {
        let ground = [(0.0, if self.multi > 0.0 { self.multi } else { 1.0 })];
        let levels: &[(f64, f64)] = if self.electronic_levels.is_empty() { &ground } else { &self.electronic_levels };
        let e = electronic_contributions(levels, temp);

        self.thermo.pfelec = e.partition_function;
        self.thermo.felec = e.free_energy;
        self.thermo.uelec = e.internal_energy;
        self.thermo.helec = e.internal_energy;
        self.thermo.selec = e.entropy;
        self.thermo.gelec = e.free_energy;

        self.thermo.cvelec = e.heat_capacity;
        self.thermo.cpelec = e.heat_capacity;
    }

    // -----------------------------------------------------------------------------------------
    // Translation (your convention):
    // q0 = lambda^3 * (RT/p), S = R(ln q0 + 5/2), U = (3/2)RT
    // => F = U - TS = -RT ln(q0) - RT
    // Define PFtrans = exp(-F/RT) = e * q0   so that PF is consistent with F.
    fn all_translation(&mut self, pressure: f64, temp: f64) {
        let RT = RGAS_AU * temp;

        let mass = self.mass * AMU_TO_ELECTRON_MASS;

        // lambda_factor = sqrt(2π m RT / h^2), then cube it
        let mut lam = f64::sqrt(TWO_PI * mass * RT / HPLANCK_AU_SQ);
        lam = lam * lam * lam;

        let Vol = RT / (pressure * PASCAL_TO_AU);

        let q0 = lam * Vol;

        // Your extra "-RT" in F is equivalent to PF = e*q0
        let pf = std::f64::consts::E * q0;

        let F = -RT * pf.ln();
        let U = 1.5 * RT;
        let H = 2.5 * RT;
        let S = (U - F) / temp;
        let G = H - temp * S;

        self.thermo.pftrans = pf;
        self.thermo.ftrans = F;
        self.thermo.utrans = U;
        self.thermo.htrans = H;
        self.thermo.strans = S;
        self.thermo.gtrans = G;

        self.thermo.cvtrans = 1.5 * RGAS_AU;
        self.thermo.cptrans = 2.5 * RGAS_AU;
    }

    // -----------------------------------------------------------------------------------------
    // Rotation (rigid rotor, high-T classical):
    // PFrot = q_rot (dimensionless), F = -RT ln PF, U = (dof/2)RT, S = (U-F)/T
    // dof=2 for linear, 3 for nonlinear.
    fn all_rotations(&mut self, temp: f64) {
        let RT = RGAS_AU * temp;

        // If no rotational constants, treat as non-rotating
        if self.brot.is_empty() {
            self.thermo.pfrot = 1.0;
            self.thermo.frot = 0.0;
            self.thermo.urot = 0.0;
            self.thermo.hrot = 0.0;
            self.thermo.srot = 0.0;
            self.thermo.grot = 0.0;
            self.thermo.cvrot = 0.0;
            self.thermo.cprot = 0.0;
            return;
        }

        let sigma = if self.symnum > 0.0 { self.symnum } else { 1.0 };
        let chiral = if self.chiral > 0.0 { self.chiral } else { 1.0 };

        // Heuristic: if any B ~ 0, treat as linear (common QC output: [B,B,0])
        let eps = 1.0e-12;
        let has_zero = self.brot.iter().any(|b| *b <= eps);
        let nonzero: Vec<f64> = self.brot.iter().copied().filter(|b| *b > eps).collect();

        let (pf, dof) = if has_zero || nonzero.len() <= 1 {
            // linear: q = T / (sigma * theta_r) * chiral
            let brot_cm1 = if !nonzero.is_empty() {
                nonzero.iter().sum::<f64>() / (nonzero.len() as f64)
            } else {
                self.brot[0].max(1.0e-6)
            };
            let theta_r = brot_cm1 * CM1_TO_K; // K
            let q = (temp / theta_r) * (chiral / sigma);
            (q, 2.0)
        } else {
            // nonlinear: q = sqrt(pi) * T^(3/2) / sqrt(thetaA thetaB thetaC) * chiral/sigma
            // use first three nonzero constants (order irrelevant in product)
            let a = nonzero[0] * CM1_TO_K;
            let b = nonzero[1] * CM1_TO_K;
            let c = nonzero[2] * CM1_TO_K;
            let denom = f64::sqrt(a * b * c);
            let q = f64::sqrt(PI) * temp.powf(1.5) / denom * (chiral / sigma);
            (q, 3.0)
        };

        let F = -RT * pf.ln();
        let U = 0.5 * dof * RT;
        let H = U;
        let S = (U - F) / temp;
        let G = F;

        self.thermo.pfrot = pf;
        self.thermo.frot = F;
        self.thermo.urot = U;
        self.thermo.hrot = H;
        self.thermo.srot = S;
        self.thermo.grot = G;

        self.thermo.cvrot = 0.5 * dof * RGAS_AU;
        self.thermo.cprot = self.thermo.cvrot;
    }

    // -----------------------------------------------------------------------------------------
    // Vibrations: start from an effective PF per mode.
    // - Every mode: Grimme qRRHO mixing of the entropy (smoothly -> RRHO well above the cutoff),
    //   then F = U - TS and PF = exp(-F/RT)
    // Uses your "no ZPE in thermal vib energy" convention.
    fn all_vibrations(&mut self, temp: f64, freq_cutoff: f64) {
        let RT = RGAS_AU * temp;

        let mut Uvib = 0.0;
        let mut Hvib = 0.0;
        let mut Fvib = 0.0;
        let mut Svib = 0.0;
        let mut Cvib = 0.0;
        let mut PFvib = 1.0;

        for &omega_cm1 in &self.freq {
            if omega_cm1 <= 0.0 {
                continue;
            }

            // x = (hc/kB)*nu / T = (CM1_TO_K * nu_cm^-1)/T  (dimensionless)
            let x = (CM1_TO_K * omega_cm1) / temp;
            let ex = f64::exp(x);

            // Thermal vib energy (no ZPE):
            // U = (hc*nu) / (exp(x)-1)
            let U_mode = omega_cm1 * CM1_TO_HARTREE / (ex - 1.0);
            let H_mode = U_mode;

            // RRHO heat capacity (harmonic), in Hartree/mol/K
            let Cv_mode = RGAS_AU * x * x * ex / ((ex - 1.0) * (ex - 1.0));

            // Entropy: Grimme qRRHO mixing (Chem. Eur. J. 18, 9955 (2012)) for EVERY mode;
            // the damping w -> 1 makes high modes pure RRHO smoothly, and freq_cutoff = 0 gives RRHO.
            // (A hard switch at freq_cutoff belongs to Truhlar's method and makes S jump there.)
            let S_mode = Self::grimme_entropy_qrrho(omega_cm1, freq_cutoff, temp);

            // Define free energy from U and S (this is the qRRHO practice)
            let F_mode = U_mode - temp * S_mode;

            // Effective PF from F
            let PF_mode = f64::exp(-F_mode / RT);

            Uvib += U_mode;
            Hvib += H_mode;
            Svib += S_mode;
            Fvib += F_mode;
            Cvib += Cv_mode;
            PFvib *= PF_mode;
        }

        self.thermo.uvib = Uvib;
        self.thermo.hvib = Hvib;
        self.thermo.svib = Svib;
        self.thermo.fvib = Fvib;
        self.thermo.gvib = Fvib;

        self.thermo.cvvib = Cvib;
        self.thermo.cpvib = Cvib;
        self.thermo.pfvib = PFvib;
    }

    // =========================================================================================
    // Correct Grimme qRRHO entropy mixing (per your reference)

    // RRHO vibrational entropy (Hartree/mol/K), matches:
    // S = R * [ x/(e^x - 1) - ln(1 - e^{-x}) ], x = (hc nu)/(kT)
    fn entropy_vib_rrho(omega_cm1: f64, temp: f64) -> f64 {
        if omega_cm1 <= 0.0 {
            return 0.0;
        }
        let x = (CM1_TO_K * omega_cm1) / temp;
        if x < 1.0e-12 {
            return 0.0;
        }
        let ex = f64::exp(x);
        let term = x / (ex - 1.0) - f64::ln(1.0 - f64::exp(-x));
        RGAS_AU * term
    }

    // Free rotor entropy (Hartree/mol/K), matches reference:
    // mu = h/(8*pi^2*c*nu_tilde)
    // mu' = mu*Bav/(mu+Bav)
    // S = R*(1/2 + 1/2 ln( 8*pi^3*mu'*kT/h^2 ))
    fn entropy_free_rotor(omega_cm1: f64, temp: f64, bav_si: f64) -> f64 {
        if omega_cm1 <= 0.0 {
            return 0.0;
        }

        // omega (cm^-1) -> (m^-1)
        let omega_m1 = omega_cm1 * 100.0;

        // nu (s^-1) = c * omega_m^-1
        let nu_s1 = CLIGHT_SI * omega_m1;

        // mu (kg m^2) = h / (8*pi^2*nu)
        let mu = PLANCK_SI / (8.0 * PI_SQ * nu_s1);

        // mu'
        let mu_prime = mu * bav_si / (mu + bav_si);

        // dimensionless factor inside ln
        let factor = 8.0 * PI.powi(3) * mu_prime * BOLTZMANN_SI * temp / (PLANCK_SI * PLANCK_SI);

        // S in J/mol/K then convert to Hartree/mol/K
        let s_si = RGAS_SI * (0.5 + 0.5 * factor.ln());
        s_si / J_PER_HARTREE_PER_MOL
    }

    // Damping w = 1/(1+(nu_cut/nu)^alpha), alpha=4
    fn grimme_damp(omega_cm1: f64, freq_cutoff_cm1: f64) -> f64 {
        let ratio = freq_cutoff_cm1 / omega_cm1;
        1.0 / (1.0 + ratio.powf(GRIMME_ALPHA))
    }

    // Grimme qRRHO mixed entropy: S = w*S_RRHO + (1-w)*S_free_rotor
    fn grimme_entropy_qrrho(omega_cm1: f64, freq_cutoff_cm1: f64, temp: f64) -> f64 {
        if omega_cm1 <= 0.0 {
            return 0.0;
        }
        let w = Self::grimme_damp(omega_cm1, freq_cutoff_cm1);
        let s_rrho = Self::entropy_vib_rrho(omega_cm1, temp);
        let s_fr = Self::entropy_free_rotor(omega_cm1, temp, GRIMME_BAV_SI);
        w * s_rrho + (1.0 - w) * s_fr
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::{MolType, MoleculeBuilder};

    #[test]
    fn test_eval_all_therm_func_grimme() {
        // - uses a NONZERO Grimme cutoff
        // - also prints per-mode RRHO vs free-rotor vs mixed entropies for inspection

        let name = "Water".to_string();
        let moltype = MolType::mol;

        let temp = 298.15;
        let pressure = 101_325.0;

        // Typical Grimme cutoff ~ 100 cm^-1 (common default in qRRHO practice)
        let freq_cutoff = 100.0;

        let mut water = MoleculeBuilder::new(name, moltype)
            .freq(vec![1626.92, 3761.93, 3876.98]) // cm^-1
            .brot(vec![26.513921, 14.346808, 9.309431]) // cm^-1
            .mass(18.02) // amu
            .ene(-76.37226823 / CM1_TO_HARTREE) // stored in cm^-1 in your struct (per your original test)
            .multi(1.0)
            .chiral(1.0)
            .symnum(1.0)
            .build();

        water.eval_all_therm_func(temp, pressure, freq_cutoff);

        println!("\n============================== INPUT ==============================");
        println!(
            "T = {:.2} K   p = {:.1} Pa   Grimme cutoff = {:.1} cm^-1",
            temp, pressure, freq_cutoff
        );
        println!("symnum:  {:?}", water.symnum);
        println!("multi:   {:?}", water.multi);
        println!("chiral:  {:?}", water.chiral);
        println!("mass:    {:?}", water.mass);
        println!("freq:    {:?}", water.freq);
        println!("brot:    {:?}", water.brot);

        println!("\n==================== ELECTRONIC / ZPE (as stored) ===================");
        println!("E0:  {:15.8} Eh", water.ene * CM1_TO_HARTREE);
        println!("H0:  {:15.8} Eh", water.dh0 * CM1_TO_HARTREE);
        println!(
            "ZPE: {:15.8} Eh   {:12.3} kcal/mol",
            water.zpe * CM1_TO_HARTREE,
            water.zpe * CM1_TO_KCAL
        );

        println!("\n========================= PARTITION FUNCTIONS =======================");
        println!("Q_elec : {:15.6e}", water.thermo.pfelec);
        println!("Q_trans: {:15.6e}", water.thermo.pftrans);
        println!("Q_rot  : {:15.6e}", water.thermo.pfrot);
        println!("Q_vib  : {:15.6e}", water.thermo.pfvib);
        println!("Q_tot  : {:15.6e}", water.thermo.pftot);

        println!("\n====================== CONTRIBUTIONS (Hartree) ======================");
        println!(
            "U_elec : {:15.8}   H_elec : {:15.8}   F_elec : {:15.8}   G_elec : {:15.8}",
            water.thermo.uelec, water.thermo.helec, water.thermo.felec, water.thermo.gelec
        );
        println!(
            "U_trans: {:15.8}   H_trans: {:15.8}   F_trans: {:15.8}   G_trans: {:15.8}",
            water.thermo.utrans, water.thermo.htrans, water.thermo.ftrans, water.thermo.gtrans
        );
        println!(
            "U_rot  : {:15.8}   H_rot  : {:15.8}   F_rot  : {:15.8}   G_rot  : {:15.8}",
            water.thermo.urot, water.thermo.hrot, water.thermo.frot, water.thermo.grot
        );
        println!(
            "U_vib  : {:15.8}   H_vib  : {:15.8}   F_vib  : {:15.8}   G_vib  : {:15.8}",
            water.thermo.uvib, water.thermo.hvib, water.thermo.fvib, water.thermo.gvib
        );

        println!("\n====================== TOTAL THERMAL (Hartree) ======================");
        println!(
            "U_therm: {:15.8}   H_therm: {:15.8}   F_therm: {:15.8}   G_therm: {:15.8}",
            water.thermo.utherm, water.thermo.htherm, water.thermo.ftherm, water.thermo.gtherm
        );
        println!(
            "U_tot  : {:15.8}   H_tot  : {:15.8}   F_tot  : {:15.8}   G_tot  : {:15.8}",
            water.thermo.utot, water.thermo.htot, water.thermo.ftot, water.thermo.gtot
        );

        println!("\n=================== ENTROPY (J/mol/K and S*T) =======================");
        println!(
            "S_elec : {:12.3} J/mol/K   (S*T = {:10.3} kcal/mol)",
            water.thermo.selec * 1000.0 * AU_TO_KJ,
            water.thermo.selec * temp * AU_TO_KCAL
        );
        println!(
            "S_trans: {:12.3} J/mol/K   (S*T = {:10.3} kcal/mol)",
            water.thermo.strans * 1000.0 * AU_TO_KJ,
            water.thermo.strans * temp * AU_TO_KCAL
        );
        println!(
            "S_rot  : {:12.3} J/mol/K   (S*T = {:10.3} kcal/mol)",
            water.thermo.srot * 1000.0 * AU_TO_KJ,
            water.thermo.srot * temp * AU_TO_KCAL
        );
        println!(
            "S_vib  : {:12.3} J/mol/K   (S*T = {:10.3} kcal/mol)",
            water.thermo.svib * 1000.0 * AU_TO_KJ,
            water.thermo.svib * temp * AU_TO_KCAL
        );
        println!(
            "S_tot  : {:12.3} J/mol/K   (S*T = {:10.3} kcal/mol)",
            water.thermo.stot * 1000.0 * AU_TO_KJ,
            water.thermo.stot * temp * AU_TO_KCAL
        );

        println!("\n=================== GRIMME CHECK (per-mode) =========================");
        println!("Mode    nu/cm^-1   w(damp)     S_RRHO(J/mol/K)   S_FR(J/mol/K)   S_mix(J/mol/K)");
        for (i, &nu) in water.freq.iter().enumerate() {
            let w = MoleculeStruct::grimme_damp(nu, freq_cutoff);

            let s_rrho_au = MoleculeStruct::entropy_vib_rrho(nu, temp);
            let s_fr_au = MoleculeStruct::entropy_free_rotor(nu, temp, GRIMME_BAV_SI);
            let s_mix_au = MoleculeStruct::grimme_entropy_qrrho(nu, freq_cutoff, temp);

            let s_rrho = s_rrho_au * 1000.0 * AU_TO_KJ;
            let s_fr = s_fr_au * 1000.0 * AU_TO_KJ;
            let s_mix = s_mix_au * 1000.0 * AU_TO_KJ;

            println!(
                "{:>3}  {:10.2}  {:8.5}      {:12.3}        {:12.3}      {:12.3}",
                i + 1,
                nu,
                w,
                s_rrho,
                s_fr,
                s_mix
            );
        }

        println!("\n============================== SANITY ===============================");
        let kbt = RGAS_AU * temp;
        println!(
            "kB*T = {:12.6} Eh  ({:8.3} kcal/mol)",
            kbt,
            kbt * AU_TO_KCAL
        );

        // Lightweight sanity assertions (won’t overconstrain your conventions)
        assert!(water.thermo.pftot.is_finite() && water.thermo.pftot > 0.0);
        assert!(water.thermo.stot.is_finite());
        assert!(water.thermo.gtot.is_finite());
    }

    #[test]
    fn one_electronic_level_gives_the_degeneracy_factor() {
        // q = g0, F = -RT ln g0, U = H = Cv = 0, S = R ln g0: the former multiplicity model.
        let t = 298.15;
        let e = electronic_contributions(&[(0.0, 3.0)], t);
        assert_eq!(e.partition_function, 3.0);
        assert!((e.free_energy + RGAS_AU * t * 3.0_f64.ln()).abs() < 1e-18);
        assert_eq!((e.internal_energy, e.heat_capacity), (0.0, 0.0));
        assert!((e.entropy - RGAS_AU * 3.0_f64.ln()).abs() < 1e-18);
    }

    #[test]
    fn two_electronic_levels_follow_the_level_sums() {
        // OH X 2Pi: 2Pi_3/2 (g = 2) and 2Pi_1/2 at 139.7 cm-1 (g = 2), at 300 K; x = CM1_TO_K eps / T
        let (t, eps) = (300.0, 139.7);
        let x = CM1_TO_K * eps / t;
        let q = 2.0 + 2.0 * (-x).exp();
        let mean = 2.0 * x * (-x).exp() / q;
        let mean_sq = 2.0 * x * x * (-x).exp() / q;
        let rt = RGAS_AU * t;
        let e = electronic_contributions(&[(0.0, 2.0), (eps, 2.0)], t);
        assert!((e.partition_function - q).abs() < 1e-14);
        assert!((e.internal_energy - rt * mean).abs() < 1e-18);
        assert!((e.free_energy + rt * q.ln()).abs() < 1e-18);
        assert!((e.entropy - (rt * mean + rt * q.ln()) / t).abs() < 1e-18);
        assert!((e.heat_capacity - RGAS_AU * (mean_sq - mean * mean)).abs() < 1e-18);
        // limits: q -> sum g at high T, q -> g0 and U, Cv -> 0 at low T
        let hot = electronic_contributions(&[(0.0, 2.0), (eps, 2.0)], 1.0e7);
        assert!((hot.partition_function - 4.0).abs() < 1e-4);
        let cold = electronic_contributions(&[(0.0, 2.0), (eps, 2.0)], 2.0);
        assert!((cold.partition_function - 2.0).abs() < 1e-20 && cold.internal_energy.abs() < 1e-25 && cold.heat_capacity.abs() < 1e-25);
    }

    #[test]
    fn molecules_with_excited_electronic_levels_sum_them_into_the_thermochemistry() {
        let levels = vec![(0.0, 2.0), (139.7, 2.0)];
        let mut oh = MoleculeBuilder::new("OH".to_string(), MolType::mol)
            .freq(vec![3737.8])
            .brot(vec![18.91])
            .mass(17.0027)
            .ene(0.0)
            .multi(2.0)
            .electronic_levels(levels.clone())
            .build();
        oh.eval_all_therm_func(300.0, 101_325.0, 0.0);
        let e = electronic_contributions(&levels, 300.0);
        assert_eq!((oh.thermo.pfelec, oh.thermo.uelec, oh.thermo.selec, oh.thermo.cvelec), (e.partition_function, e.internal_energy, e.entropy, e.heat_capacity));
        assert_eq!((oh.thermo.helec, oh.thermo.felec, oh.thermo.gelec, oh.thermo.cpelec), (e.internal_energy, e.free_energy, e.free_energy, e.heat_capacity));
        // without levels: the multiplicity alone
        let mut plain = MoleculeBuilder::new("OH".to_string(), MolType::mol).freq(vec![3737.8]).brot(vec![18.91]).mass(17.0027).ene(0.0).multi(2.0).build();
        plain.eval_all_therm_func(300.0, 101_325.0, 0.0);
        assert_eq!(plain.thermo.pfelec, 2.0);
    }

    #[test]
    fn qrrho_entropy_mixes_every_mode() {
        // Grimme, Chem. Eur. J. 18, 9955 (2012): S = w S_RRHO + (1-w) S_free-rotor for every
        // mode, w = 1/(1+(nu0/nu)^4), mu' = mu Bav/(mu + Bav). Reference values evaluated
        // independently from these formulas at T = 298.15 K, nu0 = 100 cm-1, Bav = 1e-44 kg m^2.
        // Just above the cutoff the pure RRHO value would be 14.452781 (100.001) and
        // 11.180604 (150) J/mol/K.
        for (nu, s_ref) in [(100.001, 13.198430), (150.0, 11.028562)] {
            let mut mol = MoleculeBuilder::new("mode".to_string(), MolType::mol)
                .freq(vec![nu])
                .brot(vec![1.0, 0.5, 0.25])
                .mass(18.0)
                .ene(0.0)
                .build();
            mol.eval_all_therm_func(298.15, 101_325.0, 100.0);
            let s_vib = mol.thermo.svib * 1000.0 * AU_TO_KJ; // J/mol/K
            assert!(
                (s_vib - s_ref).abs() < 2.0e-3,
                "nu = {nu}: S_vib = {s_vib:.6} J/mol/K, reference = {s_ref:.6} J/mol/K"
            );
        }
    }

    #[test]
    fn translational_entropy_matches_sackur_tetrode() {
        // Independent reference: Sackur-Tetrode equation in SI units (CODATA 2018 exact h, kB;
        // u = 1.66053906660e-27 kg):  S = R [ ln( (2 pi m kB T / h^2)^{3/2} kB T / p ) + 5/2 ]
        let temp = 298.15;
        let pressure = 101_325.0;
        let mass_amu = 18.02;

        let mut water = MoleculeBuilder::new("Water".to_string(), MolType::mol)
            .freq(vec![1626.92, 3761.93, 3876.98])
            .brot(vec![26.513921, 14.346808, 9.309431])
            .mass(mass_amu)
            .ene(0.0)
            .build();
        water.eval_all_therm_func(temp, pressure, 100.0);
        let s_trans_code = water.thermo.strans * 1000.0 * AU_TO_KJ; // J/mol/K

        let amu_kg = 1.660_539_066_60e-27;
        let m = mass_amu * amu_kg;
        let lambda_inv3 = (TWO_PI * m * BOLTZMANN_SI * temp / (PLANCK_SI * PLANCK_SI)).powf(1.5);
        let volume = BOLTZMANN_SI * temp / pressure;
        let s_trans_ref = RGAS_SI * ((lambda_inv3 * volume).ln() + 2.5);

        assert!(
            (s_trans_code - s_trans_ref).abs() < 0.01,
            "S_trans = {s_trans_code:.4} J/mol/K, Sackur-Tetrode = {s_trans_ref:.4} J/mol/K"
        );
    }
}
