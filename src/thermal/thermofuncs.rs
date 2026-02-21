// thermofuncs.rs
use crate::molecule::MoleculeStruct;

const PI: f64 = std::f64::consts::PI;
const PI_SQ: f64 = PI * PI;
const TWOPI: f64 = 2.0 * PI;

const CLIGHT: f64 = 137.035999074;
const AMU_TO_ELECMASS: f64 = 1836.15267343;

const HPLANCK_AU: f64 = TWOPI;
const HPLANCK_AU_SQ: f64 = TWOPI * TWOPI;

const AU_TO_KJ: f64 = 2625.5;
const AU_TO_KCAL: f64 = 627.51;

const RGAS_AU: f64 = 8.31446261815324 / 1000.0 / AU_TO_KJ; // Hartree/mol/K
const CM1_TO_K: f64 = 1.43877;
const CM1_TO_HARTREE: f64 = 4.55635e-6;
const CM1_TO_KCAL: f64 = 2.85914e-3;
const PASCAL_TO_AU: f64 = 1.0e-13 / 2.9421912;

// -----------------------------
// SI constants needed ONLY for Grimme free-rotor entropy (dimensionless inside ln)
// -----------------------------
const PLANCK_SI: f64 = 6.62607015e-34; // J*s
const BOLTZMANN_SI: f64 = 1.380649e-23; // J/K
const RGAS_SI: f64 = 8.31446261815324; // J/mol/K
const CLIGHT_SI: f64 = 2.99792458e8; // m/s

// 1 Hartree/mol = AU_TO_KJ kJ/mol = AU_TO_KJ*1000 J/mol
const J_PER_HARTREE_PER_MOL: f64 = AU_TO_KJ * 1000.0;

// Grimme default average moment of inertia (kg*m^2)
const GRIMME_BAV_SI: f64 = 1.0e-44;

// Damping exponent in Grimme qRRHO
const GRIMME_ALPHA: f64 = 4.0;

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
    // Electronic: PF = multiplicity, F = -RT ln PF, U=H=0, S=(U-F)/T
    fn all_electronic(&mut self, temp: f64) {
        let RT = RGAS_AU * temp;

        let pf = if self.multi > 0.0 { self.multi } else { 1.0 };
        let F = -RT * pf.ln();
        let U = 0.0;
        let H = 0.0;
        let S = (U - F) / temp;

        self.thermo.pfelec = pf;
        self.thermo.felec = F;
        self.thermo.uelec = U;
        self.thermo.helec = H;
        self.thermo.selec = S;
        self.thermo.gelec = F;

        self.thermo.cvelec = 0.0;
        self.thermo.cpelec = 0.0;
    }

    // -----------------------------------------------------------------------------------------
    // Translation (your convention):
    // q0 = lambda^3 * (RT/p), S = R(ln q0 + 5/2), U = (3/2)RT
    // => F = U - TS = -RT ln(q0) - RT
    // Define PFtrans = exp(-F/RT) = e * q0   so that PF is consistent with F.
    fn all_translation(&mut self, pressure: f64, temp: f64) {
        let RT = RGAS_AU * temp;

        let mass = self.mass * AMU_TO_ELECMASS;

        // lambda_factor = sqrt(2π m RT / h^2), then cube it
        let mut lam = f64::sqrt(TWOPI * mass * RT / HPLANCK_AU_SQ);
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
    // - For omega > cutoff: pure RRHO -> PF_mode = 1/(1-exp(-x))
    // - For omega <= cutoff: Grimme qRRHO mixing affects entropy -> define F = U - TS, PF=exp(-F/RT)
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

            // Entropy: RRHO or Grimme-mixed
            let S_mode = if omega_cm1 > freq_cutoff {
                Self::entropy_vib_rrho(omega_cm1, temp)
            } else {
                Self::grimme_entropy_qrrho(omega_cm1, freq_cutoff, temp)
            };

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
}
