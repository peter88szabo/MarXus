use crate::molecule::MoleculeStruct;

const HPLANCK: f64 = 3.3356E-11;
const HPLANCK_SQ: f64 = HPLANCK * HPLANCK;
const PI: f64 = std::f64::consts::PI;
const TWOPI: f64 = 2.0 * PI;
const PI_SQ: f64 = PI * PI;
const RGAS: f64 = 8.314;

#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
#[allow(unused_variables)]
#[allow(dead_code)]
impl MoleculeStruct {
    pub fn eval_all_thermo_func(&mut self, temp: f64, pressure: f64, freq_cutoff: f64) {
        self.all_electronic(temp);
        self.all_translation(pressure, temp);
        self.all_rotations(temp);
        self.all_vibrations(temp, freq_cutoff);
        //all_hinderedrotors(temp);

        self.thermo.utherm =
            self.thermo.uelec + self.thermo.utrans + self.thermo.urot + self.thermo.uvib; // + self.termo.uhindrot;
                                                                                          //
        self.thermo.htherm =
            self.thermo.helec + self.thermo.htrans + self.thermo.hrot + self.thermo.hvib; // + self.termo.hhindrot;
                                                                                          //
        self.thermo.stherm =
            self.thermo.selec + self.thermo.strans + self.thermo.srot + self.thermo.svib; // + self.termo.shindrot;
                                                                                          //
        self.thermo.ftherm =
            self.thermo.felec + self.thermo.ftrans + self.thermo.frot + self.thermo.fvib; // + self.termo.fhindrot;
                                                                                          //
        self.thermo.gtherm =
            self.thermo.gelec + self.thermo.gtrans + self.thermo.grot + self.thermo.gvib; // + self.termo.ghindrot;
                                                                                          //
        self.thermo.utot = self.thermo.utherm + self.dh0;
        self.thermo.htot = self.thermo.htherm + self.dh0;
        self.thermo.ftot = self.thermo.ftherm + self.dh0;
        self.thermo.gtot = self.thermo.gtherm + self.dh0;
        self.thermo.stot = self.thermo.stherm;

        self.thermo.cvtot =
            self.thermo.cvelec + self.thermo.cvtrans + self.thermo.cvrot + self.thermo.cvvib; // + self.termo.cvhindrot;
        self.thermo.cptot =
            self.thermo.cpelec + self.thermo.cptrans + self.thermo.cprot + self.thermo.cpvib; // + self.termo.cphindrot;

        self.thermo.pftot =
            self.thermo.pfelec * self.thermo.pftrans * self.thermo.pfrot * self.thermo.pfvib;
        // * self.termo.pfhindrot;
    }

    fn all_electronic(&mut self, temp: f64) {
        //When only one electronic state is available
        //at a given temp, then the partition function
        //approx Qele = multiplicity

        self.thermo.uelec = 0.0e0;
        self.thermo.helec = 0.0e0;
        self.thermo.selec = RGAS * f64::ln(self.multi);
        self.thermo.felec = -RGAS * temp * f64::ln(self.multi);
        self.thermo.gelec = self.thermo.felec;
        self.thermo.pfelec = self.multi;
        self.thermo.cvelec = 0.0e0;
        self.thermo.cpelec = 0.0e0;
    }

    fn all_translation(&mut self, pressure: f64, temp: f64) {
        let RT = RGAS * temp;
        let mut lam = f64::sqrt(TWOPI * self.totmass * RT / (HPLANCK * HPLANCK));

        lam = lam * lam * lam;

        let Vol = RT / pressure;

        //Each of the is defined for a 3D object
        self.thermo.utrans = 1.5e0 * RT;
        self.thermo.htrans = 2.5e0 * RT;
        //entropy here is not simply single molecule entropy. N! term for N atoms are incoportated
        //into translation, that is why the F not equal with 1 particle -RT*ln(Z)
        self.thermo.strans = RGAS * f64::ln(lam * Vol) + 2.5e0 * RGAS;
        self.thermo.ftrans = -RT * f64::ln(lam * Vol) - RT; //F = U - TS
        self.thermo.gtrans = -RT * f64::ln(lam * Vol);
        self.thermo.pftrans = lam * Vol;
        self.thermo.cvtrans = 1.5e0 * RGAS;
        self.thermo.cptrans = 2.5e0 * RGAS;
    }

    fn all_rotations(&mut self, temp: f64) {
        let RT = RGAS * temp;

        let mut Urot = 0.0e0;
        let mut Hrot = 0.0e0;
        let mut Frot = 0.0e0;
        let mut Grot = 0.0e0;
        let mut Srot = 0.0e0;
        let mut Crot = 0.0e0;
        let mut PFrot = 1.0e0;

        for &brot in &self.brot {
            let intert = 1.0 / brot;
            let theta = (8.0 * PI_SQ * brot * RT / HPLANCK_SQ).sqrt();

            let Urot_mode = 0.5 * RGAS * temp;
            let Hrot_mode = Urot_mode;
            let Srot_mode = RGAS * theta.ln() + 0.5 * RGAS;
            let Frot_mode = Urot - temp * Srot; // -Rgas*Temp*log(Theta/sigma)
            let Grot_mode = Frot;
            let PFrot_mode = Frot;
            let Crot_mode = Frot;

            Urot += Urot_mode;
            Hrot += Hrot_mode;
            Srot += Srot_mode;
            Frot += Frot_mode;
            Grot += Grot_mode;
            Crot += Crot_mode; //Cv = Cp, hence we have only a single expression here
            PFrot *= PFrot_mode;
        }
    }

    fn all_vibrations(&mut self, temp: f64, freq_cutoff: f64) {
        let RT = RGAS * temp;

        let mut Uvib = 0.0e0;
        let mut Hvib = 0.0e0;
        let mut Fvib = 0.0e0;
        let mut Gvib = 0.0e0;
        let mut Svib = 0.0e0;
        let mut Cvib = 0.0e0;
        let mut PFvib = 1.0e0;

        for &omega in &self.freq {
            let theta = omega / RGAS;

            let Uvib_mode = RGAS * theta / (f64::exp(theta / temp) - 1.0e0);
            let Hvib_mode = Uvib_mode;

            let Svib_mode: f64;

            if omega < freq_cutoff {
                Svib_mode = Self::entropy_vib(omega, temp);
            } else {
                Svib_mode = Self::Grimme_entropy(omega, freq_cutoff, temp);
            }

            let Fvib_mode = Uvib_mode - temp * Svib_mode;
            let Gvib_mode = Hvib_mode - temp * Svib_mode;
            let PFvib_mode = f64::exp(-Fvib_mode / RT);

            let Cvib_mode = RGAS * theta * theta * f64::exp(theta / temp)
                / ((temp * (f64::exp(theta / temp) - 1.0e0)).powi(2));

            Uvib += Uvib_mode;
            Hvib += Hvib_mode;
            Svib += Svib_mode;
            Fvib += Fvib_mode;
            Gvib += Gvib_mode;
            Cvib += Cvib_mode; //Cv = Cp, hence we have only a single expression here
            PFvib *= PFvib_mode;
        }

        self.thermo.uvib = Uvib;
        self.thermo.hvib = Hvib;
        self.thermo.svib = Svib;
        self.thermo.fvib = Fvib;
        self.thermo.gvib = Gvib;
        self.thermo.cvvib = Cvib;
        self.thermo.cpvib = Cvib;
        self.thermo.pfvib = PFvib;
    }

    fn Grimme_entropy(omega: f64, freq_cutoff: f64, temp: f64) -> f64 {
        let Bav = 1.0;
        //Effective inertia of rotation with same period as the low-frew vibration mode
        let mu: f64 = HPLANCK / (8.0 * PI * PI * omega);

        let mu_prime: f64 = mu * Bav / (mu + Bav);

        //Weight function to smoothly switch between vib --> rot
        let dum = (freq_cutoff / omega).powi(4);
        let wgt = 1.0 / (1.0 + dum);

        let svib =
            wgt * Self::entropy_vib(omega, temp) + (1.0 - wgt) * Self::entropy_rot(mu_prime, temp);

        return svib;
    }

    fn entropy_rot(mu: f64, temp: f64) -> f64 {
        let dum = (((8.0 * PI * PI * PI * mu * RGAS * temp).sqrt()).ln()) / (HPLANCK * HPLANCK);

        let srot = RGAS * (0.5 + dum);

        return srot;
    }

    fn entropy_vib(omega: f64, temp: f64) -> f64 {
        let theta = HPLANCK * omega / (RGAS * temp);

        let dum = theta / (f64::exp(theta) - 1.0e0) - f64::ln(1.0e0 - f64::exp(-theta));

        let svib = RGAS * dum;

        return svib;
    }
}

/*

//=============================================================================================
fn partfunc_nd_rot_classic(RT: f64, nrot: usize, Brot: &[f64]) -> f64{

   const H_PLANCK: f64 = 3.3356E-11;
   const H_PLANCK_SQ: f64 = 3.3356E-11 * 3.3356E-11;
   const PI_VAL : f64 = (std::f64::consts::PI);
   const PI_SQ : f64 = (std::f64::consts::PI) * (std::f64::consts::PI);

// Product of rotational constants
   let mut prod_Brot=1.0;
   for i in 0..nrot{
      prod_Brot *= Brot[i];
   }
   prod_Brot = prod_Brot.sqrt();

   let rdim = (nrot as f64) /2.0;

   let crt = f64::powf(8.0 * PI_SQ * RT / (H_PLANCK_SQ),  nrot as f64);

   let const_rho = 0.0;//crt/gamma_func(rdim);


   res[i] = const_W*(f64::powf(Ei, rdim));


   let mut pf_rot = 0.0;
   for i in 0..nvib {
       pf_rot *= 1.0 / (1.0 - f64::exp(1.0 - H_PLANCK/RT));
   }

   return pf_rot;

}
//=============================================================================================

//=============================================================================================

fn partfunc_electronic(RT: f64, nstate: u32, ene_exc: &[f64, nstate], degeneracy: &[f64, nstate]) -> f64{

   //if ene_exc[0] != 0.0 { panic!};

   //let mut pf_elec = degeneracy[0];
   //for i in 1..nstate {
   //    pf_elec +=  degeneracy[i] * f64::exp(ene_exc[i]/RT);
   //}
   let pf_elec = 0.0;

   return pf_elec;
}
//=============================================================================================


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        assert_eq!(add(2, 3), 5);
    }
}
*/
