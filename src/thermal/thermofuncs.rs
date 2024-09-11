use crate::molecule::MoleculeStruct;
const H_PLANCK: f64 = 3.3356E-11;
const TWOPI : f64 = 2.0 * (std::f64::consts::PI);
const RGAS : f64 = 8.314;

#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
#[allow(unused_variables)]
#[allow(dead_code)]
impl MoleculeStruct {

    pub fn eval_all_thermo_func(&mut self, temp: f64, pressure: f64, freq_cutoff: f64){
        let RT = RGAS * temp;
        self.all_electronic(temp);
        self.all_translation(pressure, temp);
        //self.all_rotations(temp);
        //self.all_vibrations(temp, freq_cutoff);
        //all_hinderedrotors(temp);

        self.thermo.utot = self.thermo.uelec + self.thermo.utrans + self.thermo.urot + self.thermo.uvib;// + self.termo.uhindrot;
        self.thermo.htot = self.thermo.helec + self.thermo.htrans + self.thermo.hrot + self.thermo.hvib;// + self.termo.hhindrot;
        self.thermo.stot = self.thermo.selec + self.thermo.strans + self.thermo.srot + self.thermo.svib;// + self.termo.shindrot;
        self.thermo.ftot = self.thermo.felec + self.thermo.ftrans + self.thermo.frot + self.thermo.fvib;// + self.termo.fhindrot;
        self.thermo.gtot = self.thermo.gelec + self.thermo.gtrans + self.thermo.grot + self.thermo.gvib;// + self.termo.ghindrot;

        self.thermo.cvtot = self.thermo.cvelec + self.thermo.cvtrans + self.thermo.cvrot + self.thermo.cvvib; // + self.termo.cvhindrot;
        self.thermo.cptot = self.thermo.cpelec + self.thermo.cptrans + self.thermo.cprot + self.thermo.cpvib; // + self.termo.cphindrot;
        
        self.thermo.pftot = self.thermo.pfelec * self.thermo.pftrans * self.thermo.pfrot * self.thermo.pfvib; // * self.termo.pfhindrot;
    } 

    fn all_electronic(&mut self, temp: f64){
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

    fn all_translation(&mut self, pressure: f64, temp: f64){
        let RT = RGAS * temp;

        let mut lam = f64::sqrt(TWOPI * self.totmass * RT /(H_PLANCK * H_PLANCK));

        lam = lam * lam * lam;

        let Vol = RT / pressure;

        //Each of the is defined for a 3D object
        self.thermo.utrans  = 1.5e0 * RT;
        self.thermo.htrans  = 2.5e0 * RT;
        //entropy here is not simply single molecule entropy. N! term for N atoms are incoportated
        //into translation, that is why the F not equal with 1 particle -RT*ln(Z)
        self.thermo.strans  =  RGAS * f64::ln(lam * Vol) + 2.5e0 * RGAS;
        self.thermo.ftrans  = -RT * f64::ln(lam * Vol) - RT; //F = U - TS
        self.thermo.gtrans  = -RT * f64::ln(lam * Vol); 
        self.thermo.pftrans =  lam * Vol; 
        self.thermo.cvtrans  = 1.5e0 * RGAS;
        self.thermo.cptrans  = 2.5e0 * RGAS;

    }


    fn get_full_internal_energy(&mut self) {
        // Example calculation for the internal energy
        self.thermo.utot = self.ene + self.zpe;
        println!("Internal energy calculated: {}", self.thermo.utot);
    }



}


//=============================================================================================
fn partfunc_nd_harmosc(RT: f64, nvib: usize, omega: &[f64]) -> f64 {
    const H_PLANCK: f64 = 3.3356E-11;

    let mut pf_vib = 0.0;
    for i in 0..nvib {
        pf_vib *= 1.0 / (1.0 - f64::exp(1.0 - H_PLANCK * omega[i] / RT));
    }

    return pf_vib;
}
//=============================================================================================

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
fn partfunc_3d_translation(RT: f64, masstot: f64, pressure: f64) -> f64{

   const H_PLANCK: f64 = 3.3356E-11;
   const TWOPI : f64 = 2.0 * (std::f64::consts::PI);

   let mut lam = f64::sqrt(TWOPI * masstot * RT /(H_PLANCK * H_PLANCK));

   lam = lam * lam * lam;

   let mut Vol = RT/pressure;

   let pf_trans = Vol * lam;

   return pf_trans;

}
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

*/
