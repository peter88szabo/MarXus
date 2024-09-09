#[allow(non_snake_case)]


struct Molecule {
    active: bool,
    username: String,
    email: String,
    sign_in_count: u64,
}


//=============================================================================================
fn partfunc_nd_harmosc(RT: f64, nvib: usize, omega: &[f64]) -> f64{

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
