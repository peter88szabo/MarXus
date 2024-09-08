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
       pf_vib *= 1.0 / (1.0 - f64::exp(1.0 - H_PLANCK * omega[i] / RT);
   }

   return pf_vib;
}
//=============================================================================================

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

   let crt = f64::powf(8.0 * PI_SQ * RT / (H_PLANCK_SQ),  nrot as f64) 

   let const_rho = crt/gamma_func(rdim);


   res[i] = const_W*(f64::powf(Ei, rdim)),


   let mut pf_rot = 0.0;
   for i in 0..nvib {
       pf_rot *= 1.0 / (1.0 - f64::exp(1.0 - H_PLANCK*omega[i]/RT);
   }

   return pf_rot;

}
//=============================================================================================

//=============================================================================================
fn partfunc_3d_translation(RT: f64, masstot: f64, pressure: f64) -> f64{

   const H_PLANCK: f64 = 3.3356E-11;
   const TWOPI : f64 = 2.0 * (std::f64::consts::PI); 
 
   let mut lam = f64::powf(2.0 * TWOPI*TWOPI * masstot*RT/(H_PLANCK*H_PLANCK),  nrot as f64);
   lam = lam*lam*lam

   let mut Vol = RT/pressure

   return pf_trans;

}
//=============================================================================================







//=====================================================================================================================
use std::time::{Duration, Instant};
fn format_duration(duration: Duration) -> String {
    let seconds = duration.as_secs();
    let hours = seconds / 3600;
    let minutes = (seconds % 3600) / 60;
    let seconds = seconds % 60;

    format!("{:02}:{:02}:{:02}", hours, minutes, seconds)
}
//=====================================================================================================================



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fn main() {
   let start_time = Instant::now();

 // Page 92 Holbrook Pilling book:
   let omega_t = [260.0, 260.0, 364.0, 667.0, 760.0, 760.0, 1205.0, 1205.0, 3033.0];
   let nvib_t = omega_t.len();
   let ZPE_t = 0.5 * omega_t.iter().sum::<f64>();

   let Brot_t = [];
   let nrot_t = 0 as usize; 

   let dE =  300.0; //cm-1
   let Emax = 40000.0; // cm-1
   let nebin = (Emax / dE + 0.5) as usize;


   let mut freq_bin_t = vec![0; nvib_t];
   for i in 0..nvib_t {
       freq_bin_t[i] = (omega_t[i] / dE + 0.5) as usize;
   }

   let mut WE_t = get_rovib_WE_or_rhoE("sum".to_string(), nvib_t, nebin, dE, nrot_t, &freq_bin_t, &Brot_t);

   let mut Ev = 83.5934*50.92; //59.2 kJmol -- > cm-1
   println!("{:>12.3} {:>12.3}", Ev, ZPE_t);

   use std::fs;
   use std::io::Write;
   let mut file = fs::File::create("test_WEv.dat").unwrap();
   for i in 1..=nebin {
       let ene = (i as f64 -0.5) * dE;
       writeln!(&mut file, "{:>10.2} {:>15.3e}", ene/ZPE_t,  WE_t[i]).unwrap();
   }


 // Page 111 Holbrook Pilling book:
   let omega_t = [733.0, 1018.0, 1018.0, 1355.0, 1488.0, 1488.0, 2968.0, 3044.0, 3044.0];
   let nvib_t = omega_t.len();
   let ZPE_t = 0.5 * omega_t.iter().sum::<f64>();

   let Brot_t = [5.10, 0.433, 0.433];
   let nrot_t = Brot_t.len() as usize;

   let dE =  10.0; //cm-1
   let Emax = 10000.0; // cm-1
   let nebin = (Emax / dE + 0.5) as usize;


   let mut freq_bin_t = vec![0; nvib_t];
   for i in 0..nvib_t {
       freq_bin_t[i] = (omega_t[i] / dE + 0.5) as usize;
   }

   let mut WE_t = get_rovib_WE_or_rhoE("den".to_string(), nvib_t, nebin, dE, nrot_t, &freq_bin_t, &Brot_t);

    let mut file = fs::File::create("test_rhoEVR.dat").unwrap();
    for i in 1..=nebin {
        let ene = (i as f64 -0.5) * dE;
        writeln!(&mut file, "{:>10.2} {:>15.3}", ene,  WE_t[i].log10()).unwrap();
    }






  //==================================================
  // Input for complex and TS structures to run RRKM
  //==================================================
  // Energy grain
  //--------------------------------------------------
    let dE =  20.0; //cm-1
    let Emax = 20000.0; // cm-1
    let nebin = (Emax / dE + 0.5) as usize; 

    println!();          
    println!("ΔE (energy grain):             {:>10.1} cm-1", dE);
    println!("Emax:                          {:>10.1} cm-1 \n", Emax);
  //--------------------------------------------------
  // Vibrations
  //--------------------------------------------------
    let omega_cpx = [600.0, 1000.0, 1500.0];
    let omega_ts = [600.0, 1200.0, 1500.0];
    let nvib_cpx = omega_cpx.len();
    let nvib_ts = omega_ts.len();

    //use std::iter::Sum;
    let ZPE_cpx = 0.5 * omega_cpx.iter().sum::<f64>(); //0.5*(omega_cpx.iter().sum());
    let ZPE_ts = 0.5 * omega_ts.iter().sum::<f64>(); //0.5*(omega_cpx.iter().sum());

    println!("ZPE of complex:                {:>10.1} cm-1", ZPE_cpx);
    println!("ZPE of TS:                     {:>10.1} cm-1 \n", ZPE_ts );
  //--------------------------------------------------
  // Rotations
  //--------------------------------------------------
    let Brot_cpx = [600.0, 1200.0, 1500.0];
    let Brot_ts = [600.0, 1200.0, 1500.0];
    let nrot_cpx = Brot_cpx.len();
    let nrot_ts = Brot_ts.len();
 //--------------------------------------------------
  // E0 = reaction energy (only electronic no ZPE)
  // ΔH0 = 0K heat of formation (including ZPEs)
  //
  // ΔH0 = E0 + ZPE_ts - ZPE_cpx
  //--------------------------------------------------
    let Ezero = 4100.0; // in cm-1
    let dH0 = Ezero + ZPE_ts - ZPE_cpx;
    println!("E0  (reaction energy at 0K):   {:>10.1} cm-1  (no ZPE only pure electronic)", Ezero);
    println!("ΔH0 (reaction enthalpy at 0K): {:>10.1} cm-1  (ΔH0 = E0 + ZPE_ts - ZPE_cpx)\n", dH0);

    if Ezero >= Emax {
        println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        println!{"Code stoped beacuse Emax is too small comapred to reaction energy (E0)"};
        println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        panic!{"!!!! Increase Emax !!!!"};
    }
 //--------------------------------------------------
 // Symmetry numbers
 //--------------------------------------------------
    let sigma_cpx = 1.0;
    let sigma_ts = 1.0;
    println!("symm_num_cpx:                   {:>10.0}", sigma_cpx);
    println!("symm_num_ts:                    {:>10.0} \n", sigma_ts);
  //==================================================
  // End of Input
  //==================================================

  //Minimum energy (including ZPE) of the reaction as integer energy bin
   let nbin_dH0 = (dH0/dE + 0.5) as usize;

  // Compute RRKM rate constant
    let kRRKM = get_kE(nebin, dE, nvib_ts, nvib_cpx, &omega_ts, &omega_cpx, nrot_ts, nrot_cpx, &Brot_ts, &Brot_cpx, sigma_ts, sigma_cpx, dH0);

  //Print rate constant k(E)
    //use std::fs;
    //use std::io::Write;
    let mut file = fs::File::create("microcanonical_rate.dat").unwrap();
    //for i in 1..=nebin {
    for i in nbin_dH0..=nebin {
        let ene = i as f64 * dE;
        writeln!(&mut file, "{:>10.1} {:>15.6e}", ene, kRRKM[i]).unwrap();
    }


   let end_time = Instant::now();

   let duration = end_time - start_time;
   let duration_seconds = duration.as_secs_f64();



   println!("--------------------------------------------------------");
   println!("RRKM calculation finished \n");
   println!("Data written into file: microcanonical_rate.dat \n");

   println!("Computational time: {:>15.4} sec", duration_seconds);
   println!("Computational time:          {} hh:mm:ss\n", format_duration(duration));
   
}

