#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
#[allow(unused_variables)]
#[allow(dead_code)]

#[derive(Debug)]
enum MolType {
    mol,
    ts,
    bimol,
    reactant,
}
#[derive(Debug)]
#[allow(unused_variables)]
#[allow(dead_code)]
struct ThermoStruct {
    //partition functions
    pftot: f64,
    pfelec: f64,
    pftrans: f64,
    pfrot: f64,
    pfvib: f64,

    //entropy functions
    stot: f64,
    selec: f64,
    strans: f64,
    srot: f64,
    svib: f64,

    //internal energy functions
    utot: f64,
    uelec: f64,
    utrans: f64,
    urot: f64,
    uvib: f64,

    //enthalpy functions
    htot: f64,
    helec: f64,
    htrans: f64,
    hrot: f64,
    hvib: f64,

    //Helmholtz free energy functions
    ftot: f64,
    felec: f64,
    ftrans: f64,
    frot: f64,
    fvib: f64,

    //Gibbs free energy functions
    gtot: f64,
    gelec: f64,
    gtrans: f64,
    grot: f64,
    gvib: f64,
}

#[allow(unused_variables)]
#[allow(dead_code)]
impl Default for ThermoStruct {
    fn default() -> Self {
        ThermoStruct {
            pftot: 0.0,
            pfelec: 0.0,
            pftrans: 0.0,
            pfrot: 0.0,
            pfvib: 0.0,
            stot: 0.0,
            selec: 0.0,
            strans: 0.0,
            srot: 0.0,
            svib: 0.0,
            utot: 0.0,
            uelec: 0.0,
            utrans: 0.0,
            urot: 0.0,
            uvib: 0.0,
            htot: 0.0,
            helec: 0.0,
            htrans: 0.0,
            hrot: 0.0,
            hvib: 0.0,
            ftot: 0.0,
            felec: 0.0,
            ftrans: 0.0,
            frot: 0.0,
            fvib: 0.0,
            gtot: 0.0,
            gelec: 0.0,
            gtrans: 0.0,
            grot: 0.0,
            gvib: 0.0,
        }
    }
}

#[derive(Debug)]
#[allow(unused_variables)]
#[allow(dead_code)]
struct Tunneling {
    freq_imag: f64,   // imag freq of TS mode
    vfor: f64,        // forward barrier for Eckart
    vback: f64,       // backward barrier for Eckart
    tunprop: Vec<f64>, // energy-dependent tunneling probability
    kappa: f64,       // tunneling correction
}

#[allow(unused_variables)]
#[allow(dead_code)]
impl Default for Tunneling {
    fn default() -> Self {
        Tunneling {
            freq_imag: 0.0,
            vfor: 0.0,
            vback: 0.0,
            tunprop: Vec::new(),
            kappa: 0.0,
        }
    }
}

#[derive(Debug)]
#[allow(unused_variables)]
#[allow(dead_code)]
struct MoleculeStruct {
    name: String,          // name of the species
    nvib: u32,             // number of vibrational modes
    nrot: u32,             // number of rotational modes
    nlin: bool,            // linear (nlin=1) or not (nlin=0)
    zpe: f64,              // zero-point energy
    ene: f64,              // electronic energy (without ZPE)
    dh0: f64,              // ene + zpe
    symnum: f64,           // rotational symm num
    chiral: f64,           // number of enantiomers
    multi: f64,            // degeneracy factor, usually spin multiplicity
    totmass: f64,          // total mass
    freq: Vec<f64>,        // harmonic frequencies
    brot: Vec<f64>,        // rotational constants
    we: Vec<f64>,          // sum of states
    rhoe: Vec<f64>,        // density of states
    moltype: MolType,  // molecule type
    tunnel: Tunneling, // tunneling information
    thermo: ThermoStruct,  // thermodynamic properties
}

#[allow(unused_variables)]
#[allow(dead_code)]
impl Default for MoleculeStruct {
    fn default() -> Self {
        MoleculeStruct {
            name: String::new(),
            nvib: 0,
            nrot: 0,
            nlin: false,
            zpe: 0.0,
            ene: 0.0,
            dh0: 0.0,
            symnum: 0.0,
            chiral: 0.0,
            multi: 0.0,
            totmass: 0.0,
            freq: Vec::new(),
            brot: Vec::new(),
            we: Vec::new(),
            rhoe: Vec::new(),
            moltype: MolType::mol,
            tunnel: Tunneling::default(),
            thermo: ThermoStruct::default(),
        }
    }
}

#[allow(unused_variables)]
#[allow(dead_code)]
struct MoleculeBuilder {
    name: String,
    moltype: MolType,
    nlin: bool,
    freq: Vec<f64>,
    brot: Vec<f64>,
    qxyz: Vec<f64>,
    ene:  Option<f64>,
    dh0:  Option<f64>,
    multi: f64,
    chiral: f64,
    symnum: f64,
}

#[allow(unused_variables)]
#[allow(dead_code)]
impl MoleculeBuilder {
    fn new(name: String, moltype: MolType) -> Self {
        MoleculeBuilder {
            name,
            moltype,
            nlin: false,  // Default value for nlin
            freq: Vec::new(),
            brot: Vec::new(),
            qxyz: Vec::new(),
            ene: None,    // None means the energy is not provided yet
            dh0: None,    // None means the enthalpy is not provided yet
            multi: 1.0,   // Default multiplicity
            chiral: 1.0,   // Default multiplicity
            symnum: 1.0,  // Default symmetry number
        }
    }

    fn nlin(mut self, nlin: bool) -> Self {
        self.nlin = nlin;
        self
    }

    fn freq(mut self, freq: Vec<f64>) -> Self {
        self.freq = freq;
        self
    }

    fn brot(mut self, brot: Vec<f64>) -> Self {
        self.brot = brot;
        self
    }

    fn qxyz(mut self, qxyz: Vec<f64>) -> Self {
        self.qxyz = qxyz;
        self
    }

    fn ene(mut self, ene: f64) -> Self {
        self.ene = Some(ene);
        self
    }

    fn dh0(mut self, dh0: f64) -> Self {
        self.dh0 = Some(dh0);
        self
    }

    fn multi(mut self, multi: f64) -> Self {
        self.multi = multi;
        self
    }

    fn chiral(mut self, chiral: f64) -> Self {
        self.chiral = chiral;
        self
    }

    fn symnum(mut self, symnum: f64) -> Self {
        self.symnum = symnum;
        self
    }

    fn build(self) -> MoleculeStruct {
        // Check if `freq`, `brot` or `qxyz`, `ene` or `dh0` are `None`
        if self.freq.is_empty() {
            panic!("\n Error: Frequency vector (freq) must be provided.\n");
        }
        if self.brot.is_empty() && self.qxyz.is_empty(){
            panic!("\n Error: Rotational constants (brot) or geometry (qxyz) must be provided.\n");
        }
        if self.ene.is_none() && self.dh0.is_none() {
            panic!("\n Error: Either energy (ene) or enthalpy (dh0) must be provided.\n");
        }

        // zero-point energy:
        let zpe: f64 = self.freq.iter().sum::<f64>() / 2.0;

        // Calculate ene or dh0 based on which is provided
        let ene = match self.ene {
            Some(ene_val) => ene_val,
            None => self.dh0.unwrap() - zpe, // Calculate ene if dh0 is provided 
        };

        let dh0 = match self.dh0 {
            Some(dh0_val) => dh0_val,
            None => self.ene.unwrap() + zpe, // Calculate dh0 if ene is provided
        };

        let brot = if !self.brot.is_empty() {
            self.brot.clone() // Use the provided brot if it's not empty
        } else {
            vec![4.0, 4.0, 4.0] // inertia::get_rotconst(self.qxyz)
        };

         return MoleculeStruct {
            name: self.name,
            nvib: self.freq.len() as u32,
            nrot: self.brot.len() as u32,
            nlin: self.nlin,
            zpe,
            ene,
            dh0,
            symnum: self.symnum,
            multi: self.multi,
            chiral: self.chiral,
            freq: self.freq,
            brot: self.brot, 
            ..Default::default() // Use defaults for other fields
        };
    }
}

fn main() {

    let name = "Water".to_string();
    let moltype = MolType::mol;

    let water = MoleculeBuilder::new(name, moltype)
        .freq(vec![440.0, 1600.0, 3600.0])
        .brot(vec![10.0, 10.0, 20.0])
        .dh0(199.9)
        .multi(3.0)
        .chiral(22.0)
        .symnum(6.0)
        .build();


    println!("ene is provided");
    println!("symnum: {:?}", water.symnum);
    println!("multi: {:?}", water.multi);
    println!("chiral: {:?}", water.chiral);
    println!("zpe: {:?}", water.zpe);
    println!("ene: {:?}", water.ene);  // Output: 199.9
    println!("dh0: {:?}", water.dh0);  // Output: ene + zpe
    println!("freq: {:?}", water.freq);
    println!("brot: {:?}", water.brot);
    
    println!("\nThermo");
    println!("gtot: {:?}", water.thermo.gtot);
    println!("srot: {:?}", water.thermo.srot);
    println!("hvib: {:?}", water.thermo.hvib);
    println!("pfvib: {:?}", water.thermo.pfvib);
                                       


                                       

}

