#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
#[allow(unused_variables)]
#[allow(dead_code)]

#[derive(Debug)]
pub enum MolType {
    mol,
    ts,
    bimol,
    reactant,
}
#[derive(Debug)]
#[allow(unused_variables)]
#[allow(dead_code)]
pub struct ThermoStruct {
    //partition functions
    pub pftot: f64,
    pub pfelec: f64,
    pub pftrans: f64,
    pub pfrot: f64,
    pub pfvib: f64,
    //pub pfhindrot: f64,

    //entropy functions
    pub stot: f64,
    pub stherm: f64,
    pub selec: f64,
    pub strans: f64,
    pub srot: f64,
    pub svib: f64,
    //pub shindrot: f64,

    //internal energy functions
    pub utot: f64,
    pub utherm: f64,
    pub uelec: f64,
    pub utrans: f64,
    pub urot: f64,
    pub uvib: f64,
    //pub uhindrot: f64,

    //enthalpy functions
    pub htot: f64,
    pub htherm: f64,
    pub helec: f64,
    pub htrans: f64,
    pub hrot: f64,
    pub hvib: f64,
    //pub hhindrot: f64,

    //Helmholtz free energy functions
    pub ftot: f64,
    pub ftherm: f64,
    pub felec: f64,
    pub ftrans: f64,
    pub frot: f64,
    pub fvib: f64,
    //pub fhindrot: f64,

    //Gibbs free energy functions
    pub gtot: f64,
    pub gtherm: f64,
    pub gelec: f64,
    pub gtrans: f64,
    pub grot: f64,
    pub gvib: f64,
    //pub ghindrot: f64,

    //Cv heat capacity
    pub cvtot: f64,
    pub cvelec: f64,
    pub cvtrans: f64,
    pub cvrot: f64,
    pub cvvib: f64,
    //pub cvhindrot: f64,

    //Cp heat capacity
    pub cptot: f64,
    pub cpelec: f64,
    pub cptrans: f64,
    pub cprot: f64,
    pub cpvib: f64,
    //pub cphindrot: f64,

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
            stherm: 0.0,
            selec: 0.0,
            strans: 0.0,
            srot: 0.0,
            svib: 0.0,
            utot: 0.0,
            utherm: 0.0,
            uelec: 0.0,
            utrans: 0.0,
            urot: 0.0,
            uvib: 0.0,
            htot: 0.0,
            htherm: 0.0,
            helec: 0.0,
            htrans: 0.0,
            hrot: 0.0,
            hvib: 0.0,
            ftot: 0.0,
            ftherm: 0.0,
            felec: 0.0,
            ftrans: 0.0,
            frot: 0.0,
            fvib: 0.0,
            gtot: 0.0,
            gtherm: 0.0,
            gelec: 0.0,
            gtrans: 0.0,
            grot: 0.0,
            gvib: 0.0,
            cvtot: 0.0,
            cvelec: 0.0,
            cvtrans: 0.0,
            cvrot: 0.0,
            cvvib: 0.0,
            cptot: 0.0,
            cpelec: 0.0,
            cptrans: 0.0,
            cprot: 0.0,
            cpvib: 0.0,

        }
    }
}

#[derive(Debug)]
#[allow(unused_variables)]
#[allow(dead_code)]
pub struct Tunneling {
    pub freq_imag: f64,   // imag freq of TS mode
    pub vfor: f64,        // forward barrier for Eckart
    pub vback: f64,       // backward barrier for Eckart
    pub tunprop: Vec<f64>, // energy-dependent tunneling probability
    pub kappa: f64,       // tunneling correction
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
pub struct MoleculeStruct {
    pub name: String,          // name of the species
    pub nvib: u32,             // number of vibrational modes
    pub nrot: u32,             // number of rotational modes
    pub nlin: bool,            // linear (nlin=1) or not (nlin=0)
    pub zpe: f64,              // zero-point energy
    pub ene: f64,              // electronic energy (without ZPE)
    pub dh0: f64,              // ene + zpe
    pub symnum: f64,           // rotational symm num
    pub chiral: f64,           // number of enantiomers
    pub multi: f64,            // degeneracy factor, usually spin multiplicity
    pub totmass: f64,          // total mass
    pub freq: Vec<f64>,        // harmonic frequencies
    pub brot: Vec<f64>,        // rotational constants
    pub we: Vec<f64>,          // sum of states
    pub rhoe: Vec<f64>,        // density of states
    pub moltype: MolType,  // molecule type
    pub tunnel: Tunneling, // tunneling information
    pub thermo: ThermoStruct,  // thermodynamic properties
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
pub struct MoleculeBuilder {
    pub name: String,
    pub moltype: MolType,
    pub nlin: bool,
    pub mass: Vec<f64>,
    pub freq: Vec<f64>,
    pub brot: Vec<f64>,
    pub qxyz: Vec<f64>,
    pub ene:  Option<f64>,
    pub dh0:  Option<f64>,
    pub multi: f64,
    pub chiral: f64,
    pub symnum: f64,
}

//use crate::inertia::inertia::get_brot;
#[allow(unused_variables)]
#[allow(dead_code)]
impl MoleculeBuilder {
    pub fn new(name: String, moltype: MolType) -> Self {
        MoleculeBuilder {
            name,
            moltype,
            nlin: false,  // Default value for nlin
            mass: Vec::new(),
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

    pub fn nlin(mut self, nlin: bool) -> Self {
        self.nlin = nlin;
        self
    }

    pub fn freq(mut self, freq: Vec<f64>) -> Self {
        self.freq = freq;
        self
    }

    pub fn brot(mut self, brot: Vec<f64>) -> Self {
        self.brot = brot;
        self
    }

    pub fn qxyz(mut self, qxyz: Vec<f64>) -> Self {
        self.qxyz = qxyz;
        self
    }

    pub fn mass(mut self, mass: Vec<f64>) -> Self {
        self.mass = mass;
        self
    }

    pub fn ene(mut self, ene: f64) -> Self {
        self.ene = Some(ene);
        self
    }

    pub fn dh0(mut self, dh0: f64) -> Self {
        self.dh0 = Some(dh0);
        self
    }

    pub fn multi(mut self, multi: f64) -> Self {
        self.multi = multi;
        self
    }

    pub fn chiral(mut self, chiral: f64) -> Self {
        self.chiral = chiral;
        self
    }

    pub fn symnum(mut self, symnum: f64) -> Self {
        self.symnum = symnum;
        self
    }

    pub fn build(self) -> MoleculeStruct {
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
            vec![0.0, 0.0, 0.0] //get_brot(&self.qxyz, &self.mass)
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

/*
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
*/

