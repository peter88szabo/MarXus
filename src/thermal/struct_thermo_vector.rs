#[allow(non_snake_case)]
#[derive(Debug)]
enum MolTypeEnum{
    cpx,
    ts,
    bimol
}


#[derive(Debug)]
struct MoleculeStruct{
    natom:     i32, // number of atoms
    nvib:      i32, // number of vibrational modes
    nrot:      i32, // number of rotational modes
    nlin:      i32, // linear (nlin=1) or not (nlin=0)
    zpe:       f64, // zero-point energy
    ezero:     f64, // electronic energy (without ZPE)
    dh0:       f64, // ezero + zpe
    symnum:    f64, // rotational symm num
    chiral:    f64, // number of eniantiomers
    spinmulti: f64, // spin multiplicity
    totmass:   f64, // total mass 
    name:      String, // name of the species
    freq:      Vec<f64>, // harmonic frequencies
    brot:      Vec<f64>, // rotational constants
    mass:      Vec<f64>, // rotational constants
    qxyz:      Vec<f64>, // Cartesian coords
    atom:      Vec<String>, // atomic symbols
    we:        Vec<f64>, // sum of states, micorcanonical W(E)
    rhoe:      Vec<f64>, // density of states, microcanonical rho(E)
    moltype:   MolTypeEnum, // what kind of species (well, TS, react, prod)
    freq_imag: f64, // imag freq of TS mode
    vfor:      f64, // forward barrier for Eckart
    vback:     f64, // backward barrier for Eckart 
    tunprop:   Vec<f64>, // energy-dependent tunneling probability
    kappa:     f64, // tunneling corr
}


impl MoleculeStruct {
    fn new(name: String, moltype: MolTypeEnum) -> MoleculeStruct {
        MoleculeStruct {
            natom: 0, nvib: 0, nrot: 0, nlin: 0,
            zpe: 0.0, ezero: 0.0, dh0: 0.0,
            symnum: 0.0, chiral: 0.0, multi: 0.0, totmass: 0.0, name,
            freq: Vec::new(), brot: Vec::new(), qxyz: Vec::new(), atom: Vec::new(),
            we: Vec::new(), rhoe: Vec::new(),
            moltype,
            freq_imag: 0.0, vfor: 0.0, vback: 0.0, kappa: 0.0,
        }
    }
}


fn main() {
    // Define the properties of the water molecule
    let natom = 3; // H2O has 3 atoms
    let nvib = 3; // 3 vibrational modes for H2O
    let nrot = 3; // 3 rotational modes for a non-linear molecule like H2O
    let nlin = 0; // Water is not a linear molecule
    let zpe = 0.5; // Hypothetical value for zero-point energy
    let ezero = 10.0; // Hypothetical value for electronic energy
    let dh0 = zpe + ezero;
    let symnum = 2.0; // Hypothetical value for symmetry number
    let chiral = 1.0; // Water has 1 enantiomer
    let multi = 1.0; // Hypothetical value for spin multiplicity
    let totmass = 18.01528; // Approximate mass of H2O in atomic mass units
    let name = String::from("Water");
    let moltype = MolTypeEnum::Complex; // Assuming 'Molecule' is a variant of MolType_enum

    // Initialize the water molecule
    let oxygen_molecule = MoleculeStruct::new("oxygen".to_string(), MolTypeEnum::Complex);

    let water_molecule = MoleculeStruct {
        natom, nvib, nrot, nlin,
        zpe, ezero, dh0,
        symnum, chiral, multi, totmass, name,
        freq: vec![440.0, 1600.0, 3600.0], // Hypothetical vibrational frequencies for H2O
        brot: vec![20.0, 20.0, 20.0], // Hypothetical rotational constants for H2O
        qxyz: vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0], // Hypothetical coordinates for H2O
        atom: vec![String::from("O"), String::from("H"), String::from("H")], // Atomic symbols
        we: vec![], 
        rhoe: vec![],
        moltype,
        freq_imag: 0.0, vfor: 0.0, vback: 0.0, kappa: 0.0, 
    };

    // Printing the water molecule details
    println!("{:?}", water_molecule);
    println!();
    println!("{:?}", oxygen_molecule);
}
