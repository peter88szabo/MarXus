struct Molecule {
    natom: usize,         // number of atoms
    //nvib: usize           // number vibrational modes
    ene: f64,             // electronic energy (no ZPE)
    //zpe: f64,             // zero-point energy
    //dh0: f64,             // ene + ZPE
    //sigma: f64,           // symmetry number
    //chiral: f64,          // number of enantiomers
    //imag_freq: f64,       // imaginary frequency for TS
    brot: Vec<f64>,       // rotational constants
    //freq: Vec<f64>,       // vibrational freqs
    //omega: Vec<f64>,      // vibrational freqs in energy unit (together with hbar)
    //we: Vec<f64>,         // sum of states W(E)
    //rhoe: Vec<f64>,       // denisty of states rho(E)
    //coord: Vec<Vec<f64>>, // coordinates of the equilbrium geometry 
}

/*impl Molecule {
    fn valami(&mut self){
        self.dh = self.ene + self.zpe
    }

    fn valami() --> Self {
    }
}
*/

fn main() {
    let nmol = 2;
    let mut molec: Vec<Molecule> = Vec::with_capacity(nmol);
    let vector_length = molec.len();
    println!("Length of molec: {}", vector_length);

    for i in 0..nmol {
        let g = i as f64;
        let this = Molecule {
            natom: 2 + i,             
            ene: 1.2 + g,              
            brot: vec![2.0 + g, 1.0 + g, 4.0 + g],    
        };
        molec.push(this);
    }

    let vector_length = molec.len();
    println!("Length of molec: {}", vector_length);

    // Accessing fields of the first instance in the vector
    println!("natom: {}", molec[0].natom);
    println!("ene: {}", molec[0].ene);
    println!("brot: {:?}", molec[0].brot);

    println!("natom: {}", molec[1].natom);
    println!("ene: {}", molec[1].ene);
    println!("brot: {:?}", molec[1].brot);


}


