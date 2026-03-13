/// Simpified version of phase space theory (PST):
/// Ref: Troe & Ushakov, J. Phys. Chem A, Vol 110, 6732-6741, (2006)
/// and the other refereinces therein
///
/// W0(E,J) is the number of open channels for a primitive PST where
/// we ignore the presence of the centrifugal barrier and we
/// consider only angular momentum conservation rules for
/// the speratred reactants based on their composition:
///    - atom + linear
///    - atom + spherical top
///    - linear + linear
///    - linear + spherical top
///    - spherical top + spherical top
/// For such cases the value of W0(E,J) is known in analytical form
///
/// Then the total number of open channels can be approximated
/// in presence of centrifugal barrier as:
///
/// W_PST(E,J) = W0(E,J) * wcorr(E,J)
///
/// where the wcorr(E,J) encodes the approximate correction for the presence
/// of the centrifugal barrier at angular momentum J.
/// wcorr(E,J) also known in simple analytical form for the above-mentioned
/// simple collisions partners.


use crate::barrierless::pst::collision_type::{CollisionType, ReactantType};

//----------------------------------------------------------------------------------
// Correction factor for primitive W0(E,J) to get W_PST(E,J) = W0(E,J) * wcorr(E,J)
// see formula 2.10 and 2.11 in Ref. Troe & Ushakov 2006
//----------------------------------------------------------------------------------
pub fn wcorr(colltype: CollisionType, ene: f64, j: usize, vmax: &[f64]) -> f64 {
    let ntype = match colltype {
        CollisionType::AtomLinear => 1.0,
        CollisionType::AtomSphericalTop => 1.5,
        CollisionType::LinearLinear => 2.0,
        CollisionType::LinearSphericalTop => 2.5,
        CollisionType::SphericalTopSphericalTop => 3.0,
         _ => panic!("unsupported collision type"),
    };

    let base = 1.0 - vmax[j] / ene;

    let mut wcorr = 0.0;

    if base > 0.0 {
        wcorr = base.powf(ntype);
    }

    wcorr
}
//----------------------------------------------------------------------------------





#[cfg(test)]
mod tests {

    use super::*;

    use crate::barrierless::pst::locate_centrifugal_barrier::centrifugal_barrier_bisect;
    use crate::barrierless::pst::morsepot::morse_value_and_derivative;

    // ---------------------------------------------------------------
    // Build vmax(J) vector automatically
    // ---------------------------------------------------------------
    fn build_vmax_vector(
        j_max: i32,
        mu: f64,
        de: f64,
        beta: f64,
        re: f64,
        r_lo: f64,
        r_hi: f64,
    ) -> (Vec<i32>, Vec<f64>) {

        let mut j_valid = Vec::new();
        let mut vmax = Vec::new();

        for j in 0..=j_max {

            let pot = |r: f64| morse_value_and_derivative(r, de, beta, re);

            if let Some((_rmax, vmax_j)) =
                centrifugal_barrier_bisect(j, mu, pot, de, r_lo, r_hi)
            {
                j_valid.push(j);
                vmax.push(vmax_j);
            }
            else {
                break;
            }
        }

        (j_valid, vmax)
    }

    // ---------------------------------------------------------------
    // Test wcorr(E,J)
    // ---------------------------------------------------------------
    #[test]
    fn test_wcorr() {

        let colltype = CollisionType::LinearLinear;

        // Example parameters
        let de = 0.02;
        let beta = 1.2;
        let re = 3.0;
        let mu = 1000.0;

        let r_lo = 2.0;
        let r_hi = 20.0;

        let j_max = 200;

        let (j_valid, vmax) =
            build_vmax_vector(j_max, mu, de, beta, re, r_lo, r_hi);

        println!("Valid J values: {:?}", j_valid);

        let energies = vec![
            0.005,
            0.010,
            0.050,
            0.100,
            0.500,
            1.200
        ];

        println!("\n==== wcorr(E,J) table ====");

        for e in energies {

            println!("\nE = {:.6}", e);

            for (i, j) in j_valid.iter().enumerate() {

                let w = wcorr(colltype, e, i, &vmax);

                println!(
                    "J = {:3}   vmax = {:12.6e}   wcorr = {:12.6e}",
                    j,
                    vmax[i],
                    w
                );
            }
        }
    }
}

