//=====================================================================================
// Lanczos recursive approximation of Gamma function
// Tested agains Fortran intrinsic DGAMMA function
//=====================================================================================
pub fn gamma_func(a: f64) -> f64 {
    const PI: f64 = 3.14159265358979324;
    const CG: usize = 7;

    // These numbers are taken from the sample code in Wikipedia
    // and the sample itself takes them from the GNU Scientific Library
    const P: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];

    let mut x = a;
    let gamm: f64;

    if x < 0.5 {
        gamm = PI / (f64::sin(x*PI) * gamma_func(1.0 - x));
    } else {
        x = x - 1.0;
        let mut t = P[0];
        for i in 1..CG + 2 {
            t = t + P[i] / (x + i as f64);
        }
        let w = x + (CG as f64) + 0.5;
        gamm = f64::sqrt(2.0 * PI) * f64::powf(w,x + 0.5) * (f64::exp(-w)) * t;
    }

    gamm
}
//=============================================================================================

