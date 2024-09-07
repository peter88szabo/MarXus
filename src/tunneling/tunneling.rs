use std::f64::consts::PI;

pub fn wigner(beta: f64, omega: f64) -> f64 {
    let u = omega * beta;

    return 1.0 + u * u / 24.0;
}

pub fn bell(beta: f64, omega: f64) -> f64 {
    let u = omega * beta;

    return (0.5 * u / f64::sin(0.5 * u)).abs();
}

pub fn skodje_truhlar(beta: f64, omega: f64, v0: f64) -> f64 {
    let alpha = 2.0 * PI / omega;

    if beta > alpha {
        // When omega is large and low T
        let kappa = (((beta - alpha) * v0).exp() - 1.0) * beta / (beta - alpha);

        //println!("Skodje, beta*omega: {}", beta * omega);
        //println!("beta*omega > twopi");

        return kappa;
    } else if beta < alpha {
        // When omega is small and high T
        let dum0 = PI * beta / alpha;
        let dum1 = dum0 / dum0.sin();
        let dum2 = ((beta - alpha) * v0).exp() * beta / (beta - alpha);
        let kappa = dum1 + dum2;

        //println!("Skodje, beta*omega: {}", beta * omega);
        //println!("beta*omega < twopi");

        return kappa;
    } else {
        panic!("Wrong mode-argument in Skodje_Truhlar routine");
    }
}

pub fn skodje_truhlar_exact(beta: f64, omega: f64, v0: f64) -> f64 {
    let alpha = 2.0 * PI / omega;

    const NMAX: usize = 100;

    let mut res = 0.0;

    for n in 0..=NMAX {
        let numerator = 1.0 - ((beta - (n + 1) as f64 * alpha).exp() * v0);
        let denom = (n + 1) as f64 * alpha - beta;
        let dum = 1.0 / (n as f64 * alpha + beta);
        let mut alter = 1.0;

        if n % 2 == 1 {
            alter = -1.0;
        }

        res += alter * beta * (numerator / denom + dum);
    }

    return res;
}

pub fn eckart(beta: f64, omega: f64, vf: f64, vb: f64, de: f64, emax: f64) -> (f64, f64) {
    let n_emax = (emax / de) as usize;
    let mut to_int1 = vec![0.0; n_emax];
    let mut to_int2 = vec![0.0; n_emax];

    let alpha1 = 2.0 * PI * vf / omega;
    let alpha2 = 2.0 * PI * vb / omega;

    let a = (vf.sqrt() + vb.sqrt()).powi(2);
    let b = vf - vb;
    let d = ( (b*b - a*a).powi(2) / (a*a*a) / 8.0).sqrt() /  omega;

    let e0 = 0.0;
    for i in 0..n_emax {
        let e = e0 + (i as f64) * de;
        let csi = e / vf;

        let ptun1 = tunprop1(alpha1, alpha2, csi);
        let ptun2 = tunprop2(a, b, d, e);

        to_int1[i] = (-e * beta).exp() * ptun1;
        to_int2[i] = (-e * beta).exp() * ptun2;

    }
    let n_emin = 0;

    //let mut res1 = 0.0;
    //for i in n_emin+1..n_emax {
    //    res1 += 0.5*(to_int1[i]+to_int1[i-1])*de;
    //}

    let res1 = simpson_integrate(&to_int1, n_emin, n_emax, de);
    let kappa1 = res1 * (beta * vf).exp() * beta;

    let res2 = simpson_integrate(&to_int2, n_emin, n_emax, de);
    let kappa2 = res2 * (beta * vf).exp() * beta;

    return (kappa1, kappa2);
}

pub fn tunprop1(alpha1: f64, alpha2: f64, csi: f64) -> f64 {
    let denom = 1.0 / alpha1.sqrt() + 1.0 / alpha2.sqrt();
    let twopi_a = 2.0 * (alpha1 * csi).sqrt() / denom;
    let twopi_b = 2.0 * ((alpha1 * csi) + (alpha1 - alpha2).abs()).sqrt() / denom;
    let twopi_d = 2.0 * ((alpha1 * alpha2) - (PI * PI / 4.0)).sqrt();


    let csh_amb = (twopi_a - twopi_b).cosh();
    let csh_apb = (twopi_a + twopi_b).cosh();
    let csh_d = twopi_d.cosh();

    let res = 1.0 - (csh_amb + csh_d) / (csh_apb + csh_d);
    return res;
}

pub fn tunprop2(a: f64, b: f64, d: f64, e: f64) -> f64 {
    let twopi_d = 2.0 * PI * d;

    let alpha = twopi_d * (2.0 * e).sqrt();
    let beta  = twopi_d * (2.0 * (e - b)).sqrt();
    let delta = twopi_d * ((2.0 * a) - (1.0 / (2.0 * d * d))).sqrt();

    let csh_amb = (alpha - beta).cosh();
    let csh_apb = (alpha + beta).cosh();
    let csh_d = delta.cosh();

    let res = 1.0 - (csh_amb + csh_d) / (csh_apb + csh_d);
    return res;
}

fn simpson_integrate(func: &[f64], nmin: usize, nmax: usize, step: f64) -> f64 {
    let mut s0 = 0.0;
    let mut s1 = 0.0;
    let mut s2 = 0.0;

    if nmin > nmax {
        panic!("Wrong boundaries in Simpson Integration: nmin or nmax");
    }

    let ndata: usize = nmax - nmin;

    for i in (nmin..nmax-2).step_by(2) {
         s1 += func[i];
         s0 += func[i+1];
         s2 += func[i + 2];
     }

     //println!("sm1 =  {}, s0 = {}, sp1 = {}", s1, s0, s2);
     let mut res = step * (s1 + 4.0 * s0 + s2) / 3.0;
     //println!("step =  {}, res = {}", step, res);

     // If n is even, add the last slice separately
     if ndata % 2 == 0 {
        res += step * (5.0 * func[nmax-1] + 8.0 * func[nmax - 2] - func[nmax - 3]) / 12.0;
     }

    return res;

}
