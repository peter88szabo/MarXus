/// Locate the centrifugal barrier (local maximum) of
///     V_eff(r) = V(r) + J(J+1) / (2 * mu * r^2)
/// on r in [r_lo, r_hi] (Bohr). Atomic units throughout.
///
/// Inputs:
/// - J: angular momentum quantum number
/// - mu: reduced mass (a.u.)
/// - V: callable that returns (V(r), dV/dr) in (Hartree, Hartree/Bohr) for r in Bohr
/// - D_e: dissociation energy (Hartree). Used only when J == 0.
/// - r_lo, r_hi: bracketing interval (Bohr)
///
/// Output:
/// - Some((r_max, V_max)) if a centrifugal barrier exists in the interval
/// - None if V_eff is barrierless (no local maximum) in the interval
///
/// Behavior:
/// - If J == 0: returns Some((r_hi, D_e))
/// - No panics for missing barrier; returns None
pub fn centrifugal_barrier_bisect<F>(
    J: i32,
    mu: f64,
    V: F,
    D_e: f64,
    r_lo: f64,
    r_hi: f64,
) -> Option<(f64, f64)>
where
    F: Fn(f64) -> (f64, f64),
{
    assert!(mu > 0.0, "mu must be positive");
    assert!(r_lo > 0.0 && r_hi > r_lo, "invalid r interval");

    if J == 0 {
        return Some((r_hi, D_e));
    }

    let jj = (J as f64) * ((J + 1) as f64);

    // d/dr [ jj/(2*mu*r^2) ] = - jj/(mu*r^3)
    let dVeff_dr = |r: f64| -> f64 {
        let (_v, dv) = V(r);
        dv - jj / (mu * r * r * r)
    };

    let veff = |r: f64| -> f64 {
        let (v, _dv) = V(r);
        v + jj / (2.0 * mu * r * r)
    };

    // Scan to bracket a maximum: derivative crosses from + to -
    let n_scan = 4000usize;
    let mut r_prev = r_lo;
    let mut f_prev = dVeff_dr(r_prev);

    let mut a = 0.0;
    let mut b = 0.0;
    let mut found = false;

    for i in 1..=n_scan {
        let t = i as f64 / n_scan as f64;
        let r = r_lo + t * (r_hi - r_lo);
        let f = dVeff_dr(r);

        if f_prev.is_finite() && f.is_finite() && f_prev > 0.0 && f < 0.0 {
            a = r_prev;
            b = r;
            found = true;
            break;
        }

        r_prev = r;
        f_prev = f;
    }

    if !found {
        return None;
    }

    // Bisection on dVeff/dr = 0 within [a,b]
    let fa0 = dVeff_dr(a);
    let fb0 = dVeff_dr(b);
    if !(fa0.is_finite() && fb0.is_finite() && fa0 > 0.0 && fb0 < 0.0) {
        return None;
    }

    let tol_r = 1e-8;
    let tol_f = 1e-10;
    let max_iter = 200;

    for _ in 0..max_iter {
        let m = 0.5 * (a + b);
        let fm = dVeff_dr(m);

        if (b - a) < tol_r || fm.abs() < tol_f {
            return Some((m, veff(m)));
        }

        if fm < 0.0 {
            b = m;
        } else {
            a = m;
        }
    }

    None
}
