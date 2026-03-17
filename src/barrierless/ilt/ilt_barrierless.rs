// All energies are in cm-1.
// Output w_e is also in cm-1 
//
// Assumptions:
// - rho(E) or rho_r(E) returns density of states in states / cm-1
// - For the bimolecular case, the translational/spin prefactor must already be
//   supplied in a form consistent with the same cm-1 convention.

pub mod ilt_barrierless {
    /// Boltzmann constant in cm^-1 K^-1
    pub const KB_CM: f64 = 0.695_034_76;

    /// Unimolecular barrierless ILT expression
    ///
    /// w_e(E) = k(E) * rho(E)
    ///
    /// Formula:
    /// w_e(E) = (A * beta0^n / Gamma(n))
    ///          * ∫_{Ea}^{E} rho(E - tau) * (tau - Ea)^(n-1) d tau
    ///
    /// Units:
    /// - e, ea, t0 -> cm-1, cm-1, K
    /// - rho(x) -> states / cm-1
    /// - result w_e -> cm-1 
    pub fn w_uni<F, G>(
        e: f64,
        a: f64,
        ea: f64,
        t0: f64,
        n: f64) -> f64
    {
        if e <= ea {
            return 0.0;
        }

        let beta0 = 1.0 / (KB_CM * t0);
        let prefactor = a * beta0.powf(n) / gamma_fn(n);

        let n_steps = 4000usize;
        let h = (e - ea) / n_steps as f64;

        let mut sum = 0.0;

        for i in 0..=n_steps {
            let x = i as f64 * h; // x = tau - Ea
            let tau = ea + x;

            let weight = if i == 0 || i == n_steps { 0.5 } else { 1.0 };

            let kernel = if x == 0.0 {
                if n > 1.0 { 0.0 } else { 0.0 }
            } else {
                x.powf(n - 1.0)
            };

            let val = rho(e - tau) * kernel;
            sum += weight * val;
        }

        prefactor * sum * h
    }

    /// Bimolecular barrierless ILT expression
    ///
    /// w_e(E) = k(E) * rho(E)
    ///
    /// Formula:
    /// w_e(E) = bimol_prefactor_cm
    ///          * (A * beta0^n / Gamma(n + 1.5))
    ///          * ∫_{Ea + delta_h0}^{E}
    ///              rho_r(E - tau) * (tau - Ea - delta_h0)^(n + 0.5) d tau
    ///
    /// Important:
    /// bimol_prefactor_cm must already include the translational and degeneracy factor
    /// in a unit system consistent with cm-1 energies:
    ///
    ///     bimol_prefactor_cm = ((2πμ/h^2)^(3/2)) * (g_a g_b / g_c)
    ///
    /// but converted beforehand to the same wavenumber-based convention we use elsewhere.
    ///
    /// Units:
    /// - e, ea, delta_h0 -> cm-1
    /// - t0 -> K
    /// - rho_r(x) -> states / cm-1
    /// - result w_e -> cm-1
    pub fn w_bimol<F, G>(
        e: f64,
        a: f64,
        ea: f64,
        delta_h0: f64,
        t0: f64,
        n: f64,
        bimol_prefactor_cm: f64) -> f64
    {
        let eth = ea + delta_h0;

        if e <= eth {
            return 0.0;
        }

        let beta0 = 1.0 / (KB_CM * t0);
        let prefactor = bimol_prefactor_cm * a * beta0.powf(n) / gamma_fn(n + 1.5);

        let n_steps = 4000usize;
        let h = (e - eth) / n_steps as f64;

        let mut sum = 0.0;

        for i in 0..=n_steps {
            let x = i as f64 * h; // x = tau - Ea - delta_h0
            let tau = eth + x;

            let weight = if i == 0 || i == n_steps { 0.5 } else { 1.0 };

            let kernel = x.powf(n + 0.5);
            let val = rho_r(e - tau) * kernel;

            sum += weight * val;
        }

        prefactor * sum * h
    }
}
