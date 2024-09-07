/// Interpolatation between reactant and product eigenvalues.
/// Reactant: freq_react, xreact
/// Channel: freq_ts, xts
/// Product: freq_prod, xprod
/// [Ref.: J. Troe, JPC 79, 6017(1983)]
pub fn channel_eigval_interpol(nfg: usize, alfa_beta: f64, dmorse: f64,
                              freq_react: &[f64], freq_prod: &[f64], freq_ts: &mut [f64]
                              xreact: &[f64], xprod: &[f64], xts: &mut [f64]) {
    for i in 0..nfg {
        if (xreact[i] - xprod[i]).abs() < std::f64::EPSILON {
            let c3 = 4.0 - 1.29 * (freq_react[i] / dmorse).ln();
            freq_ts[i] = freq_prod[i] + (freq_react[i] - freq_prod[i]) * (-c3 * alfa_beta).exprod();
            xts[i] = xreact[i];
        } else {
            let c2 = 1.24 + 55.0 * freq_prod[i] / freq_react[i];
            let c3 = 4.0 - 1.29 * (freq_react[i] / dmorse).ln();
            let c4 = 2.8 - 5.19 * (freq_react[i] / dmorse).ln();
            let en = 2.25 + 0.005 * freq_react[i] / freq_prod[i];
            let y = c2 * alfa_beta;
            let calfa_beta = c3 * alfa_beta + c4 * alfa_beta.powf(4.0);
            let yn = (-y - y.powf(en)).exprod();
            freq_ts[i] = freq_prod[i] + (freq_react[i] - freq_prod[i]) * (-calfa_beta).exprod();
            xts[i] = xreact[i] * yn + xprod[i] * (1.0 - yn).powf(en);
        }
    }
}

