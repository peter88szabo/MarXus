/// Interpolation between reactant and product eigenvalues.
/// Reactant: freq_react, x_react
/// Channel: freq_ts, x_ts
/// Product: freq_prod, x_prod
/// [Ref.: J. Troe, JPC 79, 6017 (1983)]
#[derive(Debug, Clone)]
pub struct ChannelEigenvalues {
    pub freq_ts: Vec<f64>,
    pub x_ts: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct SacmInterpolationInput {
    pub alfa_beta: f64,
    pub dmorse: f64,
    pub freq_react: Vec<f64>,
    pub freq_prod: Vec<f64>,
    pub x_react: Vec<f64>,
    pub x_prod: Vec<f64>,
}

pub fn channel_eigval_interpol(
    alfa_beta: f64,
    dmorse: f64,
    freq_react: &[f64],
    freq_prod: &[f64],
    x_react: &[f64],
    x_prod: &[f64],
) -> ChannelEigenvalues {
    let nfg = freq_react
        .len()
        .min(freq_prod.len())
        .min(x_react.len())
        .min(x_prod.len());
    let mut freq_ts = vec![0.0; nfg];
    let mut x_ts = vec![0.0; nfg];

    for i in 0..nfg {
        if (x_react[i] - x_prod[i]).abs() < std::f64::EPSILON {
            let c3 = 4.0 - 1.29 * (freq_react[i] / dmorse).ln();
            freq_ts[i] = freq_prod[i] + (freq_react[i] - freq_prod[i]) * (-c3 * alfa_beta).exp();
            x_ts[i] = x_react[i];
        } else {
            let c2 = 1.24 + 55.0 * freq_prod[i] / freq_react[i];
            let c3 = 4.0 - 1.29 * (freq_react[i] / dmorse).ln();
            let c4 = 2.8 - 5.19 * (freq_react[i] / dmorse).ln();
            let en = 2.25 + 0.005 * freq_react[i] / freq_prod[i];
            let y = c2 * alfa_beta;
            let calfa_beta = c3 * alfa_beta + c4 * alfa_beta.powf(4.0);
            let yn = (-y - y.powf(en)).exp();
            freq_ts[i] = freq_prod[i] + (freq_react[i] - freq_prod[i]) * (-calfa_beta).exp();
            x_ts[i] = x_react[i] * yn + x_prod[i] * (1.0 - yn).powf(en);
        }
    }

    ChannelEigenvalues { freq_ts, x_ts }
}

pub fn channel_eigval_interpol_from_input(input: &SacmInterpolationInput) -> ChannelEigenvalues {
    channel_eigval_interpol(
        input.alfa_beta,
        input.dmorse,
        &input.freq_react,
        &input.freq_prod,
        &input.x_react,
        &input.x_prod,
    )
}

pub fn format_correlation_table(
    freq_react: &[f64],
    x_react: &[f64],
    freq_ts: &[f64],
    x_ts: &[f64],
    freq_prod: &[f64],
    x_prod: &[f64],
) -> String {
    let n = freq_react
        .len()
        .min(freq_ts.len())
        .min(freq_prod.len())
        .min(x_react.len())
        .min(x_ts.len())
        .min(x_prod.len());
    let mut out = String::new();
    out.push_str("         CORRELATIONS (VIB, ROT):\n\n");
    out.push_str("              REACTANT               CHANNEL                 PRODUCT\n\n");
    for i in 0..n {
        let tag_react = if x_react[i] <= 0.6 { "rot" } else { "vib" };
        let tag_ts = if x_ts[i] <= 0.6 { "rot" } else { "vib" };
        let tag_prod = if x_prod[i] <= 0.6 { "rot" } else { "vib" };
        out.push_str(&format!(
            "          {:>7.2}   ({:>3})  ---->  {:>7.2}   ({:>3})  ---->  {:>7.2}   ({:>3})\n",
            freq_react[i], tag_react, freq_ts[i], tag_ts, freq_prod[i], tag_prod,
        ));
    }
    out
}
