#[derive(Debug, Clone, Copy)]
pub struct SacmCentrifugal {
    pub eq_rot_const: f64,
    pub diss_energy: f64,
    pub corr_a1: f64,
    pub corr_a2: f64,
}

impl SacmCentrifugal {
    /// J-dependent centrifugal rotational constant for the lowest channel.
    pub fn bcentrifugal(self, j: usize) -> f64 {
        if j == 0 {
            return self.eq_rot_const;
        }

        let jval = j as f64;
        let ecent = self.eq_rot_const * jval * (jval + 1.0) / self.diss_energy;
        let mut z = 1.0e-10;
        let mut last = z;
        let mut iter = 0;

        loop {
            iter += 1;
            last = z;
            let f = lowest_channel_derivative(z, ecent, self.corr_a1, self.corr_a2);
            let df = lowest_channel_second_derivative(z, ecent, self.corr_a1, self.corr_a2);
            z = last - f / df;

            if (z - last).abs() <= 1.0e-4 {
                break;
            }

            if iter > 1000 {
                break;
            }
        }

        let bracket = 1.0 + self.corr_a1 * z + self.corr_a2 * z * z;
        self.eq_rot_const / bracket
    }
}

fn lowest_channel_derivative(z: f64, ecent: f64, a1: f64, a2: f64) -> f64 {
    let bracket = 1.0 + a1 * z + a2 * z * z;
    let exp_neg = (-z).exp();
    2.0 * (1.0 - exp_neg) * exp_neg - ecent * (a1 + 2.0 * a2 * z) / (bracket * bracket)
}

fn lowest_channel_second_derivative(z: f64, ecent: f64, a1: f64, a2: f64) -> f64 {
    let bracket = 1.0 + a1 * z + a2 * z * z;
    let exp_neg = (-z).exp();
    let numerator = (a1 + 2.0 * a2 * z).powi(2) - a2 * bracket;
    let factor = numerator / (bracket * bracket * bracket);
    4.0 * (-2.0 * z).exp() - 2.0 * exp_neg + 2.0 * ecent * factor
}
