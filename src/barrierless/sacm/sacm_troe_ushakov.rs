use super::pst_channels::{PstChannels, KB_CM};
use super::types::SacmEnergyGrid;

#[derive(Debug, Clone, Copy)]
pub struct RigidityTerm {
    pub weight: f64,
    pub temperature: f64,
}

#[derive(Debug, Clone)]
pub enum RigidityModel {
    Exponential { temperature: f64 },
    SumExponential { base: f64, terms: Vec<RigidityTerm> },
}

/// Apply Troe-Ushakov rigidity factors f_rigid(z) to PST W(E,J).
pub fn apply_troe_ushakov_rigidity(
    pst: &PstChannels,
    grid: SacmEnergyGrid,
    model: RigidityModel,
) -> Vec<f64> {
    let mut out = vec![0.0; pst.w_e.len()];
    for (i, &we) in pst.w_e.iter().enumerate() {
        let z = i as f64 * grid.dE;
        let factor = rigidity_factor(z, &model);
        out[i] = we * factor;
    }
    out
}

fn rigidity_factor(z: f64, model: &RigidityModel) -> f64 {
    if z <= 0.0 {
        return 0.0;
    }
    match *model {
        RigidityModel::Exponential { temperature } => {
            if temperature <= 0.0 {
                0.0
            } else {
                (-(z / (KB_CM * temperature))).exp()
            }
        }
        RigidityModel::SumExponential { base, ref terms } => {
            let mut out = base;
            for term in terms {
                if term.temperature > 0.0 {
                    out += term.weight * (-(z / (KB_CM * term.temperature))).exp();
                }
            }
            out
        }
    }
}
