//! Collisional energy-transfer kernels of the energy-grained master equation.
//!
//! Both kernels return transition probabilities per collision, P(t|j) for a collision that takes a
//! molecule from grain j to grain t, on the grid of one well (grain 0 = well bottom). Both obey
//! detailed balance exactly, P(t|j) f_j = P(j|t) f_t with f_i = rho_i exp(-E_i/kT):
//!
//! - Exponential down: P(t|j) ∝ exp(-(E_j - E_t)/<dE_down>) for t <= j, activating probabilities from
//!   detailed balance, normalization by back substitution from the top grain (Robertson, Comprehensive
//!   Chemical Kinetics 43 (2019), eqs. 4.4, 4.7, 4.11, 4.16). Where the density of states rises steeply at
//!   low energy, all transition probabilities involving a low grain are reduced by a factor attached to
//!   that (lower) grain, which keeps both normalization and detailed balance exact.
//! - Stepladder: steps of size dE_SL between grains i and i+n, P(i+n|i) = A/(1+A),
//!   P(i|i+n) = 1 - P(i+n|i), A = (rho_{i+n}/rho_i) exp(-dE_SL/kT)
//!   (Olzmann, Gebhardt, Scherzer, Int. J. Chem. Kinet. 23, 825 (1991), eqs. 13-18).

/// Transition probabilities per collision on the grid of one well.
#[derive(Debug, Clone)]
pub struct CollisionKernel {
    /// For each source grain j: (target grain t, P(t|j)) for all t != j with P(t|j) > 0.
    pub transitions: Vec<Vec<(usize, f64)>>,
    /// Elastic (no change) probability P(j|j) for each source grain.
    pub elastic: Vec<f64>,
    /// Grains below this index carry a low-energy reduction factor (exponential down only).
    pub low_energy_cut_grain: usize,
    /// Exponent of the low-energy reduction factor that made the normalization possible.
    pub reduction_exponent: f64,
    /// Stepladder step in grains (stepladder only, 0 otherwise).
    pub step_grains: usize,
}

/// Density-of-states gradient above which the low-energy reduction applies:
/// rho(i + n_ref)/rho(i) > threshold with n_ref = floor(1.5 <dE_down>/dE) + 1.
pub const LOW_ENERGY_RHO_GRADIENT_THRESHOLD: f64 = 3.0;
/// Range and step of the reduction exponent m in redfac_i = max(1, (rho(i + n_ref)/rho(i))^m).
pub const REDUCTION_EXPONENT_MIN: f64 = 1.0;
pub const REDUCTION_EXPONENT_MAX: f64 = 3.05;
pub const REDUCTION_EXPONENT_STEP: f64 = 0.1;

/// Exponential-down kernel on a well grid with densities `rho` (states per cm-1, grain i at E = i dE),
/// mean energy transferred in deactivating collisions `mean_down_cm1`, thermal energy `kt_cm1` and
/// collision band half width `band` (grains).
pub fn exponential_down_kernel(
    rho: &[f64],
    grain_width: f64,
    mean_down_cm1: f64,
    kt_cm1: f64,
    band: usize,
) -> Result<CollisionKernel, String> {
    let n = rho.len();
    if n == 0 {
        return Err("Exponential-down kernel: empty energy grid.".into());
    }
    if !(grain_width > 0.0) || !(mean_down_cm1 > 0.0) || !(kt_cm1 > 0.0) {
        return Err(format!(
            "Exponential-down kernel: grain width ({grain_width}), <dE_down> ({mean_down_cm1}) and kT \
             ({kt_cm1}) must be positive."
        ));
    }
    if rho.iter().any(|r| !(*r > 0.0) || !r.is_finite()) {
        return Err("Exponential-down kernel: the density of states must be positive in every grain.".into());
    }

    // Exponential down for deactivating collisions, Robertson (2019) eq. 4.7:
    //   P(t|j) = A_j exp(-(E_j - E_t)/<dE_down>),  E_j >= E_t,
    // on the grain grid one step of dE contributes the factor
    //   beta = exp(-dE/<dE_down>).
    // Activating collisions follow from detailed balance, Robertson eqs. 4.6, 4.10-4.11:
    //   P(t|j) = A_t (rho_t/rho_j) exp(-(E_t - E_j)(1/<dE_down> + 1/kT)),  E_t > E_j,
    // i.e. per grain step the factor beta * gamma with gamma = exp(-dE/kT).
    let beta = (-grain_width / mean_down_cm1).exp();
    let gamma = (-grain_width / kt_cm1).exp();
    let band = band.max(1);

    // Low-energy reduction. For sparse states near the well bottom the activating probabilities of
    // eq. 4.11 can exceed the available normalization and the back substitution of eq. 4.16 returns
    // non-positive coefficients; Robertson (2019, p. 294) traces this to the sparsity of the states,
    // not to numerical error. Below the grain `low_cut`, where the density of states still rises by
    // more than LOW_ENERGY_RHO_GRADIENT_THRESHOLD over n_ref = floor(1.5 <dE_down>/dE) + 1 grains,
    // every transition probability involving a grain i < low_cut is divided by
    //   redfac_i = max(1, (rho_{i+n_ref}/rho_i)^m),
    // attached to the LOWER grain of the pair in both directions, so that detailed balance remains
    // exact. The exponent m is increased from REDUCTION_EXPONENT_MIN in steps of
    // REDUCTION_EXPONENT_STEP until the back substitution succeeds.
    let n_ref = (1.5 * mean_down_cm1 / grain_width) as usize + 1;
    let mut low_cut = 0usize;
    if n > n_ref {
        for i in (0..n - n_ref).rev() {
            if rho[i + n_ref] / rho[i] > LOW_ENERGY_RHO_GRADIENT_THRESHOLD {
                low_cut = i + 1;
                break;
            }
        }
    }

    let mut exponent = REDUCTION_EXPONENT_MIN;
    loop {
        let redfac: Vec<f64> = (0..low_cut)
            .map(|i| (rho[i + n_ref] / rho[i]).powf(exponent).max(1.0))
            .collect();
        match normalize_exponential_down(rho, beta, gamma, band, low_cut, &redfac) {
            Ok((transitions, elastic)) => {
                return Ok(CollisionKernel {
                    transitions,
                    elastic,
                    low_energy_cut_grain: low_cut,
                    reduction_exponent: if low_cut > 0 { exponent } else { 0.0 },
                    step_grains: 0,
                });
            }
            Err(grain) => {
                exponent += REDUCTION_EXPONENT_STEP;
                if exponent > REDUCTION_EXPONENT_MAX {
                    return Err(format!(
                        "Exponential-down normalization failed at grain {grain}: the activating probabilities \
                         exceed 1 even with the low-energy reduction exponent {REDUCTION_EXPONENT_MAX}."
                    ));
                }
            }
        }
    }
}

/// Back substitution of the normalization conditions (Robertson 2019, eq. 4.16) with the low-energy
/// reduction factors `redfac` (grains below `low_cut`). Returns the transitions and elastic
/// probabilities, or the grain at which the activating probabilities alone reach 1.
fn normalize_exponential_down(
    rho: &[f64],
    beta: f64,
    gamma: f64,
    band: usize,
    low_cut: usize,
    redfac: &[f64],
) -> Result<(Vec<Vec<(usize, f64)>>, Vec<f64>), usize> {
    let n = rho.len();
    // Reduction factor of grain i (1 above the low-energy cut).
    let reduction = |i: usize| if i < low_cut { redfac[i] } else { 1.0 };

    // norm_j = 1/A_j of Robertson eq. 4.7. Eq. 4.16 for grain j reads
    //   (1/norm_j) sum_{t<=j} beta^(j-t)/r_t  +  sum_{t>j} (1/norm_t)(rho_t/rho_j)(beta gamma)^(t-j)/r_j = 1,
    // with r the reduction factor of the lower grain of each pair (r = 1 for the elastic term t = j).
    // It is upper triangular: the activating sum needs norm_t of HIGHER grains only, so it is solved
    // from the top grain downward (the top grain has only deactivating collisions, i.e. a reflecting
    // upper boundary; Robertson 2019, p. 278).
    let mut norm = vec![0.0; n];
    let mut transitions: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut elastic = vec![0.0; n];

    for j in (0..n).rev() {
        let lowest_target = j.saturating_sub(band);
        let highest_target = (j + band).min(n - 1);

        // Deactivating part including the elastic term (t = j): sum_{t<=j} beta^(j-t)/r_t.
        let mut down_sum = 1.0;
        for t in lowest_target..j {
            down_sum += beta.powi((j - t) as i32) / reduction(t);
        }

        // Activating probabilities from detailed balance (Robertson eq. 4.11), using the already
        // known normalization of the higher target grain t and the reduction of the lower grain j.
        let mut up = Vec::with_capacity(highest_target.saturating_sub(j));
        let mut up_sum = 0.0;
        for t in (j + 1)..=highest_target {
            let p = (beta * gamma).powi((t - j) as i32) * (rho[t] / rho[j]) / norm[t] / reduction(j);
            up_sum += p;
            up.push((t, p));
        }
        if up_sum >= 1.0 {
            return Err(j);
        }

        // Remaining probability is shared by the deactivating collisions: norm_j = down_sum/(1 - up_sum).
        norm[j] = down_sum / (1.0 - up_sum);

        let mut list = Vec::with_capacity(highest_target - lowest_target);
        for t in lowest_target..j {
            list.push((t, beta.powi((j - t) as i32) / reduction(t) / norm[j]));
        }
        list.extend(up);
        transitions[j] = list;
        elastic[j] = 1.0 / norm[j];
    }

    Ok((transitions, elastic))
}

/// Stepladder kernel with step size `step_cm1` (rounded to a whole number of grains, at least one).
pub fn stepladder_kernel(
    rho: &[f64],
    grain_width: f64,
    step_cm1: f64,
    kt_cm1: f64,
) -> Result<CollisionKernel, String> {
    let n = rho.len();
    if n == 0 {
        return Err("Stepladder kernel: empty energy grid.".into());
    }
    if !(grain_width > 0.0) || !(kt_cm1 > 0.0) {
        return Err("Stepladder kernel: grain width and kT must be positive.".into());
    }
    if rho.iter().any(|r| !(*r > 0.0) || !r.is_finite()) {
        return Err("Stepladder kernel: the density of states must be positive in every grain.".into());
    }

    // The step dE_SL "represents the average amount of energy transferred in down collisions"
    // (Gonzalez-Garcia, Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010), text before eq. 16). The master
    // equation keeps its fine grain dE; a step of dE_SL couples grain i with grain i + n,
    // n = dE_SL/dE (Olzmann, Phys. Chem. Chem. Phys. 4, 3614 (2002), p. 3616: the master equation is
    // split into dE_SL/dE energetically shifted sub-equations). The step is rounded to whole grains.
    let step_grains = (step_cm1 / grain_width).round();
    if !(step_grains >= 1.0) {
        return Err(format!(
            "Stepladder kernel: the step {step_cm1} cm-1 is smaller than one grain ({grain_width} cm-1)."
        ));
    }
    let step = step_grains as usize;
    let step_energy = step as f64 * grain_width;

    let mut transitions: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut out_probability = vec![0.0; n];

    // Olzmann et al. (1991), Int. J. Chem. Kinet. 23, 825:
    //   completeness   P_{i+1,i} + P_{i-1,i} = 1                                   (eq. 13)
    //   detailed bal.  P_{i+1,i}/P_{i,i+1} = (rho_{i+1}/rho_i) exp(-(E_{i+1}-E_i)/kT) (eq. 14)
    // with eq. 13 replaced, "in a very good approximation", by
    //                  P_{i+1,i} + P_{i,i+1} = 1                                   (eq. 15)
    // which gives for a stepladder with step size dE_SL
    //                  P_{i+1,i} = A/(1 + A)                                       (eq. 16)
    //                  P_{i,i+1} = 1 - P_{i+1,i}                                   (eq. 17)
    //                  A = (rho_{i+1}/rho_i) exp(-dE_SL/kT)                        (eq. 18)
    // Index i here counts stepladder steps; on the fine grain grid step i -> i+1 is grain i -> i + n.
    // Detailed balance (eq. 14) holds exactly: P_up/P_down = A. The completeness of eq. 13 is
    // approximate ("with the exception of the first and the last step, all transitions are nearly
    // exactly detailed balanced", p. 831); the remaining probability is elastic.
    for i in 0..n.saturating_sub(step) {
        let a = (rho[i + step] / rho[i]) * (-step_energy / kt_cm1).exp();
        let p_up = a / (1.0 + a); // eq. 16: grain i -> i + n
        let p_down = 1.0 - p_up; // eq. 17: grain i + n -> i
        transitions[i].push((i + step, p_up));
        transitions[i + step].push((i, p_down));
        out_probability[i] += p_up;
        out_probability[i + step] += p_down;
    }

    let elastic = out_probability.iter().map(|p| 1.0 - p).collect();
    Ok(CollisionKernel {
        transitions,
        elastic,
        low_energy_cut_grain: 0,
        reduction_exponent: 0.0,
        step_grains: step,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const KT: f64 = 0.695_034_76 * 500.0;

    fn boltzmann(rho: &[f64], d_e: f64) -> Vec<f64> {
        (0..rho.len()).map(|i| rho[i] * (-(i as f64) * d_e / KT).exp()).collect()
    }

    fn probability(kernel: &CollisionKernel, target: usize, source: usize) -> f64 {
        if target == source {
            return kernel.elastic[source];
        }
        kernel.transitions[source]
            .iter()
            .find(|(t, _)| *t == target)
            .map(|(_, p)| *p)
            .unwrap_or(0.0)
    }

    fn assert_detailed_balance(kernel: &CollisionKernel, f: &[f64]) {
        for j in 0..f.len() {
            for &(t, p) in &kernel.transitions[j] {
                let lhs = p * f[j];
                let rhs = probability(kernel, j, t) * f[t];
                assert!(
                    ((lhs - rhs) / lhs.max(rhs)).abs() < 1e-12,
                    "detailed balance violated for {j}->{t}: {lhs:e} vs {rhs:e}"
                );
            }
        }
    }

    #[test]
    fn exponential_down_is_normalized_and_detailed_balanced_for_sparse_low_energy_states() {
        // Steeply rising density of states at the well bottom, where unscaled back substitution fails.
        let d_e = 20.0;
        let rho: Vec<f64> = (0..60).map(|i| (1.0 + 0.05 * i as f64).powi(10)).collect();
        let kernel = exponential_down_kernel(&rho, d_e, 166.4, KT, 20).unwrap();
        assert!(kernel.low_energy_cut_grain > 0, "the low-energy reduction should be active");
        for j in 0..rho.len() {
            let total = kernel.elastic[j] + kernel.transitions[j].iter().map(|(_, p)| p).sum::<f64>();
            assert!((total - 1.0).abs() < 1e-12, "grain {j}: sum_t P(t|j) = {total}");
            assert!(kernel.elastic[j] > 0.0);
        }
        assert_detailed_balance(&kernel, &boltzmann(&rho, d_e));
    }

    #[test]
    fn exponential_down_reduces_to_unscaled_back_substitution_for_smooth_densities() {
        // Slowly rising density: no reduction; probabilities must equal Robertson eq. 4.16 solved
        // independently here: A_j from back substitution, P(t|j) = A_j b^(j-t) (t <= j) and
        // A_t (rho_t/rho_j) b^(t-j) g^(t-j) (t > j), b = exp(-dE/<dE_down>), g = exp(-dE/kT).
        let d_e = 20.0;
        let alpha = 200.0;
        let band = 30;
        let n = 80;
        // rho(i + n_ref)/rho(i) = exp(0.002 dE n_ref) = exp(0.64) = 1.9 < 3 (n_ref = 16): no reduction.
        let rho: Vec<f64> = (0..n).map(|i| (0.002 * i as f64 * d_e).exp()).collect();
        let kernel = exponential_down_kernel(&rho, d_e, alpha, KT, band).unwrap();
        assert_eq!(kernel.low_energy_cut_grain, 0);

        let b = (-d_e / alpha).exp();
        let g = (-d_e / KT).exp();
        let mut a = vec![0.0; n];
        for j in (0..n).rev() {
            let down: f64 = (j.saturating_sub(band)..=j).map(|t| b.powi((j - t) as i32)).sum();
            let up: f64 = (j + 1..(j + band + 1).min(n))
                .map(|t| a[t] * rho[t] / rho[j] * (b * g).powi((t - j) as i32))
                .sum();
            a[j] = (1.0 - up) / down;
        }
        for j in 0..n {
            for t in j.saturating_sub(band)..(j + band + 1).min(n) {
                let expected = if t <= j {
                    a[j] * b.powi((j - t) as i32)
                } else {
                    a[t] * rho[t] / rho[j] * (b * g).powi((t - j) as i32)
                };
                let got = probability(&kernel, t, j);
                assert!(((got - expected) / expected).abs() < 1e-12, "P({t}|{j}) = {got}, expected {expected}");
            }
        }
    }

    #[test]
    fn stepladder_follows_olzmann_eqs_16_to_18() {
        let d_e = 10.0;
        let step = 190.0; // 19 grains
        let rho: Vec<f64> = (0..200).map(|i| (1.0 + 0.02 * i as f64).powi(12)).collect();
        let kernel = stepladder_kernel(&rho, d_e, step, KT).unwrap();
        assert_eq!(kernel.step_grains, 19);
        let n = 19;
        for i in 0..rho.len() - n {
            let a = rho[i + n] / rho[i] * (-step / KT).exp();
            let up = probability(&kernel, i + n, i);
            let down = probability(&kernel, i, i + n);
            assert!((up - a / (1.0 + a)).abs() < 1e-14, "P({}|{i}) = {up}", i + n);
            assert!((down - (1.0 - up)).abs() < 1e-14, "P({i}|{}) = {down}", i + n);
        }
        // Only steps of exactly n grains occur.
        for j in 0..rho.len() {
            for &(t, _) in &kernel.transitions[j] {
                assert_eq!(t.abs_diff(j), n);
            }
        }
        assert_detailed_balance(&kernel, &boltzmann(&rho, d_e));
    }

    #[test]
    fn stepladder_step_must_cover_at_least_one_grain() {
        let rho = vec![1.0; 10];
        assert!(stepladder_kernel(&rho, 20.0, 5.0, KT).is_err());
    }
}
