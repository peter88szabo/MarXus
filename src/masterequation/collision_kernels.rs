//! Collisional energy-transfer kernels of the energy-grained master equation.
//!
//! Both kernels return transition probabilities per collision, P(t|j) for a collision that takes a
//! molecule from grain j to grain t, on the grid of one well (grain 0 = well bottom). Both obey
//! detailed balance exactly, P(t|j) f_j = P(j|t) f_t with f_i = rho_i exp(-E_i/kT):
//!
//! - Exponential down: P(t|j) ∝ exp(-(E_j - E_t)/<dE_down>) for t <= j, activating probabilities from
//!   detailed balance, normalization by back substitution from the top grain (Robertson, Comprehensive
//!   Chemical Kinetics 43 (2019), eqs. 4.4, 4.7, 4.11, 4.16), from the top grain down. Where the back
//!   substitution fails (sparse states at low energy, Robertson 2019, p. 294), the failing grain and all
//!   grains below it form a reservoir: one thermalized state in the master equation, as the reservoir
//!   state of MESMER (manual, Sec. 14.2.1). The kernel keeps the normalization of every grain above, and
//!   their transitions into reservoir grains; it gives no transitions out of reservoir grains (activation
//!   from the reservoir follows by detailed balance in the operator). Their number is `reservoir_grains`.
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
    /// Stepladder step in grains (stepladder only, 0 otherwise).
    pub step_grains: usize,
    /// Exponential down: the grains 0 .. reservoir_grains form the reservoir, because eq. 4.16 cannot be
    /// satisfied at the highest of them; they have no transitions here (and P(j|j) = 0). 0 when eq. 4.16
    /// holds at every grain.
    pub reservoir_grains: usize,
}

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

    // For sparse states near the well bottom the activating probabilities of eq. 4.11 can exceed the
    // available normalization, and the back substitution of eq. 4.16 has no positive solution; Robertson
    // (2019, p. 294) traces this to the sparsity of the states, not to numerical error. The back substitution
    // runs from the top; the grains above the first failure keep their normalization (it needs the higher
    // grains only), and the failing grain and all grains below form the reservoir (MESMER manual,
    // Sec. 14.2.1: "a collection of grains that are represented with one grain because we assume that these
    // grains are always thermalized").
    let (transitions, elastic, failed) = normalize_exponential_down(rho, beta, gamma, band);
    let reservoir_grains = failed.map_or(0, |j| j + 1);
    if reservoir_grains >= n {
        return Err(format!(
            "Exponential-down normalization (Robertson 2019, eq. 4.16) fails at the top grain {}: the activating \
             probabilities exceed 1 at every grain of the well.",
            n - 1
        ));
    }
    Ok(CollisionKernel { transitions, elastic, step_grains: 0, reservoir_grains })
}

/// Back substitution of the normalization conditions (Robertson 2019, eq. 4.16) from the top grain down.
/// Returns the transitions and elastic probabilities of every grain above the first grain at which the
/// activating probabilities alone reach 1, and that grain (None if eq. 4.16 holds everywhere); the grains
/// from it down have no transitions.
fn normalize_exponential_down(
    rho: &[f64],
    beta: f64,
    gamma: f64,
    band: usize,
) -> (Vec<Vec<(usize, f64)>>, Vec<f64>, Option<usize>) {
    let n = rho.len();

    // norm_j = 1/A_j of Robertson eq. 4.7. Eq. 4.16 for grain j reads
    //   (1/norm_j) sum_{t<=j} beta^(j-t)  +  sum_{t>j} (1/norm_t)(rho_t/rho_j)(beta gamma)^(t-j) = 1.
    // It is upper triangular: the activating sum needs norm_t of HIGHER grains only, so it is solved
    // from the top grain downward (the top grain has only deactivating collisions, i.e. a reflecting
    // upper boundary; Robertson 2019, p. 278).
    let mut norm = vec![0.0; n];
    let mut transitions: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut elastic = vec![0.0; n];

    for j in (0..n).rev() {
        let lowest_target = j.saturating_sub(band);
        let highest_target = (j + band).min(n - 1);

        // Deactivating part including the elastic term (t = j): sum_{t<=j} beta^(j-t).
        let mut down_sum = 1.0;
        for t in lowest_target..j {
            down_sum += beta.powi((j - t) as i32);
        }

        // Activating probabilities from detailed balance (Robertson eq. 4.11), using the already
        // known normalization of the higher target grain t.
        let mut up = Vec::with_capacity(highest_target.saturating_sub(j));
        let mut up_sum = 0.0;
        for t in (j + 1)..=highest_target {
            let p = (beta * gamma).powi((t - j) as i32) * (rho[t] / rho[j]) / norm[t];
            up_sum += p;
            up.push((t, p));
        }
        if up_sum >= 1.0 {
            return (transitions, elastic, Some(j));
        }

        // Remaining probability is shared by the deactivating collisions: norm_j = down_sum/(1 - up_sum).
        norm[j] = down_sum / (1.0 - up_sum);

        let mut list = Vec::with_capacity(highest_target - lowest_target);
        for t in lowest_target..j {
            list.push((t, beta.powi((j - t) as i32) / norm[j]));
        }
        list.extend(up);
        transitions[j] = list;
        elastic[j] = 1.0 / norm[j];
    }

    (transitions, elastic, None)
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
        step_grains: step,
        reservoir_grains: 0,
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
    fn exponential_down_lumps_the_grains_from_the_first_failure_down_into_a_reservoir() {
        // Steeply rising density of states at the well bottom, where the back substitution of eq. 4.16 fails
        // (at grain 18). As the reservoir state of MESMER (manual, Sec. 14.2.1): the normalization from the
        // top is kept for every grain above the failure, and the failing grain and all grains below it form
        // the reservoir, which the operator treats as one thermalized state. The kernel gives no transitions
        // out of reservoir grains; the transitions of the grains above, including those into reservoir
        // grains, are those of eq. 4.16 on the complete grid.
        let (d_e, alpha, band) = (20.0, 166.4, 20);
        let rho: Vec<f64> = (0..60).map(|i| (1.0 + 0.05 * i as f64).powi(10)).collect();
        let n = rho.len();
        let kernel = exponential_down_kernel(&rho, d_e, alpha, KT, band).unwrap();

        // Independent: eq. 4.16 from the top grain down to the first failure.
        let b = (-d_e / alpha).exp();
        let g = (-d_e / KT).exp();
        let mut a = vec![0.0; n];
        let mut failed = None;
        for j in (0..n).rev() {
            let down: f64 = (j.saturating_sub(band)..=j).map(|t| b.powi((j - t) as i32)).sum();
            let up: f64 = (j + 1..(j + band + 1).min(n)).map(|t| a[t] * rho[t] / rho[j] * (b * g).powi((t - j) as i32)).sum();
            if up >= 1.0 {
                failed = Some(j);
                break;
            }
            a[j] = (1.0 - up) / down;
        }
        let first = failed.expect("eq. 4.16 fails for this density") + 1;
        assert_eq!(first, 19);
        assert_eq!(kernel.reservoir_grains, first);

        for j in 0..first {
            assert!(kernel.transitions[j].is_empty() && kernel.elastic[j] == 0.0, "grain {j} is in the reservoir");
        }
        for j in first..n {
            let total = kernel.elastic[j] + kernel.transitions[j].iter().map(|(_, p)| p).sum::<f64>();
            assert!((total - 1.0).abs() < 1e-12, "grain {j}: sum_t P(t|j) = {total}");
            for t in j.saturating_sub(band)..(j + band + 1).min(n) {
                let expected = if t <= j { a[j] * b.powi((j - t) as i32) } else { a[t] * rho[t] / rho[j] * (b * g).powi((t - j) as i32) };
                let got = probability(&kernel, t, j);
                assert!(((got - expected) / expected).abs() < 1e-12, "P({t}|{j}) = {got}, expected {expected}");
            }
        }
        // Detailed balance between the grains above the reservoir.
        let f = boltzmann(&rho, d_e);
        for j in first..n {
            for &(t, p) in kernel.transitions[j].iter().filter(|&&(t, _)| t >= first) {
                let (lhs, rhs) = (p * f[j], probability(&kernel, j, t) * f[t]);
                assert!(((lhs - rhs) / lhs.max(rhs)).abs() < 1e-12, "detailed balance {j}->{t}");
            }
        }
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
        assert_eq!(kernel.reservoir_grains, 0);

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

    // Density of states rising steeply from the well bottom (rho(i + 9)/rho(i) > 3 for the first grains),
    // for which the back substitution of eq. 4.16 nevertheless succeeds at every grain.
    fn steep_but_normalizable_density() -> Vec<f64> {
        (0..120).map(|i| (1.0 + i as f64 * 38.0 / 2000.0).powi(10)).collect()
    }

    #[test]
    fn exponential_down_uses_plain_back_substitution_wherever_it_succeeds() {
        // No low-energy reduction where eq. 4.16 can be solved, however steep the density of states:
        // the probabilities equal eq. 4.16 solved independently here.
        let (d_e, alpha, band) = (38.0, 202.6, 80);
        let rho = steep_but_normalizable_density();
        let n = rho.len();
        let kernel = exponential_down_kernel(&rho, d_e, alpha, KT, band).unwrap();
        assert_eq!(kernel.reservoir_grains, 0);

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
    fn exponential_down_kernel_is_continuous_in_the_mean_energy_transfer() {
        // <dE_down> = 202.6666 and 202.6668 cm-1 on 38 cm-1 grains: 1.5 <dE_down>/dE crosses 8 in between
        // (the integer window of the low-energy reduction would change from 8 to 9 grains there).
        // The smooth dependence changes P(t|j) by about |t - j| dE d<dE_down>/<dE_down>^2 <= 1.2e-5
        // (|t - j| <= 80); no probability may change by more than 1e-4.
        let (d_e, band) = (38.0, 80);
        let rho = steep_but_normalizable_density();
        let below = exponential_down_kernel(&rho, d_e, 202.6666, KT, band).unwrap();
        let above = exponential_down_kernel(&rho, d_e, 202.6668, KT, band).unwrap();
        for j in 0..rho.len() {
            for t in j.saturating_sub(band)..(j + band + 1).min(rho.len()) {
                let (p, q) = (probability(&below, t, j), probability(&above, t, j));
                assert!(((p - q) / p).abs() < 1e-4, "P({t}|{j}) jumps from {p:e} to {q:e}");
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
