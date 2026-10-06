//! Parallel evaluation of a run over its conditions (T, p) with rayon.
//!
//! Every condition is an independent master equation: its own operator, factorization and solution. The
//! conditions are therefore computed concurrently on a local thread pool (`ConditionPool`), and the results
//! are returned in the order of the grid (temperatures outer, pressures inner), so that the output does not
//! depend on the number of threads.
//!
//! Thread count: the `threads` argument if given, otherwise the environment variable RAYON_NUM_THREADS,
//! otherwise the number of logical cores (rayon's default).

use rayon::prelude::*;

/// A thread pool for the conditions of a run.
pub struct ConditionPool {
    pool: rayon::ThreadPool,
}

impl ConditionPool {
    /// Pool with `threads` workers (None: RAYON_NUM_THREADS, otherwise the number of logical cores).
    pub fn new(threads: Option<usize>) -> Result<Self, String> {
        if threads == Some(0) {
            return Err("ConditionPool: the number of threads must be at least 1.".into());
        }
        // num_threads(0) is rayon's default: RAYON_NUM_THREADS, otherwise the number of logical cores.
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.unwrap_or(0))
            .build()
            .map_err(|e| format!("ConditionPool: {e}"))?;
        Ok(Self { pool })
    }

    /// Number of worker threads.
    pub fn threads(&self) -> usize {
        self.pool.current_num_threads()
    }

    /// `f(T, p)` for every condition of the grid, evaluated in parallel; results in grid order
    /// (temperatures outer, pressures inner).
    pub fn map_conditions<R: Send>(
        &self,
        temperatures: &[f64],
        pressures: &[f64],
        f: impl Fn(f64, f64) -> R + Sync + Send,
    ) -> Vec<R> {
        let grid: Vec<(f64, f64)> = temperatures
            .iter()
            .flat_map(|&t| pressures.iter().map(move |&p| (t, p)))
            .collect();
        // An indexed parallel iterator collects in the order of the grid, whatever the order of completion.
        self.pool
            .install(|| grid.par_iter().map(|&(t, p)| f(t, p)).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_pool_has_the_requested_number_of_threads() {
        assert_eq!(ConditionPool::new(Some(2)).unwrap().threads(), 2);
        assert_eq!(ConditionPool::new(Some(1)).unwrap().threads(), 1);
        assert!(ConditionPool::new(Some(0)).is_err());
    }

    #[test]
    fn results_come_in_grid_order_whatever_the_order_of_completion() {
        let pool = ConditionPool::new(Some(3)).unwrap();
        // Later conditions finish first (shorter work), so completion order differs from grid order.
        let out = pool.map_conditions(&[100.0, 200.0, 300.0], &[1.0, 10.0], |t, p| {
            std::thread::sleep(std::time::Duration::from_millis((1000.0 / t) as u64));
            (t, p)
        });
        assert_eq!(
            out,
            vec![
                (100.0, 1.0),
                (100.0, 10.0),
                (200.0, 1.0),
                (200.0, 10.0),
                (300.0, 1.0),
                (300.0, 10.0)
            ]
        );
    }

    #[test]
    fn parallel_master_equations_equal_the_sequential_ones() {
        use crate::masterequation::chemical_activation_driver::{
            run_chemical_activation, ChemicalActivationRun, SourceSpecification,
        };
        use crate::masterequation::chemical_activation_network::{
            AbsorbingBarrier, ChemicalActivationOptions, CollisionModel, SteadyState,
        };
        use crate::masterequation::chemical_activation_operator::tests::two_well_network;
        use crate::masterequation::chemical_activation_steady_state::LinearSolver;
        let network = two_well_network();
        let run = |t: f64, p: f64| ChemicalActivationRun {
            temperatures_kelvin: vec![t],
            pressures_torr: vec![p],
            options: ChemicalActivationOptions {
                collision_model: CollisionModel::ExponentialDown {
                    cutoff_in_mean_down: 10.0,
                },
                steady_state: SteadyState::Intermediate {
                    barrier: AbsorbingBarrier::default(),
                },
            },
            solver: LinearSolver::BandedCholesky,
            source: SourceSpecification::Fixed(
                network
                    .wells
                    .iter()
                    .map(|w| {
                        let mut f = vec![0.0; w.grain_count()];
                        f[w.grain_count() - 1] = 1.0;
                        f
                    })
                    .collect(),
            ),
            tolerance: 1e-8,
        };
        let temperatures = [250.0, 300.0];
        let pressures = [10.0, 100.0, 760.0];
        let pool = ConditionPool::new(Some(3)).unwrap();
        let parallel = pool.map_conditions(&temperatures, &pressures, |t, p| {
            run_chemical_activation(&network, &run(t, p))
                .unwrap()
                .remove(0)
        });
        let mut k = 0;
        for &t in &temperatures {
            for &p in &pressures {
                let sequential = run_chemical_activation(&network, &run(t, p))
                    .unwrap()
                    .remove(0);
                assert_eq!(parallel[k].conditions.temperature_kelvin, t);
                assert_eq!(parallel[k].conditions.pressure_torr, p);
                for (a, b) in parallel[k]
                    .result
                    .channels
                    .iter()
                    .zip(&sequential.result.channels)
                {
                    assert_eq!(
                        a.flux.to_bits(),
                        b.flux.to_bits(),
                        "bitwise identical at {t} K, {p} Torr"
                    );
                }
                k += 1;
            }
        }
    }
}
