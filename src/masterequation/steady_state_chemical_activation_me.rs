use super::matrix_physics_assembly::{absolute_energy_cm1, assemble_raw_operator};
use super::reaction_network::MultiwellLinearSolver;
use super::reaction_network::{
    ChemicalActivationDefinition, MasterEquationSettings, MicrocanonicalProvider, SolutionResults,
    WellDefinition,
};
use super::state_index::GlobalLayout;
use crate::numeric::iterative_solvers::{solve_bicgstab_left_jacobi_dense, BiCgStabDiagnostics};
use crate::numeric::krylov::{
    solve_bicgstab_left_preconditioned, solve_gmres_left_preconditioned_restarted,
    JacobiPreconditioner, KrylovDiagnostics, LinearOperator,
};
use crate::numeric::ldlt_solvers::{
    solve_symmetric_indefinite_ldlt_bunch_kaufman, LdltDiagnostics,
};
use crate::numeric::linear_algebra::{
    cholesky_solve_spd_with_diagnostics, CholeskyDiagnostics, DenseMatrix, DiagonalScale,
};

pub enum LinearSolveMethod {
    CholeskySpd,
    LdltSymmetricIndefinite,
    BiCgStab,
    Gmres,
}

pub struct MasterEquationSolveDiagnostics {
    pub collision_conservation_max_abs: f64,
    pub transformed_symmetry_relative_frobenius: f64,
    pub solve_method: LinearSolveMethod,
    pub cholesky: Option<CholeskyDiagnostics>,
    pub ldlt: Option<LdltDiagnostics>,
    pub bicgstab: Option<BiCgStabDiagnostics>,
    pub krylov: Option<KrylovDiagnostics>,
}

/// Main engine: assembles and solves the steady-state chemical activation master equation.
pub struct MasterEquationEngine {
    wells: Vec<WellDefinition>,
    settings: MasterEquationSettings,
}

impl MasterEquationEngine {
    pub fn new(wells: Vec<WellDefinition>, settings: MasterEquationSettings) -> Self {
        Self { wells, settings }
    }

    pub fn solve_steady_state_chemical_activation_with_diagnostics(
        &self,
        micro: &dyn MicrocanonicalProvider,
        activation: ChemicalActivationDefinition,
    ) -> Result<(SolutionResults, MasterEquationSolveDiagnostics), String> {
        let global_layout = GlobalLayout::from_wells(&self.wells)?;

        // 1) Assemble the raw operator L (collisions + sinks + inter-well couplings)
        let (raw_operator, assembly_diag) =
            assemble_raw_operator(&self.wells, &self.settings, micro, &global_layout)?;

        // 2) Build chemical activation source vector s (raw, then normalized)
        let source_vector = self.build_activation_source(micro, &global_layout, &activation)?;

        // 3) Similarity transform weights W (diagonal)
        let similarity_scale = self.build_similarity_scale(micro, &global_layout)?;

        // L[target, source] obeys detailed balance L_ij f_j = L_ji f_i with W = sqrt(f), so the
        // symmetrized matrix is A = W^{-1} L W. With rhs b = W^{-1} s and p = W q,
        // (-A) q = b is equivalent to -L p = s.
        let transformed_rhs = similarity_scale.inverse_apply_to_vector(&source_vector);

        // We solve (-A) q = b, where A = W^{-1} L W.
        // For iterative modes we apply -A via matvec without explicitly forming A.
        let symmetry_rel = {
            let n = raw_operator.size();
            let mut norm2 = 0.0;
            let mut diff2 = 0.0;

            for i in 0..n {
                let aii = -raw_operator.get(i, i);
                norm2 += aii * aii;
            }

            for i in 0..n {
                for j in (i + 1)..n {
                    let aij = -raw_operator.get(i, j)
                        * (similarity_scale.diagonal[j] / similarity_scale.diagonal[i]);
                    let aji = -raw_operator.get(j, i)
                        * (similarity_scale.diagonal[i] / similarity_scale.diagonal[j]);
                    norm2 += 2.0 * aij * aij;
                    let d = aij - aji;
                    diff2 += 2.0 * d * d;
                }
            }

            (diff2.sqrt()) / (norm2.sqrt().max(1e-300))
        };

        // 4) Solve using SPD Cholesky if possible, else fall back.
        let mut diag = MasterEquationSolveDiagnostics {
            collision_conservation_max_abs: assembly_diag.collision_conservation_max_abs,
            transformed_symmetry_relative_frobenius: symmetry_rel,
            solve_method: LinearSolveMethod::CholeskySpd,
            cholesky: None,
            ldlt: None,
            bicgstab: None,
            krylov: None,
        };

        let transformed_solution = match self.settings.linear_solver {
            MultiwellLinearSolver::Direct => {
                let negative_transformed_dense = raw_operator
                    .similarity_transform(&similarity_scale)
                    .scaled(-1.0);
                // Cholesky and LDLT read only one triangle, so they are valid only when the
                // transformed matrix is symmetric (detailed balance). Otherwise go straight to
                // the non-symmetric BiCGSTAB solver.
                let sym_tol = 1e-10;
                if symmetry_rel > sym_tol {
                    let (x, bicg) = solve_bicgstab_left_jacobi_dense(
                        &negative_transformed_dense,
                        &transformed_rhs,
                        self.settings.krylov_tolerance,
                        self.settings.krylov_max_iter,
                    )
                    .map_err(|e| {
                        format!(
                            "Transformed matrix is not symmetric (relative asymmetry {symmetry_rel:e}); BiCGSTAB failed: {e}"
                        )
                    })?;
                    diag.solve_method = LinearSolveMethod::BiCgStab;
                    diag.bicgstab = Some(bicg);
                    x
                } else {
                    match cholesky_solve_spd_with_diagnostics(
                        &negative_transformed_dense,
                        &transformed_rhs,
                    ) {
                        Ok((x, chol)) => {
                            diag.solve_method = LinearSolveMethod::CholeskySpd;
                            diag.cholesky = Some(chol);
                            x
                        }
                        Err(chol_err) => match solve_symmetric_indefinite_ldlt_bunch_kaufman(
                            &negative_transformed_dense,
                            &transformed_rhs,
                        ) {
                            Ok((x, ldlt)) => {
                                diag.solve_method = LinearSolveMethod::LdltSymmetricIndefinite;
                                diag.ldlt = Some(ldlt);
                                x
                            }
                            Err(ldlt_err) => {
                                let (x, bicg) = solve_bicgstab_left_jacobi_dense(
                                    &negative_transformed_dense,
                                    &transformed_rhs,
                                    self.settings.krylov_tolerance,
                                    self.settings.krylov_max_iter,
                                )
                                .map_err(|e| {
                                    format!(
                                        "All solvers failed.\nCholesky: {chol_err}\nLDLT: {ldlt_err}\nBiCGSTAB: {e}"
                                    )
                                })?;
                                diag.solve_method = LinearSolveMethod::BiCgStab;
                                diag.bicgstab = Some(bicg);
                                x
                            }
                        },
                    }
                }
            }
            MultiwellLinearSolver::Gmres | MultiwellLinearSolver::BiCgStab => {
                struct NegativeSimilarityOp<'a> {
                    l: &'a DenseMatrix,
                    w: &'a DiagonalScale,
                    tmp: Vec<f64>,
                    tmp2: Vec<f64>,
                }

                impl<'a> NegativeSimilarityOp<'a> {
                    fn new(l: &'a DenseMatrix, w: &'a DiagonalScale) -> Self {
                        let n = l.size();
                        Self {
                            l,
                            w,
                            tmp: vec![0.0; n],
                            tmp2: vec![0.0; n],
                        }
                    }
                }

                impl LinearOperator for NegativeSimilarityOp<'_> {
                    fn dim(&self) -> usize {
                        self.l.size()
                    }

                    fn matvec(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), String> {
                        let n = self.l.size();
                        if x.len() != n || y.len() != n {
                            return Err("Dimension mismatch in NegativeSimilarityOp::matvec".into());
                        }
                        // y = -A x with A = W^{-1} L W
                        // tmp = W x
                        for i in 0..n {
                            self.tmp[i] = x[i] * self.w.diagonal[i];
                        }
                        // tmp2 = L tmp
                        self.l.matvec_into(&self.tmp, &mut self.tmp2)?;
                        // y = - W^{-1} tmp2
                        for i in 0..n {
                            y[i] = -self.tmp2[i] / self.w.diagonal[i];
                        }
                        Ok(())
                    }
                }

                // Jacobi preconditioner on diag(-A) = -diag(L).
                let mut diag_a = vec![0.0; raw_operator.size()];
                for i in 0..raw_operator.size() {
                    diag_a[i] = -raw_operator.get(i, i);
                }
                let jacobi = JacobiPreconditioner::from_diagonal(&diag_a)?;

                let mut op = NegativeSimilarityOp::new(&raw_operator, &similarity_scale);
                let (x, kdiag) = match self.settings.linear_solver {
                    MultiwellLinearSolver::Gmres => {
                        let (x, kd) = solve_gmres_left_preconditioned_restarted(
                            &mut op,
                            &jacobi,
                            &transformed_rhs,
                            self.settings.krylov_tolerance,
                            self.settings.krylov_max_iter,
                            self.settings.gmres_restart,
                        )?;
                        diag.solve_method = LinearSolveMethod::Gmres;
                        (x, kd)
                    }
                    MultiwellLinearSolver::BiCgStab => {
                        let (x, kd) = solve_bicgstab_left_preconditioned(
                            &mut op,
                            &jacobi,
                            &transformed_rhs,
                            self.settings.krylov_tolerance,
                            self.settings.krylov_max_iter,
                        )?;
                        diag.solve_method = LinearSolveMethod::BiCgStab;
                        (x, kd)
                    }
                    MultiwellLinearSolver::Direct => unreachable!(),
                };
                diag.krylov = Some(kdiag);
                x
            }
        };

        // 5) Recover p = W q
        let steady_state_population = similarity_scale.apply_to_vector(&transformed_solution);

        // 6) Compute macroscopic (averaged) rate constants
        let macro_rates =
            self.compute_macroscopic_rates(micro, &global_layout, &steady_state_population)?;

        Ok((macro_rates, diag))
    }

    /// Solve the steady-state chemical activation problem.
    ///
    /// Returns:
    /// - steady_state_population: p (stacked over all wells/grains)
    /// - macroscopic_rate_constants: per (well, channel) averaged rates k(T,p)
    /// - total_outgoing_rate_constant: sum of outgoing channels (+ stabilization leakage if truncation is used)
    pub fn solve_steady_state_chemical_activation(
        &self,
        micro: &dyn MicrocanonicalProvider,
        activation: ChemicalActivationDefinition,
    ) -> Result<SolutionResults, String> {
        let (res, _diag) =
            self.solve_steady_state_chemical_activation_with_diagnostics(micro, activation)?;
        Ok(res)
    }

    /// Build the chemical activation source vector (normalized).
    ///
    /// Implemented injection distribution (per your doc’s CA concept):
    ///
    ///   s_raw(w*, i) = ρ(i) * exp( - i ΔE / (k_B T) ) * k_recomb(i)
    ///   s(w*, i)     = s_raw / Σ_i s_raw
    ///
    /// All other wells: s = 0.
    fn build_activation_source(
        &self,
        micro: &dyn MicrocanonicalProvider,
        layout: &GlobalLayout,
        activation: &ChemicalActivationDefinition,
    ) -> Result<Vec<f64>, String> {
        let mut source = vec![0.0; layout.total_state_count];

        let w_star = activation.activated_well_index;
        if w_star >= self.wells.len() {
            return Err("Activated well index out of range".into());
        }
        let well = &self.wells[w_star];
        let c_star = activation.recombination_channel_index;
        if c_star >= well.channels.len() {
            return Err("Recombination channel index out of range".into());
        }

        let temperature = self.settings.temperature_kelvin;
        let boltzmann = self.settings.boltzmann_constant_wavenumber_per_kelvin;

        let mut normalization_sum = 0.0;

        for local_grain in
            well.lowest_included_grain_index..well.one_past_highest_included_grain_index
        {
            let rho = micro.density_of_states(w_star, local_grain);
            if rho <= 0.0 {
                continue;
            }

            let k_recomb = if local_grain >= well.nonreactive_grain_count {
                micro
                    .microcanonical_rate(w_star, c_star, local_grain)
                    .max(0.0)
            } else {
                0.0
            };

            let energy_cm1 = absolute_energy_cm1(well, local_grain);
            let boltzmann_weight = (-energy_cm1 / (boltzmann * temperature)).exp();

            let raw = rho * boltzmann_weight * k_recomb;

            let global_idx = layout.global_index_of(w_star, local_grain)?;
            source[global_idx] = raw;
            normalization_sum += raw;
        }

        if normalization_sum <= 0.0 {
            return Err(
                "Chemical activation source normalization is zero; check k(E) and rho(E) inputs."
                    .into(),
            );
        }

        for x in &mut source {
            *x /= normalization_sum;
        }

        Ok(source)
    }

    /// Build similarity transform scale W (diagonal).
    ///
    /// We use an equilibrium-like weight:
    ///
    ///   W(state) = sqrt( ρ(i) * exp( -E_abs / (k_B T) ) )
    ///
    /// where E_abs uses the alignment offset:
    ///
    ///   E_abs = (i + offset_w) * ΔE
    fn build_similarity_scale(
        &self,
        micro: &dyn MicrocanonicalProvider,
        layout: &GlobalLayout,
    ) -> Result<DiagonalScale, String> {
        let mut diag = vec![0.0; layout.total_state_count];
        let temperature = self.settings.temperature_kelvin;
        let boltzmann = self.settings.boltzmann_constant_wavenumber_per_kelvin;

        for (well_index, well) in self.wells.iter().enumerate() {
            let grain_width = well.energy_grain_width_cm1;
            for local_grain in
                well.lowest_included_grain_index..well.one_past_highest_included_grain_index
            {
                let rho = micro.density_of_states(well_index, local_grain);
                if rho <= 0.0 {
                    return Err(format!(
                        "Non-positive density of states in similarity scale at well={}, grain={}",
                        well.well_name, local_grain
                    ));
                }

                let abs_grain = (local_grain as isize) + well.alignment_offset_in_grains;
                let abs_energy_cm1 = (abs_grain as f64) * grain_width;
                let abs_boltzmann = (-abs_energy_cm1 / (boltzmann * temperature)).exp();

                let global_idx = layout.global_index_of(well_index, local_grain)?;
                diag[global_idx] = (rho * abs_boltzmann).sqrt();
            }
        }

        Ok(DiagonalScale { diagonal: diag })
    }

    /// Compute macroscopic rates by population-weighted microcanonical averaging:
    ///
    /// Let population vector be p(state). Normalize:
    ///   P_total = Σ_state p(state)
    ///
    /// For each well w and channel c:
    ///   k_w,c(T,p) = (1/P_total) Σ_{i in included grains} p(w,i) * k_w,c(i)
    ///
    /// Also compute total outgoing:
    ///   k_out_total = Σ_{outgoing channels} k_w,c + (optional truncation leakage, not included in this simplified kernel)
    fn compute_macroscopic_rates(
        &self,
        micro: &dyn MicrocanonicalProvider,
        layout: &GlobalLayout,
        population: &[f64],
    ) -> Result<SolutionResults, String> {
        if population.len() != layout.total_state_count {
            return Err("Population vector length mismatch".into());
        }

        let total_population: f64 = population.iter().sum();
        if total_population <= 0.0 {
            return Err("Total population is non-positive after solve".into());
        }

        let mut per_well_per_channel_rates: Vec<Vec<f64>> = Vec::new();
        let mut total_outgoing = 0.0;

        for (well_index, well) in self.wells.iter().enumerate() {
            let mut channel_rates = vec![0.0; well.channels.len()];

            for local_grain in
                well.lowest_included_grain_index..well.one_past_highest_included_grain_index
            {
                let global_idx = layout.global_index_of(well_index, local_grain)?;
                let weight = population[global_idx] / total_population;

                if local_grain < well.nonreactive_grain_count {
                    continue;
                }

                for (channel_index, channel) in well.channels.iter().enumerate() {
                    let k = micro
                        .microcanonical_rate(well_index, channel_index, local_grain)
                        .max(0.0);
                    channel_rates[channel_index] += weight * k;

                    if channel.connected_well_index.is_none() {
                        total_outgoing += weight * k;
                    }
                }
            }

            per_well_per_channel_rates.push(channel_rates);
        }

        Ok(SolutionResults {
            steady_state_population: population.to_vec(),
            per_well_per_channel_rates,
            total_outgoing_rate_constant: total_outgoing,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::KB_CM;
    use crate::masterequation::microcanonical_builder::{
        ArrayMicrocanonicalProvider, MicrocanonicalNetworkData,
    };
    use crate::masterequation::reaction_network::{
        CollisionKernelImplementation, CollisionModelParams, ReactionChannel,
    };

    fn one_well_problem(
        kernel: CollisionKernelImplementation,
        solver: MultiwellLinearSolver,
    ) -> (Vec<WellDefinition>, MasterEquationSettings, MicrocanonicalNetworkData) {
        let n_grains = 80;
        let channels = vec![ReactionChannel { name: "P".to_string(), connected_well_index: None }];
        let well = WellDefinition {
            well_name: "W".to_string(),
            energy_grain_width_cm1: 20.0,
            lowest_included_grain_index: 0,
            one_past_highest_included_grain_index: n_grains,
            alignment_offset_in_grains: 0,
            nonreactive_grain_count: 0,
            collision_params: CollisionModelParams {
                lennard_jones_sigma_angstrom: 5.0,
                lennard_jones_epsilon_kelvin: 300.0,
                reduced_mass_amu: 25.0,
                alpha_at_1000K_cm1: 300.0,
                alpha_temperature_exponent: 0.85,
            },
            channels: channels.clone(),
        };
        let settings = MasterEquationSettings {
            temperature_kelvin: 500.0,
            pressure_torr: 10.0,
            boltzmann_constant_wavenumber_per_kelvin: KB_CM,
            collision_band_half_width: 20,
            collision_kernel_implementation: kernel,
            outgoing_rate_threshold: 0.0,
            internal_rate_threshold: 0.0,
            enforce_interwell_detailed_balance: false,
            linear_solver: solver,
            krylov_tolerance: 1e-12,
            krylov_max_iter: 5000,
            gmres_restart: 80,
        };
        let rho: Vec<f64> = (0..n_grains).map(|i| (1.0 + 0.05 * i as f64).powi(10)).collect();
        let k: Vec<f64> = (0..n_grains)
            .map(|i| if i < 40 { 0.0 } else { 1.0e7 * (i - 39) as f64 })
            .collect();
        let data = MicrocanonicalNetworkData {
            rho_by_well: vec![rho],
            k_by_well_by_channel: vec![vec![k]],
            channels_by_well: vec![channels],
        };
        (vec![well], settings, data)
    }

    #[test]
    fn steady_state_populations_solve_the_master_equation() {
        // Whatever the internal similarity scaling and solver, the returned populations must
        // satisfy the untransformed steady-state equation  L p + s = 0.
        for kernel in [CollisionKernelImplementation::Spd, CollisionKernelImplementation::Mess] {
            for solver in [
                MultiwellLinearSolver::Direct,
                MultiwellLinearSolver::Gmres,
                MultiwellLinearSolver::BiCgStab,
            ] {
                let (wells, settings, data) = one_well_problem(kernel, solver);
                let micro = ArrayMicrocanonicalProvider::new(&data);
                let activation = ChemicalActivationDefinition {
                    activated_well_index: 0,
                    recombination_channel_index: 0,
                };
                let engine = MasterEquationEngine::new(wells.clone(), settings.clone());
                let (results, _) = engine
                    .solve_steady_state_chemical_activation_with_diagnostics(&micro, activation.clone())
                    .unwrap();
                let p = results.steady_state_population;

                let layout = GlobalLayout::from_wells(&wells).unwrap();
                let (l, _) = assemble_raw_operator(&wells, &settings, &micro, &layout).unwrap();
                let s = engine.build_activation_source(&micro, &layout, &activation).unwrap();
                let n = l.size();
                let mut res2 = 0.0;
                let mut s2 = 0.0;
                for i in 0..n {
                    let mut lp = 0.0;
                    for j in 0..n {
                        lp += l.get(i, j) * p[j];
                    }
                    res2 += (lp + s[i]).powi(2);
                    s2 += s[i] * s[i];
                }
                let rel = (res2 / s2).sqrt();
                assert!(rel < 1e-6, "{kernel:?} / {solver:?}: |L p + s| / |s| = {rel:e}");
            }
        }
    }
}
