# Parallel master-equation runs over the conditions (T, p)

**MarXus, 2026-10-06.** Peter's requests:
- "do you use Rayon for parallelization? because for example different temperature and/or pressures can be independently, embarrassingly parallel"
- "build first the parallelization with Rayon as you find the best way, then rerun the validation test systems with the parallel code (also in the output display how many cores used for the run)"

## 1. Design

**The conditions are independent.** Every condition (T, p) is an independent master equation: its own operator J, factorization, solution or eigendecomposition. The conditions are therefore computed concurrently.

**`src/masterequation/parallel_conditions.rs`: `ConditionPool`.**
- **A local rayon thread pool**, so there is no global state and different pool sizes can be tested in one process.
- **Thread count:** `--threads N` (now `NCores` in the deck or `--ncore N`, Section 5), otherwise RAYON_NUM_THREADS, otherwise all logical cores (rayon's default).
- **`map_conditions(temperatures, pressures, f)`:** f(T, p) for every condition; the results come in grid order (temperatures outer, pressures inner). An indexed parallel iterator collects in order, so the output does not depend on the thread count or on the order of completion.

**BLAS threads** (`numeric/lapack_interface.rs::set_blas_threads`, wrapping `openblas_set_num_threads` / `openblas_get_num_threads` of the system OpenBLAS):
- The conditions leave $\max(1,\ \mathrm{threads}/\min(\mathrm{threads}, \mathrm{conditions}))$ threads to OpenBLAS: 1 when there are at least as many conditions as threads.
- Workers × BLAS threads therefore stay within the requested core count. Without this, 4 workers × 4 OpenBLAS threads would use 16 cores.

**Where the parallelism is.**
- **Parallel:** the program (`examples/chemical_activation_from_deck.rs`) runs each of its four per-condition loops through `map_conditions`: steady states, thermal eigenpair, thermal fates of the wells, CSE.
- **Unchanged:** a condition without a result is still reported as "not available" without stopping the others. Warnings and the machine-readable blocks are written afterwards, in grid order.
- **Sequential:** the library run functions (`run_chemical_activation`, …). Otherwise `cargo test --test-threads=4` plus a rayon pool would exceed the 4-core limit.

**Output.** RUN SETTINGS shows "Parallel: N threads (rayon) over the M conditions (T, p); K conditions at a time / BLAS threads per LAPACK call: B".

**Dependency.** `rayon = "1"` (1.12 with rayon-core 1.13, crossbeam and either; pure Rust, from the local cargo cache). This is within the easy-install policy: `cargo build` fetches and compiles it.

## 2. Tests (TDD, each seen failing on a stub)

| test | checks |
|---|---|
| `parallel_conditions::the_pool_has_the_requested_number_of_threads` | 1 and 2 threads; 0 is refused |
| `…::results_come_in_grid_order_whatever_the_order_of_completion` | later conditions finish first (sleep 1000/T ms), results still in grid order |
| `…::parallel_master_equations_equal_the_sequential_ones` | two-well network, 2 T × 3 p, intermediate steady state: every channel flux bitwise identical to the sequential run |
| `lapack_interface::the_blas_thread_count_can_be_set_and_lapack_still_agrees` | after setting 1 BLAS thread, DSYEVD agrees with Householder/QL to 10⁻¹⁰ |

**Program check.** Case 2 with 1 and with 4 threads, steady state (final + thermal + fates) and CSE: the machine-readable files and the tables files are byte-identical (`cmp`).

## 3. Timings

Wall clock, one run of the example including the setup (state counting, k(E)).

| system | method | conditions | sequential | 4 threads |
|---|---|---|---|---|
| ZZ-allyl + O₂ (4 wells, 1709 grains) | intermediate steady state | 21 | 3.0 s | |
| | final + thermal eigenpair (inverse) + well fates | 21 | 18.3 s | 5.7 s (3.2×) |
| | CSE, LAPACK | 21 | 10.1 s (4 BLAS threads), 18.7 s (1 BLAS thread) | 8.2 s |
| C₂H₃ full deck (2634 grains) | intermediate steady state | 40 | 19.9 s | |
| | final + thermal (inverse) + well fates | 40 | 95.8 s | 26.4 s (3.6×) |
| | final + thermal (LAPACK) + well fates | 40 | 123.9 s | |
| | CSE, LAPACK | 40 | 72.6 s (4 BLAS threads) | 50.9 s (1.4×) |

The dense LAPACK runs gain less, because DSYEVD on 2634² matrices already used 4 BLAS threads in the sequential run.

## 4. Validation re-run with 4 threads

All three validation directories were re-run with `--threads 4` in the run scripts (`RAYON_NUM_THREADS=4` and `OPENBLAS_NUM_THREADS=4` stay set).

**Time** (sequential before, 4 threads after):

| directory | sequential | 4 threads |
|---|---|---|
| Case 2 | 76 s | 31 s |
| C₂H₃ eigen | 3 min 40 s | 1 min 32 s |
| C₂H₃ intermediate | 2 min 38 s | 46 s |

**Results:**
- **Identical:** all comparison CSVs of the steady-state runs, and every machine-readable file not involving LAPACK.
- **The two LAPACK-based outputs changed at the rounding level,** because the BLAS thread count per call changed from 4 to 1, which changes DSYEVD's summation order:
  - **Case 2 CSE tables:** entries above 1% of their row maximum are unchanged at the printed precision; entries above 10⁻⁴ of the row maximum change by at most 4.2·10⁻⁶; entries above 10⁻⁶ by at most 2.2·10⁻⁵ (R → P7, about 10⁻¹⁸ cm³/s). Smaller entries are the documented rounding noise (10⁻⁹ … 10⁻²³) and change by large factors.
  - **C₂H₃ LAPACK eigen run:** only the λ₁ values far below the double-precision floor changed (warnings). `comparison_table.csv` is identical.

## 5. Number of cores in the deck: `NCores` (2026-10-06)

**Request (Peter):** "we also need in the input a ncore = 8 variable to tell the code how many processors to use for the run and batch accordingly the given ncore the T,P jobs".

**Keyword.** `NCores N` in the `MarXus ... End` block of the deck header, in the CamelCase style of the other keywords there.
- **Storage.** It goes into `SolutionSettings::cores`. It is a run setting of every method, so it is never listed as an unused setting.
- **Command line.** `--ncore N` overrides it, like every command-line option overrides its deck keyword (`SolutionSettings::overridden_by`). The deck may also have no `NCores`. (`--ncore` replaced the earlier `--threads` on Peter's request, 2026-10-06; `--threads` no longer exists.)
- **Errors.** 0, a non-integer or a negative number is an error that names `NCores`.

**Batching.**
- `ConditionPool::new(cores)` makes a rayon pool of N worker threads. The conditions (T, p) are computed in batches of up to N at a time, and the results are kept in grid order.
- LAPACK calls get the cores the conditions leave free: BLAS threads = N / (conditions at a time), at least 1. Workers × BLAS threads therefore stay within N.
- **Without `NCores` and `--ncore`:** RAYON_NUM_THREADS, otherwise all logical cores.

**Report.** The `Parallel:` line of RUN SETTINGS gives the number of cores and where it came from (`NCores in the deck`, `--ncore`, or the default), then the batch size.

**Tests:**
- `mess_input::tests::the_marxus_header_block_gives_the_number_of_cores`: `NCores 8` is read; 0, 2.5 and −1 are refused.
- `solution_method::tests::the_number_of_cores_is_a_run_setting_of_every_method`: the deck value, the command-line override, never unused for any of the four methods, and 0 refused.
