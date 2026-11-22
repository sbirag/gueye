# gueye 0.7.2 release notes

## Correctness fixes
- `Solver.resonances(n)` now includes the screening/mass term `mu**2`: the discrete eigen-wavenumbers satisfy `k_n**2 = lambda_n + mu**2`.
- Complex prescribed Dirichlet values are preserved by 1D field expansion, 2D face evaluation/expansion, and 3D face evaluation, right-hand-side assembly, and expansion.
- Point extraction with `at()` now returns complex Dirichlet face values without discarding their imaginary parts.

## Regression coverage
- Added tests for nonzero `mu` in 1D, 2D, and 3D, including a check that the shifted resonance is detected as a singular operator.
- Added tests for constant complex Dirichlet data in 1D–3D and coordinate-dependent complex data in 2D–3D.

## Reproducible performance measurements
- Added `validation/benchmark_reproducible.py`, which reports modal setup, modal solve, sparse matrix assembly, sparse factorization, and sparse solves for one and multiple right-hand sides separately, for both homogeneous and mixed Dirichlet–Neumann cases in 2D and 3D.
- The five-repeat benchmark results are recorded in `validation/benchmark_results_0.7.2.csv`, including the measured environment and numerical relative errors. Timings are machine-specific and are not performance guarantees.

## Independent reference validation added for the SoftwareX review
- `validation/test_independent_references.py` checks complete Green columns against independently assembled sparse systems in 1D–3D (maximum relative errors: 2.42e-14, 1.45e-14, 9.88e-16).
- It checks three-frequency sweeps against separate sparse solves at each complex wavenumber (3.28e-14, 1.02e-14, 6.44e-16).
- It checks discrete resonance wavenumbers against eigenvalues of independently assembled matrices (5.33e-15, 1.94e-15, 1.91e-15).
- It checks the 2D Dirichlet-obstacle capacitance formulation against direct sparse elimination for circular and rectangular obstacles on two grids (maximum field deviation 1.58e-14; maximum obstacle-node field 7.22e-16). These are discrete-algebraic checks, not continuum obstacle validation.
- Open-boundary tests compare PML fields with the discrete transparent condition over selected 1D, 2D and 3D cases. The reported deviations are configuration-dependent field errors, not reflection coefficients. Only the 1D explicit-layer elimination comparison is included as a reproducible exact algebraic-elimination check.

## Validation status
- Full automated test suite: 166 tests passed in the local Linux environment (Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0). The CI workflow is configured for the 10 Linux/Windows × Python 3.9–3.13 matrix, the dedicated independent scientific-validation job, manuscript compilation, source-distribution contents checks, wheel installation and `python -m gueye` smoke test. The final-candidate CI run must pass after all dossier changes before the candidate is considered frozen; its run record is available in the repository's GitHub Actions history. macOS is not yet validated.
- Package smoke test: `python -m gueye`.
- Benchmark: `python validation/benchmark_reproducible.py --repeats 5`.

This release does not change the documented limitation that exact discrete transparent conditions are available only on the last axis. PML remains an approximate open-boundary treatment physically, despite exact algebraic elimination of the discrete layer.
