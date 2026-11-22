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

## Validation status
- Full automated test suite: 166 tests passed in the local Linux environment (Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0). GitHub Actions also passed the full test and installation-smoke-test matrix on Python 3.9–3.13 (Ubuntu 24.04).
- Package smoke test: `python -m gueye`.
- Benchmark: `python validation/benchmark_reproducible.py --repeats 5`.

This release does not change the documented limitation that exact discrete transparent conditions are available only on the last axis. PML remains an approximate open-boundary treatment physically, despite exact algebraic elimination of the discrete layer.
