# SoftwareX review-response map — gueye 0.7.2

**Manuscript:** “gueye: semi-analytic inversion of discrete Laplace/Helmholtz operators in 1D–2D–3D”  
**Manuscript number:** SOFTX-D-26-01518  
**Purpose:** technical transfer dossier for a future submission. This document records the changes and evidence; it is not a claim that the rejected submission has been formally revised or accepted.

## Reviewer 1

### R1.1 — Resonance wavenumbers and the screening term

**Comment.** The reported resonance wavenumber did not match the singularity of the discrete Helmholtz matrix when the screening parameter `mu` was nonzero.

**Change.** Corrected `Solver.resonances()` to include the shift `mu**2`. The discrete relation is now stated as
`k_n^2 = lambda_n + mu^2`,
where `lambda_n` is an eigenvalue of the discrete negative Laplacian.

**Evidence.** Regression tests cover nonzero `mu` in 1D, 2D and 3D, and verify singularity detection at the shifted resonance. The independent-reference validation compares API resonance wavenumbers with eigenvalues of separately assembled discrete matrices; maximum relative deviations are $5.33\times10^{-15}$, $3.51\times10^{-15}$, and $2.38\times10^{-15}$ in 1D, 2D and 3D, respectively.

### R1.2 — Transparent boundary scope and complex boundary data

**Comment.** The exact transparent condition was only implemented on the last axis; the handling of complex prescribed boundary values needed clarification.

**Change.** The API, README and manuscript now explicitly state that the exact discrete transparent condition `A` is supported on the last axis only. PML can be used on all axes. Complex prescribed Dirichlet data are preserved in 1D–3D boundary evaluation, right-hand-side assembly, field expansion and selected point extraction.

**Evidence.** Regression tests cover constant complex Dirichlet data in 1D–3D and coordinate-dependent complex data in 2D–3D. Open-boundary scripts compare selected 1D, 2D and 3D PML cases against discrete transparent references. The manuscript reports field deviations and does not label them as reflection coefficients.

### R1.3 — Reproducibility of the performance table

**Comment.** The benchmark table needed reproducible scripts, explicit error norms and separate timing of setup/factorization and solve stages.

**Change.** Added `validation/benchmark_reproducible.py` and `validation/benchmark_results_0.7.2.csv`. The script separates modal setup, modal solve, sparse assembly, sparse LU factorization and sparse solves for one and four right-hand sides. The manuscript identifies the environment and does not claim universal speed superiority.

**Evidence.** Five-repeat measurements and relative errors are recorded in the CSV. The manuscript notes that a sparse solve can be faster for one right-hand side on the tested small grids and identifies the discrete sine transform as a specialised comparator for the all-Dirichlet Poisson case.

## Reviewer 2

### R2.1 — Manuscript completeness and readability

**Change.** Reworked the SoftwareX manuscript into a complete source document with abstract, required software metadata, motivation/significance, architecture, functionalities, examples, validation, performance, impact, limitations, conclusions, declarations and references. The PDF is compiled from the tracked TeX source in CI.

**Evidence.** The build workflow compiles the manuscript with LaTeX, checks that the PDF has at least eight pages, extracts its text and checks for the Validation and Limitations sections.

### R2.2 — Contribution, background, impact, reuse and use cases

**Change.** The manuscript distinguishes the package contribution from standard tensor-product and transform methods: a scaled closed-form tridiagonal inverse, transverse diagonalisation, mixed face conditions, selected Green values, frequency sweeps and discrete resonance safeguards. It includes a software architecture description, usage examples, a local GUI example and explicit impact and limitation statements.

**Qualification.** The package is intended for rectangular finite-difference Laplace/Helmholtz problems and the stated 2D embedded-obstacle use case. It is not presented as a general-purpose solver for arbitrary geometries.

### R2.3 — Scientific validation, comparisons and reproducibility

**Changes and evidence.**

- **Discrete algebraic verification:** 1D kernel checks against dense references; 2D and 3D modal solutions against independently assembled sparse systems.
- **Continuous-solution convergence:** manufactured-solution tests show second-order convergence in the tested 1D, 2D and 3D configurations.
- **Green functions:** complete Green columns are checked against independent sparse direct solves with the source scaled by the cell volume.
- **Frequency sweeps:** three complex wavenumbers in each dimension are checked against separate sparse direct solves at each frequency.
- **Discrete resonances:** API values are checked against eigenvalues of independently assembled matrices in 1D–3D.
- **Dirichlet obstacles:** the 2D capacitance formulation is checked against direct sparse elimination for circular and rectangular obstacles on two grids.
- **Open boundaries:** selected PML field solutions are compared with discrete transparent references in 1D–3D; 1D time responses are compared with a larger-domain Newmark reference.
- **Time domain:** Newmark modal solutions are compared with sparse LU time stepping in 2D and 3D.

The independent-reference script is `validation/test_independent_references.py`; the validation scripts are run by the dedicated GitHub Actions job and their text outputs are uploaded as a workflow artifact.

**Scope limitation.** These checks establish agreement with independently assembled discrete systems and convergence for the tested manufactured solutions. They do not constitute exhaustive frequency, angle, PML-parameter or geometry studies, nor do they validate a continuum obstacle model for arbitrary shapes.

### R2.4 — Algorithmic detail, stability and computational cost

**Change.** The manuscript now gives the tridiagonal inverse formula, describes scaled evaluation regimes, explains transverse diagonalisation in 2D and 3D, states the role of the modal stack, documents the singularity threshold and distinguishes setup from solve costs.

**Qualification.** The manuscript explicitly notes that a general eigensolver is used for nonsymmetric transverse operators, that its setup can dominate at large transverse dimensions, and that modal arrays require storage proportional to the grid size. The reported benchmark is representative rather than an asymptotic performance claim.

### R2.5 — Software identity, installation and cross-platform reproducibility

**Change.** The repository metadata, version, license, citation metadata, README, installation instructions, source distribution and wheel checks have been aligned to version 0.7.2. The workflow tests Python 3.9–3.13 on Ubuntu 24.04 and Windows Server 2022, builds the wheel and source distribution, installs the wheel with dependencies and runs `python -m gueye`.

**Evidence.** The source distribution is checked for the manuscript, tests and validation scripts. The dedicated scientific-validation and distribution-build jobs are part of CI.

**Remaining archival step.** A public version tag/release and a permanent archival identifier (for example, a Zenodo DOI) have not been created. They must be added before a submission that requires a DOI or immutable release record. macOS has not been tested.

### R2.6 — Feature limitations and boundary-condition scope

**Change.** The README and manuscript explicitly state the current limits:

- exact discrete transparent condition `A`: last axis only;
- PML: approximate physical absorption with configuration-dependent residual field error, despite exact algebraic layer elimination in the tested 1D comparison;
- Newmark: no step-by-step PML; open-domain time responses use frequency synthesis;
- first-order operator term: 1D only;
- embedded Dirichlet obstacle: 2D only;
- rectangular computational box, with the stated 2D embedded-obstacle extension.

The manuscript distinguishes algebraic elimination error from physical PML absorption error and reports the latter as field deviations against a discrete transparent reference, not as reflection coefficients.

## Final transfer checklist

- [x] Resonance correction and regression tests.
- [x] Complex Dirichlet data regression tests.
- [x] Independent Green, sweep and resonance sparse-reference checks.
- [x] Independent 2D obstacle sparse-elimination checks for two shapes and two grids.
- [x] Reproducible benchmark script and CSV.
- [x] PML/transparent and time-domain validation scripts with stated scope.
- [x] Python 3.9–3.13 Linux/Windows test matrix.
- [x] Source distribution and wheel build/install smoke test in CI.
- [x] Manuscript source compiled and checked in CI.
- [ ] Create an immutable public release tag after author approval.
- [ ] Deposit the frozen source archive and manuscript in a permanent repository and record its DOI.
- [ ] Final author proofread of the compiled PDF, references, figures and journal formatting.
- [ ] Prepare a journal-specific cover letter and final transfer/submission package.

