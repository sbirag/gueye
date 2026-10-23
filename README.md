# gueye 0.7.1 — Semi-analytic inversion of discrete Laplace/Helmholtz operators in 1D–2D–3D

`gueye` is a Python package for solving finite-difference Laplace/Helmholtz problems on rectangular grids in one, two, and three dimensions. It combines a scaled closed-form tridiagonal inverse with transverse diagonalization and supports mixed boundary conditions, selected Green-function values, frequency sweeps, discrete resonances, and time-domain responses.

## Main features

- 1D, 2D, and 3D rectangular domains using second-order finite differences.
- Face-by-face Dirichlet (`D`), Neumann (`N`), Robin (`R`), discrete transparent (`A`), and perfectly matched layer (`P`) conditions, subject to the limitations below.
- Full-field solves, Green-function evaluation, values at selected points, frequency sweeps, and discrete resonances.
- A 2D embedded Dirichlet obstacle and a Newmark time-domain solver.
- Open-domain time responses through frequency synthesis.
- Optional browser-based local GUI for field exploration, source/receiver placement, frequency sweeps, resonance display, a 2D circular obstacle, and wave animation.
- Explicit singularity checks and validation of coordinates and input parameters.

## Installation

Python 3.9 or later is required. The numerical core depends on NumPy and SciPy.

```bash
python -m pip install .
```

For development and tests:

```bash
python -m pip install -e ".[test]"
python -m pytest -q tests
python -m gueye
```

The smoke test should finish with `TOUT EST OK`. Matplotlib is optional and is only needed for plotting workflows that use it.

## Launch the graphical interface

```bash
python -m gueye.gui
# or, after installation:
gueye-gui
```

The interface opens in the default browser. The HTTP server binds to `127.0.0.1` and is designed for local interactive use only; do not expose it to a network, tunnel, or public proxy. The GUI is an exploratory visualization tool and is not a substitute for numerical validation.

## Quick example

```python
import numpy as np
import gueye as gy

# Rectangle of size 3 x 2, with 300 x 200 cells and mixed boundary conditions.
pb = gy.problem(
    size=(3.0, 2.0),
    cells=(300, 200),
    k=9.0,
    bc=[("D", ("R", 3.0)), ("N", "D")],
)

u = pb.solve(lambda x, y: np.exp(-((x - 1)**2 + (y - 1)**2) / 0.01))
g = pb.green((0.9, 0.7))
resonances = pb.resonances(6)
response = pb.sweep(np.linspace(1, 12, 500), src=(0.9, 0.7), rec=[(2.0, 1.2)])
```

A pair of boundary conditions is specified for each axis, in lower/upper order. The supported forms are:

| Form | Meaning |
|---|---|
| `"D"` or `("D", g)` | Dirichlet: `u = g` |
| `"N"` or `("N", g)` | Neumann: outward normal derivative equals `g` |
| `("R", alpha)` or `("R", alpha, g)` | Robin: outward normal derivative + `alpha*u = g` |
| `"A"` | Exact discrete transparent condition; last axis only |
| `("P", m, R0)` | PML with `m` cells and target reflection `R0` |

## Numerical safeguards and limitations

- Singular or near-singular discrete operators are detected by default; frequency sweeps return `NaN` at singular frequencies.
- Sources on Dirichlet nodes and points outside the domain are rejected. A receiver on a Dirichlet face returns the prescribed value through the full-field API and zero for the homogeneous Green response.
- The exact transparent condition is supported on the last axis only. PML can be used on all axes, but its residual reflection depends on thickness, target reflection, angle, and grid resolution; it is not a perfectly reflection-free physical boundary.
- PML is not available as a step-by-step boundary treatment in the Newmark time marcher. Open time responses are instead synthesized over complex frequencies.
- First-order operator terms are supported in 1D only; 2D and 3D use the Laplace/Helmholtz operator described here. The embedded obstacle capability is 2D.
- The domain is rectangular except for the supported 2D embedded obstacle.
- For a 3D all-Dirichlet Poisson problem requiring the full field, a discrete sine transform (DST) can be faster. `gueye` is most useful when mixed/open boundary conditions, Green values, frequency sweeps, or discrete resonances are needed.
- A practical wave-resolution guideline is about 10 grid cells per wavelength; the second-order grid-dispersion error scales approximately as `(k h)^2/24`.

## Validation

The repository includes validation scripts and reference outputs. The release-preparation environment ran 158 automated tests successfully and the package smoke test ended with `TOUT EST OK`. Validation compares the tridiagonal kernel with dense references, multidimensional results with sparse assembly, PML elimination with explicit-layer models, and selected open-domain results with enlarged domains. See `validation/` for scripts and logs. Performance figures are environment-specific and should not be treated as universal speed guarantees.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). If you use `gueye`, cite the software release and, when available, the associated SoftwareX article. Do not cite an unpublished article as published.

## Project information

- **Author:** Serigne Bira Gueye, Department of Physics, Faculty of Sciences and Techniques, Cheikh Anta Diop University, Dakar, Senegal.
- **Contact:** bira.gueye@ucad.edu.sn
- **License:** MIT; see [`LICENSE`](LICENSE) and [`Licence.txt`](Licence.txt).
- **Repository:** https://github.com/sbirag/gueye
- **Release:** 0.7.1

## AI-assisted development disclosure

ChatGPT (OpenAI) and Claude (Anthropic) were used as aids for drafting, code review and verification, troubleshooting, and technical advice during software development and manuscript preparation. The maintainer remains responsible for the software, test execution and interpretation, and all scientific claims.

A French-language version of the documentation is available in [`README_FR.md`](README_FR.md).
