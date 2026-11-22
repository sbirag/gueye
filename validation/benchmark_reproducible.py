"""Reproducible setup/factorization/solve benchmark for gueye 0.7.2.

Run from the repository root after installing the package:
    python validation/benchmark_reproducible.py --repeats 7

The direct reference is split into sparse-matrix assembly, SuperLU
factorization, and triangular solves. The gueye setup (including transverse
modal decomposition) is timed separately from repeated field solves.
All reported times are medians; this is a reproducibility aid, not a
hardware-independent performance claim.
"""
from __future__ import annotations
import argparse
import csv
import platform
import statistics
import time
from pathlib import Path

import numpy as np
import scipy
import scipy.sparse.linalg as spla
import gueye


def median_seconds(fn, repeats: int):
    samples = []
    result = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = fn()
        samples.append(time.perf_counter() - t0)
    return statistics.median(samples), result


def benchmark(label, size, cells, repeats):
    rng = np.random.default_rng(20261122)
    full_shape = tuple(n + 1 for n in cells[::-1])
    F = rng.standard_normal(full_shape)
    bc = [("D", "D")] * len(size)

    setup_s, pb = median_seconds(
        lambda: gueye.problem(size, cells, bc, k=0.0, mu=0.0, check=False),
        repeats,
    )
    modal_s, u = median_seconds(lambda: pb.solve(F), repeats)
    interior = pb.backend.rhs(F)
    if len(size) == 1:
        axes = [pb.backend]
    elif len(size) == 2:
        axes = [pb.backend.Py, pb.backend.Px]
    else:
        axes = [pb.backend.Pz, pb.backend.Py, pb.backend.Px]
    slices = tuple(slice(P.i0, P.i1 + 1) for P in axes)
    u_int = u[slices]

    assembly_s, A = median_seconds(lambda: pb.backend.dense().tocsc(), repeats)
    factor_s, lu = median_seconds(lambda: spla.splu(A), repeats)
    sparse_s, u_ref = median_seconds(lambda: lu.solve(interior.ravel()), repeats)
    multi_rhs = np.column_stack([interior.ravel(), 0.5 * interior.ravel(), -interior.ravel(), 2.0 * interior.ravel()])
    multi_s, _ = median_seconds(lambda: lu.solve(multi_rhs), repeats)
    u_ref = u_ref.reshape(interior.shape)
    rel = float(np.max(np.abs(u_int - u_ref)) / max(np.max(np.abs(u_ref)), 1e-300))

    return {
        "case": label,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "repeats": repeats,
        "size": "x".join(map(str, cells)),
        "unknowns": int(interior.size),
        "setup_modal_s": setup_s,
        "gueye_solve_1rhs_s": modal_s,
        "sparse_assembly_s": assembly_s,
        "sparse_factorization_s": factor_s,
        "sparse_solve_1rhs_s": sparse_s,
        "sparse_solve_4rhs_s": multi_s,
        "relative_error": rel,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--output", default="validation/benchmark_results_0.7.2.csv")
    args = parser.parse_args()
    if args.repeats < 3:
        parser.error("--repeats doit etre >= 3")

    print(f"gueye={gueye.__version__}; Python={platform.python_version()}; "
          f"NumPy={np.__version__}; SciPy={scipy.__version__}; "
          f"platform={platform.platform()}; repeats={args.repeats}")
    cases = [
        ("2D-DDDD", (1.0, 1.0), (40, 32)),
        ("3D-DDDDDD", (1.0, 1.0, 1.0), (12, 10, 8)),
    ]
    rows = [benchmark(name, size, cells, args.repeats) for name, size, cells in cases]
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print("\n" + row["case"])
        for key, value in row.items():
            if key != "case": print(f"  {key}: {value}")
    print(f"\nCSV: {out}")


if __name__ == "__main__":
    main()
