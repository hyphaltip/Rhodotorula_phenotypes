#!/usr/bin/env python3
"""Compute the GRM eigenvalue spectrum / condition number for a GEMMA kinship matrix
and flag singularity risk before trusting any downstream kinship-only LMM p-value
(spec S4, quant-genetics review fix #3 -- mandatory on ANY strain-state change).
"""
import argparse
import json
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kinship", required=True, help="GEMMA gwas.cXX.txt output")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    K = np.loadtxt(args.kinship)
    n = K.shape[0]
    assert K.shape == (n, n), f"kinship matrix not square: {K.shape}"

    eigvals = np.linalg.eigvalsh(K)
    eigvals_sorted = np.sort(eigvals)[::-1]
    max_eig, min_eig = eigvals_sorted[0], eigvals_sorted[-1]
    # Guard against a literal zero/negative min eigenvalue (numerically singular)
    # inflating the condition number to infinity uninformatively.
    eps = 1e-10
    condition_number = max_eig / max(min_eig, eps)
    n_near_zero = int(np.sum(eigvals_sorted < 1e-6))

    # Threshold informed by the original D-9 finding: singularity there arose from
    # 22/201 (~11%) near-clone strains. Flag risk at >5% near-zero eigenvalues or
    # condition number > 1e6 (GEMMA's GSL solver failed outright in the original
    # 10-PC/3-PC covariate runs at comparable conditioning).
    frac_near_zero = n_near_zero / n
    verdict = "singular_risk" if (frac_near_zero > 0.05 or condition_number > 1e6) else "well_posed"

    report = {
        "n": n,
        "max_eigenvalue": float(max_eig),
        "min_eigenvalue": float(min_eig),
        "condition_number": float(condition_number),
        "n_near_zero_eigenvalues": n_near_zero,
        "fraction_near_zero": frac_near_zero,
        "verdict": verdict,
    }
    with open(args.out, "w") as f:
        json.dump(report, f, indent=2)

    print(f"n={n}  condition_number={condition_number:.3e}  near_zero_eigvals={n_near_zero} ({frac_near_zero:.1%})")
    print(f"Verdict: {verdict}")
    if verdict == "singular_risk":
        print("WARNING: kinship-only LMM may still be reliable (D-9 ran cleanly despite "
              "singular-adjacent conditioning), but this MUST be reported in GWAS.md "
              "alongside the Tier A/B/C results, not silently passed over.")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
