"""Run White's original v0 convex program on Table 2 row A7.

A7: (h1,h2), (p1,p2), (q1,q2) = (0, 0.08), (0, 1), (-0.05, -0.025)
    N, T, R = 10000, 4000, 10
    Paper claims: 0.38
    EMOP.py (v3+v4+v6) gives: 0.37451
"""
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import sys
import time
from pathlib import Path

V0_DIR = "/home/sykim/Theory_Harness/2.Erdos minimum overlap problem/v0"
sys.path.insert(0, V0_DIR)

import cvxpy as cp
from convex_code import WhiteConvexProblem  # type: ignore

N, T, R = 10000, 4000, 10

t0 = time.time()
prob = WhiteConvexProblem(
    N=N, L=2.0/N, T=T, R=R,
    h_1=0.0, h_2=0.08,
    p_1=0.0, p_2=1.0,
    q_1=-0.050, q_2=-0.025,
)
result = prob.solve_and_verify(
    solver=cp.MOSEK, verbose=False,
    mosek_params={"MSK_IPAR_NUM_THREADS": 8},
)
elapsed = time.time() - t0

print(f"\n=== A7 under v0 (no Toeplitz) ===")
print(f"N={N}, T={T}, R={R}")
print(f"omega         = {result['primal_optimum']:.6f}")
print(f"max_violation = {result['max_violation']:.2e}")
print(f"status        = {result['status']}")
print(f"elapsed       = {elapsed:.1f}s")
print()
print(f"Reference values:")
print(f"  paper claim    : 0.38")
print(f"  EMOP.py (v6)   : 0.374513")
print(f"  v0 (this run)  : {result['primal_optimum']:.6f}")
