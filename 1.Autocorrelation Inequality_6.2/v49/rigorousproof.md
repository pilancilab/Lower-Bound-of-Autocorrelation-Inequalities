# v49 — v45 pushed to higher $K$

## Goal

Extend v45 (v44 + Bochner–Toeplitz) to $K = 128, 160$ at moderate $N$
to combine the small Bochner tightening ($\sim 10^{-5}$ per point)
with the substantial Cesàro-error reduction from larger $K$
($\sim 10^{-3}$ per $K$-step).

## Setup and validity

Identical to v45 (all constraints inherited, just larger parameters).

## Status

Running in parallel with v48 (v44 at larger $K$) to see whether
Bochner–Toeplitz is still essentially free or whether it slows the
solve prohibitively at $K=128$.
