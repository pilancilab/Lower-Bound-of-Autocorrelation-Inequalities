We want a **lower bound** on
\[
\inf_f\; \sup_{x\in[-2,2]} \int_{-1}^1 f(t)\, g(x+t)\,dt,
\]
under the assumptions
\[
f,g:[-1,1]\to[0,1],\qquad g=1-f \text{ on }[-1,1].
\]

Since \(g\) is determined by \(f\) on \([-1,1]\), this is really an optimization over a single function \(f\). The quantity inside,
\[
\int_{-1}^1 f(t)\,g(x+t)\,dt,
\]
measures the overlap between \(f\) and a shifted copy of its complement \(g\). For each shift \(x\in[-2,2]\), it gives the overlap at that displacement, and then the supremum over \(x\) picks the worst-case overlap. The infimum over \(f\) asks for the smallest possible such worst-case overlap.

To prove a **valid lower bound** \(L\), we need to show that every admissible \(f\) satisfies
\[
\sup_{x\in[-2,2]} \int_{-1}^1 f(t)\,g(x+t)\,dt \ge L.
\]

A natural way to do this is to replace the original infinite-dimensional problem by a **constrained optimization problem**. The idea is to impose a collection of constraints that every genuine admissible pair \((f,g)\) must satisfy, and then optimize over all objects satisfying those constraints.

This is useful because of the relaxation principle:

- every admissible \(f\) gives a feasible point of the constrained problem,
- the constrained problem is designed to be easier to analyze or solve,
- and minimizing over a larger feasible set can only decrease the value.

Therefore, if the constrained optimization problem has optimum \(L_{\mathrm{relax}}\), then automatically
\[
L_{\mathrm{relax}}
\le
\inf_f\; \sup_{x\in[-2,2]} \int_{-1}^1 f(t)\, g(x+t)\,dt.
\]
So \(L_{\mathrm{relax}}\) is a **valid lower bound** for the original problem.

In other words, to obtain a rigorous lower bound, we should construct a tractable optimization problem whose feasible set contains all admissible \(f\) and encodes necessary conditions such as
\[
0\le f\le 1,\qquad 0\le g\le 1,\qquad f+g=1 \text{ on }[-1,1],
\]
together with any additional structural constraints we can prove. Solving that constrained problem then gives a mathematically valid lower bound for the original inf-sup quantity.

So the overall strategy is:

1. start from the exact inf-sup problem,
2. derive necessary constraints satisfied by every admissible \(f\),
3. build a constrained optimization problem from those necessary conditions,
4. use its optimum as a provable lower bound.

The current SOTA is 0.379005.