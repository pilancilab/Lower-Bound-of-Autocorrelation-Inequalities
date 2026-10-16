Hi Claude, our mission is to find rigorous advancements of a series of Autocorrelation inequalities that has been proposed. Imagine you are an expert in functional analysis, optimization, Fourier analysis, and combinatorics.

Your job is to propose a convex-optimization like formulation that may be motivated from the form proposed in White's paper (the paper directory is in 2201.05704v1.pdf) and constantly improve it. Basically you are going to show a lower-bound or an upper-bound that is achieved by relaxing the complicated problem into a decomposed set of constrained convex problems, then optimize each problem. So your main job is to find novel constraints that were overlooked by previous versions of the program - and use it to optimize. See 2201.05704v1.pdf for an insight: Lemma 3,4,5 is used to construct a convex optimiation problem in section 5.

The problem that you need to solve is specified in Problem.md. Also in Problem.md there are tips for the problem. Note that we are trying to tackle these problems in a similar approach - but we are not trying to solve only the erdos minimum overlap problem.

Your workflow will be the following:

0. This is the initial step. See if you have folder v0. If you do, then run the code as a baseline to check that it works. If there is no v0, you can come up with the most naive version of the convex optimization problem that can be meaningful.

After this you have an iterative loop.

1. Read the most recent version folder. Then, create {nextversion} folder and copy the contents from the previous folder.
2. Come up with a novel proposal that improves upon the previous implementation. The proposal should not be something like a hyperparameter tweak; we will enhance the hyperparameters after the convex program has converged. 
3. After you come up with a proposal, you have to generate a "rigorous" proof that the constraint indeed holds. For an example see Lemma 3,4,5 in the convex optimization problem of White and how it is related to the convex optimization problem in section 5.
4. Store the proof in the current version folder. Name it as rigorousproof.md
5. Now this is the important: you should be idle until the theorist agent reads your proof and gives a verdict. It will be in the same folder, verdict.md. If the verdict is VALID, it means your proposed proof is rigorous and makes sense. If the verdict is INVALID, it means your proposed proof has a flaw.
6. IF the verdict was VALID, implement the constraint in the convex problem.
7. IF the verdict was INVALID, go to 2 and propose a different constraint that makes mathematical sense.

*** THINGS TO CONSIDER ***
- A validation should take less than 10 minutes. Don't heavily optimize the parameter space (split the original problem into too many different convex optimization problems) or increase the problem parameters too much. Our job is NOT to improve the bound for now, but to improve the PROGRAM that will eventually lead to better bounds. For this problem specifically, choose a hard square region h_1 = h_2 = 0.015, p_1 = p_2 = 0.381, -q_1=q_2=0.02 to verify and see if the bound is improving.

- MATHEMATICAL RIGOR IS THE MOST IMPORTANT!!! Do not simply assume that you are doing an approximation that is good enough or do not propose things that may or may not be a mathematically valid lower or upper bound.

- NEVER STOP: Once the experiment loop has begun (after the initial setup), do NOT pause to ask the human if you should continue. Do NOT ask "should I keep going?" or "is this a good stopping point?". The human might be asleep, or gone from a computer and expects you to continue working indefinitely until you are manually stopped. You are autonomous. If you run out of ideas, think harder — read papers referenced in the code, re-read the in-scope files for new angles, try combining previous near-misses, try more radical architectural changes. The loop runs until the human interrupts you, period.

- Especially, I will not initiate a new interaction when the verdict.md is given. So when I say idle I mean don't stop the current session and wait.