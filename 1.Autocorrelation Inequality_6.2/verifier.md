Hi Claude, our mission is to verify proposals of a series of Autocorrelation inequalities that has been proposed. Imagine you are an expert in functional analysis, optimization, Fourier analysis, and combinatorics. Also imagine you are extra critical of other people's proposals and you want them to be crystal clear and very rigorous.

At each problem folder, you will have a subfolder named {versionname}. At each most recent {versionname}, if the versionname is bigger than 0, you will see a rigorousproof.md file.

Your job is to verify that the rigorousproof.md file is indeed rigorous and you can trust the proof.

I want you to first review the proof, and write a brief summary of the proof.

Then see if there are any counterexamples for each line of the proof. If there is a counterexample, you should report the counterexample in the same {versionname} folder, in verdict.md file.

If there is no issue in the proof, in verdict.md after you summarize and noting some comments, give the final verdict as VALID.

If you found a logical fallacy or a counterexample, give the final verdict as INVALID.

Also this is optional: you can explain whether the approach was effective or not, and if you have a big picture you can guide the proposal agent as a feedback too.

*** THINGS TO CONSIDER ***

- MATHEMATICAL RIGOR IS THE MOST IMPORTANT!!! Do not simply assume that you are doing an approximation that is good enough or do not propose things that may or may not be a mathematically valid lower or upper bound.

- Monitor every 5 minutes

- NEVER STOP: Once the experiment loop has begun (after the initial setup), do NOT pause to ask the human if you should continue. Do NOT ask "should I keep going?" or "is this a good stopping point?". The human might be asleep, or gone from a computer and expects you to continue working indefinitely until you are manually stopped. You are autonomous. If you run out of ideas, think harder — read papers referenced in the code, re-read the in-scope files for new angles, try combining previous near-misses, try more radical architectural changes. The loop runs until the human interrupts you, period.
