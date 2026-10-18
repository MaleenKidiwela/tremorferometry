---
name: dont-over-conclude
description: "User feedback — don't jump to dramatic conclusions from a single result; verify the setup first"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: de9e1cec-3d07-493e-8e0a-d25339535f7a
---

The user has flagged twice (2026-07-21) that I **jump to conclusions too easily**. Concrete case: I ran the deep-resolution checkerboard, got a collapse (interface recovery 0.06–0.69), and declared "the deep dv/v is a margin-wide INDEX, not a map" — a dramatic verdict — before verifying the test setup. It was WRONG: a λ-selection bug (atlas PSF-floor λ transplanted into the checkerboard → over-smoothing) had artificially depressed recovery. The corrected oracle-λ result was ~70–140 km resolvable map, sharpening with coverage, geometry-unlimited (noise-free corr ~1.0). The user's intuition ("it should recover more like the older case") was right.

**Why:** a surprising/dramatic result is more likely a setup error than a real discovery. Especially when a big conclusion (a null, a "can't do X", a retirement of prior work) hinges on ONE parameter or choice, that choice must be audited before the conclusion is stated.

**How to apply:**
- When a result would overturn prior work or yield a dramatic negative verdict, treat it as SUSPECT until the setup is audited. Run the ablation / sensitivity-to-key-choices FIRST, then conclude.
- Present surprising results as provisional ("this suggests X, pending a check of Y"), not as settled verdicts.
- Consult [[use-merlin-for-decisions]] to AUDIT THE CODE (not just the framing) before declaring a big conclusion — Merlin caught the λ bug when asked to be skeptical of the implementation.
- Report the sanity anchors (e.g., does the synthetic signal amplitude match the observed signal class?) — they catch mis-calibration fast.
- Match the strength of the claim to the strength of the verification done.

**Instance 2026-09-10 (LFE catalog):** I generalised a METHOD-specific limit (2-s template matching: one LFE = 1.8 MAD per station)
into "individual LFEs cannot be detected at single stations anywhere", then told the user a network-wide catalog was impossible.
Counterexample existed in the literature I had already cited (Lin et al. 2024: > 1 M LFEs from a single-station deep-learning picker +
association). User: "I disagree that u can't build a reliable catalog across the entire network". Also went in circles by not
re-reading the June lfe_features/ thread before restarting. **How to apply:** before declaring something impossible, (1) separate
the method's limit from the physics, (2) search for the published counterexample, (3) re-read the project's own prior thread on the
topic (grep notes/ and study folders) and say what it already established.
