# E4 step 7 — sanity check of the prediction model on Bostock's own array

Bostock's 21 template stations, median family epicentral distances 22.6-55.6 km (median 31.1 km).
The model, with the SAME registration-vs-distance fit measured on the 2010 ring and the published-stack template quality, predicts a network excess of 7.0 MAD at the M1.8-2.1 reference and a completeness magnitude at 8 MAD of **M 2.10**.
His catalogue's median event is M 1.70 (2012 SSE) and its 10th percentile M 1.43.

**The model does NOT reproduce his array** (predicted Mc above his median magnitude M1.70) -- the sigma / gain / distance assumptions are PESSIMISTIC and the region-year numbers below are upper bounds on the completeness magnitude, not predictions. Stated rather than tuned, as the prescription requires.

Model inputs (all measured in E4):
  per-station 3-C excess, surface broadband: 2.15 MAD (< 40 km), 0.83 MAD (40-60 km)
  per-station 3-C excess, borehole:          1.53 MAD (< 40 km), 1.53 MAD (40-60 km)
  amplitude vs magnitude: 10^(0.38 dM) measured WITHIN family (moment scaling would be 10^(1.50 dM))
  network degradation vs the equal-weight prediction: 0.863
  template quality q(N) with r1 = 0.054: q(500) = 0.77 against q = 0.89 for the published stacks
  detection threshold: 8 MAD on the channel-normalised network sum
