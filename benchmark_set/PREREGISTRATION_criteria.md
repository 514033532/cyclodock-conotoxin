# Pre-registered Benchmark Evaluation Criteria

**Date fixed: 2026-11-02 (before data inspection of the crossover design).**
**Project: Benchmarking blind and site-informed docking of conotoxin peptides against human membrane targets.**

## Design

- Benchmark set: 16 conotoxins with literature-verified activity at human targets (Table 1 of the manuscript).
- Two docking modes: (A) blind (FPocket pocket-1 protocol); (B) site-informed (AIRs anchored to established binding regions; anchor definitions graded by strength, Methods 2.2).
- Each pair scored identically with HADDOCK3 v2025.11.0 (top-cluster score).
- Comparisons only within targets (HADDOCK scores are not comparable across targets).

## Pre-registered interpretation rules

1. Blind-mode failure together with site-informed success indicates pocket misassignment and licenses a site-informed confirmation stage (Stage 2).
2. **Failure of both modes requires qualification of all ranking claims** (terminal rule).
3. Any apparent success must be checked for circularity: the docked structure and the site definition must not derive from the same experimental source.
4. PTM sensitivity criterion (pre-registered): within-target Spearman correlation between modified and unmodified scores r ≥ 0.9 and median percentile shift < 10 percentage points for qualitative robustness.

## Execution record

- 2026-11-02: blind mode executed (4/4 within-target comparisons failed; rule 1 not triggered).
- 2026-11-02: site-informed mode executed; non-circular crossover (NaV1.4/NaV1.8/CaV2.2) failed 3/3 (rule 2 triggered → all ranking claims qualified/withdrawn).
- 2026-11-02: circular test (6J8E + 2LXG) separated; ruled non-evidence by rule 3.
- 2026-11-03: PTM Tier-1 sensitivity: r = 0.80 / −0.20 / −0.30 (all below threshold; rule 4 not met → rankings not robust to PTMs).

*This file is part of the submission package and is intended for deposition in the project GitHub repository.*
