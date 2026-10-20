# Progress

**Last updated:** 26-10-15

**Current milestone:** 2 — Basis set parser

## Milestone status

| # | Milestone | Status | Test passing |
|---|-----------|--------|--------------|
| 1 | Boys function | done | yes (215) |
| 2 | Basis set parser | not started | — |
| 3 | Overlap matrix S | not started | — |
| 4 | Kinetic and nuclear attraction matrices | not started | — |
| 5 | Two-electron repulsion integrals | not started | — |
| 6 | SCF for H₂ / STO-3G | not started | — |
| 7 | SCF for H₂O / STO-3G and 6-31G | not started | — |
| 8 | DIIS acceleration | not started | — |
| 9 | UHF | not started | — |
| 10 | MP2 correlation energy | not started | — |
| 11 | Dissociation curve and analysis | not started | — |

## Working and tested

- `src/boys.py` — `boys_array(n_max, x)` returns `[F_0(x), ..., F_n_max(x)]`.
  Three branches: Taylor series, upward recurrence, asymptotic formula.
  Max measured relative error: 1.5e-15.
- `tests/boys_reference.py` — two independent oracles: `boys_reference`
  (numerical integration of the definition) and `boys_hyp1f1` (confluent
  hypergeometric function).
- `tests/test_boys.py` — 215 tests: exact value at x=0, oracle agreement,
  implementation vs oracle, and properties.
- `studies/threshold_sweep.py` — parameter sweep that produced the threshold
  band recorded below.

## Partially implemented

- (nothing yet)

## Design decisions made along the way

- 26-10-20 — `SMALL_X = 2.0`, `LARGE_X = 55`, series term budget `range(40)`.
  The error is flat at machine-precision level across the whole band
  (`SMALL_X` in [1, 3.5], `LARGE_X` in [45, inf)).
  Below the band the upward recurrence suffers from catastrophic cancellation,
  above it the Taylor series runs out of terms (with max 40 terms). 
  `LARGE_X` has a cliff only on the low side (below x ≈ 43 the asymptotic form has
  not converged; above it both methods are at machine precision and the upward
  recurrence stays exact indefinitely), so it was set well away from that point.

- 26-10-20 — Two independent oracles plus a test that compares them against
  each other. `hyp1f1` shares its derivation with the Taylor series — with the
  parameters used, the hypergeometric series IS that Taylor series — so on
  its own it could not validate that branch.

- 26-10-20 — Property tests with no oracle: array length, positivity, and
  strict decrease in n. These caught an error in the recurrence without
  any reference value.

## ORCA reference values

| System | Basis | Method | Energy (Hartree) | Input file |
|--------|-------|--------|------------------|------------|
| | | | | |

## PySCF reference values

| Quantity | System | Basis | Notes |
|----------|--------|-------|-------|
| Overlap S | H₂ | STO-3G | off-diagonal ≈ 0.6593 at 0.74 Å |
