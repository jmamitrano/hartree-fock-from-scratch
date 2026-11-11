# Progress

> Update at the end of each working session. Two minutes here saves twenty
> re-establishing context next time.

**Last updated:** 26-11-11

**Current milestone:** 3 — Overlap matrix S

## Milestone status

| # | Milestone | Status | Test passing |
|---|-----------|--------|--------------|
| 1 | Boys function | done | yes (230) |
| 2 | Basis set parser | done | yes (4) |
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
  Three branches: Taylor series, `erf` + upward recurrence, asymptotic formula.
  Max measured relative error: 1.5e-15 (worst case sits at the SMALL_X
  boundary). Holds at 1.6e-13 for `n_max = 8`, so the implementation would
  survive d functions unchanged.
- `tests/boys_reference.py` — two independent oracles: `boys_reference`
  (numerical integration of the definition) and `boys_hyp1f1` (confluent
  hypergeometric function).
- `tests/test_boys.py` — 230 tests: exact value at x=0, oracle agreement,
  implementation vs oracle, and properties.
- `studies/threshold_sweep.py` — parameter sweep that produced the threshold
  band recorded below.
- `src/basis.py` — `build_basis(z_list, coords_list, basis_name)` returns a
  flat list of `BasisFunction` (frozen dataclass: `center`, `lmn`,
  `exponents`, `coefficients`).
  Coordinates in bohr, shape (n_atoms, 3). Stored coefficients are the final
  `c * d_k * N_k`. H₂O gives 7 functions in STO-3G and 13 in 6-31G.
- `tests/basis_reference.py` — `contracted_norm_reference(L, exponents,
  coefficients)`: norm of a contracted function by radial `quad` times the
  angular factor (4π for s, 4π/3 for p).
- `tests/test_basis.py` — 4 tests: function count per basis, and unit norm of
  every function against the quad reference.

## Partially implemented

- (nothing yet)

## Design decisions made along the way

- 26-10-15 — `SMALL_X = 2.0`, `LARGE_X = 55`, series term budget `range(40)`.
  The error is flat at machine-precision level across the whole band
  (`SMALL_X` in [1, 3.5], `LARGE_X` in [45, inf)), so error alone does not
  discriminate. The tie-breaker was margin to the failure modes: below the
  band the upward recurrence suffers cancellation, above it the Taylor series
  runs out of terms (cliff at x ≈ 4). `SMALL_X = 2.0` uses 25 of 30 terms with
  the old budget, so the budget was raised to 40 to restore margin. `LARGE_X`
  has a cliff only on the low side, so it was pushed above the crossover.
  Band measured with `studies/threshold_sweep.py`.

- 26-10-15 — Two independent oracles plus a test that compares them against
  each other. `hyp1f1` shares its derivation with the Taylor series — with the
  parameters used, the hypergeometric series *is* that Taylor series — so on
  its own it could not validate that branch. The bridge test settles this
  empirically instead of relying on trust.

- 26-10-15 — `epsabs=0` in `quad`. `epsabs` acts as a floor: when the integral
  is smaller than that tolerance, `quad` stops refining and the relative error
  explodes (1.9e-04 at x=500, n=8 with `epsabs=1e-15`). With `epsabs=0` and
  `epsrel=1.5e-13` (the tightest scipy accepts when `epsabs=0`, bounded by
  50*eps) accuracy holds down to F_n ≈ 1e-36.

- 26-10-15 — Property tests with no oracle: array length, positivity, and
  strict decrease in n. These caught an off-by-one in the recurrence without
  any reference value. Pattern to reuse: S symmetric with ones on the diagonal
  in milestone 3, the eight permutation symmetries in milestone 5.

- 26-10-15 — Test x values are derived from `SMALL_X` and `LARGE_X` rather than
  hardcoded, so the boundary probes follow the thresholds if they change. Note
  the sweep script must *not* do this: there the evaluation grid has to stay
  fixed across candidates, or the comparison is meaningless.

- 26-11-11 — Basis Set Exchange returns exponents and coefficients as strings.
  They are converted with `np.array(..., dtype=float)` before any arithmetic.

- 26-11-11 — SP shells (`angular_momentum = [0, 1]`) share one exponent list
  and carry one coefficient list per L, in the same order as
  `angular_momentum`.

- 26-11-11 — Two-level normalization. Primitives:
  N_k = (2α/π)^(3/4) (4α)^(L/2), valid for L ≤ 1 (the double-factorial factor
  is 1 for s and p). Contraction: c = (dᵀ S d)^(-1/2), with
  S_kk' = (4 α_k α_k' / (α_k + α_k')²)^((2L+3)/4) between *normalized*
  primitives and raw BSE d_k (which do not include N_k; passing d_k N_k here
  would count N twice).

- 26-11-11 — Contractions are normalized even though BSE data already come
  normalized to ~1e-10 (STO-3G H: 1 + 7.0e-11; largest deviation 6.6e-10 in
  the inner SP s of O in 6-31G; The deviation is due to the rounding of the 
  10-digit tabulated coefficients).

- 26-11-11 — `BasisFunction.coefficients` stores the product c · d_k · N_k,
  computed once. The reference receives these final coefficients, so it
  integrates exactly the function `basis.py` built and catches an error in
  any of the three factors.

- 26-11-11 — `rtol=1e-11` in the normalization test. Floor: quad accuracy
  (~1e-13). Ceiling: the smallest deviation that appears if c is skipped
  (~2e-11, p shell of O in STO-3G).

- 26-11-11 — Unsupported shells fail loudly. A general contraction (one L,
  several coefficient lists, e.g. H in cc-pVDZ) raises `ValueError`. 
  Shells with L > 1 also raise `ValueError`.

## ORCA reference values

| System | Basis | Method | Energy (Hartree) | Input file |
|--------|-------|--------|------------------|------------|
| | | | | |

## PySCF reference values

| Quantity | System | Basis | Notes |
|----------|--------|-------|-------|
| Overlap S | H₂ | STO-3G | off-diagonal ≈ 0.6593 at 0.74 Å |
