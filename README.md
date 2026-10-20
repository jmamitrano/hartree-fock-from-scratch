# Hartree-Fock from scratch

A Python implementation of restricted and unrestricted Hartree-Fock and second-order
Møller-Plesset perturbation theory, built from the Gaussian integrals upward. No quantum
chemistry library does any of the physics: the Boys function, the integral recursions, the
SCF cycle and the MP2 energy are all written by hand.

## Scope

**In:** contracted Cartesian Gaussian basis sets with `s` and `p` functions, overlap,
kinetic, nuclear attraction and two-electron repulsion integrals, RHF, UHF, DIIS
acceleration, MP2 correlation energy, and a dissociation curve as the closing analysis.

**Out:** `d` functions and beyond, density functional theory, geometry optimization,
integral screening, symmetry, and any performance work. A 7-function water calculation is
the target scale.

## Status

Milestone 1 of 11 complete — the Boys function, with 230 tests passing.

See [`PROGRESS.md`](PROGRESS.md) for the milestone table, the design decisions taken along
the way and their rationale.

## Running it

The environment is managed with [uv](https://docs.astral.sh/uv/):

```bash
uv sync
uv run pytest
```

## How correctness is established

Each milestone has a test that must pass before the next one starts, and the tests are
built on two principles.

**No single oracle.** The Boys function, for instance, is checked against 
both numerical quadrature of its defining integral and the confluent hypergeometric 
function, and against a test that compares those two oracles to each other, 
because one of them shares a derivation with the implementation
it is supposed to validate and so cannot certify it alone.

**Properties over values.** Every quantity has invariants that hold regardless of the
input: the Boys function is positive and strictly decreasing in its order, the overlap
matrix is symmetric with unit diagonal, the two-electron integrals obey eight permutational
symmetries. These catch structural errors that value-by-value comparison can miss, and need
no reference data at all.

Integrals will be verified element by element against PySCF, which localizes an error to a
specific matrix element rather than surfacing it as a wrong total energy. Total energies
will be verified against ORCA, with the inputs and geometries versioned alongside the
numbers so every entry in the validation table can be regenerated.

Numerical thresholds are chosen from measured error data rather than by eye;
[`studies/`](studies/) holds the sweeps that produced them.

See [`docs/development.md`](docs/development.md) for the full practices.

## Layout

```
src/        numerical core — written by hand
tests/      pytest suite, including the independent reference implementations
studies/    parameter sweeps backing the numerical choices
docs/       development practices
```

## References

- Szabo & Ostlund, *Modern Quantum Chemistry* — the derivations this follows
- Helgaker, Jørgensen & Olsen, *Molecular Electronic-Structure Theory* — integral evaluation
- McMurchie & Davidson, *J. Comput. Phys.* **26** (1978) 218 — the integral recursions used here
