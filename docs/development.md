# Development Practices

## Scope of this document

This project is an educational reimplementation of Hartree-Fock and MP2 theory. Its purpose
is pedagogical: the numerical and physical core is written by hand in order to understand
how an electronic structure program is built, rather than to compete with production
quantum chemistry codes in performance or generality.

This document records the practices followed during development.

## Hand-written core

The scientific core of this repository is implemented manually, without automated code
generation. This includes:

- Gaussian integral evaluation (Boys function, McMurchie-Davidson recursions, and the
  overlap, kinetic, nuclear attraction and electron repulsion integrals)
- Fock matrix construction and the SCF procedure
- DIIS convergence acceleration
- AO to MO integral transformation and the MP2 correlation energy

Automated assistance was restricted to test scaffolding, continuous integration
configuration, plotting code, batch reference-data scripts, and documentation. This
boundary was enforced through a local configuration file that is not versioned.

## Validation strategy

Correctness is established at two independent levels rather than by comparison against a
single reference.

**Integrals** are verified element by element against an established quantum chemistry
library, which exposes the same quantities as arrays. This localizes errors to individual
matrix elements instead of surfacing them only as a wrong total energy.

**Total energies** are verified against ORCA. Reference calculations are stored in
`reference/` together with their input files, geometries and basis specifications, so that
every number in the validation table can be regenerated independently.

Reference calculations are run without resolution-of-identity approximations, with
Cartesian rather than spherical basis functions, and with convergence thresholds
substantially tighter than the tolerance used for comparison.

## Testing

Every implementation milestone has an associated test that must pass before the next
milestone is started. Numerical comparisons use explicit tolerances. Permutational
symmetries of the two-electron integrals are verified as properties rather than through
hand-picked cases.

## Conventions

- Atomic units are used internally; unit conversion happens only at input boundaries
- Core numerical modules live in `src/`, with analysis confined to `notebooks/`
- Reference data is immutable once generated and versioned alongside its inputs
