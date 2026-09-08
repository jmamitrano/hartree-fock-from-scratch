# CLAUDE.md

> Local file. Not versioned — listed in .gitignore.
> The public counterpart of these rules is docs/development.md.

## What this repository is

A from-scratch implementation of Hartree-Fock (RHF and UHF) and MP2 with contracted
Gaussian basis sets, in Python with numpy. This is a **learning project**: the goal is not
the fastest possible code, but that I understand and write the physics and numerics by
hand.

## Primary rule

**Do not edit or write code in the scientific core.** This covers everything related to:

- Gaussian integrals (Boys function, McMurchie-Davidson, S, T, V, ERI)
- Fock matrix construction and the SCF cycle
- DIIS
- AO to MO transformation and MP2

I write those files. You may **read, explain, critique and help me debug them**, but not
modify them.

If you find a bug there, tell me where it is and why, and leave the fix to me.

## What you may write and edit freely

- `tests/` — pytest suite, fixtures, parametrization, comparisons against PySCF
- `.github/` — CI workflows
- `pyproject.toml`, `ruff` and `pre-commit` configuration
- `scripts/` — batch generation of ORCA inputs, output parsing, reference value extraction
- Plotting code for the final analysis
- `README.md` and files under `docs/`

## How to help me debug

1. Run the tests and show me what fails
2. Ask me what I expected to happen
3. Give me a hint about where to look, not the solution
4. Escalate to a fuller explanation only if I ask

## Conventions

- Python 3.11+, numpy, scipy. Environment managed with `uv`.
- Atomic units internally. Anything entering in Ångström is converted at the boundary.
- Tests use `numpy.testing.assert_allclose` with explicit tolerances.
- ORCA reference values in `reference/` are versioned with their inputs and geometries.
  Do not regenerate or modify them.
- All repository content is in English: code, docstrings, commit messages, documentation.
