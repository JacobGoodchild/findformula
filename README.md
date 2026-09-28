# findformula

Computer-assisted search for new, verified formulas.

## Papers

Two preprints written up from these results (published on Zenodo, not peer reviewed):

1. **Exact localization of Grover quantum walks on crystal lattices: transfer currents, a universal
   bound, and uniform-spanning-tree degree laws.** J. Goodchild (2026).
   Zenodo: https://zenodo.org/records/23010606 · DOI: [10.5281/zenodo.23010606](https://doi.org/10.5281/zenodo.23010606)
   · Source: `paper1_quantum_walks_and_graph_laws.tex`
2. **Exact height-one probabilities of three-dimensional Abelian sandpiles, spanning-tree entropies of
   the snub-square and Cairo lattices, and elementary resistor-network reductions.** J. Goodchild (2026).
   Zenodo: https://zenodo.org/records/23010883 · DOI: [10.5281/zenodo.23010883](https://doi.org/10.5281/zenodo.23010883)
   · Source: `paper2_sandpiles_and_lattice_entropy.tex`

**AI-use disclosure:** the research in this repository, including choosing the problems, discovering
the formulas, running all computations and writing both papers, was carried out with extensive use of an
AI system (Claude Code, by Anthropic) under the author's direction. `verification/` contains independent
numerical re-checks of the main formulas.

## Contents

- `formulas.txt` — verified results only (derivation, verification, novelty status). An index of all
  entries is at the top. Main themes:
  - Quantum walks: exact long-time trapping of Grover walks on lattices (Formulas 3-23, 29, 32) and a
    universal lower bound (Formula 24).
  - Resistor networks and spanning trees: exact resistances and spanning-tree entropies for decorated
    square, honeycomb and diamond networks (Formulas 25-28, 30, 31, 33).
  - Abelian sandpile: exact height-1 probabilities on 2D and 3D lattices (Formulas 34-36).
- `lattices/`, `lattices3d/` — exact Brillouin-zone engines (residues + quadrature), builders, checks.
- `quantum/`, `quantum2d/`, `aniso/`, `moving/`, `lazyff/`, `diamond/` — quantum-walk computations.
- `sandpile/` — Majumdar-Dhar computations, exact finite-box checks, simulations.
- `physics/` — fluid-dynamics work (added-mass images, Stokes wall-correction constants).
- `experiments/` — combinatorics scans (spanning trees, tilings, linear extensions); nothing new found there.
- `NOTES.md` — log of what was tried, including negative results and open leads.
- `verification/` — independent recomputation of key resistances, sandpile values and identities from their defining integrals.
