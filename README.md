# findformula

Computer-assisted search for new, verified formulas.

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
