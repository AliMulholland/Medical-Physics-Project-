# Carbon-ion vs photon radiotherapy in Geant4

Analysis code for my Medical Physics (PHYS11064) project at the University of Edinburgh, 2026/27.
It compares 10 MeV photons with 2750 MeV (229 MeV/u) carbon-12 ions in a water phantom
and builds a spread-out Bragg peak (SOBP) from eight carbon energies.

These scripts only analyse the simulation output. They do not run Geant4. The simulations were
run separately with the course code, [EdMedPhysics](https://github.com/bwynneHEP/EdMedPhysics),
and the scripts read the `.root` files those runs produced.

## Runs analysed

| Beam | Energy | Particles per run | Runs |
|---|---|---|---|
| Photons | 10 MeV | 10 000 | 5 seeds |
| Carbon-12 | 2750 MeV | 5000 | 5 seeds |
| Carbon-12 (SOBP) | 2300-2750 MeV, 65 MeV steps | 2000 | 8 energies |

## Scripts and their inputs

| Script | Input | Output |
|---|---|---|
| `collect_depth_dose.py` | `output.c2750_s1.1.root` ... `output.c2750_s5.1.root`, `output.g10_s1.1.root` ... `output.g10_s5.1.root` | `seed_depth_dose.npz` (depth-dose histogram of each of the ten runs) |
| `analysis.py` | `seed_depth_dose.npz`; `output.c2750_s1.1.root` and `output.g10_s1.1.root` for the lateral spread | prints the values in Sections 4 and 5 and Table 2 |
| `sobp.py` | `sobp_2300_hist.root` ... `sobp_2750_hist.root` (one per SOBP energy), `g10_hist.root` (five photon runs added together) | prints the values in Section 5 and Table 3; draws Figure 4 (`plots/sobp.pdf`) |
| `make_figures.py` | `seed_depth_dose.npz`; `output.c2750_s1.1.root` and `output.g10_s1.1.root` for Figure 3 | draws Figures 1-3 (`plots/`) |

From the `output.*.root` files, the scripts use the depth-dose histogram `Edep_vs_z`
(energy deposited in each 1 mm layer) and, for the lateral spread and Figure 3, the step table
`EdMedPh` (every energy deposit with its position).

## How to use

1. Put the `.root` files listed above in the same folder as the scripts.
2. Run `python collect_depth_dose.py` once to make `seed_depth_dose.npz`.
3. Run `python analysis.py`, `python sobp.py` and `python make_figures.py`.

Needs Python 3 with `numpy`, `matplotlib` and `uproot` libraries.
