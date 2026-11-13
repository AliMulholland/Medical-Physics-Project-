# Saves the depth-dose histogram of the ten main runs to seed_depth_dose.npz

import numpy as np
import uproot

histograms = {}
for beam in ["c2750", "g10"]:  # carbon 2750 MeV, photons 10 MeV
    for seed in [1, 2, 3, 4, 5]:
        run = f"{beam}_s{seed}"
        # energy (MeV) in each 1 mm layer
        histograms[run] = uproot.open(f"output.{run}.1.root")["Edep_vs_z"].values()

np.savez("seed_depth_dose.npz", **histograms)
print("saved seed_depth_dose.npz")
