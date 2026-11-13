# Figures 1-3 (Section 4)
# needs seed_depth_dose.npz, and output.c2750_s1.1.root / output.g10_s1.1.root for Figure 3

import os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

N_CARBON = 5000  # ions per run
N_PHOTON = 10000  # photons per run
BLUE = "#0072B2"
ORANGE = "#D55E00"

depth = np.arange(500) + 0.5  # layer centres (mm)
depth_cm = depth / 10
os.makedirs("plots", exist_ok=True)

# mean of the five runs, per particle fired
runs = np.load("seed_depth_dose.npz")
carbon = np.mean([runs[f"c2750_s{seed}"] for seed in [1, 2, 3, 4, 5]], axis=0) / N_CARBON
photon = np.mean([runs[f"g10_s{seed}"] for seed in [1, 2, 3, 4, 5]], axis=0) / N_PHOTON


# Figure 1

# exponential fit, 10-45 cm
fit_region = (depth > 100) & (depth < 450)
slope, intercept = np.polyfit(depth[fit_region], np.log(photon[fit_region]), 1)
mu_eff = -slope * 10
fit_curve = np.exp(intercept + slope * depth[fit_region])

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4))

left.step(depth_cm, 1000 * photon, where="mid", color=BLUE, lw=0.8, label="simulation (50 000 photons)")
left.plot(depth_cm[fit_region], 1000 * fit_curve, "k--", lw=1.5,
          label=f"exponential fit, $\\mu_\\mathrm{{eff}}$ = {mu_eff:.4f} cm$^{{-1}}$")
left.set_title("(a) 10 MeV photons")
left.set_xlabel("depth $z$ (cm)")
left.set_ylabel("energy deposited (keV per mm per photon)")
left.set_xlim(0, 50)
left.set_ylim(0, 17)
left.legend(frameon=False)

right.step(depth_cm, carbon, where="mid", color=ORANGE, lw=1.2)
right.set_title("(b) 2750 MeV $^{12}$C (25 000 ions)")
right.set_xlabel("depth $z$ (cm)")
right.set_ylabel("energy deposited (MeV per mm per ion)")
right.set_xlim(0, 30)
right.set_ylim(0, 88)
right.annotate("Bragg peak\n110.5 mm", xy=(10.9, 78), xytext=(4.5, 70), arrowprops={"arrowstyle": "->"}, fontsize=9)
right.annotate("plateau", xy=(4, carbon[40]), xytext=(3, 30), arrowprops={"arrowstyle": "->"}, fontsize=9)

# log inset for the fragment tail
inset = right.inset_axes([0.55, 0.33, 0.42, 0.45])
inset.step(depth_cm, carbon, where="mid", color=ORANGE, lw=1)
inset.set_yscale("log")
inset.set_xlim(0, 50)
inset.set_ylim(0.01, 200)
inset.set_title("log scale, full depth: fragment tail", fontsize=8)
inset.tick_params(labelsize=7)

plt.tight_layout()
plt.savefig("plots/fig_depth_dose.pdf")
plt.savefig("plots/fig_depth_dose.png", dpi=200)
plt.close()


# Figure 2

peak = carbon.argmax()
carbon_at_tumour = carbon[peak]
photon_at_tumour = photon[peak - 4:peak + 5].mean()  # 9 mm average

plt.figure(figsize=(7, 4))
plt.step(depth_cm, photon / photon_at_tumour, where="mid", color=BLUE, lw=0.9, label="10 MeV photons")
plt.step(depth_cm, carbon / carbon_at_tumour, where="mid", color=ORANGE, lw=1.3, label="2750 MeV $^{12}$C")
plt.axvline(depth_cm[peak], color="grey", ls=":", lw=1)
plt.text(depth_cm[peak] + 0.3, 1.22, "tumour depth\n110.5 mm", fontsize=8, color="grey")
plt.xlim(0, 50)
plt.ylim(0, 1.4)
plt.xlabel("depth $z$ (cm)")
plt.ylabel("dose relative to dose at tumour depth")
plt.legend(frameon=False, loc="upper right")
plt.tight_layout()
plt.savefig("plots/fig_overlay.pdf")
plt.savefig("plots/fig_overlay.png", dpi=200)
plt.close()


# Figure 3 (step tables, seed 1 only)

if os.path.exists("output.g10_s1.1.root") and os.path.exists("output.c2750_s1.1.root"):
    z_edges = np.linspace(0, 30, 301)  # 1 mm bins, in cm
    y_edges = np.linspace(-6, 6, 121)  # 1 mm bins, in cm

    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
    panels = [(left, "output.g10_s1.1.root", "(a) 10 MeV photons"),
              (right, "output.c2750_s1.1.root", "(b) 2750 MeV $^{12}$C")]
    for axis, filename, title in panels:
        table = uproot.open(filename)["EdMedPh"].arrays(["Edep", "Y", "Z"], library="np")
        # positions in the file are in mm
        energy_map, _, _ = np.histogram2d(table["Z"] / 10, table["Y"] / 10,
                                          bins=[z_edges, y_edges], weights=table["Edep"])
        energy_map = energy_map / energy_map.max()
        image = axis.pcolormesh(z_edges, y_edges, energy_map.T, norm=LogNorm(1e-5, 1),
                                cmap="viridis", rasterized=True)
        axis.set_title(title)
        axis.set_xlabel("depth $z$ (cm)")
    left.set_ylabel("lateral position $y$ (cm)")
    fig.colorbar(image, ax=[left, right], label="energy per bin / maximum", pad=0.02)
    plt.savefig("plots/fig_maps.pdf", bbox_inches="tight")
    plt.savefig("plots/fig_maps.png", dpi=200, bbox_inches="tight")
    plt.close()
else:
    print("full output files not found, energy maps skipped")
