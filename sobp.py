# SOBP from the eight carbon energies, compared with 10 MeV photons (Section 5)
# needs sobp_<E>_hist.root and g10_hist.root, saves plots/sobp.pdf

import os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ENERGIES = [2300, 2365, 2430, 2495, 2560, 2625, 2690, 2750]
N_IONS = 2000  # ions per energy
N_PHOTONS = 50000  # g10_hist.root = 5 runs summed

depth = np.arange(500) + 0.5  # layer centres (mm)

# per particle fired
carbon = []
for energy in ENERGIES:
    histogram = uproot.open(f"sobp_{energy}_hist.root")["Edep_vs_z"].values()
    carbon.append(histogram / N_IONS)
carbon = np.array(carbon)
photon = uproot.open("g10_hist.root")["Edep_vs_z"].values() / N_PHOTONS


def build_sobp(curves):
    # weights so the summed dose is equal at every peak (deepest beam = 1),
    # repeated because shallow beams add a little dose at the deeper peaks
    peaks = curves.argmax(axis=1)
    weights = np.zeros(len(curves))
    weights[-1] = 1.0
    target_dose = curves[-1, peaks[-1]]

    for repeat in range(50):
        for i in reversed(range(len(curves))):
            total_here = np.sum(weights * curves[:, peaks[i]])
            weights[i] += (target_dose - total_here) / curves[i, peaks[i]]
            weights[i] = max(weights[i], 0.0)

    total = np.sum(weights[:, None] * curves, axis=0)
    first, last = peaks[0], peaks[-1]  # target = first to last peak
    target_mean = total[first:last + 1].mean()
    return weights, total / target_mean, first, last


def report(name, curves):
    weights, sobp, first, last = build_sobp(curves)
    in_target = sobp[first:last + 1]
    ripple = 100 * (in_target.max() - in_target.min()) / (in_target.max() + in_target.min())
    photon_scaled = photon / photon[first:last + 1].mean()

    print(f"{name}: {len(curves)} beams, target {first}-{last + 1} mm, ripple +/- {ripple:.1f} %")
    print(f"    deepest beam: {100 * weights[-1] / weights.sum():.0f} % of the ions")
    print("                             carbon   photons")
    print(f"    entrance 0-10 mm         {sobp[0:10].mean():.3f}    {photon_scaled[0:10].mean():.2f}")
    print(f"    at 5 cm                  {sobp[45:55].mean():.3f}    {photon_scaled[45:55].mean():.2f}")
    print(f"    behind 120-200 mm        {sobp[120:200].mean():.3f}    {photon_scaled[120:200].mean():.2f}")
    print()


# single 2750 MeV peak

single = carbon[-1]
peak = single.argmax()
single = single / single[peak]
print("single peak (2750 MeV), dose at the peak = 1")
print(f"    entrance 0-10 mm         {single[0:10].mean():.3f}")
print(f"    behind 120-200 mm        {single[120:200].mean():.3f}")
print()

# 3 cm SOBP

report("3 cm target", carbon)

# Figure 4

weights, sobp, first, last = build_sobp(carbon)
scale = np.sum(weights[:, None] * carbon, axis=0)[first:last + 1].mean()
photon_scaled = photon / photon[first:last + 1].mean()

print("weights of the eight beams (3 cm target):")
for energy, weight in zip(ENERGIES, weights):
    print(f"    {energy} MeV   weight {weight:.3f}   {100 * weight / weights.sum():.1f} % of the ions")

os.makedirs("plots", exist_ok=True)
colours = plt.get_cmap("viridis")
plt.figure(figsize=(11, 5.2))
plt.axvspan(depth[first] - 0.5, depth[last] + 0.5, color="0.85", label="target")
for i, energy in enumerate(ENERGIES):
    plt.step(depth, weights[i] * carbon[i] / scale, where="mid", lw=0.9, color=colours(i / 7),
             label=f"{energy} MeV (weight {weights[i]:.2f})")
plt.step(depth, photon_scaled, where="mid", lw=1.6, ls="--", color="0.35", label="10 MeV photons")
plt.step(depth, sobp, where="mid", lw=2.4, color="black", label="carbon SOBP (sum)")
plt.xlim(0, 200)
plt.ylim(0, 1.4)
plt.xlabel("depth in water (mm)")
plt.ylabel("dose relative to mean dose in target")
plt.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
plt.tight_layout()
plt.savefig("plots/sobp.png", dpi=200)
plt.savefig("plots/sobp.pdf")

# other comparisons with photons

crossing = np.argmax(photon_scaled > sobp)
total_carbon = np.sum(weights[:, None] * carbon, axis=0)
print()
print("3 cm target, carbon against photons:")
print(f"    photons give less dose than carbon only in the first {crossing} mm")
print(f"    average dose in front of the target (0-{first} mm): carbon {sobp[:first].mean():.2f}, photons {photon_scaled[:first].mean():.2f}")
print(f"    share of all deposited energy that is inside the target: "
      f"carbon {100 * total_carbon[first:last + 1].sum() / total_carbon.sum():.0f} %, "
      f"photons {100 * photon[first:last + 1].sum() / photon.sum():.0f} %")
print(f"    highest and lowest dose inside the target: {sobp[first:last + 1].max():.2f} and {sobp[first:last + 1].min():.2f}")
