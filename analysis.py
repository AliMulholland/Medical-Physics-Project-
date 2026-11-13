# Numbers for Sections 4 and 5
# needs seed_depth_dose.npz, and output.c2750_s1.1.root / output.g10_s1.1.root for the lateral spread

import os
import numpy as np
import uproot

N_CARBON = 5000   # ions per run
N_PHOTON = 10000  # photons per run
E_CARBON = 2750.0  # MeV
E_PHOTON = 10.0
SEEDS = [1, 2, 3, 4, 5]

depth = np.arange(500) + 0.5  # layer centres (mm)


def mean_and_error(values):
    # mean and standard error over the runs
    values = np.array(values)
    return values.mean(), values.std(ddof=1) / np.sqrt(len(values))


def show(label, values, unit):
    mean, error = mean_and_error(values)
    print(f"  {label:<45} {mean:.4g} +/- {error:.2g} {unit}")


def depth_where_curve_falls_to(curve, level, start, step):
    # step from 'start' (+1 deeper, -1 shallower) until below 'level', then interpolate
    i = start
    while curve[i + step] > level:
        i += step
    fraction = (curve[i] - level) / (curve[i] - curve[i + step])
    return depth[i] + step * fraction


runs = np.load("seed_depth_dose.npz")


# carbon 2750 MeV

peak_depth, peak_height, plateau, ratio = [], [], [], []
width, beam_range, tail, beyond_peak, kept_carbon = [], [], [], [], []
tail_20cm, tail_33cm = [], []

for seed in SEEDS:
    curve = runs[f"c2750_s{seed}"] / N_CARBON  # MeV per mm per ion
    peak = curve.argmax()  # peak layer

    peak_depth.append(depth[peak])
    peak_height.append(curve[peak])
    plateau.append(curve[1:10].mean())  # plateau, 1-10 mm
    ratio.append(curve[peak] / curve[1:10].mean())

    # FWHM
    half = curve[peak] / 2
    back = depth_where_curve_falls_to(curve, half, peak, +1)
    front = depth_where_curve_falls_to(curve, half, peak, -1)
    width.append(back - front)

    # range = distal 80% point
    beam_range.append(depth_where_curve_falls_to(curve, 0.8 * curve[peak], peak, +1))

    tail.append(100 * curve[115:125].mean() / curve[peak])
    tail_20cm.append(100 * curve[195:205].mean() / curve[peak])
    tail_33cm.append(100 * curve[325:335].mean() / curve[peak])
    beyond_peak.append(100 * curve[115:].sum() / curve.sum())
    kept_carbon.append(100 * curve.sum() / E_CARBON)  # % of beam energy

print("CARBON 2750 MeV (mean of 5 runs +/- standard error)")
show("depth of Bragg peak", peak_depth, "mm")
show("energy in peak layer", peak_height, "MeV per mm per ion")
show("entrance plateau (also LET in keV/um)", plateau, "MeV per mm per ion")
show("peak / entrance", ratio, "")
show("peak width at half height", width, "mm")
show("range (80% of peak)", beam_range, "mm")
show("tail at 115-125 mm", tail, "% of peak")
show("tail at 20 cm (195-205 mm)", tail_20cm, "% of peak")
show("tail at 33 cm (325-335 mm)", tail_33cm, "% of peak")
show("energy deposited beyond 115 mm", beyond_peak, "% of deposited energy")
show("beam energy kept in phantom", kept_carbon, "%")


# photons 10 MeV

maximum, at_tumour, surface, exit_dose, mu_eff, kept_photon = [], [], [], [], [], []
tumour_layer = 110  # carbon peak layer (110.5 mm)
fit_region = (depth > 100) & (depth < 450)  # 10-45 cm

for seed in SEEDS:
    curve = runs[f"g10_s{seed}"] / N_PHOTON

    top = curve[40:90].mean()  # maximum, 4-9 cm
    maximum.append(1000 * top)  # keV
    # 9 mm average around the tumour depth
    at_tumour.append(1000 * curve[tumour_layer - 4:tumour_layer + 5].mean())
    surface.append(100 * curve[0] / top)
    exit_dose.append(100 * curve[490:].mean() / top)

    # straight-line fit to log(energy)
    slope, intercept = np.polyfit(depth[fit_region], np.log(curve[fit_region]), 1)
    mu_eff.append(-slope * 10)  # per cm

    kept_photon.append(100 * curve.sum() / E_PHOTON)

print()
print("PHOTONS 10 MeV (mean of 5 runs +/- standard error)")
show("maximum (average over 4-9 cm)", maximum, "keV per mm per photon")
show("at tumour depth (9 mm average)", at_tumour, "keV per mm per photon")
show("tumour depth / maximum", 100 * np.array(at_tumour) / np.array(maximum), "%")
show("surface (first mm) / maximum", surface, "%")
show("exit (last cm) / maximum", exit_dose, "%")
show("mu_eff from fit, 10-45 cm", mu_eff, "per cm")
show("beam energy kept in phantom", kept_photon, "%")
show("photons needed to match one carbon ion",
     1000 * np.array(peak_height) / np.array(at_tumour), "")


# energy more than 1 cm from the axis (seed 1 only)

def sideways_spread(filename, name, slice_depths_cm):
    table = uproot.open(filename)["EdMedPh"].arrays(["Edep", "X", "Y", "Z"], library="np")
    energy = table["Edep"]
    r = np.sqrt(table["X"] ** 2 + table["Y"] ** 2)  # mm
    far = r > 10

    print(f"  {name}: {100 * energy[far].sum() / energy.sum():.1f} % of all energy is more than 1 cm from the axis")
    for centre in slice_depths_cm:
        in_slice = np.abs(table["Z"] - 10 * centre) < 5  # 1 cm slice
        share = 100 * energy[in_slice & far].sum() / energy[in_slice].sum()
        print(f"      slice at {centre:4.1f} cm: {share:5.1f} %")


print()
print("SIDEWAYS SPREAD (one run of each beam)")
if os.path.exists("output.c2750_s1.1.root") and os.path.exists("output.g10_s1.1.root"):
    sideways_spread("output.c2750_s1.1.root", "carbon", [2, 5, 8, 10.5])
    sideways_spread("output.g10_s1.1.root", "photons", [])
else:
    print("  full output files not found, skipped")
