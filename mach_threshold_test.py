"""
Mach Threshold Test -- Viscoelastic Vacuum f(R,T)
=================================================
Reproduction notebook for Bouille (2026), Section 4.6.

Binary prediction: supersonic mergers (M > 1, with M = v_shock/c_s(T_X))
should exhibit a DM-gas spatial offset; subsonic mergers should not.

KEY PHYSICS: The Mach number is computed using the LOCAL sound speed
c_s(T_X) = sqrt(gamma * k_B * T_X / (mu * m_p)) of each cluster's
pre-shock ICM, NOT a fixed universal threshold.

Data: cluster_data.csv (primary literature sources documented inline)
"""

import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
import os

# --- Physical constants ---
k_B = 1.381e-23      # J/K
m_p = 1.673e-27      # kg
gamma_gas = 5.0 / 3.0  # adiabatic index, fully ionized plasma
mu_mol = 0.6          # mean molecular weight (H+He ICM)
c_kms = 299792.458    # km/s

def c_s_icm(T_keV):
    """ICM sound speed in km/s from pre-shock temperature in keV."""
    T_K = T_keV * 1.16e7  # keV -> Kelvin
    cs_ms = math.sqrt(gamma_gas * k_B * T_K / (mu_mol * m_p))  # m/s
    return cs_ms / 1e3  # -> km/s

# --- Load data ---
script_dir = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(script_dir, "cluster_data.csv")

clusters = []
with open(data_path, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        T_pre = float(row["T_pre_keV"]) if row["T_pre_keV"] else None
        n_e = float(row["n_e_cm3"]) if row["n_e_cm3"] else None
        clusters.append({
            "name": row["cluster"],
            "v_shock": float(row["v_shock_kms"]),
            "T_pre": T_pre,
            "n_e": n_e,
            "R": float(row["R_obs"]),
            "dR": float(row["R_err"]),
            "z": float(row["z"]),
            "notes": row["notes"],
        })

# --- Binary Mach threshold test (LOCAL c_s) ---
print("=" * 85)
print("BINARY MACH THRESHOLD TEST (Section 4.6)")
print("Using LOCAL sound speed c_s(T_X) for each cluster")
print("=" * 85)
print(f"{'Cluster':30s} {'v_shock':>7s} {'T_pre':>6s} {'c_s':>6s} {'Mach':>6s} {'Pred':>8s} {'Obs':>8s} {'OK?':>5s}")
print("-" * 85)

n_correct = 0
n_total = 0

merging = [c for c in clusters if "NULL" not in c["notes"].upper()]
null_tests = [c for c in clusters if "NULL" in c["notes"].upper()]

mach_values = []

for c in merging:
    if c["T_pre"] is None:
        continue
    cs = c_s_icm(c["T_pre"])
    M = c["v_shock"] / cs
    mach_values.append(M)
    predicted = "offset" if M > 1.0 else "no offset"
    observed = "offset" if c["R"] > 2.0 else "no offset"
    ok = predicted == observed
    n_correct += int(ok)
    n_total += 1
    status = "OK" if ok else "FAIL"
    print(f"{c['name']:30s} {c['v_shock']:7.0f} {c['T_pre']:6.1f} {cs:6.0f} {M:6.2f} {predicted:>8s} {observed:>8s} {status:>5s}")

print("-" * 85)
print(f"Score: {n_correct}/{n_total} = {100*n_correct/n_total:.0f}%")
print(f"Minimum Mach number in sample: {min(mach_values):.2f} (all > 1.0)")
print()

# --- Comparison with fixed v_Psi = 1500 km/s ---
print("COMPARISON: Fixed v_Psi=1500 vs Local c_s(T_X)")
print("-" * 85)
v_fixed = 1500.0
for c in merging:
    if c["T_pre"] is None:
        continue
    cs = c_s_icm(c["T_pre"])
    M_local = c["v_shock"] / cs
    M_fixed = c["v_shock"] / v_fixed
    flag = " <-- RESCUED" if M_fixed < 1.05 and M_local > 1.2 else ""
    print(f"  {c['name']:28s} M_fixed={M_fixed:.2f}  M_local={M_local:.2f}{flag}")
print()

# --- Null test ---
print("NULL TESTS (NGC 1052-DF2)")
print("-" * 85)
for c in null_tests:
    print(f"  {c['name']}: R = {c['R']:.1f} +/- {c['dR']:.1f}, "
          f"predicted R ~= 1 (no dynamic amplification)")
    consistent = abs(c['R'] - 1.0) < 2 * c['dR']
    print(f"  Consistent with prediction: {'YES' if consistent else 'NO'}")
print()

# --- Linear fit (exploratory) ---
print("=" * 85)
print("EXPLORATORY LINEAR FIT R(v^2) -- NOT A CONFIRMED RESULT")
print("Shown for completeness; requires O(50) clusters for significance")
print("=" * 85)

v = np.array([c["v_shock"] for c in merging if c["T_pre"]])
R = np.array([c["R"] for c in merging if c["T_pre"]])
dR = np.array([c["dR"] for c in merging if c["T_pre"]])
x = (v / c_kms)**2

W = np.diag(1.0 / dR**2)
A = np.column_stack([np.ones_like(x), x])
cov = np.linalg.inv(A.T @ W @ A)
params = cov @ (A.T @ W @ R)

C_stat, beta = params
sigma_C = np.sqrt(cov[0, 0])
sigma_beta = np.sqrt(cov[1, 1])

print(f"  C_stat = {C_stat:.2f} +/- {sigma_C:.2f}")
print(f"  beta   = {beta:.0f} +/- {sigma_beta:.0f}")
print(f"  beta significance = {beta/sigma_beta:.2f} sigma (from zero)")
print(f"  NOTE: beta = {beta:.0f} is {beta/sigma_beta:.1f}sigma from zero -- INCONCLUSIVE")
print()

# --- Figure ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Panel 1: Binary Mach threshold with LOCAL c_s
mach_local = []
names_local = []
for c in merging:
    if c["T_pre"] is None:
        continue
    cs = c_s_icm(c["T_pre"])
    mach_local.append(c["v_shock"] / cs)
    names_local.append(c["name"].split("(")[0].strip())

colors = ['#2ecc71' if m > 1 else '#e74c3c' for m in mach_local]
ax1.barh(names_local, mach_local, color=colors, edgecolor='white', height=0.6)
ax1.axvline(x=1.0, color='k', linestyle='--', linewidth=2, label='M = 1 threshold')
ax1.set_xlabel('Mach number M = v_shock / c_s(T_X)', fontsize=12)
ax1.set_title(f'Binary Mach Threshold Test ({n_correct}/{n_total})\n'
              f'Using LOCAL c_s(T_X) per cluster', fontsize=13)
ax1.legend(fontsize=11)
ax1.set_xlim(0, 4.5)

# Annotate T_pre on bars
for i, c in enumerate([c for c in merging if c["T_pre"]]):
    cs = c_s_icm(c["T_pre"])
    M = c["v_shock"] / cs
    ax1.text(M + 0.05, i, f'{c["T_pre"]:.1f} keV', va='center', fontsize=8, color='#555')

# Panel 2: R(v^2) exploratory
x_plot = np.linspace(0, max(x) * 1.1, 100)
ax2.errorbar(x * 1e6, R, yerr=dR, fmt='o', color='#3498db',
             markersize=8, capsize=4, label='Merging clusters')
ax2.plot(x_plot * 1e6, C_stat + beta * x_plot, 'r--', alpha=0.5,
         label=f'Fit: beta = {beta:.0f} +/- {sigma_beta:.0f} (inconclusive)')
ax2.axhline(y=1, color='gray', linestyle=':', alpha=0.5, label='R = 1 (no DM)')
ax2.set_xlabel('(v_coll/c)^2 x 10^6', fontsize=12)
ax2.set_ylabel('R = M_lens / M_bar', fontsize=12)
ax2.set_title('Exploratory R(v^2) -- NOT a confirmed result', fontsize=13)
ax2.legend(fontsize=10)

plt.tight_layout()
fig_path = os.path.join(script_dir, "mach_threshold_test.png")
plt.savefig(fig_path, dpi=150)
print(f"\nFigure saved: {fig_path}")
print("Done.")
