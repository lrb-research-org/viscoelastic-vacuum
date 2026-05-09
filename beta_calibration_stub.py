#!/usr/bin/env python3
"""
beta-coefficient calibration for the lensing amplification R(v) scaling.

STATUS: STUB -- methodology documented, awaiting Euclid/Rubin LSST ~50-cluster sample.
REVIEWER: R4 S2.2 (quantitative beta calibration).

The lensing-to-baryon mass ratio R = M_lens / M_baryon scales with the
collision velocity via the Mach-cone scalar retardation mechanism:

    R(v) = 1 + beta * (M^2 - 1)    for M = v/c_s > 1
    R(v) = 1                        for M = v/c_s <= 1

where beta encodes the efficiency of scalar momentum transfer from the
Lienard-Wiechert retarded potential to the effective lensing convergence.
The (M^2 - 1) form enforces R -> 1 at the sonic threshold.

WHAT THIS SCRIPT DOES (current state):
  - Loads the 9-cluster sample + 1 null test from cluster_data.csv
  - Propagates per-cluster R_err (weak-lensing uncertainty)
  - Fits beta from the (R-1) vs (M^2-1) relation using weighted least-squares
  - Reports the best-fit beta, its uncertainty, and chi^2/dof
  - Documents the systematic limitations of the current sample

WHAT IS NEEDED FOR PUBLICATION-GRADE CALIBRATION:
  1. Euclid / Rubin LSST weak-lensing mass estimates for ~50 merging clusters
  2. Chandra/eROSITA X-ray temperatures and shock velocities
  3. Per-cluster covariance (shear noise + photo-z + ICM T_X systematics)
  4. Selection function correction (X-ray vs optically selected mergers)

Reference: Manuscript S4, Eq. 35-36; Table 6.
"""

import numpy as np
import os
import sys


def load_cluster_data():
    """Load the 9-cluster sample + 1 null test from cluster_data.csv."""
    csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "cluster_data.csv")

    # CSV columns: cluster(0), v_shock_kms(1), v_shock_ref(2), T_pre_keV(3),
    #   T_pre_ref(4), n_e_cm3(5), n_e_ref(6), M_gas_1e14Msun(7),
    #   M_gas_ref(8), R_obs(9), R_err(10), R_method(11), ...
    clusters = []
    with open(csv_path, 'r') as f:
        header = f.readline().strip().split(',')
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 11:
                continue
            name = parts[0].strip()
            try:
                v_shock = float(parts[1])   # km/s
                R_ratio = float(parts[9])   # R_obs = M_lens / M_bar
                R_err = float(parts[10])    # 1-sigma uncertainty on R
            except (ValueError, IndexError):
                continue
            # T_pre_keV may be empty for null test (NGC 1052-DF2)
            try:
                T_X = float(parts[3])       # keV
            except (ValueError, IndexError):
                T_X = None
            clusters.append({
                'name': name,
                'v_shock': v_shock,
                'T_X': T_X,
                'R': R_ratio,
                'R_err': R_err
            })
    return clusters


def sound_speed_from_TX(T_keV):
    """ICM sound speed from X-ray pre-shock temperature.

    c_s = sqrt(gamma * k_B * T / (mu * m_p))
    with gamma = 5/3, mu = 0.6 (fully ionized ICM).
    """
    gamma = 5.0 / 3.0
    mu = 0.6
    m_p = 1.6726e-27  # kg
    k_B = 1.3807e-23  # J/K
    T_K = T_keV * 1.1605e7  # keV to Kelvin
    c_s = np.sqrt(gamma * k_B * T_K / (mu * m_p))  # m/s
    return c_s / 1e3  # km/s


def fit_beta(clusters):
    """Fit beta from (R-1) = beta * (M^2 - 1) using weighted least-squares.

    Only supersonic mergers (M > 1) with valid T_X are included.
    Weights = 1/R_err^2 (inverse variance).
    """
    M2m1_vals = []
    Rm1_vals = []
    w_vals = []
    names = []

    for c in clusters:
        if c['T_X'] is None or c['R'] <= 0 or c['R_err'] <= 0:
            continue
        c_s = sound_speed_from_TX(c['T_X'])
        M = c['v_shock'] / c_s
        if M > 1.0:
            M2m1_vals.append(M**2 - 1)  # threshold-corrected
            Rm1_vals.append(c['R'] - 1)
            w_vals.append(1.0 / c['R_err']**2)
            names.append(c['name'])

    if len(M2m1_vals) < 2:
        print("ERROR: insufficient supersonic mergers for fit")
        return None

    x = np.array(M2m1_vals)
    y = np.array(Rm1_vals)
    w = np.array(w_vals)

    # Weighted OLS through origin: y = beta * x
    # beta = sum(w*x*y) / sum(w*x^2)
    beta = np.sum(w * x * y) / np.sum(w * x**2)
    y_pred = beta * x
    residuals = y - y_pred

    # Weighted chi^2
    chi2 = np.sum(w * residuals**2)
    dof = len(x) - 1
    chi2_red = chi2 / dof if dof > 0 else float('inf')

    # Standard error on beta
    sigma_beta = np.sqrt(1.0 / np.sum(w * x**2))

    # Weighted R^2
    y_mean = np.sum(w * y) / np.sum(w)
    ss_tot = np.sum(w * (y - y_mean)**2)
    ss_res = np.sum(w * residuals**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return {
        'beta': beta,
        'sigma_beta': sigma_beta,
        'n_clusters': len(x),
        'M2m1': x,
        'Rm1': y,
        'names': names,
        'residuals': residuals,
        'chi2_red': chi2_red,
        'R2': R2,
        'dof': dof
    }


def main():
    sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None

    print("=" * 70)
    print("beta-coefficient calibration -- STUB (R4 S2.2)")
    print("=" * 70)
    print()
    print("  Model: R - 1 = beta * (M^2 - 1)   [M = v_shock / c_s]")
    print("  R = M_lens / M_bar from weak-lensing mass estimates")
    print("  Threshold: R = 1 at M = 1 (sonic barrier, by construction)")
    print()

    clusters = load_cluster_data()
    print(f"  Loaded {len(clusters)} entries from cluster_data.csv")
    print()

    # Display table with Mach numbers
    print(f"  {'Cluster':<30s} {'v_shock':>7s} {'c_s':>7s} {'M':>5s} "
          f"{'R_obs':>5s} {'R_err':>5s} {'M^2-1':>6s}")
    print("  " + "-" * 68)
    for c in clusters:
        if c['T_X'] is not None:
            c_s = sound_speed_from_TX(c['T_X'])
            M = c['v_shock'] / c_s
            M2m1 = max(M**2 - 1, 0)
            print(f"  {c['name']:<30s} {c['v_shock']:>7.0f} {c_s:>7.0f} "
                  f"{M:>5.2f} {c['R']:>5.1f} {c['R_err']:>5.1f} {M2m1:>6.2f}")
        else:
            print(f"  {c['name']:<30s} {c['v_shock']:>7.0f}     {'--':>3s} "
                  f"{'--':>5s} {c['R']:>5.1f} {c['R_err']:>5.1f}    {'--':>3s}"
                  f"  [NULL TEST]")

    print()
    result = fit_beta(clusters)
    if result:
        print("  WEIGHTED LEAST-SQUARES FIT: (R-1) = beta * (M^2 - 1)")
        print(f"  --------------------------------------------------------")
        print(f"  beta      = {result['beta']:.4f} +/- {result['sigma_beta']:.4f}")
        print(f"  chi^2/dof = {result['chi2_red']:.2f}  ({result['dof']} dof)")
        print(f"  N_clusters= {result['n_clusters']}")
        if result['chi2_red'] > 2:
            print(f"  NOTE: chi^2/dof > 2 indicates excess scatter beyond R_err.")
            print(f"        Expected with heterogeneous WL pipelines (see caveats).")
        print()

        # Per-cluster residuals
        print(f"  {'Cluster':<30s} {'M^2-1':>6s} {'R-1':>6s} {'Pred':>6s} "
              f"{'Resid':>6s}")
        print("  " + "-" * 58)
        for i, name in enumerate(result['names']):
            pred = result['beta'] * result['M2m1'][i]
            print(f"  {name:<30s} {result['M2m1'][i]:>6.2f} "
                  f"{result['Rm1'][i]:>6.2f} {pred:>6.2f} "
                  f"{result['residuals'][i]:>+6.2f}")

        print()
        print("  PHYSICAL INTERPRETATION:")
        print(f"    beta ~ {result['beta']:.2f}: each unit of M^2 above threshold")
        print(f"    amplifies the lensing-to-baryon ratio by {result['beta']:.2f}")
        print(f"    (scalar wake efficiency ~ {result['beta']*100:.0f}% per unit M^2)")

        print()
        print("  CAVEATS (current 9-cluster sample):")
        print("    * Sample too small for robust beta calibration")
        print("    * R_obs provenance heterogeneous (different WL pipelines)")
        print("    * No per-cluster covariance matrix (T_X-R correlation)")
        print("    * Selection bias: X-ray bow-shock detection requires M > 1")
        print("    * v_shock estimates span different epochs/methods")
        print()
        print("  ROADMAP:")
        print("    Euclid DR1 (~2027) + eROSITA clusters -> ~50 mergers")
        print("    with homogeneous WL masses and X-ray T_X.")
        print("    Required: per-cluster posterior P(M_lens, T_X | data)")
        print("    for joint hierarchical beta inference.")


if __name__ == "__main__":
    main()
