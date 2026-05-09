#!/usr/bin/env python3
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass
"""
β-coefficient calibration for the lensing amplification R(v) scaling.

STATUS: STUB — methodology documented, awaiting Euclid/Rubin LSST ~50-cluster sample.
REVIEWER: R4 §2.2 (quantitative β calibration).

The lensing-to-baryon mass ratio R = M_lens / M_baryon scales with the
collision velocity via the Mach-cone scalar retardation mechanism:

    R(v) = 1 + β · (v / c_s)^2      for M = v/c_s > 1
    R(v) = 1                          for M = v/c_s ≤ 1

where β encodes the efficiency of scalar momentum transfer from the
Liénard–Wiechert retarded potential to the effective lensing convergence.

WHAT THIS SCRIPT DOES (current state):
  - Loads the 9-cluster sample from cluster_data.csv
  - Fits β from the R vs M^2 relation
  - Reports the best-fit β and its uncertainty
  - Documents the systematic limitations of the current sample

WHAT IS NEEDED FOR PUBLICATION-GRADE CALIBRATION:
  1. Euclid / Rubin LSST weak-lensing mass estimates for ~50 merging clusters
  2. Chandra/eROSITA X-ray temperatures and shock velocities
  3. Per-cluster covariance (shear noise + photo-z + ICM T_X systematics)
  4. Selection function correction (X-ray vs optically selected mergers)

Reference: Manuscript §4, Eq. 35-36; Table 6.
"""

import numpy as np
import os

def load_cluster_data():
    """Load the 9-cluster sample from cluster_data.csv."""
    csv_path = os.path.join(os.path.dirname(__file__), "cluster_data.csv")
    
    # CSV columns: cluster, v_shock_kms(1), v_shock_ref(2), T_pre_keV(3),
    #   T_pre_ref(4), n_e_cm3(5), n_e_ref(6), M_gas_1e14Msun(7),
    #   M_gas_ref(8), R_obs(9), R_err(10), ...
    clusters = []
    with open(csv_path, 'r') as f:
        header = f.readline().strip().split(',')
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 10:
                continue
            name = parts[0].strip()
            try:
                v_shock = float(parts[1])   # km/s
                T_X = float(parts[3])       # keV (T_pre_keV)
                R_ratio = float(parts[9])   # R_obs
            except (ValueError, IndexError):
                continue
            clusters.append({
                'name': name,
                'v_shock': v_shock,
                'T_X': T_X,
                'R': R_ratio
            })
    return clusters

def sound_speed_from_TX(T_keV):
    """ICM sound speed from X-ray temperature.
    
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
    """Fit β from R = 1 + β · M^2 using least-squares.
    
    Only supersonic mergers (M > 1) are included.
    """
    M2_vals = []
    R_vals = []
    names = []
    
    for c in clusters:
        if np.isnan(c['R']) or c['R'] <= 0:
            continue
        c_s = sound_speed_from_TX(c['T_X'])
        M = c['v_shock'] / c_s
        if M > 1.0:
            M2_vals.append(M**2)
            R_vals.append(c['R'])
            names.append(c['name'])
    
    if len(M2_vals) < 2:
        print("ERROR: insufficient supersonic mergers for fit")
        return None
    
    M2 = np.array(M2_vals)
    R = np.array(R_vals)
    
    # R - 1 = β · M^2  →  β = <(R-1)> / <M^2>  (weighted)
    # Simple OLS: (R-1) = β * M^2
    beta = np.sum((R - 1) * M2) / np.sum(M2**2)
    residuals = (R - 1) - beta * M2
    sigma_beta = np.sqrt(np.sum(residuals**2) / ((len(M2) - 1) * np.sum(M2**2)))
    
    return {
        'beta': beta,
        'sigma_beta': sigma_beta,
        'n_clusters': len(M2),
        'M2': M2,
        'R': R,
        'names': names,
        'residuals': residuals
    }

def main():
    print("=" * 65)
    print("β-coefficient calibration — STUB (R4 §2.2)")
    print("=" * 65)
    print()
    
    clusters = load_cluster_data()
    print(f"Loaded {len(clusters)} clusters from cluster_data.csv")
    print()
    
    # Compute Mach numbers
    print(f"{'Cluster':<25} {'v_shock':>8} {'c_s':>8} {'M':>6} {'R':>6}")
    print("-" * 60)
    for c in clusters:
        c_s = sound_speed_from_TX(c['T_X'])
        M = c['v_shock'] / c_s
        R_str = f"{c['R']:.2f}" if not np.isnan(c['R']) else "N/A"
        print(f"{c['name']:<25} {c['v_shock']:>8.0f} {c_s:>8.0f} {M:>6.2f} {R_str:>6}")
    
    print()
    result = fit_beta(clusters)
    if result:
        print(f"Best-fit β = {result['beta']:.4f} ± {result['sigma_beta']:.4f}")
        print(f"  from {result['n_clusters']} supersonic mergers")
        print(f"  R² = {1 - np.var(result['residuals'])/np.var(result['R']-1):.3f}")
        print()
        print("⚠️  CAVEATS (current 9-cluster sample):")
        print("  • Sample too small for robust β calibration")
        print("  • R_lens/R_bar heterogeneous provenance (different WL pipelines)")
        print("  • No per-cluster covariance matrix")
        print("  • Selection bias: X-ray-bright mergers over-represented")
        print()
        print("ROADMAP: Euclid DR1 (~2027) + eROSITA clusters → ~50 mergers")
        print("         with homogeneous WL masses and X-ray T_X.")

if __name__ == "__main__":
    main()
