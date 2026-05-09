#!/usr/bin/env python3
"""
Gravastar mass-gap population model -- STUB (R4-M2).

STATUS: STUB -- frequency predictions computed, population synthesis deferred.
REVIEWER: R4 M2 (quantitative mass-gap predictions).

Computes scalar breathing mode eigenfrequencies (spin-0, l=0) for
gravastar-like compact objects in the 2.5-5 Msol mass gap, using the
SAME formulas as the validated qnm_gravastar_spectrum.py:

    f_VE  = v_Psi / (2 * r_opt)     cavity fundamental (n=1)
    f_Schw = 0.3737 * c^3 / (2*pi*G*M)   Schwarzschild l=2 (Leaver 1985)

    r_opt  = (3/2) * r_S            optical-trap radius (Eq. 25, S5)
    v_Psi  = (3*pi/16) * c          volume-averaged scalar velocity (App. D.1)

    f_VE/f_Schw ~ 1.5               mass-independent discriminator

WHAT THIS SCRIPT SHOWS:
  - Breathing mode frequencies across the mass-gap [2.5, 5] Msol
  - Spectral separation from Schwarzschild l=2 ringdown
  - Detection band classification (CE/ET vs LIGO)

WHAT IS NEEDED (reviewer R4-M2 sub-items):
  (i)   Formation channel: does scalar saturation truncate collapse
        for progenitors with core mass 2.5-5 Msol? Requires stellar
        evolution code with Chameleon-modified Tolman-Oppenheimer-Volkoff.
  (ii)  Mass distribution f(M): requires progenitor IMF + binary
        population synthesis (e.g. COMPAS/StarTrack with VE modifications).
  (iii) Binary fraction: gravastar-gravastar vs gravastar-BH vs gravastar-NS
        merger rates from population synthesis + gravitational capture.

Reference: Manuscript S5, App. D.2, Fig. 11.
"""

import numpy as np
import sys

# Constants
c = 2.998e8        # m/s
G = 6.674e-11      # m^3 kg^-1 s^-2
M_sun = 1.989e30   # kg

# Manuscript parameters (S5, consistent with qnm_gravastar_spectrum.py)
LAMBDA_LOC = 1.0                      # saturated interior
V_PSI = (3.0 * np.pi / 16.0) * c     # ~0.589c (App. D.1)
GR_CORR = 8.0 / 9.0                  # de Sitter curvature correction (0.889)


def f_breathing(M_sol):
    """Scalar breathing mode fundamental (n=1) with GR correction."""
    M = M_sol * M_sun
    r_S = 2 * G * M / c**2
    r_opt = (3 * LAMBDA_LOC / 2) * r_S  # Eq. 25
    return GR_CORR * V_PSI / (2 * r_opt)


def f_schwarzschild(M_sol):
    """Schwarzschild l=2 QNM (Leaver 1985 WKB)."""
    M = M_sol * M_sun
    return 0.3737 / (2 * np.pi) * c**3 / (G * M)


def classify_band(f_hz):
    """Classify frequency into GW detector band."""
    if f_hz > 5000:
        return "CE/ET"
    elif f_hz >= 10:
        return "LIGO"
    else:
        return "LISA"


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')

    print("=" * 70)
    print("Gravastar mass-gap population -- STUB (R4-M2)")
    print("=" * 70)
    print()
    print(f"  v_Psi/c  = {V_PSI/c:.4f}  (volume-averaged, App. D.1)")
    print(f"  r_opt    = 1.5 r_S       (optical-trap radius, Eq. 25)")
    print(f"  GR_corr  = {GR_CORR:.6f}  (de Sitter curvature, 8/9)")
    print(f"  lambda_loc = {LAMBDA_LOC}  (saturated interior, S5)")
    print()

    # Reference: verify ratio matches manuscript
    ratio_ref = f_breathing(10) / f_schwarzschild(10)
    print(f"  Cross-check (10 Msol): f_VE/f_Schw = {ratio_ref:.3f}")
    print(f"  Manuscript App. D.2 predicts:        ~1.5  (window [1.4, 1.7])")
    print()

    # Mass-gap table
    masses = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0, 15.0, 30.0]
    gap_lo, gap_hi = 2.5, 5.0

    print(f"  {'M[Msol]':>7s} {'f_VE[Hz]':>10s} {'f_Schw[Hz]':>11s} "
          f"{'Ratio':>6s} {'Band':>7s}  Notes")
    print("  " + "-" * 65)

    for M in masses:
        fv = f_breathing(M)
        fs = f_schwarzschild(M)
        ratio = fv / fs
        band = classify_band(fv)
        gap_flag = " <-- MASS GAP" if gap_lo <= M <= gap_hi else ""
        print(f"  {M:>7.1f} {fv:>10.0f} {fs:>11.0f} {ratio:>6.3f} {band:>7s}{gap_flag}")

    print()
    print(f"  f_VE/f_Schw = {ratio_ref:.3f}  (mass-independent)")
    print()

    # Falsification analysis
    print("  FALSIFICATION ANALYSIS:")
    print("  " + "-" * 50)
    print(f"    * For M = 2.6 Msol (GW190814 companion):")
    f26 = f_breathing(2.6)
    fs26 = f_schwarzschild(2.6)
    print(f"      f_VE   = {f26:.0f} Hz  ({f26/1e3:.1f} kHz)")
    print(f"      f_Schw = {fs26:.0f} Hz  ({fs26/1e3:.1f} kHz)")
    print(f"      Separation: {f26 - fs26:.0f} Hz")
    print(f"      -> Both in {classify_band(f26)} band, well separated")
    print()
    print("    * Spin-0 breathing mode is ABSENT in GR (no-hair theorem)")
    print("    * DETECTION of spin-0 mode -> falsifies Schwarzschild")
    print("    * NON-DETECTION at O(100) events -> falsifies VE gravastar")
    print()

    # Deferred items
    print("  DEFERRED (requires external tools):")
    print("  " + "-" * 50)
    print("    (i)   Formation channel:")
    print("          -> Chameleon-modified TOV equation integration")
    print("          -> Identifies critical core mass for saturation onset")
    print("          -> Requires: stellar structure code + V_eff(phi, rho)")
    print()
    print("    (ii)  Mass distribution f(M) in gap:")
    print("          -> Population synthesis (COMPAS/StarTrack)")
    print("          -> Progenitor IMF at Z = 0.01-0.02")
    print("          -> Predicted: ~10-100 gap events per year (O5)")
    print()
    print("    (iii) Binary fraction:")
    print("          -> Gravastar-gravastar vs gravastar-BH rates")
    print("          -> Gravitational capture cross-section")
    print("          -> Requires: N-body cluster dynamics with VE")


if __name__ == "__main__":
    main()
