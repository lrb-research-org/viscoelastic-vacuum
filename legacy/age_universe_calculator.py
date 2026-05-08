"""
age_universe_calculator.py
==========================
Age of the Universe in the Viscoelastic Vacuum Framework.

Demonstrates that the apparent 'age tension' (t0 ~ 10.3 Gyr in EdS
vs. ~12.8 Gyr from globular clusters) is an artifact of the pure-dust
(Omega_m=1) approximation.

In the VE model, the topological condition Omega_total = 1 includes:
  - Omega_bar ~ 0.316  (baryons + Chameleon-frozen condensate)
  - Omega_phi ~ 0.684  (scalar field background energy)

The frozen scalar condensate at z < z_drag has an effective equation
of state w_phi that depends on the Chameleon potential V(phi).
For the runaway potential V(phi) = M^4 exp(M/phi), the field frozen
at its effective potential minimum has w_eff ~ -1 (potential energy
dominates kinetic), giving t0 ~ 14-15 Gyr.

USAGE:
    python age_universe_calculator.py
    python age_universe_calculator.py --H0 63.2 --Omega_bar 0.316

DEPENDENCIES:
    numpy, scipy

References:
    Bouille (2026), "Viscoelastic Vacuum" preprint, Table 9
    Khoury & Weltman (2004), PRD 69, 044026 (Chameleon fields)
    Valcin et al. (2021), JCAP 12, 002 (globular cluster ages)
"""

import sys
import numpy as np
from scipy.integrate import quad

# Default parameters
H0_DEFAULT = 63.2       # km/s/Mpc (VE Pantheon+ best-fit)
OMEGA_BAR_DEFAULT = 0.316  # baryons + frozen condensate (1612/5103)
MPC_KM = 3.086e19       # 1 Mpc in km
GYR_PER_S = 1 / 3.156e16  # seconds to Gyr


def age_of_universe(H0, Omega_bar, w_phi, Omega_phi=None):
    """
    Compute the age of the universe by integrating:
        t_0 = integral_0^infty dz / [(1+z) H(z)]

    H^2(z) = H0^2 [Omega_bar (1+z)^3 + Omega_phi (1+z)^{3(1+w_phi)}]

    Parameters:
        H0 : Hubble constant (km/s/Mpc)
        Omega_bar : matter density parameter (baryons + frozen)
        w_phi : equation of state of scalar field background
        Omega_phi : scalar field density (default: 1 - Omega_bar)

    Returns:
        t0 in Gyr
    """
    if Omega_phi is None:
        Omega_phi = 1.0 - Omega_bar

    def integrand(z):
        Hz = H0 * np.sqrt(
            Omega_bar * (1 + z) ** 3
            + Omega_phi * (1 + z) ** (3 * (1 + w_phi))
        )
        return 1.0 / ((1 + z) * Hz)

    t0_inv_H, _ = quad(integrand, 0, np.inf)
    return t0_inv_H * MPC_KM * GYR_PER_S


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')

    H0 = H0_DEFAULT
    Omega_bar = OMEGA_BAR_DEFAULT

    # Parse optional CLI args
    args = sys.argv[1:]
    for i, a in enumerate(args):
        if a == '--H0' and i + 1 < len(args):
            H0 = float(args[i + 1])
        elif a == '--Omega_bar' and i + 1 < len(args):
            Omega_bar = float(args[i + 1])

    Omega_phi = 1.0 - Omega_bar

    print()
    print("=" * 72)
    print("AGE OF THE UNIVERSE — VISCOELASTIC VACUUM FRAMEWORK")
    print("=" * 72)
    print(f"  H0         = {H0} km/s/Mpc")
    print(f"  Omega_bar  = {Omega_bar:.4f} (baryons + frozen Chameleon)")
    print(f"  Omega_phi  = {Omega_phi:.4f} (scalar field background)")
    print(f"  Omega_tot  = {Omega_bar + Omega_phi:.4f}")
    print()

    # Reference values
    t_eds = age_of_universe(H0, 1.0, w_phi=0.0, Omega_phi=0.0)
    t_lcdm = age_of_universe(67.4, 0.315, w_phi=-1.0, Omega_phi=0.685)

    print("  --- Reference values ---")
    print(f"  EdS pure dust (Omega_m=1, w=0)   : t0 = {t_eds:.2f} Gyr")
    print(f"  LCDM (H0=67.4, Omega_L=0.685)    : t0 = {t_lcdm:.2f} Gyr")
    print(f"  Globular cluster ages (obs.)      : t0 >= 12.8 +/- 0.5 Gyr")
    print(f"  JWST early-galaxy tension suggests: t0 >= 14 Gyr")
    print()

    print("  --- VE model: t0 vs w_phi ---")
    print(f"  {'w_phi':>8s}  {'t0 (Gyr)':>10s}  {'Status':>10s}  Interpretation")
    print(f"  {'-------':>8s}  {'--------':>10s}  {'------':>10s}  --------------")

    w_values = [0.0, -1/6, -1/3, -1/2, -2/3, -5/6, -1.0]
    for w in w_values:
        t = age_of_universe(H0, Omega_bar, w_phi=w)
        if t >= 12.5:
            status = "OK"
        else:
            status = "TOO YOUNG"

        if w == 0.0:
            interp = "pressureless dust (EdS-like)"
        elif abs(w - (-1/3)) < 0.01:
            interp = "curvature-like"
        elif abs(w - (-2/3)) < 0.01:
            interp = "slow-roll quintessence"
        elif abs(w - (-1.0)) < 0.01:
            interp = "frozen potential (= effective Lambda)"
        else:
            interp = ""

        print(f"  {w:8.3f}  {t:10.2f}  {status:>10s}  {interp}")

    # Best estimate
    print()
    print("  --- Best estimate ---")
    print("  For a Chameleon field frozen at its effective potential minimum,")
    print("  V(phi) = M^4 exp(M/phi), kinetic energy << potential energy,")
    print("  so w_eff -> -1 at late times (Khoury & Weltman 2004).")
    print()
    t_best = age_of_universe(H0, Omega_bar, w_phi=-1.0)
    print(f"  ==> t0(VE, best) = {t_best:.1f} Gyr")
    print(f"      vs. LCDM     = {t_lcdm:.1f} Gyr")
    print(f"      Difference   = {t_best - t_lcdm:+.1f} Gyr")
    print()
    print("  The VE model predicts an OLDER universe than LCDM,")
    print("  consistent with JWST early-galaxy observations and")
    print("  recent globular-cluster age recalibrations.")
    print()

    # Sensitivity to H0
    print("  --- Sensitivity to H0 (w_phi = -1) ---")
    for h0 in [60.0, 62.0, 63.2, 65.0, 67.4, 70.0, 73.6]:
        t = age_of_universe(h0, Omega_bar, w_phi=-1.0)
        marker = " <-- VE best-fit" if abs(h0 - 63.2) < 0.1 else ""
        marker = " <-- Planck" if abs(h0 - 67.4) < 0.1 else marker
        marker = " <-- SH0ES" if abs(h0 - 73.6) < 0.1 else marker
        print(f"  H0 = {h0:5.1f}  =>  t0 = {t:.2f} Gyr{marker}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
