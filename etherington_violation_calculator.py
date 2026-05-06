"""
etherington_violation_calculator.py
====================================
Standalone calculator for the Etherington distance-duality violation
predicted by the Viscoelastic Vacuum framework (Bouille 2026).

PREDICTION:
    In Lambda-CDM, the Etherington reciprocity relation is:
        D_L / [D_A * (1+z)^2] = 1   exactly.

    In the viscoelastic vacuum framework, photon-scalar conversion along
    the line of sight produces a chromatic opacity tau(z) ~ epsilon * 2 * ln(1+z),
    which violates this relation:
        D_L^eff / [D_A * (1+z)^2] = (1+z)^epsilon

    where epsilon = 1/d_f = 0.40 is the topological opacity exponent
    (d_f = 5/2 is the fractal dimension of the cosmic web).

DISCRIMINANT TEST:
    Cross-correlate Euclid weak-lensing angular distances D_A(z) with
    LSST supernova luminosity distances D_L(z) for the SAME redshift bin.
    Lambda-CDM predicts ratio = 1; this framework predicts (1+z)^0.40.

    At z = 1: violation = (1+1)^0.40 = 1.32 = 32% deviation, well above
    instrumental precision (Euclid ~1-2% on D_A, LSST ~5-8% on D_L).

USAGE:
    python etherington_violation_calculator.py [redshift] [epsilon]

EXAMPLES:
    python etherington_violation_calculator.py 1.0
    python etherington_violation_calculator.py 0.5 0.40
    python etherington_violation_calculator.py 2.0 0.40
"""
import sys


EPSILON_DEFAULT = 0.40   # = 1/d_f, with d_f = 5/2 (3D percolation universality)


def etherington_ratio(z, epsilon=EPSILON_DEFAULT):
    """
    Predicted ratio D_L^eff / [D_A * (1+z)^2] in viscoelastic framework.

    Returns:
        float : (1+z)^epsilon, the violation factor.
    """
    return (1.0 + z) ** epsilon


def violation_percent(z, epsilon=EPSILON_DEFAULT):
    """Deviation from Lambda-CDM (which predicts ratio = 1)."""
    return 100.0 * (etherington_ratio(z, epsilon) - 1.0)


def report(z, epsilon=EPSILON_DEFAULT):
    """Print a human-readable predictions report."""
    ratio = etherington_ratio(z, epsilon)
    pct = violation_percent(z, epsilon)
    print()
    print("=" * 70)
    print(f"ETHERINGTON VIOLATION at z = {z:.3f}, epsilon = {epsilon:.3f}")
    print("=" * 70)
    print(f"  Lambda-CDM prediction:       D_L / [D_A * (1+z)^2] = 1.000")
    print(f"  Viscoelastic prediction:     D_L / [D_A * (1+z)^2] = {ratio:.4f}")
    print(f"  Violation magnitude:         {pct:+.2f}%")
    print()
    print(f"  Implied D_L^eff increase:    +{pct:.2f}% (chromatic dimming)")
    print(f"  Equivalent SNIa magnitude:   +{2.5*pct/100:.4f} mag")
    print()


def scan_table():
    """Print a scan table over z = 0.1 to 2.5."""
    print("=" * 70)
    print("DISCRIMINANT TEST TABLE  (Euclid x LSST cross-correlation)")
    print("=" * 70)
    print(f"  Topological prediction: epsilon = 1/d_f = 1/(5/2) = {EPSILON_DEFAULT}")
    print()
    print(f"   z        ratio = (1+z)^eps    violation %    SNIa Delta_mu (mag)")
    print(f"  ------    ------------------   ------------   -------------------")
    for z in [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.2, 1.5, 1.8, 2.0, 2.3, 2.5]:
        r = etherington_ratio(z)
        pct = violation_percent(z)
        dmu = 2.5 * pct / 100
        print(f"  {z:5.2f}    {r:18.4f}   {pct:+11.2f}%   {dmu:+19.4f}")
    print()
    print("  NB: at z=1 the violation reaches 32%, far above the joint")
    print("      Euclid (D_A ~1-2%) x LSST (D_L ~5-8%) precision floor.")
    print("      A null detection of (1+z)^eps signal would falsify epsilon = 0.40.")


def main():
    if len(sys.argv) >= 2:
        try:
            z = float(sys.argv[1])
        except ValueError:
            print(f"Usage: {sys.argv[0]} [redshift] [epsilon]")
            return 1
        epsilon = float(sys.argv[2]) if len(sys.argv) >= 3 else EPSILON_DEFAULT
        report(z, epsilon)
        return 0
    else:
        scan_table()
        return 0


if __name__ == "__main__":
    sys.exit(main())
