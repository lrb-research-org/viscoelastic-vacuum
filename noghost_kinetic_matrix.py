#!/usr/bin/env python3
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass
"""
No-ghost kinetic matrix D_AB — symbolic verification stub.
STATUS: STUB (R1-C4). Full post-constraint reduction in companion paper.
Reference: Manuscript App B.3-B.5, Eq. 63-66.

Verifies that with D_AB = diag(1,1) (adopted in Eq. 65),
the dispersion relation yields omega^2 > 0 for all k (no tachyons),
and that c_T = c identically (f(R,T) linear in R).
"""
import numpy as np

def dispersion_check():
    """Check det[K·omega^2 - D·k² - M²] = 0 for stability."""
    print("No-ghost kinetic matrix D_AB — STUB (R1-C4)")
    print("="*55)
    # Adopted values (Eq. 65): K_AB = diag(1, 1), D_AB = diag(1, 1)
    K = np.diag([1.0, 1.0])
    D = np.diag([1.0, 1.0])
    # Mass matrix from Chameleon + f(R,T) trace mode
    # m_eff² = d²V_eff/dφ² evaluated at φ_min(ρ)
    m_eff_sq = 1e-26  # (eV/c²)² ~ cosmological mass
    M2 = np.diag([m_eff_sq, 0.0])  # scalar massive, tensor massless
    
    print(f"\nK_AB = {np.diag(K)}")
    print(f"D_AB = {np.diag(D)}")
    print(f"M^2_AB = {np.diag(M2)}")
    
    # For each mode, omega^2(k) = (D/K)·k² + M²/K
    for i, label in enumerate(["Scalar (Psi)", "Tensor (h)"]):
        c_eff_sq = D[i,i] / K[i,i]
        m_sq = M2[i,i] / K[i,i]
        print(f"\n  {label}:")
        print(f"    omega^2(k) = {c_eff_sq:.1f}·k² + {m_sq:.1e}")
        print(f"    c_eff²/c² = {c_eff_sq:.1f} → c_T = c ✓")
        print(f"    omega^2(k=0) = {m_sq:.1e} ≥ 0 → no tachyon ✓")
        print(f"    domega^2/dk² = {c_eff_sq:.1f} > 0 → no gradient instability ✓")
    
    print(f"\nNOTE: STRUCTURAL NOTE:")
    print(f"   D_AB = I is imposed by fiat (Eq. 65), not derived.")
    print(f"   Full derivation requires post-constraint Hamiltonian")
    print(f"   reduction of the f(R,T) + Chameleon SVT decomposition.")
    print(f"   This is deferred to the companion paper (Table 14, item 5).")
    print(f"\n   c_T = c is a STRUCTURAL feature of f(R,T) being linear")
    print(f"   in R (no higher-derivative tensor sector), NOT a tuned result.")

if __name__ == "__main__":
    dispersion_check()
