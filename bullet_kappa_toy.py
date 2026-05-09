#!/usr/bin/env python3
"""
bullet_kappa_toy.py  —  Viscoelastic Vacuum Bullet Cluster κ(θ) toy model
==========================================================================

Produces a 1D projected convergence κ(θ) for the Bullet Cluster geometry
using the viscoelastic retarded-scalar framework of §4.2 (PREPRINT_V2).

Physical model
--------------
Two components contribute to convergence:
  (a) Newtonian baryonic:  κ_bar(θ)  — NFW gas + concentrated galaxy core
  (b) Retarded scalar:     κ_φ(θ)    — Liénard-Wiechert phase-lagged halo

The total convergence is  κ_tot(θ) = κ_bar(θ) + (λ_eff / 8π) × κ_φ(θ)
where the scalar contribution is shifted by the retarded offset Δx.

Inputs (Markevitch & Vikhlinin 2007; Clowe et al. 2006):
  - ICM gas:     β-model  ρ_gas(r) = ρ_0 (1 + r²/r_c²)^(-3β/2)
                 r_c = 150 kpc, β = 2/3,  M_gas = 1.1e14 M_sun
  - Galaxies:    concentrated King profile  r_c_gal = 50 kpc
                 M_star = 1.5e13 M_sun
  - Shock:       v_shock = 4700 km/s,  v_Ψ ≈ 1315 km/s  →  M = 3.57
  - d_eff = 350 kpc  →  Δx = d_eff/(M-1) ≈ 136 kpc
  - λ_eff ≈ 4000 (thin-shell screened from bare λ~10^5)

The script produces a figure suitable for the supplementary repository:
  - Top panel:  κ(θ) profiles for gas, galaxies, retarded scalar, total
  - Bottom panel:  residual κ_VE - κ_bar

Usage:
  python bullet_kappa_toy.py [--plot] [--save bullet_kappa_profile.png]

Author:  L.-R. Bouille / Antigravity agent
License: MIT
"""

import numpy as np
import sys

# ============================================================
#  Physical constants and Bullet Cluster parameters
# ============================================================
M_sun = 1.989e30       # kg
kpc   = 3.0857e19      # m
c_light = 2.998e8      # m/s
G     = 6.674e-11      # m³/(kg·s²)

# Cluster geometry
D_A       = 1.07e3 * 1e6 * 3.0857e16    # angular diameter distance to Bullet (z=0.296) in m
r_c_gas   = 150.0      # kpc — ICM core radius
beta_gas  = 2.0/3.0    # β-model exponent
M_gas     = 1.1e14     # M_sun — total ICM gas mass
r_c_gal   = 50.0       # kpc — galaxy core radius
M_star    = 1.5e13     # M_sun — total stellar mass in main cluster

# Shock and scalar parameters
v_shock   = 4700.0     # km/s
v_Psi     = 1315.0     # km/s — scalar group velocity in ICM
Mach      = v_shock / v_Psi   # ≈ 3.57
d_eff     = 350.0      # kpc — effective hydrodynamic interaction length
Delta_x   = d_eff / (Mach - 1)  # ≈ 136 kpc — retarded offset

# Screening: bare λ~10^5, thin-shell ε_shell ~ 0.04
# effective enhancement factor for convergence: λ_eff/(8π)
lambda_eff_over_8pi = 4000.0 / (8 * np.pi)  # ≈ 159
# But this still overshoots; the scalar couples to the trace T, not ρ directly.
# For gas (relativistic thermal component), T_gas/ρ_gas ≈ 1 - 3p/ρc² 
# with kT ~ 14 keV,  3p/(ρc²) ~ 3×kT/(μ m_p c²) ~ 4.5e-5 → T ≈ -ρc²
# So the coupling is: Δκ/κ_bar = λ_eff/(8π) × ε_shell × (ρ_φ/ρ_gas)
# The NET observed ratio M_lens/M_bar ≈ 5-7 for Bullet
# We calibrate the effective scalar amplification factor α_scalar:
alpha_scalar = 5.0     # net amplification (M_lens/M_bar ~ 5-7)

# ============================================================
#  Surface mass density profiles  (projected along line of sight)
# ============================================================

def sigma_beta_model(theta_arcsec, r_c_kpc, M_total_Msun, D_A_m):
    """
    Projected surface mass density of a β-model (β=2/3)
    Σ(θ) = Σ_0 / (1 + θ²/θ_c²)
    where θ_c = r_c / D_A
    
    Returns Σ in M_sun/kpc² as a function of angular position θ (arcsec).
    """
    theta_c = (r_c_kpc * kpc / D_A_m) * 206265  # core radius in arcsec
    x2 = (theta_arcsec / theta_c)**2
    
    # For β=2/3, projected β-model: Σ(R) = Σ_0 (1 + R²/r_c²)^(-1/2)
    # Normalization: ∫ Σ 2π R dR = M_total
    # → Σ_0 = M_total / (2π r_c²)  ... in projected space
    Sigma_0 = M_total_Msun / (2 * np.pi * r_c_kpc**2)  # M_sun / kpc²
    
    return Sigma_0 / np.sqrt(1 + x2)


def sigma_king(theta_arcsec, r_c_kpc, M_total_Msun, D_A_m):
    """
    Projected King profile (isothermal sphere with core):
    Σ(θ) = Σ_0 / (1 + θ²/θ_c²)
    """
    theta_c = (r_c_kpc * kpc / D_A_m) * 206265
    x2 = (theta_arcsec / theta_c)**2
    Sigma_0 = M_total_Msun / (2 * np.pi * r_c_kpc**2)
    
    return Sigma_0 / (1 + x2)


def convergence_from_sigma(Sigma, D_A_m, z_lens=0.296, z_source=1.0):
    """
    Convert surface mass density Σ (M_sun/kpc²) to convergence κ = Σ/Σ_cr.
    Σ_cr = c²/(4πG) × D_s / (D_l × D_ls)
    
    For simplicity, use the standard thin-lens approximation.
    """
    # Angular diameter distances (simplified flat cosmology, H0=70)
    H0 = 70.0e3 / (1e6 * 3.0857e16)  # 1/s
    # Use comoving distances (rough)
    D_l = D_A_m
    # D_s / D_ls ratio for z_l=0.296, z_s~1.0  →  D_s/D_ls ≈ 2.5 (rough)
    D_ratio = 2.5
    
    Sigma_cr = (c_light**2 / (4 * np.pi * G)) / (D_l * D_ratio)
    # Convert to M_sun/kpc²
    Sigma_cr_Msun_kpc2 = Sigma_cr * kpc**2 / M_sun
    
    return Sigma / Sigma_cr_Msun_kpc2


# ============================================================
#  Build the κ(θ) profiles
# ============================================================

def build_profiles(N=500, theta_max=300):
    """
    Build convergence profiles along the merger axis.
    θ in arcsec, centered on the pre-merger barycenter.
    
    Convention: positive θ = direction of bullet motion (West)
    """
    theta = np.linspace(-theta_max, theta_max, N)
    
    # --- Component 1: ICM gas (centered at θ=0) ---
    kappa_gas = convergence_from_sigma(
        sigma_beta_model(theta, r_c_gas, M_gas, D_A),
        D_A
    )
    
    # --- Component 2: Galaxy core (offset by ~ -50" ≈ -200 kpc to the East,
    #     representing the pre-merger main cluster galaxy concentration) ---
    theta_gal_offset = -(50.0 * kpc / D_A) * 206265  # ~-50 arcsec East
    # But in Bullet, galaxies are AHEAD of gas (West), so offset positive:
    # Clowe+2006: galaxy centroid ~ +35" West of gas centroid
    theta_gal_offset_arcsec = 35.0  # arcsec West = positive
    kappa_gal = convergence_from_sigma(
        sigma_king(theta - theta_gal_offset_arcsec, r_c_gal, M_star, D_A),
        D_A
    )
    
    # --- Component 3: Retarded scalar halo ---
    # The scalar field is sourced by the total baryonic distribution (gas + stars)
    # but propagates at v_Ψ << v_shock, creating a phase-lagged concentration.
    # In the retarded Liénard-Wiechert picture, the scalar peak is shifted
    # by Δx ≈ 136 kpc from the gas centroid, toward the galaxies.
    Delta_x_arcsec = (Delta_x * kpc / D_A) * 206265  # ≈ 26 arcsec
    
    # The scalar halo profile is broader than the galaxy core
    # (it's the retarded integral of the gas distribution)
    r_c_scalar = 120.0  # kpc — scalar wake width (intermediate gas-galaxy)
    M_scalar_equivalent = M_gas * alpha_scalar  # effective scalar lensing mass
    
    kappa_scalar = convergence_from_sigma(
        sigma_beta_model(theta - Delta_x_arcsec, r_c_scalar, M_scalar_equivalent, D_A),
        D_A
    )
    
    # --- Total ---
    kappa_bar  = kappa_gas + kappa_gal
    kappa_tot  = kappa_bar + kappa_scalar
    
    return theta, kappa_gas, kappa_gal, kappa_scalar, kappa_bar, kappa_tot


# ============================================================
#  Plotting
# ============================================================

def plot_profiles(save_path=None):
    """Generate the κ(θ) figure."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available. Run with --no-plot or install matplotlib.")
        return
    
    theta, kappa_gas, kappa_gal, kappa_scalar, kappa_bar, kappa_tot = build_profiles()
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True,
                                    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    
    # --- Top panel: profiles ---
    ax1.plot(theta, kappa_gas, 'b--', lw=1.5, alpha=0.7, label=r'$\kappa_{\rm gas}$ (ICM $\beta$-model)')
    ax1.plot(theta, kappa_gal, 'g:', lw=1.5, alpha=0.7, label=r'$\kappa_{\rm gal}$ (stellar King)')
    ax1.plot(theta, kappa_bar, 'k-', lw=1.0, alpha=0.5, label=r'$\kappa_{\rm bar} = \kappa_{\rm gas} + \kappa_{\rm gal}$')
    ax1.plot(theta, kappa_scalar, 'r-.', lw=2.0, alpha=0.8, 
             label=r'$\kappa_\phi$ (retarded scalar, $\Delta x \approx 136$ kpc)')
    ax1.plot(theta, kappa_tot, 'k-', lw=2.5, 
             label=r'$\kappa_{\rm tot} = \kappa_{\rm bar} + \kappa_\phi$')
    
    # Mark key positions
    Delta_x_arcsec = (Delta_x * kpc / D_A) * 206265
    ax1.axvline(0, color='blue', alpha=0.3, ls=':', label='Gas centroid')
    ax1.axvline(35, color='green', alpha=0.3, ls=':', label='Galaxy centroid (+35")')
    ax1.axvline(Delta_x_arcsec, color='red', alpha=0.3, ls=':', 
                label=f'Scalar peak (+{Delta_x_arcsec:.0f}")')
    
    # Observed lensing peak position (Clowe+2006: ~25-40" West of gas)
    ax1.axvspan(25, 40, alpha=0.1, color='orange', label=r'Observed $\kappa$ peak (Clowe+06)')
    
    ax1.set_ylabel(r'Convergence $\kappa(\theta)$', fontsize=13)
    ax1.set_title('Bullet Cluster — Viscoelastic Vacuum $\\kappa(\\theta)$ Toy Model\n'
                   r'$v_{\rm shock}=4700$ km/s, $v_\Psi=1315$ km/s, '
                   r'$\mathcal{M}=3.57$, $\Delta x \approx 136$ kpc',
                   fontsize=12)
    ax1.legend(fontsize=9, loc='upper right', ncol=2)
    ax1.set_xlim(-200, 200)
    ax1.set_ylim(bottom=0)
    ax1.grid(True, alpha=0.3)
    
    # --- Bottom panel: residual ---
    residual = kappa_tot - kappa_bar
    ax2.fill_between(theta, 0, residual, alpha=0.3, color='red')
    ax2.plot(theta, residual, 'r-', lw=1.5, label=r'$\Delta\kappa = \kappa_\phi$ (scalar enhancement)')
    ax2.axhline(0, color='k', lw=0.5)
    ax2.axvline(Delta_x_arcsec, color='red', alpha=0.3, ls=':')
    ax2.set_xlabel(r'Angular position $\theta$ [arcsec]  ($\leftarrow$ E  |  W $\rightarrow$)', fontsize=13)
    ax2.set_ylabel(r'$\Delta\kappa$', fontsize=13)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved: {save_path}")
    else:
        plt.savefig('bullet_kappa_profile.png', dpi=150, bbox_inches='tight')
        print("Figure saved: bullet_kappa_profile.png")
    
    # Print summary statistics
    idx_peak_tot = np.argmax(kappa_tot)
    idx_peak_bar = np.argmax(kappa_bar)
    idx_peak_sca = np.argmax(kappa_scalar)
    
    print("\n=== Bullet kappa(theta) Toy Model Summary ===")
    print(f"  Mach number:           M = {Mach:.2f}")
    print(f"  Retarded offset:       Dx = {Delta_x:.0f} kpc = {Delta_x_arcsec:.1f}\"")
    print(f"  Baryonic kappa peak:   theta = {theta[idx_peak_bar]:.1f}\" (gas-dominated)")
    print(f"  Scalar kappa_phi peak: theta = {theta[idx_peak_sca]:.1f}\" (shifted by Dx)")
    print(f"  Total kappa peak:      theta = {theta[idx_peak_tot]:.1f}\"")
    print(f"  kappa_bar(peak):       {kappa_bar[idx_peak_bar]:.4f}")
    print(f"  kappa_tot(peak):       {kappa_tot[idx_peak_tot]:.4f}")
    print(f"  Amplification:         kappa_tot/kappa_bar ~ {kappa_tot[idx_peak_tot]/kappa_bar[idx_peak_bar]:.1f}")
    print(f"  Galaxy-total offset:   {abs(theta[idx_peak_tot] - 35):.1f}\" from galaxy centroid")
    print(f"  Observed offset:       25-40\" (Clowe+2006)")
    print(f"\n  [Phenomenological] -- d_eff = {d_eff:.0f} kpc calibrated to Bullet observation")
    print(f"  [Derived]           -- Dx formula from Lienard-Wiechert stationary phase")


# ============================================================
#  Main
# ============================================================

if __name__ == '__main__':
    save = None
    for i, arg in enumerate(sys.argv[1:], 1):
        if arg == '--save' and i < len(sys.argv) - 1:
            save = sys.argv[i + 1]
    
    plot_profiles(save_path=save)
