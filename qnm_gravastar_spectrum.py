"""
qnm_gravastar_spectrum.py
=========================
Quasi-Normal Mode (QNM) spectrum of scalar breathing modes
for gravastar-like objects in the Viscoelastic Vacuum framework.

PHYSICS:
    In the VE framework, geometric singularities are regularised by scalar
    saturation at r_opt = (3 lambda / 2) * r_S, where lambda ~ 10^5 is the
    cosmological-scale coupling constant (= (c / sigma_v)^2). This produces
    a gravastar-like compact object with a de Sitter interior matched to
    an exterior Schwarzschild geometry via Darmois-Israel junction conditions.

    The scalar field phi propagates inside this cavity with effective
    velocity v_Psi. In the manuscript (Outlook vi), the characteristic
    breathing-mode frequency is:
        f_n ~ v_Psi / (2 L)  with L ~ r_opt

    For M = 10 M_sun:
        r_S  = 29.5 km
        r_opt = (3 * 10^5 / 2) * r_S = 4.43 * 10^6 km
        f_0 ~ v_Psi / (4 * r_opt)

    With v_Psi ~ 0.6c, this gives f_0 ~ 10^-2 Hz (LISA band).

QNM APPROACH:
    We solve the radial Klein-Gordon eigenvalue problem for a scalar
    field in the de Sitter interior of the gravastar cavity, with
    Dirichlet boundary conditions at the viscoelastic membrane (r = r_opt).

    The eigenfrequencies omega_n are extracted numerically from the
    finite-difference discretisation of:
        -v_Psi^2 * d^2 psi/dr^2 + V_eff(r) * psi = omega^2 * psi

PREDICTIONS:
    - f_0 ~ 1.7 * 10^-2 Hz * (10 M_sun / M) for v_Psi = c
    - f_0 ~ 10^-2 Hz * (10 M_sun / M) for v_Psi ~ 0.6c
    - Schwarzschild BH has NO spin-0 QNM (no-hair theorem)
    - Detection of spin-0 ringdown = falsification of GR BH

FALSIFICATION:
    f_0 != 10^-2 Hz * (10 M_sun / M)^-1 at >30% => VE gravastar falsified.

USAGE:
    python qnm_gravastar_spectrum.py [M_solar] [n_modes]
    python qnm_gravastar_spectrum.py           # default: 10 M_sun, 6 modes
    python qnm_gravastar_spectrum.py 30        # 30 M_sun, 6 modes
    python qnm_gravastar_spectrum.py 10 10     # 10 M_sun, 10 modes

DEPENDENCIES:
    numpy, scipy, matplotlib (optional for plot)

References:
    Bouille (2026), "Viscoelastic Vacuum" preprint, Sec. 5, Outlook (vi)
    Mazur & Mottola, Proc. Natl. Acad. Sci. 101, 9545 (2004)
    Visser & Wiltshire, CQG 21, 1135 (2004)
"""

import sys
import numpy as np
try:
    from scipy import linalg
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


# =============================================================================
#  PHYSICAL CONSTANTS (SI)
# =============================================================================
c_SI       = 2.998e8       # m/s
G_SI       = 6.674e-11     # m^3 kg^-1 s^-2
M_sun_SI   = 1.989e30      # kg
hbar_SI    = 1.055e-34     # J.s

# Framework parameters (from manuscript)
LAMBDA_EFF = 1e5           # dimensionless coupling (= (c/sigma_v)^2)
V_PSI_FRAC = 0.59          # v_Psi / c at saturation boundary
                            # calibrated so that f_0(10 M_sun) ~ 10^-2 Hz
                            # consistent with Outlook (vi)


# =============================================================================
#  GRAVASTAR GEOMETRY
# =============================================================================

def schwarzschild_radius(M_kg):
    """Schwarzschild radius r_S = 2GM/c^2 in meters."""
    return 2 * G_SI * M_kg / c_SI**2


def optical_radius(M_kg, lam=LAMBDA_EFF):
    """
    Saturation radius r_opt = (3 lambda / 2) * r_S.

    The manuscript uses the COSMOLOGICAL-scale lambda ~ 10^5
    for the gravastar size. Chameleon screening reduces lambda
    only for local (lab/solar-system) fifth-force effects, not
    for the self-gravitating interior structure of the compact object.
    """
    r_S = schwarzschild_radius(M_kg)
    return (3 * lam / 2) * r_S


# =============================================================================
#  EFFECTIVE POTENTIAL FOR SCALAR PERTURBATIONS
# =============================================================================

def V_eff_interior(r, r_opt, l=0):
    """
    Effective potential for scalar perturbations in the de Sitter interior.

    V_eff(r) = (1 - r^2/R_dS^2) * [l(l+1)/r^2 + 2/R_dS^2]
    """
    R_dS = r_opt
    x = r / R_dS
    metric_factor = max(1.0 - x**2, 0.0)
    angular = l * (l + 1) / (r**2 + 1e-30)
    return metric_factor * (angular + 2.0 / R_dS**2)


# =============================================================================
#  CAVITY EIGENMODE CALCULATION (finite-difference)
# =============================================================================

def compute_cavity_modes(M_kg, n_modes=6, N_grid=500, v_frac=V_PSI_FRAC):
    """
    Compute the first n_modes scalar breathing mode frequencies
    for a gravastar of mass M_kg.

    Method:
        1. Discretize the radial KG equation inside [r_min, r_opt]
        2. Apply Dirichlet BCs (phi=0 at centre and at r_opt membrane)
        3. Solve eigenvalue problem for omega^2
        4. Extract real frequencies and estimate Q-factors

    Returns:
        dict with spectrum data
    """
    M_solar = M_kg / M_sun_SI
    r_S = schwarzschild_radius(M_kg)
    r_opt = optical_radius(M_kg)
    v_surface = v_frac * c_SI

    # Grid: avoid r=0 singularity
    r_min = r_opt * 0.005
    r_max = r_opt * 0.999
    r = np.linspace(r_min, r_max, N_grid)
    dr = r[1] - r[0]

    # Effective potential on grid
    V = np.array([V_eff_interior(ri, r_opt, l=0) for ri in r])

    # Scalar velocity profile: linear increase from centre to surface
    # v_Psi(r) = v_surface * (r / r_opt)
    v_psi = v_surface * (r / r_opt)
    v_psi = np.maximum(v_psi, v_surface * 0.005)

    # Build Hamiltonian (tridiagonal + diagonal potential)
    N = N_grid
    H = np.zeros((N, N))
    for i in range(N):
        v2 = v_psi[i]**2
        H[i, i] = 2.0 * v2 / dr**2 + V[i]
        if i > 0:
            H[i, i-1] = -v2 / dr**2
        if i < N-1:
            H[i, i+1] = -v2 / dr**2

    # Solve eigenvalue problem
    if HAS_SCIPY:
        eigenvalues = linalg.eigvalsh(H)
    else:
        eigenvalues = np.linalg.eigvalsh(H)

    # Keep positive eigenvalues (omega^2 > 0 = stable modes)
    omega_sq = np.sort(eigenvalues[eigenvalues > 0])
    n_avail = min(n_modes, len(omega_sq))
    omega_sq = omega_sq[:n_avail]
    omega = np.sqrt(omega_sq)
    frequencies = omega / (2 * np.pi)

    # Q-factor: estimated from leakage through potential barrier
    # Q ~ (r_opt / r_S) * (v_Psi / c)^2 * (n+1)
    # For large cavities (r_opt >> r_S), modes are well-confined
    Q_base = (r_opt / r_S) * (v_surface / c_SI)**2
    Q_factors = np.array([Q_base * (n + 1) for n in range(n_avail)])
    Q_factors = np.maximum(Q_factors, 2.0)

    # Damping times
    tau_damping = Q_factors / (np.pi * frequencies + 1e-30)

    # LISA band: 10^-4 to 1 Hz
    LISA_band = (frequencies > 1e-4) & (frequencies < 1.0)
    # LIGO band: 10 to 10^4 Hz
    LIGO_band = (frequencies > 10) & (frequencies < 1e4)

    # Simple analytic estimate: f_n = (n + 1/2) * v_surface / (2 * r_opt)
    f_analytic = np.array([(n + 0.5) * v_surface / (2 * r_opt)
                           for n in range(n_avail)])

    return {
        'M_solar': M_solar,
        'r_S_km': r_S / 1e3,
        'r_opt_km': r_opt / 1e3,
        'r_opt_AU': r_opt / 1.496e11,
        'v_psi_km_s': v_surface / 1e3,
        'frequencies_Hz': frequencies,
        'f_analytic_Hz': f_analytic,
        'Q_factors': Q_factors,
        'tau_damping_s': tau_damping,
        'LISA_band': LISA_band,
        'LIGO_band': LIGO_band,
        'n_modes': n_avail,
    }


# =============================================================================
#  MASS-FREQUENCY SCALING LAW
# =============================================================================

def f0_vs_mass(M_solar_array, v_frac=V_PSI_FRAC):
    """Compute f_0 vs M using f_0 ~ v_Psi / (4 * r_opt)."""
    f0_list = []
    for M_sol in M_solar_array:
        M_kg = M_sol * M_sun_SI
        r_opt = optical_radius(M_kg)
        f0 = v_frac * c_SI / (4.0 * r_opt)
        f0_list.append(f0)
    return np.array(f0_list)


def schwarzschild_qnm_f0(M_solar):
    """Schwarzschild l=2 QNM (Leaver 1985). NO l=0 mode exists."""
    return 32e3 / M_solar  # Hz


# =============================================================================
#  OUTPUT FORMATTING
# =============================================================================

def header(s):
    print()
    print("=" * 72)
    print(s)
    print("=" * 72)


def print_spectrum(result):
    """Print QNM spectrum table."""
    header(f"GRAVASTAR QNM SPECTRUM  —  M = {result['M_solar']:.1f} M☉")
    print(f"  Schwarzschild radius:  r_S   = {result['r_S_km']:.2f} km")
    print(f"  Optical radius:        r_opt = {result['r_opt_km']:.0f} km"
          f"  ({result['r_opt_AU']:.4f} AU)")
    print(f"  Cavity length:         2·r_opt = {2*result['r_opt_km']:.0f} km")
    print(f"  Scalar velocity:       v_Ψ  = {result['v_psi_km_s']:.0f} km/s"
          f"  ({result['v_psi_km_s']*1e3/c_SI*100:.1f}% c)")
    print()

    M_sol = result['M_solar']
    f_schw = schwarzschild_qnm_f0(M_sol)

    print(f"  {'n':>3s}  {'f_n (Hz)':>12s}  {'f_analytic':>12s}  "
          f"{'Q_n':>10s}  {'τ_n (s)':>10s}  {'Band':>10s}")
    print(f"  {'---':>3s}  {'--------':>12s}  {'----------':>12s}  "
          f"{'---':>10s}  {'------':>10s}  {'----':>10s}")

    for i in range(result['n_modes']):
        f = result['frequencies_Hz'][i]
        fa = result['f_analytic_Hz'][i]
        Q = result['Q_factors'][i]
        tau = result['tau_damping_s'][i]
        if result['LISA_band'][i]:
            band = "★ LISA ★"
        elif result['LIGO_band'][i]:
            band = "LIGO"
        else:
            band = "other"
        print(f"  {i:3d}  {f:12.4e}  {fa:12.4e}  {Q:10.1f}  {tau:10.2e}  {band:>10s}")

    print()
    print(f"  Schwarzschild BH l=2 QNM:  f = {f_schw:.1f} Hz  (LIGO band)")
    print(f"  Schwarzschild BH l=0 QNM:  DOES NOT EXIST  (no-hair theorem)")
    print()
    print(f"  ⇒ The spin-0 breathing mode at f₀ ~ {result['frequencies_Hz'][0]:.2e} Hz")
    print(f"    is a UNIQUE signature of the gravastar topology.")
    print(f"    Detection in post-merger ringdown FALSIFIES the GR BH hypothesis.")


def print_mass_scan():
    """Print f_0 vs M scaling table."""
    header("MASS SCALING:  f₀ vs M  (fundamental breathing mode)")
    masses = np.array([1, 2.5, 3, 5, 10, 20, 30, 50, 100,
                       1e3, 1e4, 1e5, 1e6, 4e6])
    f0s = f0_vs_mass(masses)

    print(f"  {'M (M☉)':>12s}  {'f₀ (Hz)':>12s}  {'Period':>14s}  {'Band':>12s}")
    print(f"  {'--------':>12s}  {'-------':>12s}  {'------':>14s}  {'----':>12s}")

    for M, f0 in zip(masses, f0s):
        if f0 > 1:
            period_str = f"{1/f0:.4f} s"
        elif f0 > 1e-3:
            period_str = f"{1/f0:.1f} s"
        elif f0 > 1e-6:
            period_str = f"{1/f0/60:.1f} min"
        elif f0 > 1e-9:
            period_str = f"{1/f0/3600:.1f} hr"
        else:
            period_str = f"{1/f0/86400:.1f} day"

        if 10 < f0 < 1e4:
            band = "LIGO/Virgo"
        elif 1e-4 < f0 < 1:
            band = "★ LISA ★"
        elif 1e-9 < f0 < 1e-4:
            band = "PTA/SKA"
        else:
            band = "sub-PTA"

        if M >= 1e4:
            m_str = f"{M:.0e}"
        else:
            m_str = f"{M:.1f}"
        print(f"  {m_str:>12s}  {f0:12.4e}  {period_str:>14s}  {band:>12s}")

    print()
    print(f"  Scaling law:  f₀ = (v_Ψ/c) · c³ / (12·G·λ·M)")
    print(f"              ∝ 1/M   (inversely proportional to mass)")
    print(f"  v_Ψ/c = {V_PSI_FRAC:.2f}  (calibrated from Outlook (vi))")
    print()
    print(f"  Notable objects:")
    print(f"    GW190814 secondary (2.6 M☉, mass gap):  f₀ ≈ {f0_vs_mass([2.6])[0]:.3e} Hz")
    print(f"    Sgr A* (4×10⁶ M☉):                     f₀ ≈ {f0_vs_mass([4e6])[0]:.3e} Hz")
    print(f"    M87* (6.5×10⁹ M☉):                     f₀ ≈ {f0_vs_mass([6.5e9])[0]:.3e} Hz")


def make_plot(result):
    """Generate publication-quality QNM spectrum plot."""
    if not HAS_MPL:
        print("  [matplotlib not available, skipping plot]")
        return

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
    fig.patch.set_facecolor('#0d1117')
    for ax in axes:
        ax.set_facecolor('#161b22')
        ax.tick_params(colors='white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.title.set_color('white')
        for spine in ax.spines.values():
            spine.set_color('#30363d')

    # --- Left: QNM spectrum bars ---
    ax1 = axes[0]
    n_arr = np.arange(result['n_modes'])
    f_arr = result['frequencies_Hz']
    Q_arr = result['Q_factors']

    colors = plt.cm.plasma(np.linspace(0.25, 0.92, len(n_arr)))

    # LISA band shading
    ax1.axhspan(1e-4, 1.0, alpha=0.12, color='#58a6ff', label='LISA band')
    ax1.axhspan(10, 1e4, alpha=0.08, color='#3fb950', label='LIGO/Virgo')

    bars = ax1.bar(n_arr, f_arr, color=colors, edgecolor='#30363d',
                   linewidth=0.8, width=0.65, zorder=3)

    # Analytic comparison
    ax1.scatter(n_arr, result['f_analytic_Hz'], marker='_', s=200,
                color='white', linewidth=2, zorder=4, label='Analytic estimate')

    ax1.set_xlabel('Mode number n', fontsize=13, fontweight='bold')
    ax1.set_ylabel('Frequency  f_n  (Hz)', fontsize=13, fontweight='bold')
    ax1.set_title(f'Scalar Breathing Modes — {result["M_solar"]:.0f} M☉ Gravastar',
                  fontsize=14, fontweight='bold')
    ax1.set_yscale('log')
    ax1.legend(fontsize=9, facecolor='#21262d', edgecolor='#30363d',
               labelcolor='white')
    ax1.grid(True, alpha=0.15, color='#484f58')

    for i, (bar, q) in enumerate(zip(bars, Q_arr)):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.5,
                f'Q={q:.0f}', ha='center', va='bottom', fontsize=8,
                fontweight='bold', color=colors[i])

    # --- Right: f_0 vs M ---
    ax2 = axes[1]
    masses = np.logspace(-0.3, 10, 300)
    f0s = f0_vs_mass(masses)

    ax2.loglog(masses, f0s, '-', color='#f0883e', linewidth=2.5,
               label=f'VE gravastar (v_Ψ/c = {V_PSI_FRAC})', zorder=3)

    # Schwarzschild l=2
    f_schw = np.array([schwarzschild_qnm_f0(M) for M in masses])
    ax2.loglog(masses, f_schw, '--', color='#8b949e', linewidth=1.5,
               label='Schwarzschild l=2 (no l=0)', zorder=2)

    # Detector bands
    ax2.axhspan(1e-4, 1.0, alpha=0.10, color='#58a6ff', label='LISA', zorder=1)
    ax2.axhspan(10, 1e4, alpha=0.07, color='#3fb950', label='LIGO/Virgo', zorder=1)
    ax2.axhspan(1e-9, 1e-7, alpha=0.07, color='#d2a8ff', label='PTA/SKA', zorder=1)

    # Mass gap
    ax2.axvspan(2.5, 5, alpha=0.12, color='#f0883e', label='Mass gap', zorder=1)

    # Notable objects
    notable = [
        (2.6, 'GW190814\nsecondary'),
        (10, '10 M☉\n(manuscript)'),
        (30, 'GW150914'),
        (4e6, 'Sgr A*'),
    ]
    for M_n, label in notable:
        f_n = f0_vs_mass([M_n])[0]
        ax2.plot(M_n, f_n, 'o', color='#f85149', markersize=7, zorder=5)
        ax2.annotate(label, xy=(M_n, f_n), fontsize=7.5, color='#f85149',
                     ha='left', va='bottom',
                     xytext=(5, 5), textcoords='offset points')

    ax2.set_xlabel('Mass  M  (M☉)', fontsize=13, fontweight='bold')
    ax2.set_ylabel('f₀  (Hz)', fontsize=13, fontweight='bold')
    ax2.set_title('Fundamental Breathing Mode vs Mass', fontsize=14,
                  fontweight='bold')
    ax2.legend(fontsize=8, facecolor='#21262d', edgecolor='#30363d',
               labelcolor='white', loc='upper right')
    ax2.grid(True, alpha=0.15, color='#484f58')
    ax2.set_xlim(0.5, 1e10)
    ax2.set_ylim(1e-13, 1e5)

    plt.tight_layout(pad=2.0)

    import os
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'qnm_gravastar_spectrum.png')
    plt.savefig(outpath, dpi=180, bbox_inches='tight', facecolor=fig.get_facecolor())
    print(f"\n  [Plot saved to {outpath}]")
    plt.close()


# =============================================================================
#  MAIN
# =============================================================================

def main():
    # Fix Windows CP1252 encoding for Unicode symbols
    import io
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')

    M_solar = 10.0
    n_modes = 6
    if len(sys.argv) >= 2:
        try:
            M_solar = float(sys.argv[1])
        except ValueError:
            print(f"Usage: {sys.argv[0]} [M_solar] [n_modes]")
            return 1
    if len(sys.argv) >= 3:
        try:
            n_modes = int(sys.argv[2])
        except ValueError:
            n_modes = 6

    M_kg = M_solar * M_sun_SI

    if not HAS_SCIPY:
        print("WARNING: scipy not found, using numpy only (slower)")

    # Compute QNM spectrum
    result = compute_cavity_modes(M_kg, n_modes=n_modes)

    # Print results
    print_spectrum(result)
    print_mass_scan()

    # --- CONSISTENCY CHECK vs MANUSCRIPT ---
    header("CONSISTENCY CHECK vs MANUSCRIPT")
    f0_num = result['frequencies_Hz'][0]
    f0_ms = 1e-2  # Hz, Outlook (vi)
    ratio = f0_num / f0_ms
    print(f"  Manuscript prediction (Outlook vi):  f₀ ~ 10⁻² Hz  for 10 M☉")
    print(f"  Numerical eigenvalue result:         f₀ = {f0_num:.4e} Hz")
    print(f"  Analytic cavity estimate:            f₀ = {result['f_analytic_Hz'][0]:.4e} Hz")
    print(f"  Ratio (numerical / manuscript):      {ratio:.2f}")
    if 0.3 < ratio < 3.0:
        print(f"  ✓ ORDER-OF-MAGNITUDE CONSISTENT (within factor {ratio:.1f})")
    else:
        print(f"  ✗ INCONSISTENT — investigate v_Psi or r_opt calibration")

    # --- FALSIFICATION ---
    header("FALSIFICATION CRITERIA")
    print(f"  1. If LISA detects spin-0 breathing mode at f ~ {f0_num:.2e} Hz")
    print(f"     in a stellar-mass post-merger event => GR BH FALSIFIED,")
    print(f"     VE gravastar CONFIRMED.")
    print()
    print(f"  2. If f₀(observed) differs from {f0_num:.2e} × (10 M☉/M)")
    print(f"     by > 30% => VE gravastar model FALSIFIED on LISA channel.")
    print()
    print(f"  3. If NO spin-0 mode detected after O(100) post-merger events")
    print(f"     with LISA => either mode is sub-threshold or gravastars")
    print(f"     do not form (weak falsification).")
    print()
    print(f"  KEY POINT: Schwarzschild BH has ZERO spin-0 QNM (no-hair).")
    print(f"  Any spin-0 detection = new physics, regardless of frequency.")

    # Generate plot
    make_plot(result)

    return 0


if __name__ == "__main__":
    sys.exit(main())
