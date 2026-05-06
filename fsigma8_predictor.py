"""
fsigma8_predictor.py
====================
Numerical integration of the modified Jeans equation in the viscoelastic vacuum
framework (Bouille 2026, Appendix A.2).

PHYSICS:
    Topologically Einstein-de Sitter universe (Omega_m^top = 1, no Lambda).
    Frozen Chameleon condensate provides Omega_phi^frozen ~ 0.27 as effective
    cold matter at z >> 0.

    Modified Jeans equation (manuscript Eq. (jeans_modified)):
        ddot{delta} + 2*H*dot{delta} + (c_s^2 k^2/a^2 - (3/2) H^2 F(k,a)) delta = 0

    Two scale-dependent suppression mechanisms vs pure EdS (where f = 1):

    (1) Viscoelastic acoustic pressure:
        c_s^2 k^2 / (aH)^2  acts as elastic restoring force above
        k_crit^2 = (3/2) (aH/c_s)^2.  For k > k_crit, perturbations
        oscillate rather than grow.

    (2) Chameleon low-pass screening:
        F(k,a) = 1 + (2*lambda_eff/3) * S(k,a)
        S(k,a) = (k/aH)^2 / [(k/aH)^2 + m_eff^2(a)/H^2]
        At the cosmological background, m_eff^2/H^2 is large, so screening is
        active for relevant k.

    For k ~ 0.1 h/Mpc (the sigma_8 scale), only the acoustic term is significant;
    the fifth force is screened at cosmological scales by the high lambda_eff
    needed for cluster-scale dynamics.

REFERENCE OBSERVATIONS (RSD compilation):
    z = 0.067:  fsigma8 = 0.423 +/- 0.055  (6dFGS)
    z = 0.18 :  fsigma8 = 0.360 +/- 0.090  (2dFGRS)
    z = 0.38 :  fsigma8 = 0.497 +/- 0.045  (BOSS DR12)
    z = 0.51 :  fsigma8 = 0.458 +/- 0.038  (BOSS DR12)
    z = 0.61 :  fsigma8 = 0.436 +/- 0.034  (BOSS DR12)
    z = 0.85 :  fsigma8 = 0.450 +/- 0.080  (eBOSS LRG)
    z = 1.48 :  fsigma8 = 0.462 +/- 0.052  (eBOSS QSO)

    Pure EdS prediction (f=1): fsigma8(z) = sigma_8(0) / (1+z)
    At z = 0.5: pure EdS gives 0.54 (overshoot by ~18% vs observed 0.45).

USAGE:
    python fsigma8_predictor.py
    python fsigma8_predictor.py --plot
"""
import numpy as np
from scipy.integrate import solve_ivp
import sys


# ============================================================
# Cosmological parameters
# ============================================================
H0 = 73.6                 # km/s/Mpc, SH0ES geometric anchor
sigma8_today = 0.81       # observed at z = 0
c_kms = 2.998e5           # speed of light, km/s

# Viscoelastic parameters (calibrated to give fsigma8 ~ 0.45 at z = 0.5)
C_S_KMS = 850.0           # effective sound speed of viscoelastic vacuum, km/s
                          # (tuned to fsigma8 ~ 0.45 at z = 0.5; precise value in companion MCMC)
LAMBDA_EFF_COSMO = 0.0    # cosmological-scale lambda_eff after screening at recombination
                          # (the cosmological sound speed dominates linear suppression;
                          #  lambda enhancements activate only inside virialized halos)


def jeans_rhs(N, y, k_hMpc, c_s_kms=C_S_KMS, lambda_eff=LAMBDA_EFF_COSMO):
    """
    Modified Jeans equation in pure EdS background, time variable N = ln(a).

    In EdS: H^2 = H0^2 / a^3, dlnH/dN = -3/2.

    Rewriting time derivatives via d/dt = H * d/dN:
        ddot delta + 2 H dot delta = H^2 [d^2 delta/dN^2 - (3/2)(d delta/dN) + 2(d delta/dN)]
                                   = H^2 [d^2 delta/dN^2 + (1/2)(d delta/dN)]

    Source: c_s^2 k^2 / a^2 * delta - (3/2) H^2 F * delta
          = (3/2) H^2 [(c_s^2 k^2/(aH)^2) / (3/2 a^-1) - F] * ???

    Cleanly: divide whole equation by H^2:
        d^2 delta/dN^2 + (1/2)(d delta/dN) + [c_s^2 k^2 / (a H)^2 - (3/2) F] delta = 0
        d^2 delta/dN^2 + (1/2)(d delta/dN) - [(3/2) F - c_s^2 k^2 / (a H)^2] delta = 0

    In EdS, (aH)^2 = H0^2 / a, so c_s^2 k^2 / (aH)^2 = (c_s k / H0)^2 * a.
    """
    delta, delta_prime = y
    a = np.exp(N)

    # Dimensionless k * c / H0 (k in 1/Mpc, c in km/s, H0 in km/s/Mpc)
    k_dim = k_hMpc * c_kms / H0           # ~ 407 for k = 0.1 /Mpc, H0 = 73.6
    c_s_dim = c_s_kms / c_kms             # c_s in units of c

    # Acoustic term coefficient: (c_s k / aH)^2 = (c_s k / H0)^2 * a
    sigma_acc_sq = (c_s_dim * k_dim)**2 * a

    # Chameleon screening factor (negligible at cosmological scales given screening at recombination)
    if lambda_eff > 0:
        m_eff_sq_over_H_sq = 3.0 * lambda_eff
        k_over_aH_sq = k_dim**2 * a
        S_screen = k_over_aH_sq / (k_over_aH_sq + m_eff_sq_over_H_sq)
        F_eff = 1.0 + (2.0 * lambda_eff / 3.0) * S_screen
    else:
        F_eff = 1.0

    # ODE in N-time
    d2_delta = -0.5 * delta_prime + (1.5 * F_eff - sigma_acc_sq) * delta
    return [delta_prime, d2_delta]


def integrate_growth(k_hMpc, c_s_kms=C_S_KMS, lambda_eff=LAMBDA_EFF_COSMO,
                     a_init=1e-3):
    """Integrate from a_init (matter-dominated) to a = 1 (today)."""
    N_init = np.log(a_init)
    # Matter-dominated initial conditions: delta ~ a, dN(delta) = delta
    y0 = [a_init, a_init]

    sol = solve_ivp(
        lambda N, y: jeans_rhs(N, y, k_hMpc, c_s_kms, lambda_eff),
        [N_init, 0.0], y0,
        dense_output=True, rtol=1e-9, atol=1e-12, method='RK45'
    )

    # Sample at relevant redshifts
    a_arr = np.linspace(0.30, 1.0, 200)   # z in [0, 2.33]
    N_arr = np.log(a_arr)
    y_arr = sol.sol(N_arr)
    delta_arr = y_arr[0]
    dprime_arr = y_arr[1]
    f_arr = dprime_arr / delta_arr

    delta_today = sol.sol(0.0)[0]
    sigma8_arr = sigma8_today * delta_arr / delta_today
    fsigma8_arr = f_arr * sigma8_arr

    return {
        'a': a_arr,
        'z': 1.0/a_arr - 1.0,
        'delta': delta_arr,
        'f': f_arr,
        'sigma8': sigma8_arr,
        'fsigma8': fsigma8_arr,
    }


# ============================================================
# RSD observational data
# ============================================================
RSD_DATA = [
    (0.067, 0.423, 0.055, '6dFGS'),
    (0.18,  0.360, 0.090, '2dFGRS'),
    (0.38,  0.497, 0.045, 'BOSS DR12'),
    (0.51,  0.458, 0.038, 'BOSS DR12'),
    (0.61,  0.436, 0.034, 'BOSS DR12'),
    (0.85,  0.450, 0.080, 'eBOSS LRG'),
    (1.48,  0.462, 0.052, 'eBOSS QSO'),
]


def report(k_hMpc=0.1):
    print("=" * 72)
    print(f"VISCOELASTIC fsigma8(z) prediction at k = {k_hMpc} h/Mpc")
    print("=" * 72)
    print(f"  H0 = {H0} km/s/Mpc, sigma_8(0) = {sigma8_today}")
    print(f"  Topological EdS (Omega_m^top = 1, frozen Chameleon as effective CDM)")
    print(f"  c_s (vacuum sound speed) = {C_S_KMS} km/s")
    print(f"  lambda_eff (cosmological) = {LAMBDA_EFF_COSMO} (screened at recombination)")
    print()

    result = integrate_growth(k_hMpc)
    z_arr = result['z']
    fs8_arr = result['fsigma8']
    f_arr = result['f']
    sigma8_arr = result['sigma8']

    print(f"{'z':>6} {'fsigma8 NDG':>12} {'pure EdS':>10} {'f':>7} {'sigma8(z)':>11}")
    print("-" * 72)
    for z_target in [0.067, 0.18, 0.38, 0.51, 0.61, 0.85, 1.0, 1.48, 2.0]:
        idx = np.argmin(np.abs(z_arr - z_target))
        fs8_eds = sigma8_today / (1.0 + z_arr[idx])  # pure EdS reference
        print(f"{z_arr[idx]:6.3f} {fs8_arr[idx]:12.4f} {fs8_eds:10.4f} "
              f"{f_arr[idx]:7.3f} {sigma8_arr[idx]:11.4f}")
    print()

    print(f"{'Survey':>14} {'z':>6} {'fsigma8 obs':>14} {'NDG':>9} {'dev':>8}")
    print("-" * 72)
    for z, fs8, err, survey in RSD_DATA:
        idx = np.argmin(np.abs(z_arr - z))
        delta = fs8_arr[idx] - fs8
        dev_sigma = delta / err
        print(f"{survey:>14} {z:6.3f}  {fs8:.3f} +/- {err:.3f}   "
              f"{fs8_arr[idx]:6.3f}  {dev_sigma:+5.1f}s")
    print()

    # Aggregate chi^2
    chi2 = sum(((fs8_arr[np.argmin(np.abs(z_arr - z))] - fs8) / err)**2
               for z, fs8, err, _ in RSD_DATA)
    dof = len(RSD_DATA) - 1   # one effective parameter (c_s)
    print(f"Reduced chi^2 vs RSD compilation: {chi2/dof:.2f}  ({len(RSD_DATA)} pts)")
    print()


def make_plot(filename='fsigma8_predictor.png'):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping plot")
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    # NDG curves at multiple k
    for k_hMpc, color, label_suffix in [
        (0.05, 'C0', '0.05 h/Mpc (large-scale)'),
        (0.10, 'C1', '0.10 h/Mpc (sigma_8 scale)'),
        (0.20, 'C2', '0.20 h/Mpc (small-scale)'),
    ]:
        result = integrate_growth(k_hMpc)
        ax.plot(result['z'], result['fsigma8'], '-', color=color, lw=2,
                label=f'NDG VE  k = {label_suffix}')

    # Pure EdS reference: f = 1, sigma_8 ~ a
    z_ref = np.linspace(0, 2, 100)
    a_ref = 1 / (1 + z_ref)
    fs8_eds = sigma8_today * a_ref
    ax.plot(z_ref, fs8_eds, '--', color='gray', lw=1.5,
            label=r'Pure EdS $\Omega_m=1$ (no suppression): $f\sigma_8 = \sigma_{8,0} a$')

    # LCDM Planck approximation
    Omega_m = 0.315
    f_lcdm = (Omega_m * (1+z_ref)**3 / (Omega_m*(1+z_ref)**3 + (1-Omega_m)))**0.55
    sigma8_lcdm = 0.81 * a_ref / np.sqrt(Omega_m + (1-Omega_m)*a_ref**3)**0.55  # crude
    sigma8_lcdm = 0.81 / (1 + z_ref)**0.4   # phenomenological
    ax.plot(z_ref, f_lcdm * sigma8_lcdm, ':', color='C3', lw=1.5,
            label=r'$\Lambda$CDM Planck approx.')

    # RSD data
    z_obs = [d[0] for d in RSD_DATA]
    fs8_obs = [d[1] for d in RSD_DATA]
    err_obs = [d[2] for d in RSD_DATA]
    ax.errorbar(z_obs, fs8_obs, yerr=err_obs, fmt='ko', ms=6, capsize=3,
                label='RSD compilation (6dFGS, 2dFGRS, BOSS, eBOSS)', zorder=5)

    ax.set_xlabel('Redshift $z$', fontsize=12)
    ax.set_ylabel(r'$f\sigma_8(z)$', fontsize=12)
    ax.set_xlim(0, 2.0)
    ax.set_ylim(0.30, 0.65)
    ax.set_title(r'Linear growth $f\sigma_8(z)$ in viscoelastic vacuum framework',
                 fontsize=12)
    ax.legend(loc='lower right', fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(filename, dpi=150, bbox_inches='tight')
    import matplotlib.pyplot as plt
    plt.close(fig)
    print(f"Figure saved: {filename}")


def main():
    report(k_hMpc=0.1)
    if '--plot' in sys.argv:
        make_plot()
    return 0


if __name__ == "__main__":
    sys.exit(main())
