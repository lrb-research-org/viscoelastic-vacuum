"""
seebeck_bode_susceptibility.py
==============================
Reproduces the Bode-plot of the Seebeck susceptibility chi_Seebeck(omega)
in the two laboratory regimes of the preprint:

    chi_Seebeck(omega) = beta_e^2 / (m_eff^2 - omega^2 - i*Gamma*omega)

Eq. (eq:chi_seebeck) of Bouille (2026), Section 7. Gravitational Seebeck Effect.

Two laboratory regimes (Section 7.1):
    UHV  : rho_bg ~ 5e-3 kg/m^3 -> f_res ~ 1   GHz, peak Delta V ~   1 nV
    He-4 : rho_bg = 145  kg/m^3 -> f_res ~ 3   THz, peak Delta V ~ 180 nV

Falsification test F3 (Section 9.5, Pre-registration):
    Null Seebeck signal at 5 sigma in the He-4 cryogenic configuration
    (Delta V < 36 nV) would falsify the electro-scalar coupling channel.

Usage
-----
    python seebeck_bode_susceptibility.py
    -> writes seebeck_bode.png and prints the resonance peaks.
"""

import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------
# 1. Parameters from Section 7 (eq:f_res)
# ------------------------------------------------------------------------
# Resonance frequencies are set by the Chameleon mass scale m_eff(rho_bg)
# with the runaway potential V(phi) = Lambda^5 / phi:
#   m_eff(rho)    propto rho^{3/4}      (Eq. eq:meff)
#   omega_res     = m_eff c^2 / hbar    (Eq. eq:f_res)

# UHV regime  (Section 7.3, paragraph "Detector-side unscreening")
F_RES_UHV   = 1.0e9       # 1 GHz
DV_PEAK_UHV = 1.0         # nV at resonance, Section 7.4

# Liquid He-4 regime (Section 7.3, paragraph "Liquid helium configuration")
F_RES_HE4   = 3.0e12      # 3 THz
DV_PEAK_HE4 = 180.0       # nV at resonance, Section 7.4

# Visualisation Q (Section 7, Fig. 16 caption: physical Q ~ 1e3 reduced
# to ~30 here for plot legibility -- does NOT affect the peak amplitude
# which uses the canonical Q = omega_res / Gamma normalisation).
Q_VIS = 30.0

OMEGA_UHV   = 2.0 * np.pi * F_RES_UHV
OMEGA_HE4   = 2.0 * np.pi * F_RES_HE4
GAMMA_UHV   = OMEGA_UHV / Q_VIS
GAMMA_HE4   = OMEGA_HE4 / Q_VIS

# SQUID quantum noise floor (Section 7.4, Table tab:noise_budget)
SQUID_FLOOR_NV = 0.03     # nV / sqrt(Hz) effective floor in 1 mHz BW

# ------------------------------------------------------------------------
# 2. Lorentzian magnitude and phase
# ------------------------------------------------------------------------
def lorentz_amp(omega, omega_res, gamma, peak):
    """ |chi(omega)| Lorentzian normalised so that |chi(omega_res)| == peak. """
    denom    = (omega_res**2 - omega**2)**2 + (gamma * omega)**2
    raw      = 1.0 / np.sqrt(denom)
    raw_peak = 1.0 / (gamma * omega_res)        # value at omega = omega_res
    return peak * raw / raw_peak

def lorentz_phase_deg(omega, omega_res, gamma):
    """ Standard Lorentzian phase 0 -> -90 -> -180 deg. """
    return np.degrees(np.arctan2(-gamma * omega, omega_res**2 - omega**2))

# ------------------------------------------------------------------------
# 3. Frequency sweep
# ------------------------------------------------------------------------
f       = np.logspace(7, 14, 2000)        # 10 MHz -> 100 THz
omega   = 2.0 * np.pi * f

dv_uhv   = lorentz_amp(omega, OMEGA_UHV, GAMMA_UHV, DV_PEAK_UHV)
dv_he4   = lorentz_amp(omega, OMEGA_HE4, GAMMA_HE4, DV_PEAK_HE4)
ph_uhv   = lorentz_phase_deg(omega, OMEGA_UHV, GAMMA_UHV)
ph_he4   = lorentz_phase_deg(omega, OMEGA_HE4, GAMMA_HE4)

# ------------------------------------------------------------------------
# 4. SNR sanity checks vs the noise floor
# ------------------------------------------------------------------------
snr_uhv = DV_PEAK_UHV / SQUID_FLOOR_NV
snr_he4 = DV_PEAK_HE4 / SQUID_FLOOR_NV

print("=" * 64)
print(" Seebeck susceptibility -- Bode-plot reproduction")
print("=" * 64)
print(f"  UHV regime  (rho ~ 5e-3 kg/m^3):")
print(f"    f_res    = {F_RES_UHV:.2e} Hz")
print(f"    Delta V  = {DV_PEAK_UHV:.1f} nV at resonance")
print(f"    SNR      ~ {snr_uhv:.0f}                   (5 sigma threshold met)")
print()
print(f"  He-4 regime (rho = 145 kg/m^3):")
print(f"    f_res    = {F_RES_HE4:.2e} Hz")
print(f"    Delta V  = {DV_PEAK_HE4:.1f} nV at resonance")
print(f"    SNR      ~ {snr_he4:.0f}                  (>> 5 sigma)")
print()
print("  Falsification test F3 (Section 9.5):")
print("    Null Seebeck signal at 5 sigma (Delta V < 36 nV) in the He-4")
print("    configuration would falsify the electro-scalar coupling channel.")
print("=" * 64)

# ------------------------------------------------------------------------
# 5. Bode plot (magnitude + phase)
# ------------------------------------------------------------------------
fig, (ax_amp, ax_ph) = plt.subplots(
    2, 1, figsize=(8.5, 6.0), sharex=True,
    gridspec_kw={'height_ratios': [2.2, 1.0]})

# magnitude
ax_amp.loglog(f, dv_uhv, '-', color='#B2182B', linewidth=2.0,
              label=r'UHV  ($\rho \sim 5\!\times\!10^{-3}$ kg/m$^3$)'
                    r' : $\Delta V \sim 1$ nV @ 1 GHz')
ax_amp.loglog(f, dv_he4, '-', color='#2166AC', linewidth=2.0,
              label=r'He-4 ($\rho = 145$ kg/m$^3$)'
                    r' : $\Delta V \sim 180$ nV @ 3 THz')
ax_amp.axhline(SQUID_FLOOR_NV, color='#444444', linestyle='--',
               linewidth=1.3, alpha=0.8,
               label=r'SQUID noise floor ($\sim 0.03$ nV)')
ax_amp.plot(F_RES_UHV, DV_PEAK_UHV, 'D', color='#B2182B', ms=9,
            markeredgecolor='black', mew=0.6)
ax_amp.plot(F_RES_HE4, DV_PEAK_HE4, 's', color='#2166AC', ms=9,
            markeredgecolor='black', mew=0.6)
ax_amp.set_ylabel(r'Seebeck signal $|\Delta V(\omega)|$ (nV)')
ax_amp.set_xlim(1e7, 1e14)
ax_amp.set_ylim(1e-3, 1e3)
ax_amp.legend(loc='lower right', fontsize=8.5, framealpha=0.95)
ax_amp.grid(True, which='major', linestyle='--', alpha=0.3)

# phase
ax_ph.semilogx(f, ph_uhv, '-', color='#B2182B', linewidth=1.5)
ax_ph.semilogx(f, ph_he4, '-', color='#2166AC', linewidth=1.5)
ax_ph.axhline(-90, color='#888888', linestyle=':', linewidth=0.9)
ax_ph.set_xlabel(r'Driving frequency $f$ (Hz)')
ax_ph.set_ylabel(r'Phase (deg)')
ax_ph.set_ylim(-185, 5)
ax_ph.set_yticks([-180, -135, -90, -45, 0])
ax_ph.grid(True, which='major', linestyle='--', alpha=0.3)

plt.tight_layout()
out = 'seebeck_bode.png'
plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved plot: {out}")
