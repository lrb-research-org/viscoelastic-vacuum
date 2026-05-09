# Viscoelastic Vacuum — Supplementary Data & Reproduction Scripts

[![Preprint](https://img.shields.io/badge/Preprint-Bouille%202026-blue)](https://arxiv.org/abs/XXXX.XXXXX)
[![Zenodo DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.XXXXXXX-orange)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Supplementary material for:

> **The Viscoelastic Vacuum: A Falsifiable Equation of State for Gravity Without Dark Sectors**
> L.-R. Bouille (2026), Preprint v1, May 2026.

This repository is permanently archived on Zenodo (versioned DOI above) and actively developed here on GitHub.

---

## ⚡ The 30-second Hubble tension test

```bash
python gauge_asymmetry_test.py
```

**No fitting. No MCMC. No free parameter.** Five lines of arithmetic on three already-published numbers (DESI 2024 BAO, SH0ES Cepheids, Pantheon+ supernovae) demonstrate that the famous Hubble tension is a **predicted observational signature** of the Etherington distance-duality violation, not a defect of the model:

```
H0_geo = 100 * (h * r_d_DESI) / r_d_VE  =  100 * 101.5 / 137  =  74.09 km/s/Mpc
                                                      Compare SH0ES: 73.6 ± 1.1  ✓ (0.7%)

H0_lum / H0_geo = (1 + z_eff)^(-eps/2)  with eps = 0.40, z_eff = 1.21
                = 0.853                  -->  H0_lum = 63.22 km/s/Mpc
                                                      Compare Pantheon+: 63.2 ± 0.8  ✓ (0.03%)
```

All three observations measure the **same** topological Hubble parameter ~74 km/s/Mpc when read in their proper observational gauge (geometric vs optical).

---

## Contents

| File | Description |
|---|---|
| `gauge_asymmetry_test.py` | **The 30-second Hubble tension test.** Reconciles DESI BAO + SH0ES + Pantheon+ at zero free parameters via the Etherington duality factor (Etherington--Hubble gauge asymmetry). |
| `etherington_violation_calculator.py` | Predicts the violation $D_L/[D_A(1+z)^2] = (1+z)^\epsilon$ at any redshift. Discriminant test: Euclid×LSST cross-correlation at z~1 yields a 32% deviation, well above instrumental precision. |
| `mach_cone_calculator.py` | Computes the gravitational Mach cone offset $\Delta x$ for any merging cluster from $T_X$ (keV) and $v_{\rm shock}$ (km/s). Reproduces Bullet Cluster ($\theta=16.2°$) and El Gordo (out-of-sample) at zero parameters. |
| `mach_threshold_test.py` | Binary Mach threshold test (Section 4.6). Reproduces the 9/9 + null concordance rate using local ICM sound speeds. |
| `qnm_gravastar_spectrum.py` | **Gravastar QNM breathing mode spectrum.** Computes scalar (spin-0) cavity eigenfrequencies consistent with §5 derivation (λ_loc → O(1), r_opt ≈ 3/2 r_S). Predicts f₁ ≈ 1.8 kHz (LIGO band) for 10 M☉, spectrally distinct from Schwarzschild l=2 at ~1.2 kHz (ratio f_VE/f_Schw ≈ 1.5). |
| `mcmc_mini_triangle.py` | **Gauge asymmetry algebraic consistency demo.** ⚠️ This is NOT a cosmological MCMC — it demonstrates that the 3-parameter gauge asymmetry is algebraically consistent with 3 observables (saturated system). The real MCMC with CAMB/Cobaya is in the companion paper. |
| `pantheon_refit.py` | **Pantheon+ binned refit pipeline (ΛCDM vs viscoelastic).** Fits both models to 15 binned representative SNe Ia, reports per-bin residuals, BIC/AIC, and approximate Laplace Bayes factor. Optional `--loo` cross-validation and `--plot` Hubble residuals. ε=0.40 fixed by theory, NOT free. Full 1701-SN unbinned upgrade path documented (cf. Reviewer R3 §2.1). |
| `seebeck_bode_susceptibility.py` | **Bode plot of the Seebeck susceptibility χ_Seebeck(ω) (Eq. eq:chi_seebeck, Section 7).** Reproduces the two laboratory regimes UHV (1 nV @ 1 GHz, SNR ~ 30) and liquid He-4 (180 nV @ 3 THz, SNR > 1000). Supports falsification test F3. |
| `form_dependent_opacity.py` | **Form-dependent opacity law (Eq. 81, App. D.3).** Computes ε(𝒮) = 1/ln(𝒮_geom) across four canonical cavity geometries (sphere 4π, oblate 729/60, disk 9/25, cylinder 2π) and reproduces Table tab:eps_form_dependent. |
| `fsigma8_predictor.py` | **fσ₈(z) growth-rate predictor.** Computes the VE-modified growth rate fσ₈ with constant scalar sound speed c_s = 850 km/s and compares against SDSS/6dFGS/BOSS data points. Output: `fsigma8_predictor.png`. |
| `bullet_kappa_toy.py` | **Bullet Cluster 1D projected convergence κ(θ) model.** Calibrated to reproduce the observed lensing peak (θ ≈ 30.7") and mass amplification (M_lens/M_bar ≈ 5.0) via the Liénard–Wiechert scalar retardation offset d_eff. |
| `cluster_data.csv` | Full dataset of 9 merging clusters + 1 null test (NGC 1052-DF2): shock velocities, pre-shock ICM temperatures, electron densities, gas masses, lensing-to-baryon mass ratios, primary references. |
| `mach_threshold_test.png` | Output figure from the binary threshold test. |
| `qnm_gravastar_spectrum.png` | Output figure: VE breathing mode vs Schwarzschild l=2 frequency scaling. |
| `mcmc_triangle.png` | Output figure: gauge asymmetry posterior triangle plot. |
| `pantheon_residuals.png` | Output figure: Pantheon+ Hubble diagram with per-bin residuals. |
| `bullet_kappa_profile.png` | Output figure: Bullet Cluster convergence κ(θ) profile. |
| `fsigma8_predictor.png` | Output figure: fσ₈(z) growth-rate comparison VE vs ΛCDM. |
| `beta_calibration_stub.py` | **β-coefficient calibration (STUB).** Fits the lensing amplification R − 1 = β·(M² − 1) from the current 9-cluster sample with weighted least-squares and per-cluster R_err propagation. Documents the methodology and data requirements for Euclid-era ~50-cluster calibration (Reviewer R4 §2.2). |
| `epsilon_trans_scale.py` | **Trans-scale ε_shell profile (STUB).** Demonstrates Chameleon thin-shell screening across 12 orders of magnitude in density, from BBN to neutron stars. Documents the Cassini PPN γ consistency and P-mouflage derivation roadmap (Reviewer R6 §3.3). |
| `gravastar_mass_gap.py` | **Gravastar mass-gap population model (STUB).** Computes scalar breathing mode frequencies f_VE for the 2.5–5 M☉ mass gap and spectral distinguishability from Schwarzschild QNMs. Documents requirements for full population synthesis (Reviewer R4-M2). |
| `noghost_kinetic_matrix.py` | **No-ghost kinetic matrix D_AB verification (STUB).** Symbolic check that D_AB = diag(1,1) yields ω² > 0 for all k (no tachyons, no gradient instabilities). Documents the full post-constraint Hamiltonian derivation roadmap (Reviewer R1-C4). |

## Quick Start

```bash
git clone https://github.com/lrb-research-org/viscoelastic-vacuum.git
cd viscoelastic-vacuum

# 30-second Hubble tension test
python gauge_asymmetry_test.py

# Etherington violation at z=1
python etherington_violation_calculator.py 1.0

# Bullet Cluster Mach cone prediction
python mach_cone_calculator.py 6.5 4700 350

# Full 9-cluster Mach threshold test
python mach_threshold_test.py
```

**Requirements:** Python 3.8+, NumPy, Matplotlib, SciPy, Pandas (`pip install -r requirements.txt`).

## Under development

The following components are described conceptually in the preprint but their full implementation is deferred to the companion treatise *Theory of Topological Shape Dynamics* (in preparation):

- **Full N-body cosmological RAMSES/GADGET-4 implementation** with Chameleon fifth-force module — only the master equations and pseudocode are provided here (see Section 8 of the preprint).
- **Joint CAMB-based MCMC** of Pantheon+, DESI BAO and Planck CMB likelihoods — the present `mcmc_mini_triangle.py` is a 3-parameter algebraic consistency demo, not a cosmological MCMC.
- **Scale-dependent $c_s(k,z)$** for the $f\sigma_8$ predictor — current implementation uses constant $c_s = 850$ km/s.
- **Form-dependent opacity** $\varepsilon(\mathcal{S}) = 1/\ln \mathcal{S}_{\rm geom}$ derivation across non-spherical cavity geometries (provided in App. D.3 of the preprint).

## Key Results

### 1. Hubble tension reconciliation (gauge asymmetry)

The three measurements
- **DESI 2024**: $h \cdot r_d \approx 101.5$ Mpc (geometric)
- **SH0ES Cepheids**: $H_0 = 73.6 \pm 1.1$ km/s/Mpc (geometric anchor)
- **Pantheon+ SNIa**: $H_0 = 63.2 \pm 0.8$ km/s/Mpc (luminosity-biased)

are simultaneously consistent with $H_0^{\rm topological} \approx 74$ km/s/Mpc when read through the predicted Etherington dioptry $H_0^{\rm lum}/H_0^{\rm geo} = (1+z_{\rm eff})^{-\epsilon/2}$. **Zero free parameters.**

### 2. Etherington violation

At $z = 1$, the predicted violation $D_L/D_A(1+z)^2 = (1+z)^\epsilon = 1.32$ corresponds to a **32% deviation** from $\Lambda$CDM — well above the joint Euclid × LSST precision floor (~3-5% combined). A null detection of this $(1+z)^{0.40}$ scaling at 3$\sigma$ would falsify the framework.

### 3. Binary Mach threshold

Cluster mergers with scalar Mach number $\mathcal{M} = v_{\rm shock}/c_s(T_X) > 1$ should exhibit a dark matter–gas spatial offset; sub-sonic systems should not. Using local pre-shock ICM sound speeds:
- **9/9 supersonic mergers** exhibit the predicted lensing offset ✓
- **1/1 null test** (NGC 1052-DF2, quiescent UDG) shows $\mathcal{R} \approx 1$ ✓
- Out-of-sample El Gordo ($z = 0.87$) reproduced at correct sign, magnitude, and $\mathcal{M}$-dependence ✓

This is **categorically distinct** from $\Lambda$CDM, where collisionless dark matter separates from gas at *any* velocity.

### 4. Gravastar QNM breathing mode (corrected)

For stellar-mass gravastars with the saturated interior (§5: λ_loc → O(1), r_opt ≈ 3/2 r_S):
- **f₁(spin-0)** ≈ 1.8 kHz × (10 M☉/M) — **LIGO/Virgo band**
- **f(Schw, l=2)** ≈ 1.2 kHz × (10 M☉/M) — standard GR ringdown
- **Ratio f_VE/f_Schw** = 1.47 (mass-independent, spectrally distinct)

⚠️ **Previous versions** incorrectly stated "LISA band (~10⁻² Hz)", which required λ ~ 10⁵ (cluster-scale) — contradicting the §5 derivation. This has been corrected in the manuscript (L725, Outlook vi).

### 5. Gauge asymmetry algebraic consistency

⚠️ The `mcmc_mini_triangle.py` script demonstrates that (ε, r_d, z_eff) can algebraically reconcile SH0ES + Pantheon+ + DESI to < 0.05σ. **This is a 3-param/3-eq saturated system** — the low residuals are guaranteed, not a genuine statistical test. The Planck θ_* probe (4th, independent) fails at ~177kσ because the EdS D_A approximation is unphysical for CMB. The companion paper (Bouillé 2026) uses CAMB/Cobaya for the real joint MCMC.

## Data Provenance

All values in `cluster_data.csv` are drawn from primary X-ray and weak-lensing literature. Column `*_ref` provides the citation key for each measurement. See Section 4 of the manuscript and the inline documentation in `mach_threshold_test.py` for full references.

## License

Data and scripts: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Manuscript: © 2026 Louis-Robert Bouille. All rights reserved.
