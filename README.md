# Viscoelastic Vacuum — Supplementary Data & Reproduction Scripts

[![DOI](https://img.shields.io/badge/Preprint-Bouille%202026-blue)](https://arxiv.org/abs/XXXX.XXXXX)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Supplementary material for:

> **The Viscoelastic Vacuum: A Falsifiable Equation of State for Gravity Without Dark Sectors**
> Louis-Robert Bouille (2026)

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
| `gauge_asymmetry_test.py` | **The 30-second Hubble tension test.** Reconciles DESI BAO + SH0ES + Pantheon+ at zero free parameters via the Etherington dioptry. |
| `etherington_violation_calculator.py` | Predicts the violation $D_L/[D_A(1+z)^2] = (1+z)^\epsilon$ at any redshift. Discriminant test: Euclid×LSST cross-correlation at z~1 yields a 32% deviation, well above instrumental precision. |
| `mach_cone_calculator.py` | Computes the gravitational Mach cone offset $\Delta x$ for any merging cluster from $T_X$ (keV) and $v_{\rm shock}$ (km/s). Reproduces Bullet Cluster ($\theta=16.2°$) and El Gordo (out-of-sample) at zero parameters. |
| `mach_threshold_test.py` | Binary Mach threshold test (Section 4.6). Reproduces the 9/9 + null concordance rate using local ICM sound speeds. |
| `cluster_data.csv` | Full dataset of 9 merging clusters + 1 null test (NGC 1052-DF2): shock velocities, pre-shock ICM temperatures, electron densities, gas masses, lensing-to-baryon mass ratios, primary references. |
| `mach_threshold_test.png` | Output figure from the binary threshold test. |

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

**Requirements:** Python 3.8+, NumPy, Matplotlib (standard Anaconda/pip).

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

## Data Provenance

All values in `cluster_data.csv` are drawn from primary X-ray and weak-lensing literature. Column `*_ref` provides the citation key for each measurement. See Section 4 of the manuscript and the inline documentation in `mach_threshold_test.py` for full references.

## License

Data and scripts: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Manuscript: © 2026 Louis-Robert Bouille. All rights reserved.
