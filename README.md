# Viscoelastic Vacuum — Supplementary Data & Reproduction Scripts

[![DOI](https://img.shields.io/badge/Preprint-Bouille%202026-blue)](https://arxiv.org/abs/XXXX.XXXXX)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Supplementary material for:

> **The Viscoelastic Vacuum: A Falsifiable Equation of State for Gravity Without Dark Sectors**  
> Louis-Robert Bouille (2026)

## Contents

| File | Description |
|---|---|
| `cluster_data.csv` | Full dataset of 9 merging clusters + 1 null test (NGC 1052-DF2), with shock velocities, pre-shock ICM temperatures, electron densities, gas masses, lensing-to-baryon mass ratios, and primary literature references. |
| `mach_threshold_test.py` | Binary Mach threshold test (Section 4.6). Computes environment-dependent Mach numbers using local ICM sound speed $c_s(T_X)$ for each cluster. Reproduces the 9/9 + null concordance rate. |
| `mach_threshold_test.png` | Output figure from the test script. |

## Quick Start

```bash
git clone https://github.com/lrb-research/viscoelastic-vacuum.git
cd viscoelastic-vacuum
python mach_threshold_test.py
```

**Requirements:** Python 3.8+, NumPy, Matplotlib (standard Anaconda/pip).

## Key Result

The binary Mach threshold test predicts that cluster mergers with scalar Mach number $\mathcal{M} = v_\text{shock}/c_s(T_X) > 1$ should exhibit a dark matter–gas spatial offset, while subsonic systems should not. Using the **local** pre-shock ICM sound speed $c_s(T_X) = \sqrt{\gamma k_B T_X / \mu m_p}$ for each cluster:

- **9/9 supersonic mergers** exhibit the predicted lensing offset ✓
- **1/1 null test** (NGC 1052-DF2, quiescent UDG) shows $R \approx 1$ ✓
- **Minimum Mach number in sample: 1.30** (Abell 754) — all clearly supersonic

This is **categorically distinct** from ΛCDM, where collisionless dark matter separates from gas at *any* velocity.

## Data Provenance

All values in `cluster_data.csv` are drawn from primary X-ray and weak-lensing literature. Column `*_ref` provides the citation key for each measurement. See Section 4 of the manuscript and the inline documentation in `mach_threshold_test.py` for full references.

## License

Data and scripts: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).  
Manuscript: © 2026 Louis-Robert Bouille. All rights reserved.
