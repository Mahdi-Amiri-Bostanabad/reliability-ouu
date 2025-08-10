# reliability-ouu
Final Project of the Reliability Analysis and Risk Assessment (AE45-788) course - Spring 2024

# Quadrotor Reliability & Optimization Under Uncertainty (UAV/VTOL)

> **Deterministic sizing → Reliability (FORM/SORM, DS, IS) → SORA (OUU) → Epistemic UQ (Sobol)**  
> Reproducible modeling pipeline for a quadrotor delivering a 2 kg payload over 50 km.

[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](#)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-orange.svg)](#)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## Abstract

This repository contains a complete, reproducible workflow for sizing, reliability analysis, and optimization under uncertainty of a quadrotor mission. We start with a physics-based **deterministic design** of rotor radius and cruise speed under disk-loading, blade-loading, and energy constraints. We then endow the limits with aleatory uncertainty (lognormal \(C_{d0}\) and \(\sigma\)) to define **probabilistic failure boundaries**, and quantify reliability with **FORM/SORM** (fast indices/design points) plus **Directional Sampling (DS)** and **Importance Sampling (IS)** as variance-controlled checks. Building on this, we solve **Optimization Under Uncertainty (OUU)** with **SORA**, alternating between reliability targets in standard normal space and a deterministic subproblem. The final design satisfies \(p_{f,\mathrm{struct}}\le 10^{-4}\) and \(p_{f,\mathrm{energy}}\le 10^{-3}\) with modest mass growth (e.g., \(r=0.150\) m, \(V_\infty\approx 77\) m/s, and an added energy reserve moving mass from ≈ 5.0 kg to ≈ 6.6 kg at reliability). Lastly, we perform **epistemic UQ** via **Sobol indices** over drag-mean bias and Weibull fatigue shape, showing drag bias dominates mission-risk variance with negligible interactions—clear guidance to prioritize **drag data** over redesign.

**Reproducibility.** All scripts and data to reproduce the figures and tables are available here.  
If you mirror this work, please cite (see **Citation** below).

---

## What’s inside

- **`matlab/`**
  - `deterministic_sizing_and_opt.m` — Task 1–2: deterministic baseline sizing & fmincon optimization (your provided code).
  - `runReliability2.m` — Task 3: reliability (FORM/SORM + optional DS and IS) (your provided code).
  - Outputs `baseline_constants.mat` for downstream tasks.

- **`python/`**
  - `task4_sora.py` — Task 4: SORA for OUU (final v3; chain-rule gradients, unclipped reliability, DS-weighted & IS checkers).
  - `task5_sobol.py` — Task 5: Epistemic UQ via Sobol S1/ST for mission failure probability.
  - `notebooks/` — optional step-by-step Jupyter notebooks.

- **`docs/`**
  - `report/` — LaTeX source of the full project report (Tasks 1–5) + compiled PDF.

- **`results/`** — figures, tables, and logs created when you run the scripts.

---

## Key modeling choices (short)

- **Rotor speed** fixed at `Ω = 1110 rad/s`, benchmarked so total cruise power falls in ~0.5–1.0 kW for plausible designs.
- **Energy reserve**: base battery keeps **30 % SOC unused**; design knob `extra_m ∈ [0, 0.35]` adds reserve (heavier battery → larger usable energy).
- **Uncertainties** (aleatory): `C_d0 ~ Lognormal(mean=0.012, COV=0.20)`, `σ ~ Lognormal(mean=0.13, COV=0.12)`.
- **Reliability targets**: structural \(p_f \le 10^{-4}\) ⇒ \(\beta \approx 3.719\); energy \(p_f \le 10^{-3}\) ⇒ \(\beta \approx 3.090\).

---

## Quickstart

### 1) Python environment

```bash
# Option A: pip + venv
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# Option B: conda
conda env create -f environment.yml
conda activate uav-ouu
