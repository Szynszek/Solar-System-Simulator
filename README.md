<h1 align="center">High-Precision Symplectic N-Body Orbital Propagator</h1>

<h3 align="center">
A fully vectorized, coordinate-free N-body numerical integrator designed for long-term orbital propagation, nodal regression analysis, and deep-space trajectory modeling.
</h3>

<p align="center">
<img src="assets/orbital_simulation.gif" alt="" width="400">
</p>

<p align="center">
Example visualization of J2-induced Nodal Regression (RAAN precession) over long-term integration.
</p>

<hr>

## Executive Summary

Built from scratch to maintain strict symplectic energy bounds and conserve system linear momentum (Barycenter stability). This engine bypasses standard test-mass approximations by implementing two-way reflex accelerations for gravitational harmonics ($J_2$).

## Core Architecture & Physics Engine

### 1. Symplectic Time Evolution (PEFRL)

Utilizing a 4th-order **Position Extended Forest-Ruth Like (PEFRL)** symplectic algorithm for phase-space volume conservation. It guarantees that the system's Shadow Hamiltonian remains bounded over multi-decade integrations.

### 2. **$J_2$** Harmonics

Implements a fully algebraic formulation of planetary oblateness ($J_2$). By utilizing projections of rotational poles ($z_{local}$), the engine entirely avoids computationally expensive and mathematically stiff rotation matrices (DCMs).

### 3. Strict Momentum Conservation

Eliminates the standard non-conservation of total momentum in hierarchical N-body formulations. The engine calculates and applies reflex/recoil accelerations to oblate bodies, guaranteeing no artificial Barycenter drift. The system calculates the mutual perturbation symmetrically:

$$
\vec{a}_{reflex} = - \frac{\mu_j}{\mu_i} \cdot \vec{a}_{direct}(\text{parameters}_i, \vec{r}_{ji})
$$

<hr>

## Validation & Error Analysis

The integrity of the physics engine is validated against theoretical limits and operational ephemerides:

### Shadow Hamiltonian Stability

Mathematical proof of engine stability. The relative pseudo-energy error remains bounded at the **O(10^-13)** level over a 30-year simulation of the Solar System (dt = 1800s).

> **Proof of Concept:** The high-frequency, zero-trend noise confirms that symplectic bounds are perfectly maintained despite the sharp $1/r^4$ gradients introduced by the $J_2$ potential.

<p align="center">
<img src="assets/energy_error.jpg" alt="Shadow Hamiltonian Stability" width="800">
</p>

<p align="center">
Strict energy conservation verified: Relative pseudo-energy error bounded at O(10^-13).
</p>

### JPL Horizons Benchmarking

System state outputs are actively benchmarked against NASA JPL Horizons data. Known discrepancies (e.g., Earth's secular along-track error of ~1881 km over 30 years) are strictly isolated, quantified, and attributed to the theoretical limits of classical Newtonian mechanics:

* Absence of Einstein-Infeld-Hoffmann (EIH) general relativity corrections.
* Lack of Barycentric Dynamical Time (TDB) dilation.

<img src="assets/position_error.jpg" alt="Position error" width="800">

<hr>

## Repository Structure

```
├── ephemeris.json          # Mission config, Epoch-of-Date vectors, J2 data
├── solar_system.m          # Main execution, preallocation and visualization.
├── jpl_scraper.py          # Python script for generating JSON ephemeris
└── README.md               # Architecture documentation
```

## Configuration & Data Traceability

Mission parameters and initial state vectors are decoupled from the physics engine via `ephemeris.json`.

**Data Traceability Standard:**

* **All** planetary constants (e.g., standard gravitational parameters $\mu$) and initial state vectors (Epoch-of-Date) are sourced uniformly from the NASA JPL Horizons system.

* The only exception is Jupiter, whose planetary constants are traced to NASA/TM-20210022058 (Jupiter Global Reference Atmospheric Model: User Guide; H. L. Justh et al.). 
