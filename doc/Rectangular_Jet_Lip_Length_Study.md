# Experimental Investigation of Mean and Turbulent Characteristics of a Free Air Rectangular Jet with/without Lip Length

## Paper Summary

**Authors:** Youssef Gamal Nasr Ahmed, Reza Alidoost Dafsari, Jeekeun Lee
**Affiliation:** Jeonbuk National University, Korea; Minia University, Egypt
**Contact:** leejk@jbnu.ac.kr

---

## 1. Research Objective

This study experimentally investigates how **lip length extension** on a rectangular nozzle affects the flow structure of a free turbulent jet. The goal is to understand how this simple geometric modification can be used as a **passive flow control** method to enhance mixing, air entrainment, and jet spreading.

---

## 2. Experimental Setup

### Nozzle Configuration
- **Nozzle type:** Rectangular with adaptable upper/lower plates
- **Dimensions:** Width (w) = 190 mm, Height (h) = 10 mm
- **Aspect Ratio (AR):** 19 (width/height)
- **Lip length configurations tested:**
  - L₀ = 0 mm (equal length, reference case)
  - L₁₀ = 10 mm lower lip extension
  - L₂₀ = 20 mm lower lip extension

### Flow Conditions
- **Exit velocity (U₀):** ~30 m/s
- **Reynolds number (Re_h):** ~20,100
- **Ambient temperature:** 19°C ± 0.5°C

### Measurement Techniques
1. **Constant Temperature Anemometry (CTA):** 2D hot-wire sensor for detailed turbulence measurements
   - Sampling frequency: 30 kHz
   - 200,000 instantaneous samples per measurement point

2. **Particle Image Velocimetry (PIV):** For flow visualization and POD analysis
   - Double pulse laser (200 mJ/pulse, 15 Hz)
   - 12 MP CMOS camera
   - 4000 snapshots for POD analysis

### Measurement Domain
- Axial range: 0 ≤ x/h ≤ 30
- Radial range: -2.5 ≤ y/y₀.₅ ≤ 2.5

---

## 3. Key Findings

### 3.1 Mean Velocity Characteristics

#### Axial Velocity (U/U_c)
- **L₀ nozzle:** Symmetric velocity profile about jet center
- **L₁₀ and L₂₀ nozzles:** Asymmetric profiles with broader distribution on the free side
- Maximum velocity position shifts toward the lip side as lip length increases
- Self-similarity achieved at different downstream positions:
  - L₂₀: x/h = 1 (fastest)
  - L₁₀: x/h = 3
  - L₀: x/h = 5

#### Radial Velocity (V/U_c)
- Extended lip length reduces peak radial velocity but broadens distribution
- Enhanced outward momentum and radial spreading on free side
- Important for applications like **pre-filming air-blast atomizers**

### 3.2 Potential Core Length

The potential core (where U_c ≈ 0.98 U₀) shortens significantly with lip extension:

| Nozzle | Potential Core Length (P_c) |
|--------|----------------------------|
| L₀     | 3.85h                      |
| L₁₀    | 2.82h                      |
| L₂₀    | 2.37h                      |

**Implication:** Shorter potential core indicates enhanced mixing and flow redistribution.

### 3.3 Jet Spreading Rate

Half-width correlation: y₀.₅/h = SR_y × (x/h + V_or/h)

| Nozzle | Spreading Rate (SR_y) | Virtual Origin (V_or/h) |
|--------|----------------------|------------------------|
| L₀     | 0.1166               | -3.17                  |
| L₁₀    | 0.1192               | -2.27                  |
| L₂₀    | 0.1265               | -1.49                  |

**Result:** Greater lip length → higher spreading rate → virtual origin closer to nozzle exit.

### 3.4 Air Entrainment

- 4% increase in air entrainment for L₂₀ compared to L₀ at nozzle exit
- 2% increase for L₁₀ compared to L₀
- Effect most pronounced up to x/h = 10
- Beyond x/h = 10, entrainment becomes independent of initial geometry

### 3.5 Turbulence Characteristics

#### Axial Turbulence Intensity (U_rms/U_c)
- Two separate peaks observed (dual shear layer structure)
- Peak intensities and widths increase on free side with longer lip lengths
- Asymptotic value: ~0.25

#### Radial Turbulence Intensity (V_rms/U_c)
- Transitions from two peaks (near field) to single peak (~0.18) downstream
- Shear layer broadens on free side with greater lip length

#### Reynolds Shear Stress (u'v'/U_c²)
- Largest stresses in shear layers
- Broader profiles on free side for lip length nozzles
- Indicates enhanced turbulence and momentum transfer

#### Triple Velocity Products
- **u'u'v'/U_c³:** Turbulent transport of Reynolds normal stress
- **u'v'v'/U_c³:** Turbulent transport of Reynolds shear stress
- L₂₀ produces larger inward/outward transport regions
- Results in increased centerline velocity decay and greater jet spread

---

## 4. Coherent Structure Analysis (POD)

### Proper Orthogonal Decomposition Results
- First mode energy: ~3% (indicates fully turbulent, complex jet)
- Spatial structures vary significantly with nozzle geometry
- Higher modes show increasingly asymmetric vortical structures for lip length nozzles

### Gamma Criterion Analysis
- L₀: Symmetric vortex roll-up on both shear layers
- L₁₀, L₂₀: Asymmetric vortex roll-up with delayed interaction on lip side
- Confirms altered shear layer dynamics and ambient air interaction

---

## 5. Self-Similarity Correlations

### Mean Axial Velocity
Profiles converge to Gaussian, Goertler, Bradbury, and Tollmien distributions downstream.

### Mean Radial Velocity
Semi-empirical correlation:
```
V/U_c = [a(y/y₀.₅) + b(y/y₀.₅)³] × exp(-c(y/y₀.₅)²)
```
where a = 0.074, b = -0.024, c = 0.53

### Axial Turbulence Intensity
Modified Miller-Comings formulation:
```
U_rms/U_c = a × exp(-b(y/y₀.₅)²) - c × exp(-d(y/y₀.₅)²)
```
where a = 0.509, b = 0.51, c = 0.301, d = 1.13

---

## 6. Practical Applications

### Pre-filming Air-blast Atomizers
- Broader radial velocity profile allows atomizing air to interact with liquid film over wider region
- Promotes more effective primary breakup
- Yields more uniform and widely distributed spray

### Gas Turbine Combustors
- Optimal atomization essential for stable, efficient combustion
- Simultaneous widening of axial and radial velocity profiles leads to:
  - Better mixing
  - Elevated atomization efficiency
  - More uniform spray distributions

### Flow Control
- Lip length modification offers passive flow control method
- No moving parts or energy input required
- Applicable to various industrial free turbulent jet technologies

---

## 7. Key Nomenclature

| Symbol | Description |
|--------|-------------|
| AR | Aspect ratio (width/height) |
| h | Nozzle height (10 mm) |
| L | Lip length extension |
| U₀ | Mean centerline exit velocity |
| U_c | Mean axial centerline velocity |
| U_rms | RMS of axial velocity fluctuation |
| V_rms | RMS of radial velocity fluctuation |
| y₀.₅ | Half-width (y-location where U = 0.5 U_c) |
| P_c | Potential core length |
| SR_y | Jet spreading rate |
| V_or | Virtual origin |
| Re_h | Reynolds number based on h |

---

## 8. Conclusions

1. **Asymmetry:** Increasing lip length induces asymmetry between upper and lower shear layers

2. **Enhanced Mixing:** Lip length extension leads to:
   - Reduced potential core length
   - Increased centerline turbulence
   - Greater jet half-width
   - Elevated volume flow rate

3. **Turbulence Enhancement:** All turbulence metrics (intensity, shear stress, triple products) show broader distributions and higher values with longer lip lengths

4. **Earlier Self-Similarity:** L₂₀ reaches self-similar region faster than L₁₀ and L₀

5. **Flow Control Potential:** Simple geometric modification provides effective passive flow control without additional energy input

---

## 9. References to Related Work

Key prior studies cited:
- Kiwata et al. (2009): Similar geometry, lower velocity (7 m/s)
- Hirata et al. (2010): AR=300, Re=6000, lip length effects
- Bridges & Wernet (2015): Bevelled rectangular nozzles
- Deo et al. (2007, 2008): Aspect ratio and Reynolds number effects

---

## 10. Relevance to CFD Simulation

This experimental data provides valuable validation benchmarks for:
- Lattice Boltzmann Method (LBM) simulations
- Large Eddy Simulation (LES) of turbulent jets
- Reynolds-Averaged Navier-Stokes (RANS) models

Key quantities to validate:
- Centerline velocity decay
- Jet spreading rate
- Turbulence intensity profiles
- Reynolds stress distributions
- Potential core length
