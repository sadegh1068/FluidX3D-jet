# Literature Research: LES/LBM Simulation of Turbulent Rectangular Free Jets

This document summarizes the literature research for simulating the L₀ (equal lip length) rectangular jet case from the experimental study by Ahmed et al.

**Target Case:**
- Rectangular nozzle: AR = 19 (width 190mm, height h = 10mm)
- Exit velocity: U₀ = 30 m/s
- Reynolds number: Re_h ≈ 20,100
- Quiescent ambient (free jet)

---

## Table of Contents
1. [Inlet Boundary Conditions](#1-inlet-boundary-conditions)
2. [Outlet Boundary Conditions](#2-outlet-boundary-conditions)
3. [Side/Farfield Boundary Conditions](#3-sidefarfield-boundary-conditions)
4. [Turbulence Models for LES](#4-turbulence-models-for-les)
5. [Additional Considerations](#5-additional-considerations)
6. [Key Reference Papers](#6-key-reference-papers)
7. [Summary Recommendations](#7-summary-recommendations)

---

## 1. Inlet Boundary Conditions

### 1.1 Velocity Profile Specifications

#### Top-Hat (Uniform) vs Developed Profiles

| Profile Type | Characteristics | When to Use |
|--------------|-----------------|-------------|
| **Top-Hat (Uniform)** | Flat velocity profile with thin boundary layer | High contraction ratio nozzles (r/h ≈ 3.6), representative of experimental setup |
| **Hyperbolic Tangent (tanh)** | Smooth transition at shear layer edges | Models smooth velocity transition, good for stability analysis |
| **Developed (1/7th power law)** | Turbulent pipe flow profile | Long pipe/channel upstream of nozzle |

**Recommendation for this case:** The experimental nozzle has a contraction, so a **top-hat profile** with thin boundary layers is appropriate. Use hyperbolic tangent smoothing at edges:

```
U(y) = U₀/2 * [1 - tanh((|y| - h/2) / (2θ))]
```

where θ is the momentum thickness (typically θ/h ≈ 0.01-0.05).

#### Boundary Layer Treatment
- Nozzle-exit boundary layer thickness significantly affects jet development
- Thinner boundary layer → longer potential core, lower centerline turbulence
- For Re ≈ 20,000: expect δ/h ≈ 0.05-0.1 at nozzle exit

### 1.2 Turbulence Generation at Inlet (CRITICAL)

**Key Finding:** Simply adding random noise is insufficient and leads to delayed/failed transition to turbulence. Correlated, structured fluctuations are essential.

#### Method Comparison Table

| Method | Correlation | Divergence-Free | Cost | Development Length | Recommendation |
|--------|-------------|-----------------|------|-------------------|----------------|
| **White Noise** | None | No | Very Low | Fails/Very Long | ❌ Not recommended |
| **Digital Filter (DFM)** | Spatial + Temporal | Modifiable | Medium | ~10δ | ✓ Good |
| **Synthetic Eddy (SEM)** | Spatial + Temporal | No | Medium | ~10δ | ✓ Good |
| **Divergence-Free SEM (DFSEM)** | Spatial + Temporal | Yes | Medium | ~10δ | ✓✓ Recommended |
| **Vortex Method** | Spatial + Temporal | Partial | Low-Medium | ~10-20δ | ✓ Good |
| **Spectral/RFM** | Spectral | Yes | Medium | ~10δ | ✓ Good |
| **Lund Recycling** | Full NS | Yes | High | ~50δ | Best accuracy |

#### Recommended Methods for Rectangular Jet

**Option 1: Divergence-Free Synthetic Eddy Method (DFSEM)**
- Poletto, Craft & Revell (2013)
- Produces divergence-free turbulence field
- Reduces pressure fluctuations and adaptation length
- 33% shorter recovery length than standard SEM

**Option 2: Digital Filter Method (Klein et al.)**
- Klein, Sadiki, Janicka (2003)
- Generate random data, apply Gaussian filter for length scales
- Adjust to match Reynolds stress tensor

**Option 3: Vortex Ring Perturbation (Bogey & Bailly)**
- Add divergence-free vortex rings at nozzle exit
- Use 32 azimuthal modes with random phase/amplitude
- Critical: Include first 4 azimuthal modes for proper development

#### Turbulence Intensity Levels

| Flow Type | Typical Tu |
|-----------|------------|
| Wind tunnel (low turbulence) | 0.05% - 0.1% |
| Jet from smooth contraction | 1% - 3% |
| **Typical jet simulations** | **5% - 10%** |
| Complex geometries | 10% - 20% |

**For this case:** Tu ≈ 5-8% is recommended based on similar jet studies (Jaramillo: Tu = 8.8%).

#### Turbulence Length Scales
- Should not exceed problem dimension
- Typical: L_t ≈ 0.05h to 0.1h (5-10% of nozzle height)
- For pipe flow correlation: L_t ≈ 0.038D_h

### 1.3 LBM-Specific Inlet Treatments

| Method | Description | Accuracy | Stability |
|--------|-------------|----------|-----------|
| **Zou-He BC** | Bounce-back of non-equilibrium part | 2nd order | Good |
| **Non-Equilibrium Extrapolation (NEQ)** | f = f_eq + f_neq (extrapolated) | 2nd order | Better |
| **Regularized BC** | Replace all f components via regularization | 2nd order | Highest |
| **Equilibrium BC** | f = f_eq only | 1st order | Simple |

**Recommendation:** Use **Zou-He** or **NEQ** for velocity inlet with synthetic turbulence superimposed.

**LBM Synthetic Turbulence:**
- Generate velocity fluctuations using DFSEM or similar
- Add to mean velocity at inlet nodes each timestep
- Regularized LBM reconstructs distribution functions

---

## 2. Outlet Boundary Conditions

### 2.1 Convective (Orlanski) Boundary Condition

**Formulation:**
```
∂φ/∂t + U_c * ∂φ/∂n = 0
```

**Convection Velocity Selection:**
- **Local calculation:** Can cause instabilities in turbulent flows
- **Averaged velocity:** More stable, recommended approach
- Stability requirement: 0 ≤ U_c ≤ Δx/Δt

**Issues:**
- Can modify jet spreading rates
- Results may become domain-dependent
- May require larger domains than area of interest

### 2.2 Zero-Gradient (Neumann) Conditions

```
∂φ/∂n = 0
```

**Limitations:**
- Acts as reflecting boundary for waves
- Can cause spurious reflections propagating back
- Not truly "open" boundary

### 2.3 Buffer/Sponge Zones (Recommended)

**Formulation:**
```
∂q/∂t = RHS - σ(x)(q - q_ref)
```

**Damping Profile Selection:**

| Profile | Formula | Application |
|---------|---------|-------------|
| Linear | σ ~ x | Basic |
| **Quadratic** | σ ~ x² | **Optimal for 20-60 dB damping** |
| Cubic | σ ~ x³ | >60 dB damping |

**Sponge Length:**
- Ineffective if < 30% of wavelength
- Near-ideal if > 10 wavelengths
- **Practical: 5-10h length recommended**

### 2.4 Non-Reflecting BC (NSCBC)

- Poinsot & Lele (1992) formulation
- Decomposes solution into characteristic waves
- Partially reflecting at low frequencies
- Use relaxation parameter σ ≈ 0.25-0.30 for jets

### 2.5 LBM-Specific Outlet Treatments

| Method | Description | Notes |
|--------|-------------|-------|
| **Population Extrapolation** | Extrapolate distribution functions | Simple, stable |
| **Convective BC for LBM** | Lou et al. (2013) | Good for jets |
| **NEQ Extrapolation** | Extrapolate non-equilibrium part | 2nd order |
| **Pressure Outlet** | Fixed pressure with extrapolation | Common choice |

### 2.6 Domain Length Requirements

| Study Type | Outlet Distance from Nozzle |
|------------|----------------------------|
| Basic jet LES | 25-30h |
| High-fidelity DNS | 40-50h |
| Aeroacoustic studies | 50-75h |
| Self-similarity validation | 40h minimum |

**Recommendation for this case:** Outlet at **x = 40-50h** (400-500mm) to capture the experimental measurement domain (0 ≤ x/h ≤ 30) plus buffer zone.

---

## 3. Side/Farfield Boundary Conditions

### 3.1 Entrainment Boundary Conditions

Free jets entrain ambient fluid - boundaries must allow this inflow.

**Options:**
1. **Pressure inlet/outlet at atmospheric pressure** - allows bidirectional flow
2. **Traction-free BC** - Lagrangian-based, stable for turbulent flows
3. **Opening condition** - permits inflow/outflow based on local pressure

### 3.2 Co-flow vs Quiescent Ambient

| Approach | Description | Notes |
|----------|-------------|-------|
| **Quiescent (U_∞ = 0)** | True free jet | Requires entrainment BC |
| **Small co-flow (U_∞/U₀ ≈ 0.05-0.1)** | Numerical stability | Common practice |

**Recommendation:** Use quiescent ambient with proper entrainment BC, or very small co-flow (< 3% of U₀) if stability issues arise.

### 3.3 Domain Size Requirements

**Lateral Extent (y-direction):**

| Configuration | Lateral Domain |
|---------------|---------------|
| Minimum | ±10h from centerline |
| Recommended | ±15-20h from centerline |
| With sponge zones | ±10h + 5h sponge |
| For self-similarity | ±30-60 jet radii |

**Spanwise Extent (z-direction):**

For AR = 19 rectangular jet:
- Full width needed: 190mm (19h)
- **Periodic BC acceptable** if flow is statistically homogeneous in z
- Minimum spanwise extent if periodic: 6h

### 3.4 LBM-Specific Farfield Treatments

| Method | Description | Recommendation |
|--------|-------------|----------------|
| **Equilibrium BC** | f = f_eq at boundaries | Simple, but imposes wrong gradients |
| **NEQ Extrapolation** | f = f_eq + f_neq_extrapolated | Better, 2nd order |
| **Characteristic BC** | Thompson's CBC adapted to LBM | Best for waves |
| **PML (Absorbing)** | Perfectly matched layer | Best for acoustics |

**Recommendation:** Use **equilibrium BC with low ambient velocity** or **NEQ extrapolation** for lateral boundaries.

### 3.5 Avoiding Artificial Confinement

**Checklist:**
- [ ] Lateral boundaries > 10h from jet edge
- [ ] Spreading rate matches experimental free jet values
- [ ] No artificial acceleration near boundaries
- [ ] Domain independence study performed

---

## 4. Turbulence Models for LES

### 4.1 Subgrid-Scale (SGS) Model Comparison

| Model | Strengths | Weaknesses | Jets Performance |
|-------|-----------|------------|-----------------|
| **Static Smagorinsky** | Simple, robust | Overly dissipative, no backscatter | Fair |
| **Dynamic Smagorinsky** | Self-calibrating, allows backscatter | Computational overhead | Excellent |
| **WALE** | Good wall behavior, transition prediction | Grid sensitive | Very Good |
| **Vreman** | Simple as static, accurate as dynamic | Slightly more dissipative | Very Good |
| **Sigma** | Accurate near-wall | Less validated | Good |
| **Mixed Models** | Highest correlation with true SGS stress | Complex | Excellent |

### 4.2 Recommended SGS Models for This Case

**Primary Recommendation: Smagorinsky-Lilly** (available in FluidX3D)
- Simple implementation
- Well-validated for free shear flows
- Use C_s ≈ 0.12-0.15 for free shear layers

**If available: Dynamic Smagorinsky or WALE**
- Better for transition and outer regions
- Self-calibrating coefficients

### 4.3 Smagorinsky Constant Values

| Flow Region | C_s Value |
|-------------|-----------|
| Isotropic turbulence | 0.17 - 0.23 |
| Shear flows | ~0.10 |
| **Free shear layers** | **0.12 - 0.15** |
| Near walls | 0.10 - 0.12 |

### 4.4 LBM-Specific Turbulence Modeling

**Implementation in LBM:**
```
ν_eff = ν₀ + ν_t
ν_t = (C_s * Δ)² * |S̄|
τ_eff = 3ν_eff + 0.5
```

The strain rate tensor |S̄| can be computed directly from non-equilibrium distributions in LBM:
```
S_ij = -3ω_s/(2ρ₀) * Π_ij^(neq)
```

### 4.5 Collision Operators for High-Re Flows

| Operator | Stability at High Re | Recommendation |
|----------|---------------------|----------------|
| BGK (SRT) | Poor | ❌ Not for Re > 10,000 |
| TRT | Moderate | ✓ Acceptable |
| MRT | Good | ✓✓ Recommended |
| **Cumulant** | Excellent | ✓✓✓ Best for high Re |

**For Re ≈ 20,000:** Use MRT or Cumulant collision operator if available.

### 4.6 Grid Resolution Requirements

**Pope's Criterion:** Resolve > 80% of turbulent kinetic energy

**Resolution Guidelines:**

| Region | Resolution |
|--------|------------|
| Nozzle height | 20-50 cells across h |
| Shear layer | Resolve momentum thickness (3-5 cells in δ_θ) |
| Downstream | Can coarsen gradually |
| Acoustic | ~12 cells per wavelength |

**Typical Grid Sizes:**

| Fidelity | Total Cells | Notes |
|----------|-------------|-------|
| Coarse LES | 5-10 million | Qualitative |
| Standard LES | 15-30 million | Good accuracy |
| Fine LES | 50-80 million | High accuracy |
| DNS-like | 100+ million | Research grade |

---

## 5. Additional Considerations

### 5.1 Mach Number Constraints (LBM)

| Ma Range | LBM Treatment |
|----------|---------------|
| Ma < 0.3 | Standard LBM (incompressible) |
| 0.3 < Ma < 0.5 | Extended equilibrium |
| Ma > 0.5 | Compressible LBM formulations |

**For this case:** U₀ = 30 m/s → Ma ≈ 0.09 → **Standard incompressible LBM is valid**

Density variation should be < 0.5% for incompressible treatment.

### 5.2 Time Stepping and Statistics

**CFL Requirements:**
- LBM: CFL = 1 by construction
- LES explicit schemes: CFL ≈ 0.5

**Flow-Through Times:**
- Flush transients: 2-3 domain flow-through times
- Statistics collection: **Minimum 10 large-eddy turnover times**
- Higher-order moments: 100+ turnover times

**Turnover Time Estimate:**
```
T_eddy ≈ h / U₀ = 10mm / 30m/s = 0.33 ms
```

**Sampling Frequency:**
- Sample every ~10 large-eddy times for time-averaging
- For spectra: f_sample > 2 * f_max (Nyquist)

### 5.3 Initial Conditions

**Initialization Strategy:**
1. Start with quiescent ambient (U = 0 everywhere except inlet)
2. Ramp up inlet velocity over ~10-50 timesteps
3. Run 2-3 flow-through times to flush transients
4. Begin statistics collection

**Perturbation Methods:**
- Add vortex ring perturbations at inlet
- Use synthetic turbulence from t = 0
- Small random perturbations (~0.1% U₀) can help trigger instabilities

### 5.4 Validation Metrics

**Primary Quantities to Compare:**

| Quantity | Experimental Value (L₀) | Where to Compare |
|----------|------------------------|------------------|
| Potential core length P_c | 3.85h | Centerline |
| Spreading rate SR_y | 0.1166 | Half-width vs x |
| Virtual origin V_or/h | -3.17 | Half-width fit |
| Centerline turbulence | ~20% at x/h=30 | U_rms/U_c |

**Self-Similarity Profiles:**
- U/U_c vs y/y₀.₅ (should match Gaussian/Goertler)
- U_rms/U_c vs y/y₀.₅ (two peaks at y/y₀.₅ ≈ ±1)
- u'v'/U_c² vs y/y₀.₅ (peaks in shear layers)

**Centerline Decay:**
```
U₀/U_c ~ (x/h)^0.5  (after potential core)
```

### 5.5 Common Pitfalls

| Issue | Symptom | Solution |
|-------|---------|----------|
| Insufficient inlet turbulence | Laminar jet, delayed transition | Use structured perturbations |
| Under-resolved grid | Wrong spreading rate, overly coherent structures | Refine shear layer |
| Small domain | Artificial confinement, wrong spreading | Increase lateral extent |
| Reflecting boundaries | Spurious oscillations | Add sponge zones |
| Numerical dissipation | Wrong decay rates | Use low-dissipation schemes |
| Short statistics | Noisy profiles | Longer averaging time |

### 5.6 Computational Estimates

**Domain Size for This Case:**
- Streamwise: 0 to 50h = 500mm
- Lateral: ±20h = ±200mm (400mm total)
- Spanwise: 19h = 190mm (full width, periodic)

**Grid Estimate (h = 10mm, 30 cells across h):**
- Δx = h/30 = 0.33mm
- Streamwise cells: 500/0.33 ≈ 1500
- Lateral cells: 400/0.33 ≈ 1200
- Spanwise cells: 190/0.33 ≈ 570
- **Total: ~1 billion cells** (uniform grid)

**With Variable Resolution:**
- Fine region (jet core + shear layer): ~100-200 million cells
- Can reduce to ~50 million with aggressive coarsening

**Timestep (LBM):**
```
Δt = Δx * Ma_lattice / c_s
```
For Ma_lattice = 0.1, c_s = 1/√3:
```
Δt ≈ 0.33mm * 0.1 / (343 m/s * 1/√3) ≈ 0.17 μs
```

**Physical Time for Statistics:**
- 10 flow-throughs: ~170ms
- Required timesteps: ~1 million

---

## 6. Key Reference Papers

### Inlet Boundary Conditions
1. **Klein, Sadiki, Janicka (2003)** - Digital Filter Method, J. Comput. Phys. 186:652-665
2. **Jarrin et al. (2006)** - Synthetic Eddy Method, Int. J. Heat Fluid Flow
3. **Poletto, Craft, Revell (2013)** - Divergence-Free SEM, Flow Turb. Combust.
4. **Bogey & Bailly (2010)** - Nozzle-exit boundary effects, J. Fluid Mech. 663

### Outlet/Open Boundaries
5. **Orlanski (1976)** - Convective BC, J. Comput. Phys. 21:251-269
6. **Poinsot & Lele (1992)** - NSCBC, J. Comput. Phys. 101:104-129
7. **Bodony (2006)** - Sponge zones analysis, J. Comput. Phys. 212:681-702

### LBM Boundaries
8. **Zou & He (1997)** - Velocity/pressure BC for LBM, Phys. Fluids 9:1591
9. **Guo et al. (2002)** - Non-equilibrium extrapolation, Chinese Physics 11:366
10. **Lou et al. (2013)** - Convective BC for LBM, Computers Math. Appl.

### Turbulence Models
11. **Germano et al. (1991)** - Dynamic Smagorinsky, Phys. Fluids A 3:1760
12. **Nicoud & Ducros (1999)** - WALE model, Flow Turb. Combust. 62:183
13. **Vreman (2004)** - Vreman model, Phys. Fluids 16:3670

### LBM for Jets
14. **Yu & Girimaji (2005)** - Rectangular jets LBM, Phys. Fluids 17:125106
15. **Yu, Luo, Girimaji (2006)** - MRT-LBM for square jet, Computers & Fluids
16. **Geier et al. (2017)** - Cumulant LBM, J. Comput. Phys. 348:862

### Jet Experiments (Validation)
17. **Hussein, Capp, George (1994)** - Round jet at Re~10⁵, J. Fluid Mech. 258:31
18. **Deo, Mi, Nathan (2007)** - Plane jets, Exp. Therm. Fluid Sci. 31:825
19. **Suresh et al. (2008)** - Plane jet development, Phys. Fluids 20:044105

---

## 7. Summary Recommendations

### For FluidX3D Implementation

| Aspect | Recommendation |
|--------|----------------|
| **Inlet BC** | Zou-He or Equilibrium with top-hat profile |
| **Inlet Turbulence** | Synthetic perturbations (5-8% intensity) or precursor |
| **Outlet BC** | Equilibrium BC with sponge zone (5-10h) |
| **Side BC** | Equilibrium BC at atmospheric pressure |
| **Spanwise BC** | Periodic (full width 190mm) |
| **SGS Model** | Smagorinsky-Lilly (C_s = 0.12) |
| **Collision** | SRT with SUBGRID enabled, or TRT if available |
| **Resolution** | 30-50 cells across h in shear layer |
| **Domain** | 50h × 40h × 19h (streamwise × lateral × span) |

### Critical Success Factors

1. ✓ **Structured inlet turbulence** - not just random noise
2. ✓ **Sufficient lateral domain** - avoid artificial confinement
3. ✓ **Sponge zone at outlet** - prevent reflections
4. ✓ **Adequate resolution** - especially in shear layer
5. ✓ **Long enough statistics** - 10+ eddy turnover times
6. ✓ **Validate against experiment** - P_c, spreading rate, profiles

### Next Steps

1. Review FluidX3D capabilities against requirements
2. Identify gaps (inlet turbulence generation, specific BCs)
3. Plan implementation strategy
4. Set up validation case
5. Compare with experimental data

---

*Document compiled from literature research - January 2026*
