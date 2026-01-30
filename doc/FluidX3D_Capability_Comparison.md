# FluidX3D Capability Comparison: Research Requirements vs Available Features

This document compares the research-recommended features for simulating the L₀ rectangular jet case against FluidX3D's actual capabilities.

---

## Summary Table

| Requirement | Research Recommendation | FluidX3D Status | Gap Level |
|------------|------------------------|-----------------|-----------|
| **Collision Operator** | MRT/Cumulant for Re≈20,000 | SRT, TRT | MEDIUM |
| **SGS Model** | Smagorinsky-Lilly (C_s=0.12) | Smagorinsky-Lilly (SUBGRID) | NONE |
| **Velocity Set** | D3Q19 or D3Q27 | D3Q19 (default), D3Q27 | NONE |
| **Inlet BC** | Zou-He or NEQ | TYPE_E (Equilibrium only) | LOW |
| **Inlet Turbulence** | DFSEM, Digital Filter, SEM | None (only random noise) | HIGH |
| **Outlet BC** | Convective/NSCBC with sponge | TYPE_E + Extrapolation (run loop) | MEDIUM |
| **Sponge Zone** | 5-10h damping zone | Not available | HIGH |
| **Side/Farfield BC** | Equilibrium or NEQ | TYPE_E (Equilibrium) | NONE |
| **Periodic BC** | Periodic spanwise | Available | NONE |
| **Memory Compression** | FP16 for large grids | FP16S, FP16C | NONE |
| **Multi-GPU** | Domain decomposition | Full support | NONE |

---

## Detailed Analysis

### 1. Collision Operators

**Research Recommendation:**
- MRT or Cumulant for Re > 10,000 (better stability)
- TRT as acceptable alternative
- BGK/SRT not recommended for high-Re flows

**FluidX3D Status:**
```cpp
#define SRT  // default - Single Relaxation Time (BGK)
//#define TRT // Two Relaxation Time - available
```

**Gap Assessment:** MEDIUM
- FluidX3D has SRT (default) and TRT
- Missing: MRT, Cumulant
- **Mitigation:** TRT with SUBGRID should provide adequate stability for Re≈20,000

---

### 2. Subgrid-Scale (SGS) Turbulence Model

**Research Recommendation:**
- Smagorinsky-Lilly with C_s = 0.12-0.15 for free shear layers
- Dynamic Smagorinsky preferred if available

**FluidX3D Status:**
```cpp
//#define SUBGRID // Smagorinsky-Lilly subgrid turbulence LES model
```

Implementation in `kernel.cpp` (lines 1597-1611):
```cpp
// Smagorinsky constant embedded: C_s ≈ 0.15
// Formula: 0.76421222 = 18*sqrt(2)*(C*Δ)²
const float Q = sq(Hxx)+sq(Hyy)+sq(Hzz)+2.0f*(sq(Hxy)+sq(Hxz)+sq(Hyz));
w = 2.0f/(tau0+sqrt(sq(tau0)+0.76421222f*sqrt(Q)/rhon));
```

**Gap Assessment:** NONE
- Smagorinsky-Lilly model fully available
- Constant is ~0.15 (acceptable for shear flows; 0.12 is optimal)
- **Action:** Enable `#define SUBGRID` in defines.hpp

---

### 3. Inlet Boundary Conditions

**Research Recommendation:**
- Zou-He or Non-Equilibrium Extrapolation (NEQ) for 2nd order accuracy
- Top-hat velocity profile with tanh smoothing at edges

**FluidX3D Status:**
```cpp
#define TYPE_E 0b00000010 // equilibrium boundary (inflow/outflow)
// Sets f = f_eq (1st order accurate)
```

**Gap Assessment:** LOW
- TYPE_E sets distributions to equilibrium (f = f_eq)
- This is 1st order accurate (sufficient for this application)
- Velocity profile must be specified in setup.cpp during initialization
- **Action:** Implement tanh velocity profile in setup.cpp initialization

---

### 4. Inlet Turbulence Generation (CRITICAL)

**Research Recommendation:**
- Divergence-Free Synthetic Eddy Method (DFSEM) - best
- Digital Filter Method (DFM) - good alternative
- Vortex Ring Perturbation - good for jets
- Intensity: Tu = 5-8%
- Length scale: L_t ≈ 0.05-0.1h

**FluidX3D Status:**
```cpp
// Available in utilities.hpp:
inline float random(uint& seed, const float x=1.0f);
inline float random_symmetric(uint& seed, const float x=1.0f);
class SimplexNoise { /* 2D/3D/4D simplex noise */ };
```

**Gap Assessment:** HIGH (Most Critical Gap)
- No structured turbulence generation methods
- Only random noise functions available (insufficient for realistic jet development)
- SimplexNoise provides spatial correlation but not proper turbulence characteristics

**Required Custom Implementation:**
1. **Option A: Simplified Vortex Perturbation**
   - Add random vortex rings at inlet each timestep
   - Use SimplexNoise for phase variation
   - Lower fidelity but simpler to implement

2. **Option B: Digital Filter Method**
   - Generate spatially/temporally correlated velocity fluctuations
   - Apply Gaussian filter with length/time scales
   - More complex but better turbulence

3. **Option C: Precursor Simulation**
   - Run separate periodic channel simulation
   - Store turbulent velocity field
   - Recycle at inlet boundary
   - Highest fidelity but most complex

**Recommended Approach:** Start with simplified vortex perturbation + random noise (~5% intensity), validate against experimental data, enhance if needed.

---

### 5. Outlet Boundary Conditions

**Research Recommendation:**
- Convective BC (Orlanski) or NSCBC
- Sponge zone (5-10h) to absorb reflections
- Domain length: 40-50h

**FluidX3D Status:**
```cpp
#define TYPE_E 0b00000010 // equilibrium boundary (inflow/outflow)
// Sets f = f_eq at fixed pressure/velocity
```

**CRITICAL FINDING:** TYPE_E with u=0 at outlet acts as a **wall** because it enforces equilibrium with zero velocity, blocking the jet.

**Gap Assessment:** MEDIUM (requires custom implementation)

**Solution: GPU-Side Extrapolation BC (Kernel Modification Required)**

The hydraulic jump example (setup.cpp:1001-1004) shows TYPE_E works at outlet when velocity is specified. For jets where outlet velocity is unknown, extrapolation is needed.

**CRITICAL:** CPU-side implementation (read_from_device → process → write_to_device) would kill GPU performance. The solution must be implemented in **OpenCL kernel** (`kernel.cpp`).

**Required Kernel Modification Approach:**

Option A: Modify `stream_collide()` kernel to handle outlet extrapolation:
```cpp
// In kernel.cpp stream_collide() - pseudocode
if(flagsn_bo==TYPE_E && x==def_Nx-1) {
    // Extrapolate velocity from interior cell (x-1)
    const uxx n_interior = index(x-1, y, z);
    uxn = u[n_interior];
    uyn = u[def_N + n_interior];
    uzn = u[2ul*def_N + n_interior];
    // Then compute f_eq with extrapolated velocity
}
```

Option B: Create new kernel for outlet BC that runs after stream_collide()

Option C: Use a new flag (TYPE_X or TYPE_Y available) to mark extrapolation boundaries

**This requires modifying kernel.cpp** - not a simple setup.cpp change.

---

### 6. Sponge/Buffer Zone (IMPORTANT)

**Research Recommendation:**
- Damping zone near outlet: 5-10h length
- Quadratic damping profile: σ ~ x²
- Gradually relax solution toward ambient state

**FluidX3D Status:**
- **Not available** in standard implementation
- Must be implemented in run() loop

**Gap Assessment:** MEDIUM (requires kernel modification)

**CRITICAL:** Sponge zone must be implemented in **OpenCL kernel** to maintain GPU performance. CPU-side implementation would kill simulation speed.

**Required Kernel Modification Approach:**

In `kernel.cpp` stream_collide(), add sponge damping before collision:

```cpp
// Pseudocode for kernel modification
// After calculating rhon, uxn, uyn, uzn but before collision

// Sponge zone parameters (could be passed as kernel arguments or defines)
const uint sponge_start = def_Nx - def_sponge_length;

if(x >= sponge_start) {
    float xi = (float)(x - sponge_start) / (float)def_sponge_length;
    float sigma = xi * xi; // quadratic damping profile

    // Relax toward ambient state
    uxn *= (1.0f - sigma);
    uyn *= (1.0f - sigma);
    uzn *= (1.0f - sigma);
    rhon = (1.0f - sigma) * rhon + sigma * 1.0f;
}

// Then proceed with collision using damped values
```

**Implementation requires:**
1. Add sponge zone defines in `defines.hpp` (start position, length)
2. Modify `stream_collide()` kernel in `kernel.cpp`
3. Apply damping to macroscopic variables before computing f_eq

---

### 7. Side/Farfield Boundary Conditions

**Research Recommendation:**
- Equilibrium BC at atmospheric pressure
- Allow entrainment (inflow of ambient air)
- Domain: ±15-20h from centerline

**FluidX3D Status:**
- TYPE_E can be used with ambient conditions
- Default boundary behavior depends on domain setup

**Gap Assessment:** NONE
- TYPE_E with ρ = 1.0, u = 0 provides appropriate farfield
- **Action:** Set TYPE_E on y-boundaries with ambient conditions

---

### 8. Spanwise (Periodic) Boundary

**Research Recommendation:**
- Periodic BC acceptable for full-width simulation
- Full width: 190mm (19h)

**FluidX3D Status:**
- Periodic boundaries are default behavior in FluidX3D for edges not explicitly set
- Full-width domain is feasible

**Gap Assessment:** NONE

---

## Implementation Strategy

### Phase 1: Kernel Modifications (REQUIRED FIRST)
1. **Modify `kernel.cpp`** to add:
   - Outlet extrapolation BC (for TYPE_E at outlet)
   - Sponge zone damping before collision
   - (Optional) Inlet turbulence perturbation
2. Add necessary defines in `defines.hpp` for sponge zone parameters
3. All modifications must be GPU-side (OpenCL kernel code)

### Phase 2: Simulation Setup
1. Enable `SUBGRID` and `EQUILIBRIUM_BOUNDARIES` in defines.hpp
2. Create rectangular nozzle geometry with domain size ~50h × 40h × 19h
3. Set TYPE_E inlet with top-hat velocity profile (jet exit)
4. Set TYPE_E outlet (kernel handles extrapolation)
5. Set TYPE_E farfield boundaries (y-direction) with ambient conditions
6. Leave z-direction as periodic (default)

### Phase 3: Systematic Validation
1. Monitor key parameters during simulation:
   - Centerline velocity decay U_c(x)
   - Mass conservation (inlet vs outlet flux)
   - Turbulence intensity profiles
2. Compare potential core length with experiment (P_c = 3.85h)
3. Compare spreading rate (target SR_y = 0.1166)
4. If mismatch: systematic debugging, not trial-and-error

### Phase 4: Refinement Based on Physics
1. If potential core too long → check inlet turbulence level
2. If spreading rate wrong → check SGS model, grid resolution
3. Document root cause of any discrepancies

---

## Required Code Changes Summary

### defines.hpp
```cpp
// Uncomment these lines:
#define SUBGRID                    // Enable Smagorinsky-Lilly LES
#define EQUILIBRIUM_BOUNDARIES     // Enable TYPE_E for inlet/outlet

// Comment out:
//#define BENCHMARK                // Disable benchmark mode
```

### kernel.cpp (GPU-Side Modifications - CRITICAL)
Must implement in OpenCL kernel code:
1. **Outlet Extrapolation BC**: In `stream_collide()`, extrapolate velocity from interior for TYPE_E cells at x=Nx-1
2. **Sponge Zone**: Apply damping to macroscopic variables before collision in sponge region
3. **(Optional) Inlet Turbulence**: GPU-side random perturbation at inlet

### setup.cpp (New simulation function)
Key components:
1. Domain creation with appropriate size (~50h × 40h × 19h)
2. Nozzle geometry definition (rectangular slot at inlet)
3. TYPE_E boundary assignment:
   - Inlet (x=0): TYPE_E with jet velocity profile
   - Outlet (x=Nx-1): TYPE_E (kernel handles extrapolation)
   - Farfield (y=0, y=Ny-1): TYPE_E with u=0, rho=1
   - Spanwise (z): Periodic (default, no flag needed)
4. Initial velocity profile (top-hat with tanh smoothing at shear layer edges)
5. **NO per-timestep CPU operations** - all BC handling in kernel

---

## Risk Assessment

| Risk | Severity | Mitigation |
|------|----------|------------|
| **Outlet acts as wall** | HIGH | Implement extrapolation BC in kernel.cpp - CRITICAL |
| **Performance degradation** | HIGH | All custom BCs must be GPU-side (kernel), never CPU-side |
| Insufficient inlet turbulence | HIGH | GPU-side perturbation in kernel; validate against P_c |
| Outlet reflections | MEDIUM | Implement sponge zone in kernel + extrapolation BC |
| Numerical instability at Re=20,000 | MEDIUM | Use TRT + SUBGRID, monitor divergence |
| Memory limitation (large grid) | LOW | Use FP16S, multi-GPU if available |
| Inaccurate spreading rate | MEDIUM | Validate against experiment, adjust parameters |

---

## Conclusion

FluidX3D provides a solid foundation for this simulation with:
- Appropriate SGS model (Smagorinsky-Lilly)
- Adequate collision operators (SRT/TRT)
- Equilibrium boundary conditions
- Memory optimization (FP16)
- Multi-GPU capability

**Primary gaps requiring custom implementation:**
1. **Outlet BC with extrapolation** - CRITICAL: TYPE_E with u=0 acts as wall; must implement extrapolation in `kernel.cpp`
2. **Sponge zone** - Important for preventing outlet reflections; must implement in `kernel.cpp`
3. **Inlet turbulence generation** - GPU-side implementation required for proper jet development

**Estimated implementation complexity:** MODERATE-HIGH
- Requires OpenCL kernel modifications in `kernel.cpp` (not just setup.cpp)
- All implementations must be GPU-side to maintain performance
- Scientific validation methodology: systematic debugging, no trial-and-error

---

*Document created: January 2026*
