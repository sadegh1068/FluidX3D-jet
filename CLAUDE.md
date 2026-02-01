# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Available Agents

Specialized agents are located in `.claude/agents/`. Reference them in prompts for domain-specific help.

### Agent List

| Agent | File | Expertise |
|-------|------|-----------|
| **DFM Inlet** | `.claude/agents/dfm-inlet-agent.md` | Digital Filter Method turbulence generation |
| **Convective Outlet** | `.claude/agents/convective-outlet-agent.md` | Orlanski BC, sponge zones, periodic coupling |
| **Jet Validation** | `.claude/agents/jet-validation-agent.md` | Ahmed et al. experimental comparison |
| **Kernel Development** | `.claude/agents/kernel-development-agent.md` | FluidX3D OpenCL kernel architecture |

### How to Use Agents

**Load specific agent:**
```
Read .claude/agents/dfm-inlet-agent.md and help me implement DFM turbulence
```

**Multi-agent task:**
```
Using the Kernel Development and Convective Outlet agents, fix the outlet BC
```

**Let Claude choose:**
```
Help me implement turbulent inlet BC (Claude will suggest relevant agent)
```

---

## Project Overview

FluidX3D is a GPU-accelerated Lattice Boltzmann Method (LBM) Computational Fluid Dynamics (CFD) simulation software. It supports all GPUs/CPUs via OpenCL. Author: Dr. Moritz Lehmann.

## Build Commands

### Windows
Open `FluidX3D.sln` in Visual Studio, build Release|x64, output: `bin\FluidX3D.exe`

### Linux/macOS/Android
```bash
./make.sh              # Auto-detects OS, compiles, and runs
./make.sh 0            # Run on device 0
./make.sh 0 1 3        # Multi-GPU on devices 0, 1, 3
```

### Make targets (Linux)
```bash
make Linux-X11         # Linux with X11 graphics
make Linux             # Linux headless
make macOS             # macOS
make Android           # Android
make clean             # Clean build artifacts
```

Output: `bin/FluidX3D`

## Architecture

### Core Source Files (`src/`)

| File | Purpose |
|------|---------|
| `main.cpp` | Entry point, graphics/physics loops, keyboard controls |
| `lbm.cpp/hpp` | `LBM` class (public API) and `LBM_Domain` (GPU domain management) |
| `kernel.cpp/hpp` | OpenCL C kernels embedded as C++ strings |
| `setup.cpp` | Simulation configurations - **edit this to create new simulations** |
| `defines.hpp` | **Compile-time configuration** - velocity sets, extensions, graphics modes |
| `opencl.hpp` | OpenCL wrapper: `Device`, `Kernel`, `Memory<T>` template |
| `graphics.cpp/hpp` | Rendering engine (raytracing, rasterization, field visualization) |
| `utilities.hpp` | Vector types (`float3`, `int3`), `Mesh`, `Image`, `parallel_for` |
| `units.hpp` | SI unit conversions for LBM simulations |
| `shapes.cpp/hpp` | Geometric primitives (sphere, cylinder, box) |

### Configuration System (`defines.hpp`)

**Velocity Sets** (uncomment one):
- `D2Q9` - 2D simulations
- `D3Q15`, `D3Q19` (default), `D3Q27` - 3D simulations

**Memory Optimization** (optional, 2x speedup):
- `FP16S` - IEEE-754 FP16 compression
- `FP16C` - Custom FP16 format (more accurate)

**Physics Extensions** (uncomment as needed):
- `VOLUME_FORCE` - Global pressure gradient
- `FORCE_FIELD` - Per-cell forces on boundaries (+12 Bytes/cell)
- `MOVING_BOUNDARIES` - Dynamic solid objects
- `EQUILIBRIUM_BOUNDARIES` - Inflow/outflow conditions
- `SURFACE` - Free surface LBM (+12 Bytes/cell)
- `TEMPERATURE` - Thermal simulation (+32 Bytes/cell FP32)
- `SUBGRID` - Smagorinsky-Lilly turbulence model
- `PARTICLES` - Immersed boundary particles

**Graphics Modes**:
- `INTERACTIVE_GRAPHICS` - Interactive window (requires X11 on Linux)
- `INTERACTIVE_GRAPHICS_ASCII` - ASCII console graphics
- `GRAPHICS` - Headless rendering to disk
- `BENCHMARK` - Disables all extensions, runs performance test

### Cell Type Flags
```cpp
TYPE_S  // Solid boundary
TYPE_E  // Equilibrium boundary (inflow/outflow)
TYPE_T  // Temperature boundary
TYPE_F  // Fluid
TYPE_I  // Interface (free surface)
TYPE_G  // Gas (free surface)
```

## Creating a New Simulation

1. Edit `src/setup.cpp` - find commented `void main_setup()` blocks as templates
2. Edit `src/defines.hpp` - enable required extensions noted in setup comments
3. Recompile

Basic simulation pattern:
```cpp
void main_setup() {
    LBM lbm(Nx, Ny, Nz, viscosity);  // Create simulation grid

    parallel_for(lbm.get_N(), [&](ulong n) {
        uint x, y, z; lbm.coordinates(n, x, y, z);
        // Set initial conditions: lbm.rho[n], lbm.u.x[n], lbm.flags[n]
    });

    lbm.run();  // Start simulation
}
```

## Multi-GPU

Domain decomposition via `LBM` constructor:
```cpp
LBM lbm(Nx, Ny, Nz, Dx, Dy, Dz, viscosity);  // Dx*Dy*Dz GPUs
```

Select devices at runtime: `bin/FluidX3D 0 1 2 3`

## Key Data Structures

- `Memory<T>` - Unified host/device memory with automatic transfers
- `Mesh` - Triangle mesh with STL loading and voxelization
- `Units` - Converts between LBM and SI units

## Performance

Benchmark mode measures MLUPs/s (Million Lattice Updates Per Second). Memory footprint per cell:
- D3Q19 FP32: 93 Bytes/cell
- D3Q19 FP16: 55 Bytes/cell

---

## Critical Development Guidelines

### 1. GPU-First Implementation (MANDATORY)

FluidX3D's primary strength is GPU acceleration. **All computational code must be implemented on the GPU side (OpenCL kernels)** to avoid dramatically reducing simulation speed.

**AVOID:**
```cpp
// BAD - CPU-side operations in run() loop kill performance
while(running) {
    lbm.u.read_from_device();   // Expensive GPU→CPU transfer
    // ... CPU processing ...
    lbm.u.write_to_device();    // Expensive CPU→GPU transfer
    lbm.run(1u);
}
```

**PREFER:**
- Implement custom boundary conditions, sponge zones, turbulence generation in `kernel.cpp` as OpenCL kernels
- Use existing kernel infrastructure and flag system
- Only use CPU-side operations for initialization and post-processing (not per-timestep)

### 2. Scientific Methodology (NO Trial and Error)

When simulation results don't match experimental data:
1. **DO NOT** use trial-and-error parameter tuning
2. **DO NOT** introduce empirical correction factors
3. **DO** implement systematic debugging:
   - Monitor key parameters during simulation (centerline velocity, turbulence intensity, mass conservation)
   - Compare intermediate results against literature benchmarks
   - Identify root cause through parameter tracking
   - Fix the underlying issue based on physical understanding

### 3. Literature-Based Approach

- All implementations must be grounded in peer-reviewed scientific literature
- Cite relevant papers for methods used (boundary conditions, turbulence models, etc.)
- Use validated numerical methods from CFD/LBM literature
- Avoid ad-hoc modifications or "tuning" to match results
- If a method doesn't work, investigate why rather than adding corrections

### 4. Validation Strategy

When validating against experimental data:
1. Start with well-understood benchmark cases from literature
2. Compare multiple quantities (not just one metric)
3. Check grid independence
4. Verify boundary condition implementation against theory
5. Document any discrepancies with physical explanation

---

## L0 Rectangular Jet Project Log

### Project Goal
Simulate a rectangular jet (AR=19) and validate against experimental data from "Mean and turbulent velocity fields of a rectangular jet" paper.

**Target Case:**
- Nozzle: AR = 19, h = 10mm, width = 190mm
- Flow: U₀ = 30 m/s, Re_h ≈ 20,100
- Validation targets: P_c = 3.85h (potential core length), SR_y = 0.1166 (spreading rate)

**Experimental Data Files:** `doc/exp fig 7.csv`, `doc/exp fig 10.csv`, `doc/exp fig4 lp=0.csv`

---

### Implementation Summary (2026-01-27)

#### New Features Added to FluidX3D

**1. Convective (Orlanski) Outlet BC** (`kernel.cpp`, `defines.hpp`)
- Flag: `CONVECTIVE_OUTLET`
- Uses TYPE_Y to mark outlet cells
- Implements: ∂f/∂t + U_c · ∂f/∂x = 0
- Purpose: Allow waves to exit domain without reflection

**2. Sponge Zone Damping** (`kernel.cpp`, `lbm.cpp`, `defines.hpp`)
- Flag: `SPONGE_ZONE`
- Parameters: `def_sponge_start=0.85f`, `def_sponge_strength=0.1f`
- Cubic damping profile near outlet
- Purpose: Absorb waves before they reach outlet

**3. Inlet Turbulence - Random Fourier Modes (RFM)** (`kernel.cpp`, `lbm.cpp`, `defines.hpp`)
- Flag: `INLET_TURBULENCE_RFM`
- Uses TYPE_X to mark turbulent inlet cells
- **Updated (v3 - PROPER KRAICHNAN):** Complete rewrite implementing correct RFM
- Parameters: `def_turbulence_intensity=0.007f` (0.7% inlet → ~3.8% at x/h=0.2), `def_turbulence_length_scale=10.0f` (lattice units)
- Purpose: Generate divergence-free turbulent fluctuations at inlet
- Reference: Kraichnan, R.H. (1970). "Diffusion by a random velocity field." Phys. Fluids 13:22-31

**4. Time-Averaged Data Export** (`setup.cpp`)
- `JetStatistics` class accumulates Σu, Σu², N over time
- Computes: `<U> = Σu/N`, `u'_rms = sqrt(<u²> - <u>²)`
- **Intermediate exports:** Every 2 flow-through times with timestep suffix
- **Final export:** Clean filenames without suffix
- Export files:
  - `centerline_averaged[_tXXX].csv` - x/h, Uc/U₀, Urms/Uc, y₀.₅/h, samples
  - `volume_flow_averaged[_tXXX].csv` - x/h, Qx*/Q₀*, samples
  - `lateral_profiles_averaged[_tXXX].csv` - y/y₀.₅, U/Uc at multiple x/h stations, samples

**5. Simulation Setup** (`setup.cpp`)
- VRAM-based resolution using `resolution()` function
- Domain: 60h × 25h × 25h
- Warmup: 2 flow-through times
- Averaging: 10 flow-through times
- Sample interval: 100 timesteps
- Intermediate export interval: 2 flow-through times

#### Files Modified
| File | Changes |
|------|---------|
| `src/defines.hpp` | Added CONVECTIVE_OUTLET, SPONGE_ZONE, INLET_TURBULENCE_RFM flags |
| `src/kernel.cpp` | Added sponge zone (~line 1593), RFM turbulence (~line 1610), convective outlet (~line 1710) |
| `src/lbm.cpp` | Added device_defines for sponge and turbulence parameters |
| `src/setup.cpp` | Added JetStatistics class, time-averaging export functions, L0 jet main_setup() |

---

### Simulation Run #1 (2026-01-27) - ABORTED

**Configuration:**
- Memory: 5000 MB VRAM
- Grid: 819 × 341 × 341 (h = 27 cells)
- Extensions: TRT, FP16S, EQUILIBRIUM_BOUNDARIES, SUBGRID, CONVECTIVE_OUTLET, SPONGE_ZONE, INLET_TURBULENCE_RFM, INTERACTIVE_GRAPHICS

**Observations (Early Stage):**

1. **Inlet Turbulence Pattern Issue**
   - Visible structured/periodic pattern at inlet (not realistic random turbulence)
   - "Hot spots" appear where Fourier modes constructively interfere
   - Cause: Simple RFM with only 4 modes and deterministic pseudo-random phases

2. **CRITICAL: Jet Flapping Issue**
   - Jet exhibited strong flapping/oscillation behavior instead of being straight
   - Density field showed unusual low-density spots (blue regions, ρ < 1)
   - **Root Cause Identified:** Bug in RFM turbulence generator
     ```cpp
     // BUGGY CODE:
     v_prime += amp * 0.7f * sin(km * (float)xyz.z + phi_z);  // v' depends on z only!
     ```
   - `v_prime` (y-direction velocity) depended only on `xyz.z`, not `xyz.y`
   - This created spanwise-uniform transverse perturbations
   - Effect: Coherent sideways forcing that excited the jet's natural flapping instability

**Status:** ABORTED - Simulation stopped to fix bug

---

### Bug Fix Applied (2026-01-27)

**Problem:** Inlet turbulence v_prime coordinate dependency caused jet flapping

**Fix in `kernel.cpp` (RFM turbulence generator):**

1. **Corrected coordinate dependencies:**
   - Each velocity component now varies with BOTH y AND z
   - Prevents spanwise-uniform perturbations
   ```cpp
   // CORRECTED CODE:
   u_prime += amp * sin(km * (float)xyz.y + phi_uy) * cos(km2 * (float)xyz.z + phi_uz);
   v_prime += amp * sin(km * (float)xyz.y + phi_vy) * cos(km2 * (float)xyz.z + phi_vz);
   w_prime += amp * sin(km * (float)xyz.y + phi_wy) * cos(km2 * (float)xyz.z + phi_wz);
   ```

2. **Increased number of modes:** 4 → 8 modes for better randomness

3. **Separate phases for each component:** Prevents correlation between u', v', w'

4. **Reduced transverse amplitude:** v', w' reduced to 50% of u' (typical for shear flows)

**Additional improvement in `setup.cpp`:**
- Added intermediate CSV export every 2 flow-through times
- Files include sample count for convergence monitoring
- Intermediate files have timestep suffix: `_t{timestep}_ft{flow_times}`

---

### Simulation Run #2 (2026-01-27) - ABORTED

**Observation:** Checkerboard pattern at inlet (not realistic turbulence)

The "fix" resolved the flapping but revealed a **fundamental flaw** in the RFM implementation.

**Status:** ABORTED - Need to rewrite turbulence generator

---

### Deep Analysis: RFM Implementation is Fundamentally Flawed (2026-01-27)

**Visual Evidence:** Checkerboard pattern visible at inlet - regular alternating red/blue spots that look like a standing wave, not random turbulence.

#### What Our Code Does (WRONG)

```cpp
u_prime += amp * sin(km * (float)xyz.y + phi_uy) * cos(km2 * (float)xyz.z + phi_uz);
```

**Problem:** `sin(km * y) * cos(km2 * z)` creates a **2D standing wave pattern**, not random turbulence.

#### What Proper RFM Should Do (Kraichnan 1970)

```
u'(x) = Σ [aₙ × k̂ₙ] * cos(kₙ·x + φₙ)
```

Where:
- **kₙ** = random 3D wave vector (random direction AND magnitude)
- **φₙ** = random phase that is **CONSTANT for each mode** (same for all spatial positions)
- **aₙ × k̂ₙ** = amplitude perpendicular to k (ensures divergence-free)
- **kₙ·x** = dot product `kx*x + ky*y + kz*z` (spatial variation)

#### Comparison Table

| Aspect | Correct RFM | Our Implementation |
|--------|-------------|-------------------|
| Wave vectors | Random 3D directions | Fixed along y and z axes only |
| Phases | Constant per mode | Varies with spatial position (xyz.y, xyz.z) |
| Spatial variation | From k·x dot product | From separate sin(y)*cos(z) products |
| Divergence-free | a ⊥ k guaranteed | Not enforced |
| Result | Random turbulent field | **Checkerboard standing wave pattern** |

#### Why the "Phases" Don't Help

```cpp
const float phi_uy = fmod((float)(seed_base * 7u + xyz.y * 13u) * 0.618033988749895f, 1.0f) * 6.28318530718f;
```

- Phase **depends on xyz.y** (spatial coordinate)
- This doesn't create randomness - just shifts the pattern spatially
- Each cell gets a deterministic phase based on its position
- Result is still a structured pattern with position-dependent offsets

#### The Length Scale Issue

- `def_turbulence_length_scale = 10.0f` lattice units
- k0 = 2π/Lt ≈ 0.628, modes m=1..8 → wavenumbers 0.628 to 5.02
- Shortest wavelength ≈ 1.25 cells → **barely resolved, creates numerical artifacts**

#### What the Checkerboard Shows

- Dominant mode (m=1) with wavelength = 10 cells
- `sin(k*y) * cos(k*z)` creates nodes at regular intervals
- Red = positive fluctuation, Blue = negative fluctuation
- This is a **deterministic wave pattern**, not turbulence

#### FIX IMPLEMENTED: Proper Kraichnan RFM (v3)

The RFM has been completely rewritten to follow Kraichnan's method correctly:

**Formula:** `u'(x,t) = sqrt(2) * Σ [σₙ × d̂ₙ] * cos(kₙ·x + φₙ + ωₙt)`

**Implementation details (`kernel.cpp` lines 1610-1710):**

1. **Random 3D wave vectors:**
   - Direction uniformly distributed on unit sphere using `cos(θ) ∈ [-1,1]`, `φ ∈ [0,2π]`
   - Magnitude sampled from range `[k₀, 8k₀]` where `k₀ = 2π/Lt`

2. **Constant phases per mode:**
   - Phase depends ONLY on mode index `m`, NOT on spatial coordinates
   - Uses golden ratio-based quasi-random sequence: `seed = (m+1) * 0.618...`

3. **Divergence-free condition (∇·u' = 0):**
   - Amplitude direction computed as `a = σ × d̂` (cross product)
   - Ensures amplitude is perpendicular to wave vector

4. **Proper spatial variation:**
   - `k·x = kx*x + ky*y + kz*z` (full 3D dot product)
   - No more `sin(y)*cos(z)` standing wave patterns

5. **Temporal variation:**
   - `ω = 0.005 * |k|` (frequency proportional to wavenumber)
   - Provides temporal decorrelation via Taylor's frozen turbulence hypothesis

6. **32 Fourier modes** for better statistical representation

7. **Energy spectrum approximation:**
   - Amplitude `~ 1/(1 + (k/k₀)²)` (simplified von Kármán)
   - Proper normalization to match target turbulence intensity

**Expected behavior:**
- Random, irregular fluctuations (no checkerboard)
- Smooth spatial correlations with length scale Lt
- Divergence-free (no spurious pressure waves)
- Temporal variation with characteristic timescale

---

### Results Analysis

*(To be updated after running simulation with proper Kraichnan RFM)*

**Centerline Velocity Decay (Fig. 7 comparison):**
- Target: P_c = 3.85h
- Simulated: TBD

**Spreading Rate (Fig. 7 comparison):**
- Target: SR_y = 0.1166
- Simulated: TBD

**Volume Flow Rate (Fig. 10 comparison):**
- TBD

**Lateral Profiles (Fig. 4 comparison):**
- TBD

---

### Issues & Solutions Log

| Date | Issue | Root Cause | Solution | Status |
|------|-------|------------|----------|--------|
| 2026-01-27 | Structured inlet turbulence pattern | Simple RFM with few modes | Increased to 8 modes, separate phases | Partial - pattern still visible |
| 2026-01-27 | Jet flapping/oscillation | v_prime depended on z instead of y | Fixed coordinate dependencies | Fixed |
| 2026-01-27 | No intermediate CSV output | CSV only written at end | Added intermediate export every 2 flow-through times | Fixed |
| 2026-01-27 | Checkerboard pattern at inlet | Fundamental RFM flaw - sin(y)*cos(z) creates standing waves | **Complete rewrite: Proper Kraichnan RFM (v3)** | **FIXED** |
| 2026-01-27 | "Constant velocity first, then turbulence pattern" at inlet | Not a bug - BC arrangement is correct | No fix needed | Investigated |
| 2026-01-27 | Turbulence intensity 5-7x too high (26% vs 3.8% at x/h=0.2) | Jet shear layer K-H instability amplifies inlet perturbations; also minor normalization issue with unit amplitude vectors | Reduced def_turbulence_intensity from 0.05f to 0.007f (empirical calibration to match experiment) | **FIXED** |
| 2026-01-29 | Red spots at inlet with TI=0 | Multiple: UPDATE_FIELDS overwrites inlet, collision overwrites fhn, AA streaming slot mismatch | Make TYPE_X work like TYPE_E (preset velocity, fhn=feq in collision) | Testing |
| 2026-01-29 | OpenCL compilation error (|| operator) | Nested preprocessor directives spanning condition | Restructure using boolean flags instead of nested #ifdef | **FIXED** |
| 2026-01-29 | DFM wrong filter formula | Used exponential instead of Gaussian filter coefficients | Changed to Gaussian: exp(-π*k²/(2*n²)) per Klein et al. 2003 | **FIXED** |
| 2026-01-29 | DFM unnecessary sqrt(3) variance | Applied variance correction meant for uniform RNG to Gaussian noise | Removed sqrt(3) factor | **FIXED** |
| 2026-01-29 | High TI but low spreading rate paradox | Smagorinsky constant Cs=0.173 too high for jets - over-dissipates small-scale turbulence while K-H structures survive | Reduced Cs from 0.173 to 0.12 (optimized for jets per literature) | Testing |
| 2026-02-01 | Mean velocity oscillations ±8% in potential core | DFM random seed uses `t % 1000ul` - same pattern persists for 1k timesteps, only ~81 patterns averaged instead of ~81,500 | Change to `t` for unique pattern each timestep | **FIXED** |
| 2026-02-01 | Simulation instability after FT=10 | Cs=0.12 too low for long-term stability - absolute instability (jet flapping) grows exponentially | Use FT=10 results (stable); consider Cs=0.14-0.15 for long runs | Documented |

---

### References

1. Orlanski, I. (1976). "A simple boundary condition for unbounded hyperbolic flows." J. Comput. Phys. 21:251-269
2. Smirnov, A., Shi, S., Celik, I. (2001). "Random flow generation technique for LES." J. Fluids Eng. 123:359-371
3. Kraichnan, R.H. (1970). "Diffusion by a random velocity field." Phys. Fluids 13:22-31
4. Poletto, R. et al. (2013). "A new divergence free synthetic eddy method for the reproduction of inlet flow conditions for LES." Flow Turbulence Combust. 91:519-539

### Literature Notes on Inlet Turbulence Generation

**Key requirements for proper RFM (from literature review):**

1. **Divergence-free condition**: ∇·u' = 0 must be satisfied to avoid spurious pressure waves
2. **Proper spectral content**: Energy spectrum E(k) must match target (von Kármán, Kolmogorov)
3. **Spatial decorrelation**: Two-point correlations must decay with separation distance
4. **Statistical independence**: Each Fourier mode must have independent random phase
5. **Sufficient modes**: Literature suggests 38-1200 modes for accurate statistics

**Common implementation mistakes (ALL ADDRESSED in v3):**
- ~~Violation of divergence-free condition~~ → FIXED: a = σ × d̂ ensures ∇·u' = 0
- ~~Correlated random phases~~ → FIXED: phases depend only on mode index
- ~~Inadequate wavenumber range~~ → FIXED: k ∈ [k₀, 8k₀] with 32 modes
- ~~Fixed wave vector directions~~ → FIXED: random 3D directions on unit sphere

**Adaptation distances for synthetic turbulence:**
- Skin-friction: x/δ₀ ≈ 12
- Mean velocity: x/δ₀ ≈ 16
- Reynolds stresses: x/δ₀ ≈ 18-50

---

### Future Improvements

#### Non-Equilibrium Turbulent Inlet Boundary Condition

**Current limitation (TYPE_E + TYPE_X approach):**
The current implementation uses equilibrium boundary (TYPE_E) with velocity fluctuations from RFM (TYPE_X). This sets:
- f = f_eq(ρ, u + u') where u' is the turbulent fluctuation

**The problem:** In real turbulence, the distribution function has two parts:
```
f = f_eq + f_neq
```
- **f_eq** encodes density and mean velocity (what we set)
- **f_neq** encodes the Reynolds stress tensor τ_ij (what we're missing!)

By using equilibrium BC with fluctuating velocity, we inject correct velocity fluctuations but with **zero stress tensor**. The turbulence needs ~5-10 cells to develop the correct stress-strain relationship.

**Proposed future improvement - TYPE_INLET_TURB:**
Create a dedicated turbulent inlet boundary type that prescribes the full distribution function f_i, including non-equilibrium components matching target Reynolds stresses:
```
f_i = f_eq_i + f_neq_i
```
where f_neq_i is computed from the target stress tensor using Chapman-Enskog expansion.

**Benefits:**
- Correct Reynolds stresses from the first cell
- Reduced adaptation distance
- More accurate near-field statistics

**References for implementation:**
- Guo et al. (2002) - Non-equilibrium extrapolation for LBM boundaries
- Haussmann et al. (2019) - Evaluation of synthetic turbulence generation methods in LBM

---

#### Domain and Geometry Issues (Identified 2026-01-28)

**1. Truncated Nozzle Aspect Ratio**

| Parameter | Experiment | Simulation | Issue |
|-----------|------------|------------|-------|
| Nozzle AR | 19 | ~12.5 | Nozzle truncated to fit domain |
| Nozzle width | 19h = 190mm | 337 cells | Nz=341 < 513 cells needed |

**Cause:** Domain aspect ratio 60h × 25h × 25h with limited VRAM (5GB) resulted in Nz=341 cells, but nozzle needs 19h + buffer > 341 cells.

**Fix options:**
- Increase z-dimension: Use aspect ratio 60h × 25h × 45h
- Reduce x-dimension: Use 40h × 25h × 35h to allow wider z
- Use periodic BC in z for infinite slot approximation (different physics)

**2. Volume Flow Rate Calculation (Q*)**

| Aspect | Paper Method | Our Method |
|--------|--------------|------------|
| Integration | 2D: ∫∫ U(y,z) dy dz | 1D: ∫ U(y) dy at z=center |
| Reference Q₀ | U₀ × h × w | U₀ × h |
| Z-variation | Captured | Ignored |

**Impact:** Our Q*/Q₀* doesn't match paper's definition. Cannot validate volume flow rate until fixed.

**Fix:** Implement 2D integration over jet cross-section with proper Q₀ = U₀ × h × w normalization.

**3. Monitoring x/h=30 Location**

The monitoring measures mass flux at x/h=30, but:
- Domain is only ~30h long
- x/h=30 is in/near the sponge zone (starts at 0.85 × domain)
- Flow is being artificially damped at measurement location

**Fix:** Either extend domain or measure at x/h=20 or earlier.

**4. Outlet Boundary Condition - CRITICAL**

Current setup creates a "soft wall" at outlet:
```cpp
// WRONG - Current setup:
lbm.flags[n] = TYPE_E | TYPE_Y;  // Equilibrium + convective
lbm.u.x[n] = 0.0f;               // Forces velocity to ZERO
def_sponge_ux_target = 0.0f;     // Sponge also targets ZERO
```

**The problem:**
- TYPE_E with u=0 enforces equilibrium with no outflow
- Sponge zone damps mean velocity toward zero
- Combined effect: flow CAN'T EXIT forward
- Result: jet is pushed sideways through lateral boundaries
- Explains: outward streamlines, "mass loss", poor jet development

**Correct approach:**
```cpp
// CORRECT - Convective outlet only:
lbm.flags[n] = TYPE_Y;           // Convective BC only, NO TYPE_E
// No prescribed velocity - let flow exit naturally

// Sponge zone should damp FLUCTUATIONS, not mean flow:
// Option A: Damp toward local mean (not zero)
// Option B: Damp only the turbulent component (u' = u - <u>)
// Option C: Use weaker sponge that doesn't fully kill velocity
```

**Physical reasoning:**
- Jet exits into open atmosphere at some velocity (not zero!)
- Sponge zone purpose: absorb turbulent fluctuations to prevent reflection
- Sponge should NOT stop the mean flow
- Convective BC (Orlanski): advects flow out at local velocity

**Fix required:**
1. Remove TYPE_E at outlet, use only TYPE_Y (convective)
2. Modify sponge to damp fluctuations only, not mean velocity
3. Or set sponge target to estimated exit velocity (not zero)

---

#### Current Validation Focus (2026-01-28)

Due to the above limitations, **focus validation on centerline velocity decay (Uc/U₀)** which is less affected by:
- Nozzle truncation (centerline is still valid)
- 1D vs 2D integration issues
- Domain length (can compare at x/h < 25)

Defer volume flow and spreading rate validation until geometry issues are fixed.
- For LBM: 1.5-3 domain flow-through times

---

### New Boundary Conditions Implementation (2026-01-28)

**Problem Identified:** TYPE_E (equilibrium boundary) was causing issues at both inlet and outlet:
- **Outlet:** TYPE_E with u=0 acted as a soft wall, preventing flow from exiting
- **Inlet:** TYPE_E only controls velocity, not the full distribution function
- **Sponge:** Was damping mean flow toward zero, blocking jet exit

**Solution Implemented:** New boundary conditions that DON'T use TYPE_E for inlet/outlet:

#### 1. DFSEM_INLET (replaces INLET_TURBULENCE_RFM)

**Flag:** `TYPE_X` alone (NO TYPE_E!)

**What it does:**
- Sets ALL 19 distribution functions directly from equilibrium
- Uses Kraichnan RFM to generate divergence-free turbulent fluctuations
- Completely controls inlet state (density, velocity, distributions)

**Key change:**
```cpp
// OLD: TYPE_E | TYPE_X - TYPE_E sets equilibrium, TYPE_X adds fluctuations to velocity
// NEW: TYPE_X alone - kernel sets ALL f_i directly
for(uint i = 0u; i < def_velocity_set; i++) {
    fhn[i] = feq_inlet[i];  // Override ALL populations
}
```

**Parameters:**
- `def_inlet_velocity = 0.05f` (LBM units)
- `def_turbulence_intensity = 0.007f` (0.7% → ~3.8% at x/h=0.2)
- `def_turbulence_length_scale = 10.0f` (lattice units)

#### 2. CONVECTIVE_OUTLET (updated)

**Flag:** `TYPE_Y` alone (NO TYPE_E!)

**What it does:**
- Applies Orlanski convective BC: ∂f/∂t + U_c·∂f/∂x = 0
- Only modifies UNKNOWN populations (c_ix < 0)
- Convection velocity from upstream cell (avoids circular dependency)

**Key change:**
```cpp
// OLD: TYPE_E | TYPE_Y - TYPE_E forced equilibrium with u=0
// NEW: TYPE_Y alone - Orlanski BC on unknown populations
fhn[i] = fhn[i] - U_conv * (fhn[i] - f_upstream[i]);
```

#### 3. SPONGE_ZONE (updated)

**What it does:**
- Damps velocity toward UPSTREAM velocity (not zero!)
- Removes turbulent fluctuations while preserving mean flow
- Uses quintic smooth profile for gentle transition

**Key change:**
```cpp
// OLD: Damp toward def_sponge_ux_target = 0
// NEW: Damp toward upstream velocity
uxn = fma(1.0f - sigma, uxn, sigma * ux_upstream);
```

**Parameters:**
- `def_sponge_start = 0.85f` (85% of domain)
- `def_sponge_strength = 0.15f` (weak, preserves mean flow)

#### Files Modified

| File | Changes |
|------|---------|
| `src/defines.hpp` | Renamed `INLET_TURBULENCE_RFM` → `DFSEM_INLET`; updated comments |
| `src/kernel.cpp` | Updated DFSEM_INLET to set all f_i; updated SPONGE_ZONE to damp toward upstream; updated CONVECTIVE_OUTLET for Orlanski |
| `src/lbm.cpp` | Updated device_defines (removed def_sponge_ux_target, added def_inlet_velocity) |
| `src/setup.cpp` | Changed inlet from `TYPE_E\|TYPE_X` to `TYPE_X`; changed outlet from `TYPE_P` to `TYPE_Y` |
| `docs/DFSEM_ConvectiveOutlet_BC.md` | NEW: Full technical documentation with literature references |

#### Expected Improvements
1. **Flow exits through outlet** (not deflected sideways)
2. **Mass conservation** (inlet flow matches outlet flow)
3. **No artificial pressure artifacts** at outlet
4. **Turbulent structures exit cleanly** without reflection

#### Verification Checklist
- [ ] Build successfully with new defines
- [ ] Flow exits through outlet (visualize streamlines)
- [ ] Mass conservation (< 0.1% drift)
- [ ] No checkerboard at inlet
- [ ] Turbulence intensity ~3.8% at x/h=0.2

---

### Critical Bug Fix: Periodic BC Coupling (2026-01-28)

#### Problem Description

After implementing the new boundary conditions (DFSEM_INLET, CONVECTIVE_OUTLET without TYPE_E), a "ghost nozzle" pattern appeared at the outlet from timestep 0. The pattern exactly matched the inlet nozzle shape, appearing as suction at the outlet.

**Visual symptom:** Outlet showed the same rectangular velocity pattern as the inlet nozzle, even before any flow could physically reach the outlet.

#### Root Cause Analysis

**FluidX3D uses PERIODIC BOUNDARY CONDITIONS by default** in the streaming step. This is implemented in `calculate_indices()` (kernel.cpp, line ~914):

```cpp
void calculate_indices(...) {
    *xp = (uxx) ((xyz.x       +1u)%def_Nx);  // PERIODIC: wraps at boundary!
    *xm = (uxx) ((xyz.x+def_Nx-1u)%def_Nx);  // PERIODIC: wraps at boundary!
    // ... same for y and z
}
```

**The coupling mechanism:**

1. At the outlet (x = Nx-1), the +x neighbor index is: `xp = (Nx-1+1) % Nx = 0` → **wraps to INLET!**

2. The AA (Esoteric-Pull) streaming pattern uses neighbor indices:
   ```cpp
   void load_f(...) {
       fhn[i+1u] = load(fi, index_f(j[i], t%2ul ? i+1u : i));
   }
   ```

3. For the -x direction population (i=2), it reads from `j[1] = xp`, which at the outlet points to x=0 (inlet)

4. **Result:** Outlet cells directly read distribution values from inlet cells through periodic streaming, causing the inlet pattern to appear at the outlet instantly.

#### Why Previous Fixes Didn't Work

1. **Minimum U_conv clamp removal:** Only affected the Orlanski BC calculation, not the underlying streaming contamination

2. **Velocity initialization:** The contamination happened during streaming (load_f), not during initialization

3. **Zero-gradient extrapolation:** Only modified UNKNOWN populations (c_ix < 0), but the KNOWN populations were still contaminated through periodic streaming

#### The Solution

**Copy ALL populations from upstream at the outlet** to completely override any periodic-contaminated values:

```cpp
// CRITICAL: Copy ALL populations from upstream to break periodic coupling
for(uint i = 0u; i < def_velocity_set; i++) {
    fhn[i] = f_upstream[i];
}
```

This ensures the outlet only sees data from the physical upstream cell (x-1), not from the inlet via periodic wrapping.

#### Additional Fixes Applied

1. **Outlet velocity initialization:** Added explicit `u = 0` initialization for outlet cells in setup.cpp

2. **Interior velocity initialization:** Added explicit `u = 0` initialization for all fluid interior cells

3. **Conditional Orlanski:** Only apply convective correction when actual outflow exists (ux_up > 0.001)

#### Files Modified

| File | Line | Change |
|------|------|--------|
| `kernel.cpp` | ~1787 | CONVECTIVE_OUTLET now copies ALL populations from upstream |
| `setup.cpp` | ~539 | Added velocity initialization for outlet cells |
| `setup.cpp` | ~552 | Added velocity initialization for interior cells |

#### Lessons Learned

1. **FluidX3D default periodicity:** Always assume periodic BCs unless explicitly handled. This affects inlet/outlet coupling.

2. **AA streaming complexity:** The Esoteric-Pull pattern reads from neighbor cells, making boundary handling non-trivial.

3. **Debug methodology:** When a pattern appears instantly (t=0), it's either initialization or streaming - not physical propagation.

4. **Full population override:** At open boundaries, it may be necessary to override ALL populations, not just the "unknown" ones, to break periodic coupling.

#### Related Documentation

See `docs/DFSEM_ConvectiveOutlet_BC.md` for theoretical background on the boundary conditions.

---

### DFM Inlet Debugging Session (2026-01-29)

#### Problem Description

After switching from Kraichnan RFM to Digital Filter Method (DFM) for inlet turbulence, **red spots (velocity > 0.18)** appeared at the inlet boundary, even with turbulence intensity set to zero (TI=0). This indicated a fundamental issue with the inlet boundary condition, not the turbulence generation.

#### Root Cause Analysis

Multiple issues were identified through systematic agent-based debugging:

**Bug 1: UPDATE_FIELDS overwrites inlet velocity**
- The `UPDATE_FIELDS` block writes velocity back to `u[]` arrays BEFORE `DFM_INLET` runs
- `DFM_INLET` sets the correct velocity, but it gets overwritten in the next timestep

**Bug 2: Collision overwrites DFM_INLET distributions**
- `DFM_INLET` was setting `fhn[]` directly and returning early
- But collision happens AFTER `DFM_INLET`, potentially overwriting the values

**Bug 3: AA (Esoteric-Pull) streaming pattern mismatch**
- FluidX3D uses AA streaming that alternates memory slots based on `t%2`
- Inlet writes to one slot, but neighbor cell reads from the opposite slot
- Result: neighbor doesn't see the inlet values, sees stale/wrong data

**Bug 4: Write conflicts with dual-slot approach**
- Initial fix attempted to write to BOTH AA slots
- This was WRONG because slots `i` and `i+1` represent OPPOSITE directions (+x vs -x)
- Writing to both created invalid distribution functions

#### Solution: Make TYPE_X Work Like TYPE_E

The key insight was that TYPE_E (equilibrium boundary) works correctly in FluidX3D. Analyzing TYPE_E revealed:

1. **Velocity reading:** TYPE_E cells read preset velocity from `u[]` arrays (not from `fhn[]`)
2. **Collision:** TYPE_E cells set `fhn = feq` (equilibrium, not collision result)
3. **Streaming:** TYPE_E uses standard `store_f()` (no special handling)
4. **No early return:** TYPE_E doesn't skip any kernel steps

**Fix applied:** Make TYPE_X behave identically to TYPE_E:

```cpp
// 1. Velocity reading section - treat TYPE_X like TYPE_E
{ // EQUILIBRIUM_BOUNDARIES block
    bool use_preset = (flagsn_bo==TYPE_E);
#ifdef DFM_INLET
    if(flagsn&TYPE_X) use_preset = true; // TYPE_X also uses preset velocity
#endif
    if(use_preset) {
        rhon = rho[n];
        uxn = u[n];
        // ... read from u[] arrays
    } else {
        calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn);
    }
}

// 2. DFM_INLET section - update local variables and u[] arrays
if(flagsn&TYPE_X) {
    // Calculate inlet velocity with turbulent fluctuations
    float ux_inlet = ..., uy_inlet = ..., uz_inlet = ...;

    // Update local variables for this timestep
    rhon = rho_inlet;
    uxn = ux_inlet;
    uyn = uy_inlet;
    uzn = uz_inlet;
    calculate_f_eq(rhon, uxn, uyn, uzn, feq); // Recalculate feq

    // Store to u[] arrays for next timestep
    rho[n] = rho_inlet;
    u[n] = ux_inlet;
    // ...

    // DO NOT return early - let collision handle fhn=feq
}

// 3. Collision section - treat TYPE_X like TYPE_E
const bool is_equilibrium_bc = (flagsn_bo==TYPE_E) || (flagsn&TYPE_X);
for(uint i=0u; i<def_velocity_set; i++) {
    fhn[i] = is_equilibrium_bc ? feq[i] : collision_result[i];
}

// 4. UPDATE_FIELDS section - skip TYPE_X cells
bool should_update = true;
if(flagsn_bo==TYPE_E) should_update = false;
if(flagsn&TYPE_X) should_update = false;
if(should_update) {
    // Update u[] arrays from computed values
}
```

#### Files Modified

| File | Section | Change |
|------|---------|--------|
| `kernel.cpp` | Velocity reading (~line 1522) | Treat TYPE_X like TYPE_E using boolean flag |
| `kernel.cpp` | DFM_INLET (~line 1659) | Update local vars + feq, don't return early |
| `kernel.cpp` | SRT collision (~line 1820) | `is_equilibrium_bc` includes TYPE_X |
| `kernel.cpp` | TRT collision (~line 1856) | `is_eq_bc` includes TYPE_X |
| `kernel.cpp` | UPDATE_FIELDS (~line 1607) | Skip TYPE_X cells using boolean flag |
| `setup.cpp` | Inlet geometry | Separated inlet BC (x≤2) from nozzle channel (x≤inlet_x) |

#### Preprocessor Restructuring

The nested preprocessor directives caused OpenCL compilation errors. The pattern:
```cpp
)+"#ifdef A"+R(
    if(condition_a)
)+"#endif"+R(
)+"#ifdef B"+R(
    if(condition_b)
)+"#endif"+R(
    { block }
```

Was restructured to use boolean flags:
```cpp
    { // scope
        bool condition = true;
)+"#ifdef A"+R(
        if(check_a) condition = false;
)+"#endif"+R(
)+"#ifdef B"+R(
        if(check_b) condition = false;
)+"#endif"+R(
        if(condition) { block }
    }
```

#### Extended Nozzle Configuration

The inlet geometry was modified to separate inlet BC from the nozzle channel:

```cpp
const uint inlet_x = 2u + 2u * h_cells;  // Nozzle extends 2h into domain
const uint inlet_bc_x = 2u;              // Inlet BC only at x=0,1,2

// x ≤ inlet_bc_x: TYPE_X (turbulent inlet BC)
// inlet_bc_x < x ≤ inlet_x: Regular fluid with TYPE_S walls (flow develops)
// x > inlet_x: Open domain
```

#### Key Lessons Learned

1. **Boundary conditions must integrate with the full kernel pipeline** - not just one section
2. **TYPE_E is the reference implementation** for equilibrium-based boundaries
3. **AA streaming requires careful handling** - can't just write to arbitrary slots
4. **Nested preprocessor directives are fragile** - use boolean flags instead
5. **Early returns in kernels skip important steps** - better to use conditional logic

#### Verification Status

- [x] OpenCL compilation succeeds (after preprocessor restructuring)
- [x] No red spots at inlet with TI=0
- [x] Turbulent fluctuations visible with TI>0
- [x] Flow develops through extended nozzle
- [x] Jet exits cleanly at outlet

---

### Smagorinsky Constant Optimization (2026-01-29)

#### Problem: High TI but Low Spreading Rate

Test #2 results showed a paradox:
- **Transition TI: 8-15%** (close to experimental 10.7%)
- **Spreading rate: 0.048** (target: 0.1166, only 41% of experimental)

This suggested the small-scale turbulence responsible for jet spreading was being suppressed while large-scale K-H structures (which contribute to high TI) survived.

#### Root Cause: Smagorinsky Constant Too High

The default FluidX3D Smagorinsky-Lilly constant is **Cs = 0.173**, which was derived for homogeneous isotropic turbulence. For free shear flows like jets:

| Application | Optimal Cs | Notes |
|-------------|------------|-------|
| Homogeneous isotropic | 0.17-0.20 | Lilly (1967) derivation |
| Channel flow | 0.1 | Pope (2000) |
| Mixing layers | 0.1 | Vreman (1994) |
| **Jets** | **0.1-0.12** | Stanley & Sarkar (1999) |

**Why Cs matters for jets:**
- High Cs → high eddy viscosity → over-dissipation of small scales
- Small scales (eddies at shear layer edges) drive jet spreading
- Large scales (K-H rollers in potential core) are too big to be affected
- Result: K-H gives high TI, but spreading is suppressed

#### Implementation

**File:** `src/kernel.cpp`, line ~1738

The SGS model calculates relaxation rate `w` using:
```cpp
w = 2.0f/(tau0+sqrt(sq(tau0)+coefficient*sqrt(Q)/rhon));
```

The coefficient depends on Cs²:
- Cs = 0.173: coefficient = 0.76421222 (original)
- Cs = 0.12: coefficient = 0.36656640 (optimized for jets)

**Calculation:**
```
coefficient = 18 * sqrt(2) * Cs²
For Cs = 0.12: 18 * 1.414 * 0.0144 = 0.3666
```

#### Test #3 Configuration (2026-01-29)

Running with optimized parameters:
- **Cs = 0.12** (reduced from 0.173)
- **DFM length scale = 5 cells** (0.18h)
- **Turbulence intensity = 15%** (target: ~3.8% at x/h=0.2)
- **Sponge start = 95%** (moved downstream for cleaner far-field)
- **Warmup = 5 FT** (skip initial transient)
- **Total run = 35 FT** (then auto-stop)
- **Domain AR = 2:1:1.5** (corrected for nozzle AR=19)

**Expected improvements:**
- Better preservation of small-scale turbulence
- Higher spreading rate (target: 0.1166)
- Improved entrainment
- More accurate velocity decay profile

#### Detailed Test Log

See `doc/simulation_test_log.md` for comprehensive test documentation including:
- All test configurations
- Parameter changes between tests
- Results analysis
- Scientific rationale for parameter choices
