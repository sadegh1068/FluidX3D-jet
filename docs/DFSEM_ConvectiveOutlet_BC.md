# DFSEM Inlet and Convective Outlet Boundary Conditions for FluidX3D

## Overview

This document details the implementation of two critical boundary conditions for turbulent jet simulations in FluidX3D:
1. **DFSEM (Divergence-Free Synthetic Eddy Method)** - Turbulent inlet BC
2. **Convective Outlet BC with Fluctuation-Only Sponge** - Non-reflective outlet BC

**Key Requirement:** Neither BC uses TYPE_E (equilibrium boundary). Both define distributions directly from velocity/density fields.

---

## Part 1: DFSEM Turbulent Inlet Boundary Condition

### 1.1 Background and Literature

**Primary Reference:** Poletto, R., Craft, T., Revell, A. (2013). "A New Divergence Free Synthetic Eddy Method for the Reproduction of Inlet Flow Conditions for LES." Flow, Turbulence and Combustion, 91(3), 519-539.

**Why DFSEM over standard SEM:**
- Standard SEM generates velocity fluctuations that are NOT divergence-free
- Non-zero divergence creates spurious pressure waves in incompressible flow
- DFSEM applies SEM to vorticity field, then converts to velocity via curl
- Curl of any vector field is automatically divergence-free: ∇·(∇×ω) = 0

### 1.2 DFSEM Algorithm

#### Step 1: Define Eddy Box
The eddy box surrounds the inlet plane with extent σ (eddy length scale) in all directions:
```
x ∈ [-σ, +σ]  (streamwise, centered on inlet)
y ∈ [y_min - σ, y_max + σ]
z ∈ [z_min - σ, z_max + σ]
```

#### Step 2: Initialize Eddies
For N_eddy eddies, each eddy k has:
- Position: x_k = (x_k, y_k, z_k) - random within eddy box
- Intensity: ε_k = (ε_k^x, ε_k^y, ε_k^z) - random ±1 for each component
- Convection velocity: U_c (bulk inlet velocity)

#### Step 3: Shape Function
Tent function (linear decay):
```
f(r) = sqrt(3/2) * max(0, 1 - |r|)  for |r| ≤ 1
```
Or Gaussian for smoother results:
```
f(r) = sqrt(2) * exp(-π * r²)
```

#### Step 4: Compute Vorticity Fluctuations
At each inlet point x:
```
ω'(x) = (1/sqrt(N_eddy)) * Σ_k ε_k * f_σ(x - x_k)
```
where f_σ is the shape function scaled by eddy size σ.

#### Step 5: Convert to Velocity via Stream Function
The velocity fluctuation is obtained from a vector potential (stream function):
```
ψ(x) = (σ/sqrt(N_eddy)) * Σ_k ε_k × (x - x_k) * F_σ(x - x_k)
```
where F_σ is the integral of the shape function.

Then:
```
u'(x) = ∇ × ψ(x)
```

This guarantees ∇·u' = 0.

#### Step 6: Scale to Target Reynolds Stress
```
u'_scaled = A · u'_raw
```
where A is the Cholesky decomposition of the target Reynolds stress tensor R:
```
R = <u'u'> = A · A^T
```

For isotropic turbulence with intensity TI:
```
A = TI * U_0 * I  (identity matrix)
```

#### Step 7: Advect Eddies
Each timestep:
```
x_k(t+Δt) = x_k(t) + U_c * Δt
```
When an eddy exits the box (x_k > x_inlet + σ), recycle it:
- New random position at upstream face (x = x_inlet - σ)
- New random intensity ε_k

### 1.3 LBM Implementation for Inlet

**Critical:** Do NOT use TYPE_E. Instead, reconstruct the full distribution function.

#### Non-Equilibrium Distribution Approach

The distribution function at the inlet is:
```
f_i = f_i^eq(ρ, u) + f_i^neq
```

**Equilibrium part** (known from target ρ and u = U_mean + u'):
```
f_i^eq = w_i * ρ * [1 + (c_i·u)/c_s² + (c_i·u)²/(2c_s⁴) - u²/(2c_s²)]
```

**Non-equilibrium part** from Chapman-Enskog (regularized approach):
```
f_i^neq = (w_i / 2c_s⁴) * Q_i:Π^neq
```
where:
```
Q_iαβ = c_iα*c_iβ - c_s²*δ_αβ
Π^neq_αβ = -ρ*ν*(∂u_α/∂x_β + ∂u_β/∂x_α)  (from target velocity gradient)
```

For inlet with uniform inflow, Π^neq ≈ 0, so:
```
f_i ≈ f_i^eq(ρ_inlet, u_inlet + u'_DFSEM)
```

#### Unknown Populations at Inlet (x = 0, +x direction into domain)

For D3Q19 at x=0 inlet, unknown populations are those with c_ix > 0:
- i = 1 (c = [+1, 0, 0])
- i = 7 (c = [+1, +1, 0])
- i = 9 (c = [+1, -1, 0])
- i = 13 (c = [+1, 0, +1])
- i = 15 (c = [+1, 0, -1])

These 5 populations are SET (not streamed) from the DFSEM-computed equilibrium.

#### Handling Known Populations
Populations with c_ix ≤ 0 stream normally from interior. We only override the unknowns.

### 1.4 Simplified DFSEM for GPU

Full DFSEM requires tracking N_eddy positions (memory overhead). Simplified GPU version:

**Pseudo-Random Fourier Mode Approximation:**
```cpp
// At each inlet cell (y, z) and time t:
float3 u_prime = float3(0.0f);

for(uint m = 0; m < N_modes; m++) {
    // Deterministic pseudo-random phases (no memory needed)
    float phase_y = fract(sin(m * 12.9898f + y * 78.233f) * 43758.5453f) * 2π;
    float phase_z = fract(sin(m * 93.9898f + z * 47.233f) * 24958.5453f) * 2π;
    float phase_t = fract(sin(m * 43.2918f + t * 0.001f) * 93821.2381f) * 2π;

    // Wave numbers for integral length scale Lt
    float ky = (m + 1) * 2π / Lt;
    float kz = (m + 1) * 2π / Lt;

    // Divergence-free construction: u' = curl(ψ)
    // ψ = (0, 0, ψ_z) → u' = (∂ψ_z/∂y, -∂ψ_z/∂x, 0) for 2D slice
    // Extended to 3D with multiple components

    float amp = 1.0f / sqrt((float)(m + 1)); // Energy decay
    u_prime.x += amp * cos(ky * y + phase_y) * cos(kz * z + phase_z + phase_t);
    u_prime.y += amp * 0.5f * sin(ky * y + phase_y) * cos(kz * z + phase_z + phase_t);
    u_prime.z += amp * 0.5f * cos(ky * y + phase_y) * sin(kz * z + phase_z + phase_t);
}

// Scale to target intensity
u_prime *= TI * U_0 / normalization_factor;
```

This is NOT true DFSEM but provides correlated fluctuations without memory overhead.

---

## Part 2: Convective Outlet Boundary Condition

### 2.1 Background and Literature

**Primary References:**
- Orlanski, I. (1976). "A simple boundary condition for unbounded hyperbolic flows." J. Comput. Phys. 21:251-269.
- Lou, Z., Luo, L.-S., Shi, B. (2013). "Evaluation of outflow boundary conditions for two-phase LBM." Phys. Rev. E, 87(6), 063301.

**Why Convective BC:**
- Allows structures to exit domain without reflection
- No artificial pressure specification needed
- Maintains mass conservation better than extrapolation alone
- Works for both laminar and turbulent flows

### 2.2 Orlanski Convective Boundary Condition

The convective BC advects all quantities out of the domain:
```
∂φ/∂t + U_c * ∂φ/∂n = 0
```
where:
- φ = any field variable (or distribution function f_i)
- U_c = convection velocity
- n = outward normal direction

**Discretization (first-order upwind):**
```
φ(t+Δt) = φ(t) - U_c * Δt/Δx * [φ(t) - φ_upstream(t)]
```

For LBM with Δt = Δx = 1:
```
f_i^new = f_i^old - U_c * (f_i^old - f_i,upstream^old)
```

### 2.3 Convection Velocity Estimation

**Option 1: Local velocity**
```
U_c = u_x(n)  // x-component at outlet cell
```
Problem: Can become negative, causing instability.

**Option 2: Averaged velocity**
```
U_c = <u_x>_outlet  // Average over outlet plane
```
Problem: Requires reduction operation (expensive on GPU).

**Option 3: Fixed characteristic velocity**
```
U_c = U_bulk * 0.5  // Half of inlet bulk velocity
```
Stable but less accurate for non-uniform flows.

**Recommended: Clamped local velocity**
```
U_c = clamp(u_x(n), 0.1 * U_bulk, 0.9 * c_s)
```
This ensures:
- U_c > 0 (always outward)
- U_c < c_s (CFL stability)

### 2.4 Unknown Populations at Outlet

For D3Q19 at x = Nx-1 outlet (+x boundary), unknown populations have c_ix < 0:
- i = 2 (c = [-1, 0, 0])
- i = 8 (c = [-1, +1, 0])
- i = 10 (c = [-1, -1, 0])
- i = 14 (c = [-1, 0, +1])
- i = 16 (c = [-1, 0, -1])

These arrive from outside the domain and must be specified.

### 2.5 Convective BC for Distribution Functions

Apply Orlanski to each unknown population:
```cpp
// At outlet cell n, with upstream cell n_up (at x-1)
for each unknown i in {2, 8, 10, 14, 16}:
    f_i[n] = f_i[n] - U_c * (f_i[n] - f_i[n_up])
```

**Important timing:** This must be applied AFTER collision but BEFORE streaming, operating on post-collision distributions.

### 2.6 Alternative: Non-Equilibrium Extrapolation

If convective BC alone is insufficient:
```cpp
// Extrapolate density and velocity from interior
rho_out = 2*rho[n-1] - rho[n-2]  // Linear extrapolation
u_out = 2*u[n-1] - u[n-2]

// Or zero-gradient:
rho_out = rho[n-1]
u_out = u[n-1]

// Compute equilibrium
f_eq = calculate_f_eq(rho_out, u_out)

// Extrapolate non-equilibrium part
f_neq[n-1] = f[n-1] - f_eq[n-1]
f_neq_out = f_neq[n-1]  // Zero-gradient for neq

// Reconstruct
f_out = f_eq + f_neq_out
```

---

## Part 3: Fluctuation-Only Sponge Zone

### 3.1 Problem with Standard Sponge

Standard sponge damps velocity toward zero or ambient:
```
u_damped = (1 - σ) * u + σ * u_target
```

If u_target = 0 or u_ambient, this also damps the mean jet velocity, creating artificial deceleration and backflow.

### 3.2 Fluctuation-Only Damping

**Key idea:** Decompose velocity into mean and fluctuation, damp only fluctuation.

```
u = <u> + u'
u_damped = <u> + (1 - σ) * u'
```

**Problem:** Computing <u> requires temporal averaging (memory + history).

### 3.3 Practical Implementation: Damp Toward Local Mean

**Approximation:** Use spatial filtering instead of temporal averaging.

At each sponge cell, estimate local mean from neighbors:
```cpp
// 7-point stencil average (excluding current cell)
float3 u_local_mean = (u[n-1] + u[n+1] + u[n-Nx] + u[n+Nx] +
                       u[n-Nx*Ny] + u[n+Nx*Ny]) / 6.0f;

// Fluctuation
float3 u_prime = u[n] - u_local_mean;

// Damp fluctuation only
float sigma = sponge_strength(x);  // 0 at start, increases toward outlet
u[n] = u_local_mean + (1.0f - sigma) * u_prime;
```

### 3.4 Simpler Alternative: Weak Damping to Extrapolated Value

**Approach:** Damp toward zero-gradient extrapolated velocity (not zero).

```cpp
// In sponge zone
float3 u_upstream = u[n - stride_x];  // Upstream cell
float sigma = sponge_strength(x);

// Damp toward upstream value (removes fluctuations, preserves mean transport)
u[n] = (1.0f - sigma) * u[n] + sigma * u_upstream;
```

This effectively smooths the flow, removing small-scale fluctuations while preserving the mean jet profile.

### 3.5 Sponge Profile

**Polynomial ramp (recommended: cubic or quintic):**
```cpp
float sponge_strength(uint x, uint Nx, uint sponge_start) {
    if (x < sponge_start) return 0.0f;

    float xi = (float)(x - sponge_start) / (float)(Nx - 1 - sponge_start);

    // Cubic profile (smooth start)
    return xi * xi * xi;

    // Or quintic (smoother)
    // return xi * xi * xi * (6*xi*xi - 15*xi + 10);
}
```

**Parameters:**
- Sponge start: 85-90% of domain length
- Maximum strength: 0.05-0.2 (weak!)
- Thickness: 10-15% of domain length

---

## Part 4: Combined Implementation Strategy

### 4.1 New Flags Required

```cpp
#define TYPE_I_INLET  0b01000000  // DFSEM turbulent inlet (replaces TYPE_X)
#define TYPE_C_OUTLET 0b10000000  // Convective outlet (replaces TYPE_Y)
```

### 4.2 Kernel Order of Operations

```
1. Load f_i from memory
2. Compute macroscopic (ρ, u)
3. [SPONGE ZONE] Apply fluctuation damping to u
4. Compute f_eq
5. [DFSEM INLET] Override unknown populations with DFSEM distributions
6. Collision (TRT/SRT)
7. [CONVECTIVE OUTLET] Apply Orlanski to unknown populations
8. Store f_i to memory (with implicit streaming via AA pattern)
```

### 4.3 Memory Considerations

**DFSEM:**
- Simplified RFM version: No extra memory
- Full DFSEM: N_eddy × 6 floats (position + intensity) = ~1-10 KB

**Convective Outlet:**
- Need previous timestep f_i at outlet: Already available (no extra memory)
- Upstream f_i: Already in memory

**Sponge Zone:**
- No extra memory if using spatial filtering

---

## Part 5: Validation Checklist

### 5.1 Inlet BC Validation
- [ ] Turbulence intensity at inlet matches target (5% of U_0)
- [ ] Fluctuations are spatially correlated (check autocorrelation)
- [ ] No spurious pressure waves from inlet
- [ ] Mean velocity profile is correct (top-hat for jet)

### 5.2 Outlet BC Validation
- [ ] No reflection of vortices at outlet (visualize vorticity)
- [ ] Jet exits smoothly without distortion
- [ ] Mass conservation (total mass remains constant)
- [ ] No artificial acceleration/deceleration at outlet

### 5.3 Overall Simulation
- [ ] Potential core length P_c ≈ 3.85h (experimental)
- [ ] Spreading rate SR_y ≈ 0.1166 (experimental)
- [ ] Centerline velocity decay follows U_c/U_0 ~ (x/P_c)^(-0.5)

---

## Part 6: FluidX3D-Specific Implementation Notes

### 6.1 CRITICAL: Periodic Boundary Condition Coupling

**FluidX3D uses PERIODIC BOUNDARY CONDITIONS by default** in the streaming step. This creates a hidden coupling between inlet and outlet that must be handled explicitly.

#### The Problem

In `kernel.cpp`, the `calculate_indices()` function uses modular arithmetic:
```cpp
void calculate_indices(...) {
    *xp = (uxx) ((xyz.x       +1u)%def_Nx);  // Wraps at x=Nx-1 → x=0
    *xm = (uxx) ((xyz.x+def_Nx-1u)%def_Nx);  // Wraps at x=0 → x=Nx-1
}
```

This means:
- At x = Nx-1 (outlet), the +x neighbor is x = 0 (inlet)
- At x = 0 (inlet), the -x neighbor is x = Nx-1 (outlet)

#### Impact on AA Streaming

The Esoteric-Pull (AA) streaming pattern reads certain populations from neighbor cells:
```cpp
void load_f(...) {
    fhn[i+1u] = load(fi, index_f(j[i], t%2ul ? i+1u : i));
}
```

At the outlet (x = Nx-1), when loading the -x population (i=2), it reads from `j[1] = xp = 0` (inlet!). This causes **instant contamination** of the outlet with inlet data.

#### The Solution

**Override ALL populations at the outlet** from the upstream cell (x-1):
```cpp
// CRITICAL: Copy ALL populations to break periodic coupling
for(uint i = 0u; i < def_velocity_set; i++) {
    fhn[i] = f_upstream[i];
}
```

This ensures the outlet only sees physically valid data from the interior, not contaminated data from the inlet via periodic wrapping.

### 6.2 Velocity Field Initialization

When using TYPE_X (inlet) or TYPE_Y (outlet) without TYPE_E, you MUST explicitly initialize the velocity field:

```cpp
// Outlet cells - MUST initialize velocity!
else if(x >= Nx - 3u) {
    lbm.flags[n] = TYPE_Y;
    lbm.u.x[n] = 0.0f;  // IMPORTANT: Initialize to zero
    lbm.u.y[n] = 0.0f;
    lbm.u.z[n] = 0.0f;
    lbm.rho[n] = 1.0f;
}

// Interior cells - also initialize!
else {
    lbm.u.x[n] = 0.0f;
    lbm.u.y[n] = 0.0f;
    lbm.u.z[n] = 0.0f;
    lbm.rho[n] = 1.0f;
}
```

Without TYPE_E, the equilibrium boundary code doesn't run, so initial velocities won't be overwritten. Uninitialized memory can cause artifacts.

### 6.3 Kernel Execution Order in FluidX3D

The actual order in `stream_collide()`:

```
1. Load f from memory (with AA streaming)
2. Compute (ρ, u) from f
3. [TEMPERATURE] Handle temperature if enabled
4. [VOLUME_FORCE] Apply volume forces
5. [UPDATE_FIELDS] Write (ρ, u) to memory
6. [SPONGE_ZONE] Damp velocity toward upstream
7. Calculate f_eq from (ρ, u)
8. [DFSEM_INLET] Override all f for TYPE_X cells
9. [SUBGRID] Compute turbulent viscosity
10. [COLLISION] TRT or SRT collision
11. [CONVECTIVE_OUTLET] Override all f for TYPE_Y cells from upstream
12. Store f to memory
```

**Key insight:** The DFSEM_INLET and CONVECTIVE_OUTLET both happen AFTER collision, modifying the post-collision distributions before storage.

### 6.4 Flag Bit Definitions

```cpp
#define TYPE_S 0b00000001  // Solid boundary
#define TYPE_E 0b00000010  // Equilibrium boundary (avoid for inlet/outlet!)
#define TYPE_T 0b00000100  // Temperature boundary
#define TYPE_F 0b00001000  // Fluid (free surface)
#define TYPE_I 0b00010000  // Interface (free surface)
#define TYPE_G 0b00100000  // Gas (free surface)
#define TYPE_X 0b01000000  // DFSEM inlet (use alone, no TYPE_E)
#define TYPE_Y 0b10000000  // Convective outlet (use alone, no TYPE_E)
```

**Recommended usage for jet simulations:**
- Inlet nozzle: `TYPE_X` alone
- Inlet walls: `TYPE_S`
- Outlet: `TYPE_Y` alone
- Lateral far-field: `TYPE_E` (equilibrium with u=0 is OK here)
- Interior: No flag (or 0)

---

## References

1. Poletto, R., Craft, T., Revell, A. (2013). "A New Divergence Free Synthetic Eddy Method." Flow Turb. Combust. 91(3):519-539.

2. Jarrin, N., Benhamadouche, S., Laurence, D., Prosser, R. (2006). "A synthetic-eddy-method for generating inflow conditions for LES." Int. J. Heat Fluid Flow 27(4):585-593.

3. Orlanski, I. (1976). "A simple boundary condition for unbounded hyperbolic flows." J. Comput. Phys. 21:251-269.

4. Latt, J., Chopard, B., Malaspinas, O., Deville, M., Michler, A. (2008). "Straight velocity boundaries in the lattice Boltzmann method." Phys. Rev. E 77:056703.

5. Lou, Z., Luo, L.-S., Shi, B. (2013). "Evaluation of outflow boundary conditions for two-phase LBM." Phys. Rev. E 87:063301.

6. Skordos, P.A. (1993). "Initial and boundary conditions for the lattice Boltzmann method." Phys. Rev. E 48:4823.
