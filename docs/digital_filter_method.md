# Digital Filter Method (DFM) for Inlet Turbulence Generation

## Overview

The Digital Filter Method (DFM) generates synthetic turbulent inflow conditions for LES/DNS simulations by filtering spatially and temporally uncorrelated random data to produce fluctuations with prescribed integral length scales and Reynolds stress tensors.

## Key References

1. **Klein, M., Sadiki, A., & Janicka, J. (2003)**. "A digital filter based generation of inflow data for spatially developing direct numerical or large eddy simulations." *Journal of Computational Physics*, 186(2), 652-665.

2. **Xie, Z.-T., & Castro, I.P. (2008)**. "Efficient generation of inflow conditions for large eddy simulation of street-level flows." *Flow, Turbulence and Combustion*, 81(3), 449-470.

3. **Veloudis, I., Yang, Z., McGuirk, J.J., Page, G.J., & Spencer, A. (2007)**. "Novel implementation and assessment of a digital filter based approach for the generation of LES inlet conditions." *Flow, Turbulence and Combustion*, 79(1), 1-24.

## Mathematical Formulation

### Core Principle

DFM generates correlated turbulent fluctuations by convolving uncorrelated white noise with a filter kernel that encodes the desired spatial correlation:

```
u'(y,z) = Σ_m Σ_n  b_m * b_n * r(y+m, z+n)
```

where:
- `r(y,z)` = uncorrelated Gaussian random field (white noise)
- `b_k` = filter coefficients encoding spatial correlation
- `m, n` = filter indices from `-N` to `+N`

### Filter Coefficients (Exponential Correlation)

For an exponential two-point correlation `R(r) = exp(-π|r|/L)`, the filter coefficients are:

```
b_k = b̃_k / √(Σ b̃_k²)

where:
b̃_k = exp(-π|k|/n)
n = L/Δx  (integral length scale in grid cells)
```

The filter half-width `N` should be at least `N ≥ 2n` to capture the correlation decay.

### 2D Filtering via Sequential 1D Operations

For a 2D inlet plane (y-z), the 2D convolution can be decomposed into two sequential 1D convolutions:

```
ψ(y,z) = Σ_n b_n^(z) * Σ_m b_m^(y) * r(y+m, z+n)
```

This reduces complexity from O(N²) to O(2N) per point.

### Temporal Correlation (Forward-Stepwise Method)

The Xie & Castro (2008) Forward-Stepwise Method introduces temporal correlation:

```
u'(t+Δt) = α * u'(t) + √(1-α²) * u'_new

where:
α = exp(-Δt/T)
T = L/U_mean  (integral time scale)
```

This creates an exponentially decaying temporal autocorrelation while preserving variance.

### Reynolds Stress Transformation (Cholesky)

For anisotropic turbulence, the fluctuations are transformed using the Cholesky decomposition of the target Reynolds stress tensor:

```
[u']     [a₁₁  0    0  ] [ψ₁]
[v'] = A [a₂₁ a₂₂  0  ] [ψ₂]
[w']     [a₃₁ a₃₂ a₃₃] [ψ₃]
```

where `A * Aᵀ = R_ij` (Reynolds stress tensor).

For **isotropic turbulence**, this simplifies to:
```
u' = TI * U₀ * ψ
v' = TI * U₀ * ψ
w' = TI * U₀ * ψ
```

where `TI` is the turbulence intensity and `U₀` is the mean inlet velocity.

## Algorithm Steps

### 1. Setup Phase (computed once)

```
1. Set filter half-width: N = ceil(2 * L/Δx)
2. Compute filter coefficients:
   For k = -N to N:
     b̃[k] = exp(-π|k|/n)

3. Normalize: b[k] = b̃[k] / sqrt(Σ b̃[k]²)
4. Store as constants for GPU kernel
```

### 2. Runtime Phase (each timestep, each inlet cell)

```
For each inlet cell at (y, z):

  // Step 1: Generate random field at neighbor positions
  For m = -N to N:
    For n = -N to N:
      r[m,n] = hash_random(y+m, z+n, t)  // 3 components

  // Step 2: Apply spatial filter
  ψ_x = Σ_m Σ_n b[m] * b[n] * r_x[m,n]
  ψ_y = Σ_m Σ_n b[m] * b[n] * r_y[m,n]
  ψ_z = Σ_m Σ_n b[m] * b[n] * r_z[m,n]

  // Step 3: Apply temporal correlation (requires storing previous value)
  // Note: For simplicity, can skip if temporal decorrelation is fast enough
  // u'(t) = α * u'(t-1) + sqrt(1-α²) * ψ

  // Step 4: Scale to target intensity
  u' = TI * U₀ * ψ_x
  v' = TI * U₀ * ψ_y
  w' = TI * U₀ * ψ_z

  // Step 5: Add to mean flow
  u_total = U_mean + u'
  v_total = v'
  w_total = w'
```

## Advantages Over Kraichnan RFM

| Aspect | Kraichnan RFM | DFM |
|--------|---------------|-----|
| **Amplitude Control** | Complex normalization via energy spectrum sum | Direct: `TI × U₀` |
| **Normalization** | Sum over N modes with varying amplitudes | Simple filter norm `√(Σb²)` |
| **Spatial Correlation** | Phase relationships in Fourier space | Explicit convolution kernel |
| **Implementation** | Many random numbers per mode (direction, phase, amplitude) | One random number per neighbor |
| **Failure Mode** | Normalization errors cause 100x+ amplification | Bounded by filter structure |
| **Divergence-Free** | Enforced via `a ⊥ k` | Not enforced (acceptable for LBM) |
| **GPU Suitability** | Good (embarrassingly parallel) | Good (local stencil operation) |

### Why RFM Failed in Our Case

The Kraichnan RFM implementation had issues with:

1. **Energy sum normalization**: The `energy_sum` accumulated amplitude² values, but the relationship to RMS velocity is indirect
2. **Mode amplitude spectrum**: The von Kármán simplification `1/(1 + k²)` wasn't properly calibrated
3. **The sqrt(2) factor**: Present in Kraichnan's formula but interacts with other factors
4. **Result**: 130-370x larger fluctuations than intended (0.18 instead of 0.05 velocity)

DFM avoids these issues because:
- The filter coefficients are normalized to unity variance (`Σb² = 1`)
- The random input has known variance (uniform [-1,1] → var = 1/3, or Gaussian → var = 1)
- Scaling is direct: `output_rms = TI * U₀`

## GPU/OpenCL Implementation Notes

### Filter Coefficient Storage

Since the filter is small (typically N=10-20), coefficients can be:
- Stored in constant memory (fast access, shared across all work items)
- Or computed on-the-fly using `exp(-π|k|/n)` (saves memory, adds computation)

### Random Number Generation

For each inlet cell and filter neighbor, we need random numbers. Options:

1. **Hash-based RNG** (recommended for GPU):
   ```c
   uint hash(uint y, uint z, uint t) {
     uint seed = y * 73856093u ^ z * 19349663u ^ t * 83492791u;
     seed = (seed ^ (seed >> 16)) * 0x45d9f3bu;
     seed = (seed ^ (seed >> 16)) * 0x45d9f3bu;
     return seed;
   }
   float random = (hash & 0xFFFF) / 32767.5f - 1.0f;  // [-1, 1]
   ```

2. **Precomputed random texture**: Store random values in 2D texture, index by (y,z,t)

### Temporal Correlation Storage

For the Forward-Stepwise Method, we need to store `u'(t-1)` for each inlet cell. Options:

1. **Store in separate buffer**: Extra memory, but allows any α value
2. **Seed-based regeneration**: Use `t-1` in hash to regenerate previous value (but then no true temporal correlation)
3. **Skip temporal correlation**: If advection time scale is small compared to turbulent time scale

For simplicity in the initial implementation, we skip explicit temporal correlation. The random seed changes each timestep, providing temporal variation without explicit correlation.

### Convolution Optimization

The double loop over filter indices is expensive. Optimizations:

1. **Separable filtering**: Apply 1D filter in y, then in z (reduces O(N²) to O(2N))
2. **Shared memory**: Load random values for a tile into shared memory
3. **Unrolled loops**: For small N, unroll the filter loops

### Boundary Handling

At domain boundaries, filter neighbors may be outside the domain. Options:
- **Periodic wrap**: Use modulo arithmetic (natural for periodic BCs)
- **Clamp**: Repeat boundary value
- **Zero-padding**: Treat outside as zero (reduces variance near boundaries)

For our jet simulation with periodic y-z boundaries, periodic wrapping is appropriate.

## Parameter Selection

### Filter Width (N)

```
N = ceil(2 * L / Δx)
```

where `L` is the integral length scale in physical units and `Δx` is the grid spacing.

For `L = 10 cells`: N = 20, filter width = 41 points (2N+1)

### Integral Length Scale (L)

Typical values for jet inflows:
- L ≈ 0.1 × nozzle height (h) for shear layer turbulence
- L ≈ 0.5h for core turbulence
- For our case: h ≈ 27 cells, so L ≈ 3-14 cells

### Turbulence Intensity (TI)

```
TI = u'_rms / U₀
```

For jet experiments:
- Typical inlet TI: 0.5-2%
- Our target: 0.7% at inlet → ~3.8% at x/h = 0.2 (after shear layer amplification)

## Validation

After implementation, verify:

1. **Velocity bounds**: No red dots at inlet (velocity < 0.18, should be ~0.05)
2. **Randomness**: No structured patterns (checkerboard, stripes)
3. **Spatial correlation**: Fluctuations should vary smoothly over ~L cells
4. **Amplitude**: RMS velocity fluctuation ≈ TI × U₀
5. **Jet development**: Smooth spreading without spurious oscillations

## References for Further Reading

- Di Mare, L., & Jones, W.P. (2003). "LES of turbulent flow past a swept fence." *Int. J. Heat Fluid Flow*, 24, 606-615.
- Tabor, G.R., & Baba-Ahmadi, M.H. (2010). "Inlet conditions for large eddy simulation: A review." *Computers & Fluids*, 39(4), 553-567.
- Poletto, R., Craft, T., & Revell, A. (2013). "A new divergence free synthetic eddy method for the reproduction of inlet flow conditions for LES." *Flow Turbulence Combust.*, 91, 519-539.
