# Debug Session: FluidX3D Jet Simulation - 2026-01-29

## Initial Problem Statement

The simulation produced incorrect results compared to experiment:

| Metric | Expected | Actual | Error |
|--------|----------|--------|-------|
| Turbulence Intensity at x/h=0.2 | 3.8% | 35.2% | +826% |
| Centerline velocity at x/h=3 | 0.988 | 0.756 | -23% |
| Volume flow at x/h=20 | +1.12 | -0.24 | Reverse flow |

---

## Attempted Fixes (From Original Debug Plan)

### Fix 1: DFSEM Normalization (kernel.cpp:1731)
**Change:** Removed `* 0.5f` from rms_factor calculation
```cpp
// Before:
const float rms_factor = sqrt(fmax(energy_sum * 0.5f, 1e-10f));
// After:
const float rms_factor = sqrt(fmax(energy_sum, 1e-10f));
```
**Status:** Applied, but didn't solve the problem

### Fix 2: Smagorinsky Constant (kernel.cpp:1785)
**Change:** Reduced from 0.76421222f to 0.25f
```cpp
// Before:
w = 2.0f/(tau0+sqrt(sq(tau0)+0.76421222f*sqrt(Q)/rhon));
// After:
w = 2.0f/(tau0+sqrt(sq(tau0)+0.25f*sqrt(Q)/rhon));
```
**Status:** Applied

### Fix 3: Orlanski Convection Velocity (kernel.cpp:1871)
**Original Plan:** Change to `def_c` (speed of sound)
```cpp
// Original plan suggested:
const float U_conv = def_c;  // WRONG - caused wave amplification
```
**Problem:** Using speed of sound (0.577) caused over-extrapolation and wave amplification
**Corrected to:**
```cpp
const float U_conv = clamp(ux_up, 0.01f, 0.15f);  // Local velocity
```
**Status:** Corrected after first failed attempt

### Fix 4: Sponge Zone Parameters (lbm.cpp:480-481)
**Original Plan:**
```cpp
def_sponge_start = 0.95f   // Only 5% of domain
def_sponge_strength = 0.05f // Weak
```
**Problem:** Too weak, didn't absorb waves
**Corrected to:**
```cpp
def_sponge_start = 0.80f   // 20% of domain
def_sponge_strength = 0.10f // Moderate
```
**Status:** Corrected after first failed attempt

### Fix 5: Outlet Velocity Initialization (setup.cpp:540)
**Original Plan:** Remove velocity initialization
**Problem:** Created pressure discontinuity at t=0
**Corrected to:**
```cpp
lbm.u.x[n] = 0.01f;  // Small positive velocity
lbm.u.y[n] = 0.0f;
lbm.u.z[n] = 0.0f;
```
**Status:** Corrected after first failed attempt

---

## Additional Fixes Applied During Debug Session

### Fix 6: Distribution Function Damping in Sponge Zone
**Problem:** Original sponge only modified velocity, not distribution functions. Acoustic waves carried in non-equilibrium part weren't damped.
**Solution:** Added second sponge pass after collision (kernel.cpp after line 1830):
```cpp
// SECOND PASS: Direct distribution function damping
if(x_rel > def_sponge_start && xyz.x > 1u) {
    const float sigma_f = 0.5f * xi * xi * xi;  // Up to 50% relaxation
    for(uint i = 0u; i < def_velocity_set; i++) {
        fhn[i] = fma(1.0f - sigma_f, fhn[i], sigma_f * feq[i]);
    }
}
```
**Status:** Applied - successfully eliminated wave reflection

### Fix 7: Velocity Ramp-Up at Inlet
**Problem:** Sudden velocity jump at t=0 created strong pressure wave
**Solution:** Added tanh ramp over first 5000 timesteps (kernel.cpp:1662-1664):
```cpp
const float ramp_time = 5000.0f;
const float ramp = (t < (ulong)ramp_time) ? tanh(3.0f * (float)t / ramp_time) : 1.0f;
float ux_inlet = def_inlet_velocity * ramp;
```
**Status:** Applied

### Fix 8: Turbulence Intensity Ramp-Up
**Solution:** TI also scales with ramp (kernel.cpp:1670):
```cpp
const float TI = def_turbulence_intensity * ramp;
```
**Status:** Applied

---

## RNG Fix Attempt

### Problem Identified
The original Kraichnan RFM used a broken pseudo-random number generator:
```cpp
// BROKEN - all r1-r6 linearly correlated:
const float seed = (float)(m + 1u) * golden;
const float r1 = fmod(seed * 1.0f, 1.0f);
const float r2 = fmod(seed * 7.0f, 1.0f);
// ... etc
```
This caused visible "donut/ring" interference patterns at the inlet instead of random turbulence.

### Fix Applied
Replaced with hash-based RNG (kernel.cpp:1682-1717):
```cpp
// Hash-based RNG with independent sequences:
uint h1 = m * 1664525u + 1013904223u;
h1 = ((h1 >> 16u) ^ h1) * 0x45d9f3bu;
h1 = ((h1 >> 16u) ^ h1) * 0x45d9f3bu;
const float r1 = (float)(h1 & 0xFFFFu) / 65535.0f;
// ... separate hash for each r1-r7
```

### Additional Changes in RNG Fix
- Increased N_modes from 32 to 64
- Added fallback for zero cross-product in amplitude direction
- Added sqrt(2) factor to mode_val (Kraichnan formula)
- Changed k sampling from linear to log-uniform
- Increased temporal decorrelation (omega: 0.005 -> 0.01)
- Added independent r7 for phase

**Status:** Applied - eliminated donut pattern, BUT created new problem

---

## Current Problem: Turbulence Amplitude Too High

### Observation
After RNG fix, the inlet shows:
- Mean velocity: 0.05 (correct)
- Peak velocity: 0.18+ (red dots visible at inlet)
- Actual fluctuation: ~0.13 (260% TI)
- Expected fluctuation: ~0.001 (0.7% TI)

**The fluctuations are ~130-370x larger than target!**

### Root Cause Analysis (from Code Review)

#### Identified Issue: Missing 1/sqrt(N) Factor
Standard Kraichnan formula:
```
u'(x) = sqrt(2/N) × Σ σₙ × cos(kₙ·x + φₙ)
```

Current implementation:
```cpp
mode_val = sqrt(2) × amp × cos(...)  // Missing 1/sqrt(N) = 1/8
```

This alone would cause 8x amplification.

#### Other Contributing Factors
1. Vector component (ax²) averaging to 1/3 - but this reduces amplitude
2. Peak vs RMS at specific spatial points - could add 2-4x
3. Constructive interference of modes - could add 2-4x

#### Unexplained Factor
Even accounting for all identified issues, there's still ~16x unexplained amplification.

### Mathematical Analysis

| Quantity | Calculated | Observed |
|----------|------------|----------|
| energy_sum | ~2.5 | - |
| rms_factor | ~1.6 | - |
| scale | ~0.00022 | - |
| u_prime (typical) | ~2-6 | - |
| Expected fluctuation | ~0.001 | - |
| Observed fluctuation | - | ~0.13 |
| Unexplained ratio | - | ~130x |

---

## Simulation Results Summary

### Run 1: Original Plan Fixes
- **Result:** FAILED - Wave reflection (magenta throughout domain), chaotic inlet pattern
- **Cause:** Orlanski U_conv=def_c too aggressive, sponge too weak, outlet init removed

### Run 2: Corrected Sponge and Outlet
- **Result:** FAILED - Still had wave reflection and chaotic pattern
- **Cause:** Sponge only damped velocity, not distribution functions

### Run 3: Added Distribution Function Damping + Ramp
- **Result:** PARTIAL SUCCESS - Wave reflection eliminated
- **Cause:** But inlet showed donut/ring interference pattern

### Run 4: Fixed RNG (Hash-Based)
- **Result:** PARTIAL SUCCESS - Donut pattern eliminated, jet developing
- **Cause:** But turbulence amplitude ~130-370x too high (red dots at inlet)

---

## Files Modified

| File | Lines | Changes |
|------|-------|---------|
| kernel.cpp | 1731 | Removed `* 0.5f` from rms_factor |
| kernel.cpp | 1785 | Reduced Smagorinsky constant |
| kernel.cpp | 1871 | Changed U_conv from def_c to clamped local velocity |
| kernel.cpp | 1660-1670 | Added velocity and TI ramp-up |
| kernel.cpp | 1677-1770 | Replaced RNG with hash-based, added sqrt(2), increased modes |
| kernel.cpp | 1832-1846 | Added distribution function damping (second sponge pass) |
| lbm.cpp | 480-481 | Changed sponge start to 0.80, strength to 0.10 |
| setup.cpp | 538-544 | Added small positive outlet velocity init |

---

## Next Steps

The current Kraichnan RFM approach has fundamental issues with amplitude normalization that are difficult to debug. A different approach for inlet turbulence generation will be explored.

---

## Lessons Learned

1. **Orlanski BC:** Use local flow velocity, NOT speed of sound
2. **Sponge zones:** Must damp distribution functions, not just velocity
3. **Startup transients:** Use velocity ramp to avoid pressure spikes
4. **RNG quality:** Golden ratio sequences are correlated; use proper hash functions
5. **Kraichnan RFM:** Requires careful normalization including 1/sqrt(N) factor
6. **Debugging CFD:** Visual inspection of results is essential for identifying issues
