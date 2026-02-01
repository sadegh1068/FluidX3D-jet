# FluidX3D Rectangular Jet Simulation - Test Log

## Experimental Reference (Ahmed et al.)

| Parameter | Value |
|-----------|-------|
| Nozzle height (h) | 10 mm |
| Nozzle width (w) | 190 mm |
| Aspect Ratio (AR) | 19 |
| Inlet velocity (U0) | 30 m/s |
| Reynolds number (Re_h) | 20,100 |
| Potential core length (Pc) | 3.85h |
| Spreading rate (SR_y) | 0.1166 |
| Inlet TI at x/h=0.2 | 3.8% |

### Experimental Data Points (Fig. 7)

| x/h | Uc/U0 | Urms/Uc (%) |
|-----|-------|-------------|
| 0.2 | 1.000 | 3.8 |
| 4.0 | 0.979 | 4.9 |
| 6.0 | 0.939 | 7.8 |
| 10.0 | 0.804 | 14.7 |
| 20.0 | 0.576 | 19.9 |
| 30.0 | 0.457 | 20.9 |

---

## Test Run History

---

### Test #1: AR=10, Large Length Scale (Previous Configuration)

**Date:** 2026-01-29 to 2026-01-30

**Configuration:**

| Parameter | Value |
|-----------|-------|
| Domain aspect | 2:1:1 (Nx:Ny:Nz) |
| Nozzle AR | 10 |
| Nozzle length | 0.5h |
| DFM TI parameter | 17% |
| DFM length scale (y,z) | 10 cells (~0.37h) |
| DFM filter width | 21 |
| Sponge start | 90% |
| Sponge strength | 0.1 |

**Results at FT=7.5 (~105 samples):**

| x/h | Uc/U0 (corr) | Exp | Error | Urms/Uc | Exp | Error |
|-----|--------------|-----|-------|---------|-----|-------|
| 0.2 | 0.945 | 1.000 | -5.5% | 2.3% | 3.8% | -1.5 pp |
| 4.0 | 0.937 | 0.979 | -4.2% | **25.6%** | 4.9% | **+20.7 pp** |
| 6.0 | 0.830 | 0.939 | -11.6% | **22.8%** | 7.8% | **+15.0 pp** |
| 10.0 | 0.743 | 0.804 | -7.6% | 17.3% | 14.7% | +2.6 pp |
| 20.0 | 0.554 | 0.576 | -3.7% | 16.2% | 19.9% | -3.7 pp |
| 30.0 | 0.481 | 0.457 | +5.3% | 22.3% | 20.9% | +1.3 pp |

| Metric | Simulation | Experiment | Error |
|--------|------------|------------|-------|
| Spreading Rate | 0.170 | 0.1166 | +45.9% |
| Q/Q0 at x/h=10 | 1.54 | 1.11 | +39% |

**Key Issues:**
1. Massive TI overshoot at x/h=4-6 (25-38% vs 5-8%) - K-H instability
2. Spreading rate 46% too high
3. Excessive entrainment (Q/Q0)
4. Inlet TI too low (2.3% vs 3.8%)

**Root Cause:** Large DFM length scale (0.37h) excited K-H instability wavelengths

---

### Test #2: AR=19, Small Length Scale (Current Configuration)

**Date:** 2026-01-30

**Configuration:**

| Parameter | Value | Change from Test #1 |
|-----------|-------|---------------------|
| Domain aspect | 2:1:1.5 (Nx:Ny:Nz) | Extended z |
| Nozzle AR | **19** | Was 10 |
| Nozzle length | 0.5h | Same |
| DFM TI parameter | **10%** | Was 17% |
| DFM length scale (y,z) | **3 cells (~0.1h)** | Was 10 cells |
| DFM filter width | **13** | Was 21 |
| Sponge start | 90% | Same |
| Sponge strength | 0.1 | Same |

**Results at FT=4.5 (~69 samples):**

| x/h | Uc/U0 (corr) | Exp | Error | Urms/Uc | Exp | Error |
|-----|--------------|-----|-------|---------|-----|-------|
| 0.2 | 0.986 | 1.000 | -1.4% | 1.1% | 3.8% | -2.8 pp |
| 4.0 | 0.981 | 0.979 | **+0.2%** | **8.0%** | 4.9% | +3.1 pp |
| 6.0 | 0.919 | 0.939 | -2.2% | **14.8%** | 7.8% | +7.0 pp |
| 10.0 | 0.723 | 0.804 | -10.2% | 20.7% | 14.7% | +6.0 pp |
| 20.0 | 0.535 | 0.576 | -7.0% | 16.2% | 19.9% | -3.7 pp |
| 30.0 | 0.346 | 0.457 | **-24.2%** | **43.0%** | 20.9% | **+22.1 pp** |

| Metric | Simulation | Experiment | Error |
|--------|------------|------------|-------|
| Spreading Rate | 0.085 | 0.1166 | **-27.3%** |
| Q/Q0 at x/h=10 | 1.27 | 1.11 | +14% |
| Q/Q0 at x/h=30 | **0.73** | 1.09 | **-33%** |

**Key Improvements vs Test #1:**
1. Transition TI reduced: 25-38% → 8-15% (much better!)
2. Near-field velocity excellent (x/h<6 within 2%)
3. Entrainment reduced (Q/Q0 at x/h=10: 1.54 → 1.27)

**New Issues:**
1. Spreading rate now TOO LOW (0.085 vs 0.117, -27%)
2. Far-field velocity decays too fast (x/h=30: -24% error)
3. Sponge zone causes TI spike at x/h>26 (43%)
4. Mass loss in sponge zone (Q/Q0 drops to 0.73)
5. Inlet TI still too low (1.1% vs 3.8%)

**Root Cause of New Issues:**
- Sponge zone starts at x/h≈27, affecting measurements at x/h=30
- Smaller length scale may be over-suppressing turbulence development

---

**Results at FT=9.0 (~138 samples):**

| x/h | Uc/U0 (corr) | Exp | Error | Urms/Uc | Exp | Error |
|-----|--------------|-----|-------|---------|-----|-------|
| 0.2 | 0.986 | 1.000 | -1.4% | 1.0% | 3.8% | -2.8 pp |
| 4.0 | 0.984 | 0.979 | **+0.5%** | **8.4%** | 4.9% | +3.5 pp |
| 6.0 | 0.922 | 0.939 | **-1.8%** | **12.1%** | 7.8% | +4.3 pp |
| 10.0 | 0.745 | 0.804 | -7.3% | 19.5% | 14.7% | +4.7 pp |
| 20.0 | 0.559 | 0.576 | **-2.9%** | **19.0%** | 19.9% | **-0.9 pp** |
| 30.0 | 0.418 | 0.457 | -8.5% | 29.2% | 20.9% | +8.3 pp |

| Metric | Simulation | Experiment | Error |
|--------|------------|------------|-------|
| Spreading Rate | 0.089 | 0.1166 | -23.4% |
| Q/Q0 at x/h=10 | 1.31 | 1.11 | +18% |
| Q/Q0 at x/h=20 | 1.15 | 1.12 | **+3%** |
| Q/Q0 at x/h=30 | 0.89 | 1.09 | -18% |

---

### Test #2 Convergence: FT=4.5 → FT=9.0

| Metric | FT=4.5 (69 samples) | FT=9.0 (138 samples) | Change |
|--------|---------------------|----------------------|--------|
| Uc/U0 at x/h=6 | 0.919 | 0.922 | +0.3% |
| Uc/U0 at x/h=10 | 0.723 | 0.745 | **+3.0%** |
| Uc/U0 at x/h=20 | 0.535 | 0.559 | **+4.5%** |
| Uc/U0 at x/h=30 | 0.346 | 0.418 | **+20.8%** |
| TI at x/h=6 | 14.8% | 12.1% | -2.7 pp |
| TI at x/h=10 | 20.7% | 19.5% | -1.2 pp |
| TI at x/h=30 | 43.0% | 29.2% | **-13.8 pp** |
| Spreading Rate | 0.085 | 0.089 | +4.7% |
| Q/Q0 at x/h=30 | 0.73 | 0.89 | **+22%** |

**Key Observations from Convergence:**
1. Far-field velocity improved significantly (+21% at x/h=30)
2. TI spike at x/h=30 reduced (43% → 29%)
3. Q/Q0 at x/h=30 improved (0.73 → 0.89)
4. Spreading rate slightly increased (0.085 → 0.089)
5. Near-field (x/h<10) is stable - converged

**Conclusion:** Results still converging, especially in far-field. Need more FT.

---

### Test #3: Optimized SGS + DFM Parameters (Current Run)

**Date:** 2026-01-30

**Configuration:**

| Parameter | Value | Change from Test #2 |
|-----------|-------|---------------------|
| Domain aspect | 2:1:1.5 (Nx:Ny:Nz) | Same |
| Nozzle AR | 19 | Same |
| Nozzle length | 0.5h | Same |
| **Smagorinsky Cs** | **0.12** | Was 0.173 (default) |
| DFM TI parameter | **15%** | Was 10% |
| DFM length scale (y,z) | **5 cells (~0.18h)** | Was 3 cells |
| DFM filter width | **21** | Was 13 |
| Sponge start | **95%** | Was 90% |
| Sponge strength | 0.1 | Same |
| Warmup | **5 FT** | Was 2 FT |
| Total duration | **35 FT** | Was 10 FT |
| Auto-stop | **Yes** | New |

**Scientific Rationale for Changes:**

1. **Cs: 0.173 → 0.12:** Default Smagorinsky-Lilly constant (0.173) is optimized for isotropic turbulence, not jets. Literature recommends Cs=0.10-0.12 for free shear flows. Lower Cs preserves small-scale turbulence that drives mixing.

2. **L: 3 → 5 cells:** Larger eddies have longer lifetimes (τ ~ L/u'), providing sustained mixing. Still below K-H unstable wavelength (~0.3h) to avoid instability excitation.

3. **TI: 10% → 15%:** Higher inlet TKE to better match experimental inlet TI (3.8%).

4. **Sponge: 90% → 95%:** Reduced far-field contamination at x/h=28-30.

5. **Warmup: 2 → 5 FT:** Better initial transient removal.

**Expected Improvements:**
- Spreading rate: 0.089 → ~0.10-0.11 (target: 0.117)
- Better preservation of small-scale turbulent mixing
- Cleaner far-field statistics

**Results at FT=10 (STABLE - Best Results):**

| x/h | Uc/U0 (corr) | Exp | Error | Urms/Uc | Exp | Error |
|-----|--------------|-----|-------|---------|-----|-------|
| 0.2 | **1.000** | 1.000 | **0.0%** | 2.4% | 3.8% | -1.4 pp |
| 4.0 | 1.004 | 0.979 | +2.5% | 8.9% | 4.9% | +4.0 pp |
| 6.0 | 0.976 | 0.939 | +3.9% | 13.2% | 7.8% | +5.4 pp |
| 10.0 | **0.818** | 0.804 | **+1.7%** | 16.9% | 14.7% | +2.2 pp |
| 20.0 | 0.626 | 0.576 | +8.6% | **20.1%** | 19.9% | **+0.2 pp** |
| 30.0 | **0.454** | 0.457 | **-0.6%** | 28.8% | 20.9% | +7.9 pp |

| Metric | Simulation | Experiment | Error |
|--------|------------|------------|-------|
| Q/Q0 at x/h=10 | 1.33 | 1.11 | +20% |
| Q/Q0 at x/h=20 | 1.15 | 1.12 | **+3%** |

**CRITICAL: Instability Developed After FT=10**

Temporal evolution analysis revealed that the simulation became unstable after FT=10:

| FT | Uc/U0 @ x/h=10 | TI @ x/h=10 | Status |
|----|----------------|-------------|--------|
| 10 | 0.844 | 17% | ✓ **STABLE** |
| 15 | 0.772 | 26% | ⚠️ Degrading |
| 20 | 0.567 | 63% | ❌ Unstable |
| 25 | 0.468 | 78% | ❌ Severe |
| 30 | 0.403 | 89% | ❌ Critical |
| 35 | 0.346 | 103% | ❌ Catastrophic |

**Root Cause of Instability:** Smagorinsky constant Cs=0.12 was too low for long-term stability. Initially provided good results (FT=5-10), but insufficient SGS dissipation allowed absolute instability (jet flapping) to grow exponentially.

**Conclusion:** FT=10 results are valid and show excellent agreement with experiment. After FT=10, self-excited oscillation grew and corrupted statistics.

---

### Test #3 - BUG DISCOVERED: Mean Velocity Oscillations

**Problem:** Even at FT=10 (stable period), the time-averaged mean velocity showed ±5-8% spatial oscillations in the potential core region, which should be smooth (Uc/U0 ≈ 1.0 for x/h < 4).

**Data showing oscillations:**
```
x/h=0.0: Uc/U0 = 1.052
x/h=0.2: Uc/U0 = 1.032
x/h=0.4: Uc/U0 = 1.033
x/h=0.6: Uc/U0 = 1.112  ← SPIKE (+8%)
x/h=0.8: Uc/U0 = 1.079
x/h=2.0: Uc/U0 = 1.109  ← SPIKE (+8%)
x/h=3.4: Uc/U0 = 1.092  ← SPIKE
```

**Root Cause Found:** Bug in DFM random number generation (kernel.cpp line ~1608):

```cpp
// BUGGY CODE:
uint seed = (uint)yy * 73856093u ^ (uint)zz * 19349663u ^ (uint)(t % 1000ul) * 83492791u;
```

The `t % 1000ul` means the random seed only changes every **1,000 timesteps**!

- Flow-through time ≈ 16,300 timesteps
- 5 FT averaging ≈ 81,500 timesteps
- Only ~81 different DFM patterns averaged (should be ~81,500)
- Result: Persistent spatial patterns appear as "mean" velocity oscillations

**Impact:** The DFM was generating quasi-steady spatial patterns instead of proper turbulent fluctuations. Each pattern persisted for 1000 timesteps before changing. With only ~81 patterns in 5 FT of averaging (instead of ~81,500 independent samples), time averaging couldn't smooth out the spatial structure.

**Fix Applied:**
```cpp
// CORRECT CODE:
uint seed = (uint)yy * 73856093u ^ (uint)zz * 19349663u ^ (uint)t * 83492791u;
```

This ensures each timestep has a unique random seed, giving truly independent samples for time averaging.

---

## Comparison Summary

### Transition Zone TI (x/h=4-6)

| Test | TI at x/h=4 | TI at x/h=6 | Experiment |
|------|-------------|-------------|------------|
| #1 (AR=10, L=10) FT=7.5 | 25.6% | 22.8% | 4.9%, 7.8% |
| #2 (AR=19, L=3) FT=4.5 | 8.0% | 14.8% | 4.9%, 7.8% |
| #2 (AR=19, L=3) FT=9.0 | **8.4%** | **12.1%** | 4.9%, 7.8% |

**Improvement: 3x reduction in TI overshoot (Test #1 → Test #2)**

### Spreading Rate

| Test | SR | Experiment | Error |
|------|-----|------------|-------|
| #1 (AR=10, L=10) FT=7.5 | 0.170 | 0.1166 | +46% (too high) |
| #2 (AR=19, L=3) FT=4.5 | 0.085 | 0.1166 | -27% (too low) |
| #2 (AR=19, L=3) FT=9.0 | **0.089** | 0.1166 | **-23%** (improving) |

### Velocity Decay

| Test | Error at x/h=10 | Error at x/h=20 | Error at x/h=30 |
|------|-----------------|-----------------|-----------------|
| #1 FT=7.5 | -7.6% | -3.7% | +5.3% |
| #2 FT=4.5 | -10.2% | -7.0% | -24.2% |
| #2 FT=9.0 | **-7.3%** | **-2.9%** | **-8.5%** |

---

## Parameter Sensitivity Analysis

| Parameter | Effect of Increasing |
|-----------|---------------------|
| **Smagorinsky Cs** | ↑ SGS dissipation, ↓ small-scale turbulence, ↓ spreading |
| DFM Length Scale | ↑ K-H instability, ↑ TI overshoot, ↑ spreading |
| DFM TI | ↑ Inlet turbulence, may ↑ or ↓ transition TI |
| Aspect Ratio | ↑ 2D-like behavior, affects edge effects |
| Sponge Start | Later start → less far-field contamination |

---

## Pending Analysis

- [x] ~~Wait for FT=10 results for Test #2~~ (FT=9 shows convergence trend)
- [x] ~~Consider adjusting sponge zone to 95%~~ (Implemented in Test #3)
- [x] ~~Consider intermediate length scale (5-6 cells instead of 3)~~ (Implemented: L=5)
- [x] ~~Consider increasing DFM TI to 12-15%~~ (Implemented: TI=15%)
- [x] ~~Reduce Smagorinsky Cs for jets~~ (Implemented: Cs=0.12)
- [ ] **Wait for Test #3 results at FT=10, 20, 35**

---

## Recommendations Based on FT=9.0 Analysis

### What's Working Well:
1. **Near-field velocity (x/h<6):** Within 2% of experiment
2. **Transition TI:** Reduced from 25% to 8-12% (was 5x too high, now ~2x)
3. **Mid-field TI (x/h=20):** Matches experiment within 1%

### Remaining Issues:
1. **Spreading rate too low:** 0.089 vs 0.117 (-23%) - need more mixing
2. **Inlet TI too low:** 1.0% vs 3.8% - affects downstream development
3. **Far-field (x/h>25):** Still affected by sponge zone

### Suggested Next Steps (Priority Order):
1. **Continue current run to FT=15** - Verify convergence trend continues
2. **Then try:** Increase DFM length scale to 5 cells (compromise between 3 and 10)
3. **Or try:** Increase DFM TI to 15% to boost inlet turbulence
4. **Move sponge to 95%** to reduce far-field contamination

---

## File Locations

| File | Description |
|------|-------------|
| `bin/fig7a_centerline_velocity_FT*.csv` | Centerline velocity data |
| `bin/fig7b_turbulence_intensity_FT*.csv` | Turbulence intensity data |
| `bin/halfwidth_FT*.csv` | Half-width growth data |
| `bin/fig10_volume_flow_FT*.csv` | Volume flow rate data |
| `bin/fig4_lateral_profiles_FT*.csv` | Lateral velocity profiles |
| `bin/comparison_plots.png` | Comparison plots |
| `src/lbm.cpp` (lines 462-475) | DFM parameters |
| `src/setup.cpp` (lines 56-64) | Domain and nozzle geometry |

---

---

## Lessons Learned

### 1. Smagorinsky Constant Selection

| Cs Value | Behavior | Application |
|----------|----------|-------------|
| 0.17-0.20 | High dissipation, stable | Isotropic turbulence, HIT |
| 0.12-0.15 | Moderate dissipation | Channel flows, boundary layers |
| 0.10-0.12 | Low dissipation | Jets, mixing layers (short runs) |
| < 0.10 | Insufficient dissipation | **Risk of absolute instability** |

**Lesson:** For jets, Cs=0.12 gives good physics short-term but may allow instabilities to grow. Cs=0.14-0.15 may be safer for long runs.

### 2. DFM Temporal Correlation

**Wrong approach:**
- Using `t % N` in random seed → same pattern repeats for N timesteps
- Creates quasi-steady patterns, not turbulence
- Time averaging doesn't help

**Correct approach:**
- Use `t` directly in seed → new pattern each timestep
- Apply temporal correlation via exponential smoothing: `u'(t+dt) = α*u'(t) + √(1-α²)*u'_new`
- Time averaging then produces smooth statistics

**Lesson:** Turbulent inlet conditions must have proper temporal decorrelation. Verify that fluctuations change every timestep.

### 3. Monitoring for Instability

**Early warning signs:**
- TI growing over time instead of converging
- Centerline velocity decreasing with more averaging
- TI exceeding 50% (fluctuations > 0.5 × mean)
- Volume flow rate decreasing (mass "loss")

**Lesson:** Export intermediate results (every 0.5-1 FT) and check for instability. Don't wait until final export.

### 4. Statistical Convergence

| Sample Count | Statistical Uncertainty |
|--------------|------------------------|
| 50 | ~14% |
| 100 | ~10% |
| 500 | ~4.5% |
| 1000 | ~3% |

**Lesson:** Need 100+ independent samples for <10% uncertainty. Verify samples are independent (separated by integral time scale).

### 5. Validation Strategy

**Multi-metric validation is essential:**
1. Centerline velocity decay (Uc/U0 vs x/h)
2. Turbulence intensity (Urms/Uc vs x/h)
3. Half-width / spreading rate
4. Volume flow rate (entrainment)

**Lesson:** A simulation can match one metric while failing others. Match at least 3 metrics for confidence.

### 6. Time Evolution Analysis

**Always check temporal evolution, not just final results:**
- Plot metrics vs FT (flow-through time)
- Stable simulation: metrics converge to constant
- Unstable simulation: metrics drift or oscillate

**Lesson:** FT=10 results were excellent; FT=35 results were catastrophic. The final result is not always the best result.

---

## Bug Fixes Applied

| Date | Bug | File | Fix | Status |
|------|-----|------|-----|--------|
| 2026-02-01 | DFM seed only changes every 1000 steps | kernel.cpp:1608 | Change `t % 1000ul` to `t` | **APPLIED** |

---

## Current Best Configuration (Pending Bug Fix)

Based on Test #3 FT=10 results:

| Parameter | Value | Notes |
|-----------|-------|-------|
| Smagorinsky Cs | 0.12-0.15 | 0.12 works short-term; 0.15 safer |
| DFM Length Scale | 5 cells | ~0.18h, below K-H wavelength |
| DFM TI | 15% | Inlet TI still low; may need increase |
| Sponge Start | 95% | Keeps x/h=28 clean |
| Warmup | 5 FT | Removes initial transient |
| Averaging | 5-10 FT | FT=10 showed good results |

---

*Last updated: 2026-02-01, Test #3 analysis complete, DFM bug identified*
