# Rectangular Jet Simulation Analysis - FT=3.0

**Date:** 2026-01-30
**Samples:** 135
**Flow-through times:** 3.0 (1.0 FT of averaging after 2.0 FT warmup)

---

## 1. Simulation Parameters

| Parameter | Value | Unit |
|-----------|-------|------|
| Reynolds number (Re_h) | 20,100 | - |
| Nozzle height (h) | ~16 cells | lattice units |
| Domain size | 60h × 25h × 25h | - |
| Inlet velocity (U0) | 0.05 | lattice units |
| Inlet turbulence intensity | 0.7% | - |
| Viscosity (nu) | calculated for Re=20,100 | - |
| Sponge zone start | 90% of domain (x/h ≈ 51) | - |
| Velocity set | D3Q19 | - |
| Collision operator | SRT with Smagorinsky SGS | - |

### Boundary Conditions
- **Inlet (TYPE_X):** DFM turbulent inlet with equilibrium BC
- **Outlet (TYPE_Y):** Convective (Orlanski) BC with Option B fix
- **Lateral:** Periodic BC
- **Nozzle walls:** No-slip (TYPE_S)

---

## 2. Data Correction

The experimental data uses the centerline velocity at x/h=0.2 as the reference velocity (Uc/U0 = 1.0 at x/h=0.2).

Our simulation shows Uc/U0 = 1.1404 at x/h=0, indicating either:
- Vena contracta effect at nozzle exit
- Flow acceleration through nozzle
- Different reference point

**Correction factor applied:** 1/1.1404 = 0.877

All simulation Uc/U0 values are multiplied by this factor for fair comparison.

---

## 3. Figure 7a: Centerline Velocity Decay (Uc/U0 vs x/h)

| x/h | Sim (raw) | Sim (corrected) | Experiment | Error |
|-----|-----------|-----------------|------------|-------|
| 0 | 1.140 | 1.000 | 1.00 | 0% |
| 1 | 1.138 | 0.998 | 1.00 | 0% |
| 2 | 1.135 | 0.995 | 0.99 | +0.5% |
| 3 | 1.128 | 0.989 | 0.99 | 0% |
| 4 | 1.131 | 0.992 | 0.98 | +1% |
| 5 | 1.129 | 0.989 | 0.96 | +3% |
| 6 | 1.121 | 0.983 | 0.94 | +5% |
| 7 | 1.115 | 0.978 | ~0.92 | +6% |
| 8 | 1.107 | 0.970 | ~0.90 | +8% |
| 9 | 1.071 | 0.939 | ~0.85 | +10% |
| 10 | 1.016 | 0.891 | 0.80 | +11% |
| 11 | 0.981 | 0.860 | ~0.76 | +13% |
| 12 | 0.950 | 0.833 | ~0.72 | +16% |
| 13 | 0.882 | 0.773 | ~0.68 | +14% |
| 14 | 0.793 | 0.695 | ~0.65 | +7% |
| 15 | 0.758 | 0.665 | ~0.64 | +4% |
| 20 | 0.638 | 0.559 | 0.58 | -4% |
| 30 | ~0.55 | ~0.48 | 0.46 | +4% |

### Observations:
1. **Near-field (x/h 0-4):** Excellent match after correction (error < 1%)
2. **Transition zone (x/h 6-12):** Simulation decays slower than experiment (error 5-16%)
3. **Far-field (x/h > 15):** Good match returns (error < 5%)

### Potential Core Length:
- **Experiment:** P_c ≈ 3.85h (Uc starts decaying at x/h ≈ 4)
- **Simulation:** P_c ≈ 10-11h (Uc starts decaying at x/h ≈ 10)
- **Discrepancy:** Simulation potential core is ~2.7x longer

---

## 4. Figure 7b: Centerline Turbulence Intensity (Urms/Uc vs x/h)

Conversion: Urms/Uc = (Uc_rms/U0) / (Uc/U0)

| x/h | Sim Uc_rms/U0 | Sim Uc/U0 | Sim Urms/Uc | Exp Urms/Uc | Error |
|-----|---------------|-----------|-------------|-------------|-------|
| 0 | 0.0016 | 1.140 | 0.14% | 3.8% | -3.7% |
| 2 | 0.0051 | 1.135 | 0.45% | 4.1% | -3.6% |
| 4 | 0.0145 | 1.131 | 1.3% | 4.9% | -3.6% |
| 6 | 0.0417 | 1.121 | 3.7% | 7.8% | -4.1% |
| 8 | 0.0871 | 1.107 | 7.9% | ~11% | -3% |
| 10 | 0.156 | 1.016 | 15.3% | 14.7% | +0.6% |
| 12 | 0.176 | 0.950 | 18.5% | ~17% | +1.5% |
| 15 | 0.153 | 0.758 | 20.2% | ~19% | +1% |
| 20 | 0.126 | 0.638 | 19.7% | 19.9% | -0.2% |

### Observations:
1. **Near-field (x/h < 6):** Simulation turbulence is much too low
   - At x/h=0: Sim 0.14% vs Exp 3.8% (simulation is 27x lower!)
   - Inlet DFM turbulence not generating enough fluctuations

2. **Far-field (x/h ≥ 10):** Excellent match
   - Turbulence develops naturally through shear layer instability
   - Self-generated turbulence reaches correct levels

3. **Peak turbulence location:**
   - Experiment: Peak ~20% around x/h=20-30
   - Simulation: Peak ~20% around x/h=15-20
   - Reasonable agreement

---

## 5. Half-width (y0.5/h vs x/h)

| x/h | Simulation y0.5/h | Expected (SR=0.1166) |
|-----|-------------------|----------------------|
| 0 | 0.46 | 0.50 |
| 5 | 0.52 | ~0.50 |
| 10 | 0.63 | ~0.72 |
| 15 | 1.23 | ~1.30 |
| 20 | 2.00 | ~1.88 |

### Spreading Rate:
- **Experiment:** SR_y = 0.1166
- **Simulation (x/h 15-20):** SR = (2.00 - 1.23) / 5 = 0.154
- **Discrepancy:** Simulation spreading rate is 32% higher

---

## 6. Root Cause Analysis

### Primary Issue: Insufficient Inlet Turbulence

| Evidence | Impact |
|----------|--------|
| Inlet Urms/Uc = 0.14% vs experiment 3.8% | 27x too low |
| Extended potential core (10h vs 3.85h) | Delayed mixing |
| Slow velocity decay in transition zone | Less entrainment |
| Far-field matches well | Self-generated turbulence eventually dominates |

### Why DFM Turbulence is Too Low:

1. **def_turbulence_intensity = 0.007 (0.7%)**
   - This targets Urms/U0 = 0.7%
   - But measured Urms/Uc at inlet is only 0.14%
   - Factor of 5 discrepancy suggests DFM amplitude scaling issue

2. **Possible causes:**
   - Filter normalization reducing amplitude
   - Spatial averaging during equilibrium conversion
   - Or the 0.7% is being applied incorrectly

### Secondary Issue: Potential Core Extension

The extended potential core is a direct consequence of low inlet turbulence:
- Without sufficient turbulent mixing, the jet core remains intact longer
- Shear layer instabilities take longer to penetrate to centerline
- This is a known LBM behavior with laminar/low-turbulence inlets

---

## 7. Summary of Discrepancies

| Metric | Simulation | Experiment | Status |
|--------|------------|------------|--------|
| Uc/U0 at x/h=0 (raw) | 1.14 | 1.00 | Needs correction factor |
| Uc/U0 near-field (corrected) | ✓ | ✓ | Good match |
| Uc/U0 transition (x/h 6-12) | Too high | - | Extended potential core |
| Uc/U0 far-field | ✓ | ✓ | Good match |
| Inlet turbulence (Urms/Uc) | 0.14% | 3.8% | 27x too low |
| Far-field turbulence | ~20% | ~20% | Good match |
| Potential core length | ~10h | 3.85h | 2.7x too long |
| Spreading rate | 0.154 | 0.117 | 32% too high |

---

## 8. Recommendations

### Immediate Fix:
1. **Increase inlet turbulence intensity** from 0.007 to ~0.05-0.06
   - Target: Urms/Uc ≈ 4% at inlet to match experiment
   - This should reduce potential core length

### Investigation Needed:
2. **Check DFM amplitude scaling**
   - Why is 0.7% TI parameter producing only 0.14% measured turbulence?
   - May need to adjust filter normalization or amplitude factor

### Future Improvements:
3. **Longer averaging time** - Current 135 samples may have statistical noise
4. **Check nozzle exit condition** - Why is Uc/U0 = 1.14 instead of 1.0?

---

## 9. Files Generated

| File | Description |
|------|-------------|
| `fig7a_centerline_velocity_FT3.0.csv` | Uc/U0 and Uc_rms/U0 vs x/h |
| `fig7b_halfwidth_FT3.0.csv` | y0.5/h vs x/h (mislabeled, not paper's Fig 7b) |
| `fig4_lateral_profiles_FT3.0.csv` | U/Uc vs y/h at various x/h stations |
| `fig10_volume_flow_FT3.0.csv` | Q/Q0 vs x/h |

---

## 10. Experimental Reference

**Paper:** "Mean and turbulent velocity fields of a rectangular jet"
**Conditions:** AR=19, Re_h=20,100, h=10mm, U0=30 m/s

**Key experimental values:**
- Potential core length: P_c = 3.85h
- Spreading rate: SR_y = 0.1166
- Inlet turbulence: Urms/Uc ≈ 3.8% at x/h=0.2
- Far-field turbulence: Urms/Uc ≈ 20% at x/h=20-30
