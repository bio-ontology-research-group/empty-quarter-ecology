# 19. Thermal performance curve + cross-desert salinity calibration

**Question.** Are there well-defined thermal/salinity niche bounds for the EQ community, and how do they compare to other arid systems?

**Method.** Fit Gaussian (or skewed) response curves to taxonomic abundance vs MAT and salinity proxies. Calibrate against Atacama/Namib/Gurbantunggut where overlapping. `scripts/run_thermal_and_calibration.py`, outputs `cache/thermal_performance_curve.tsv`, `cache/cross_desert_salinity_calibration.tsv`.

**Key results.**
- EQ thermal optima cluster around **32–38 °C** for surface/deep dominants.
- Salinity tolerances overlap heavily with Namib + Atacama salt-flats; **no unique EQ thermal niche.**

**Status.** caveated. The "thermal bound" claim from earlier referee-defense (see #54) was retracted upon finding broad cross-desert overlap.

**Cross-refs.** 47, 54.
