
using BenchmarkTools

@benchmark j_mean_H2O2(3600 * 12.0, 30.0f0, 200.0f0)

@benchmark j_mean_CH2Oa(3600 * 12.0, 30.0f0, 200.0f0)

@benchmark j_mean_CH3OOH(3600 * 12.0, 30.0f0, 200.0f0)

@benchmark j_mean_NO2(3600 * 12.0, 30.0f0, 200.0f0)

@benchmark j_mean_o31D(3600 * 12.0, 30.0f0, 200.0f0)
