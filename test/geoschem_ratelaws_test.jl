@testsnippet RateLawSetup begin
    using DynamicQuantities
    using SymbolicIndexingInterface: getsym
    using SciMLBase: ODEProblem
    using ModelingToolkit: t, @parameters, mtkcompile, System, Equation, @named, extend
    using GasChem, Test

    @parameters T=293.0 [unit = u"K"]
    N_A = 6.02214076e23 # Avogadro's number
    cm3_m3 = 1e6
    @parameters num_density = 2.7e19 / N_A * cm3_m3,
    [unit = u"mol/m^3", description = "Number density of air."]
end

@testitem "constant_k" setup=[RateLawSetup] begin
    c = 1.8e-12
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.constant_k(t, T, num_density, c), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    @test getsym(prob, sys.k)(prob) ≈ 0.0486
end

@testitem "constant_k_1" setup=[RateLawSetup] begin
    c = 1.8e-12
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.constant_k_1(t, c), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    @test getsym(prob, sys.k)(prob) ≈ 1.8e-12
end

@testitem "regress_T" setup=[RateLawSetup] begin
    a_0, b_0, T_0 = 1.33e-13, 3.82e-11, -2000.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.regress_T(t, T, num_density, a_0, b_0, T_0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    @test getsym(prob, sys.k)(prob) ≈ 0.004710333943456039
end

@testitem "arrhenius_ppb" setup=[RateLawSetup] begin
    a0, b0, c0 = 3.00e-12, 0.0, -1500.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.arrhenius_ppb(t, T, num_density, a0, b0, c0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    @test getsym(prob, sys.k)(prob) ≈ 0.00048432225861535784
end

@testitem "arrhenius_ppb_3" setup=[RateLawSetup] begin
    a0, b0, c0 = 2.90e-16, 0.0, -1000.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.arrhenius_ppb_3(t, T, num_density, a0, b0, c0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    @test getsym(prob, sys.k)(prob) ≈ 6964.52977737289
end

@testitem "arrhenius_mlc_SI" setup=[RateLawSetup] begin
    a0, b0, c0 = 3.00e-12, 0.0, -1500.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.arrhenius_mlc_SI(t, T, a0, b0, c0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    @test getsym(prob, sys.k)(prob) ≈ 1.7937861430198436e-20
end

@testitem "arrhenius_mlc_1" setup=[RateLawSetup] begin
    a0, b0, c0 = 3.00e-12, 0.0, -1500.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.arrhenius_mlc_1(t, T, a0, b0, c0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    @test getsym(prob, sys.k)(prob) ≈ 1.793786143019844e-14
end

@testitem "arr_3rdbody" setup=[RateLawSetup] begin
    a1, b1, c1 = 6.90e-31, 1.0, 0.0
    a2, b2, c2 = 2.6e-11, 0.0, 0.0
    fv = 0.6
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(
        GasChem.arr_3rdbody(
            t, T, num_density, a1, b1, c1, a2, b2, c2, fv), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.17987080444910986
end

@testitem "arr_3rdbody_1" setup=[RateLawSetup] begin
    a1, b1, c1 = 6.90e-31, 1.0, 0.0
    a2, b2, c2 = 2.6e-11, 0.0, 0.0
    fv = 0.6
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(
        GasChem.arr_3rdbody_1(
            t, T, num_density, a1, b1, c1, a2, b2, c2, fv), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 6.661881646263327e-12
end

@testitem "rate_HO2HO2" setup=[RateLawSetup] begin
    @parameters H2O=1000.0, [unit=u"1.0", description="H2O; Water vapor"]
    a0, c0 = 3.00e-13, 460.0
    a1, c1 = 2.1e-33, 920.0
    @named base_sys = System(Equation[], t, [], [T, num_density, H2O])
    sys = mtkcompile(extend(GasChem.rate_HO2HO2(t, T, num_density, H2O, a0, c0, a1, c1), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.07430493859316299
end

@testitem "rate_OHCO" setup=[RateLawSetup] begin
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_OHCO(t, T, num_density), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.006602767113875056
end

@testitem "rate_RO2NO_a1" setup=[RateLawSetup] begin
    a0, c0 = 2.80e-12, 300.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_RO2NO_a1(t, T, num_density, a0, c0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 6.314124896654882e-5
end

@testitem "rate_RO2NO_b1" setup=[RateLawSetup] begin
    a0, c0 = 2.80e-12, 300.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_RO2NO_b1(t, T, num_density, a0, c0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.21040768863952958
end

@testitem "rate_RO2NO_a2" setup=[RateLawSetup] begin
    a0, c0, a1 = 2.60e-12, 365.0, 2.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_RO2NO_a2(t, T, num_density, a0, c0, a1), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.006257279178031173
end

@testitem "rate_RO2NO_b2" setup=[RateLawSetup] begin
    a0, c0, a1 = 2.60e-12, 365.0, 2.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_RO2NO_b2(t, T, num_density, a0, c0, a1), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.2377217070157082
end

@testitem "tbranch" setup=[RateLawSetup] begin
    a0, b0, c0 = 9.50e-14, 0.0, 390.0
    a1, b1, c1 = 2.62e1, 0.0, -1130.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.tbranch(t, T, num_density, a0, b0, c0, a1, b1, c1), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.006248197805954408
end

@testitem "rate_OHHNO3" setup=[RateLawSetup] begin
    a0, c0 = 2.41e-14, 460.0
    a1, c1 = 6.51e-34, 1335.0
    a2, c2 = 2.69e-17, 2199.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_OHHNO3(t, T, num_density, a0, c0, a1, c1, a2, c2), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.0031275792345493445
end

@testitem "eq_const" setup=[RateLawSetup] begin
    a0, c0 = 4.90e-03, 12100.0
    a1, b1 = 9.70e-29, 5.6
    a2, b2 = 9.3e-12, 1.5
    fv = 0.6
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(
        GasChem.eq_const(
            t, T, num_density, a0, c0, a1, b1, a2, b2, fv), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 2.1206340759025532e-27
end

@testitem "eq_const_1" setup=[RateLawSetup] begin
    a0, c0 = 4.90e-03, 12100.0
    a1, b1 = 9.70e-29, 5.6
    a2, b2 = 9.3e-12, 1.5
    fv = 0.6
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(
        GasChem.eq_const_1(
            t, T, num_density, a0, c0, a1, b1, a2, b2, fv), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 2.1206340759025532e-27
end

@testitem "rate_RO2HO2" setup=[RateLawSetup] begin
    a0, c0, a1 = 2.91e-13, 1300.0, 4.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_RO2HO2(t, T, num_density, a0, c0, a1), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.4147862847681038
end

@testitem "rate_GLYCOH_a" setup=[RateLawSetup] begin
    a0 = 3.1e-12
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_GLYCOH_a(t, T, num_density, a0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.06695598252070009
end

@testitem "rate_GLYCOH_b" setup=[RateLawSetup] begin
    a0 = 3.1e-12
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_GLYCOH_b(t, T, num_density, a0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.016744017479299923
end

@testitem "rate_HACOH_a" setup=[RateLawSetup] begin
    a0, c0 = 3.15e-14, 920.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_HACOH_a(t, T, num_density, a0, c0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.01612257331622036
end

@testitem "rate_HACOH_b" setup=[RateLawSetup] begin
    a0, c0 = 3.15e-14, 920.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_HACOH_b(t, T, num_density, a0, c0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.0035259242069770776
end

@testitem "rate_DMSOH" setup=[RateLawSetup] begin
    a0, c0 = 1.12e-11, -250.0
    a1, c1 = 3.00e-31, 4020.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_DMSOH(t, T, num_density, a0, c0, a1, c1), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 7.287375473764786e17
end

@testitem "rate_GLYXNO3" setup=[RateLawSetup] begin
    a0, c0 = 1.4e-12, -1900.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_GLYXNO3(t, T, num_density, a0, c0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 3.567255335036053e-5
end

@testitem "arrplus_mlc_1" setup=[RateLawSetup] begin
    a0, b0, c0, d0, e0 = 1.0e-12, 298.0, 0.0, 1.0, 0.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.arrplus_mlc_1(t, T, a0, b0, c0, d0, e0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 3.6165489651031117e-13
end

@testitem "arrplus_ppb" setup=[RateLawSetup] begin
    a0, b0, c0, d0, e0 = 1.0e-12, 298.0, 0.0, 1.0, 0.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.arrplus_ppb(t, T, num_density, a0, b0, c0, d0, e0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.009764682205778403
end

@testitem "tunplus_mlc_1" setup=[RateLawSetup] begin
    a0, b0, c0, d0, e0 = 1.0e-12, 298.0, 1.0e8, 1.0, 0.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.tunplus_mlc_1(t, T, a0, b0, c0, d0, e0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 1.473210332683109e-10
end

@testitem "tunplus_ppb" setup=[RateLawSetup] begin
    a0, b0, c0, d0, e0 = 1.0e-12, 298.0, 1.0e8, 1.0, 0.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.tunplus_ppb(t, T, num_density, a0, b0, c0, d0, e0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 3.977667898244395
end

@testitem "rate_ISO1" setup=[RateLawSetup] begin
    a0, b0, c0, d0, e0, f0, g0 = 1.0e-11, 390.0, 0.3, 1.0, 380.0, 1.0, 380.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(
        GasChem.rate_ISO1(
            t, T, num_density, a0, b0, c0, d0, e0, f0, g0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.7210220267414859
end

@testitem "rate_ISO2" setup=[RateLawSetup] begin
    a0, b0, c0, d0, e0, f0, g0 = 1.0e-11, 390.0, 0.3, 1.0, 380.0, 1.0, 380.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(
        GasChem.rate_ISO2(
            t, T, num_density, a0, b0, c0, d0, e0, f0, g0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.3009403732544113
end

@testitem "rate_EPO" setup=[RateLawSetup] begin
    a1, e1, m1 = 1.0e-11, 390.0, 1.5e-4
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_EPO(t, T, num_density, a1, e1, m1), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 2.523363950607153e-16
end

@testitem "rate_PAN_abab" setup=[RateLawSetup] begin
    a0, b0 = 4.9e-3, 12100.0
    a1, b1 = 5.4e16, -13830.0
    cf = 0.3
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_PAN_abab(t, T, num_density, a0, b0, a1, b1, cf), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.00017076882173869417
end

@testitem "rate_PAN_acac" setup=[RateLawSetup] begin
    a0, c0 = 9.0e-28, 8.9
    a1, c1 = 7.7e-12, 0.2
    cf = 0.6
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_PAN_acac(t, T, num_density, a0, c0, a1, c1, cf), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.1981654191522894
end

@testitem "rate_NIT" setup=[RateLawSetup] begin
    a0, b0, c0, n, x0, y0 = 2.9e-12, 350.0, 0.08, 1.0, 8.1e-4, 0.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_NIT(t, T, num_density, a0, b0, c0, n, x0, y0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 2.5119154809815447e-5
end

@testitem "rate_ALK" setup=[RateLawSetup] begin
    a0, b0, c0, n, x0, y0 = 2.9e-12, 350.0, 0.08, 1.0, 8.1e-4, 0.0
    @named base_sys = System(Equation[], t, [], [T, num_density])
    sys = mtkcompile(extend(GasChem.rate_ALK(t, T, num_density, a0, b0, c0, n, x0, y0), base_sys))
    prob = ODEProblem(sys, [], (0.0, 3600.0))
    k_val = getsym(prob, sys.k)(prob)
    @test k_val ≈ 0.00018517572785290386
end
