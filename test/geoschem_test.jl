using GasChem
using Test
using DifferentialEquations, ModelingToolkit

tspan = (0.0,360.0)

# Unit Test 0: Base case
@test begin
    u_0 = [
        2.132595818651838
        64.6838219200333
        83.23806225285274
        17.647692925315166
        514.9056372536066
        5.861495732807038
        56.62922153797334
        ]

    rs = FullChem().sys
    @unpack O3, NO2, ISOP, O1D, OH, DMS, H2O, H2O2, N2O5 = rs
    rs_solved = solve(ODEProblem(rs, [], tspan, [], combinatoric_ratelaws=false),AutoTsit5(Rosenbrock23()))
    test0 = [rs_solved[O3][end],rs_solved[NO2][end],rs_solved[ISOP][end],rs_solved[O1D][end],rs_solved[OH][end],rs_solved[DMS][end],rs_solved[H2O][end]]

    test0 ≈ u_0
end


# Unit Test 1: O1D sensitivity to O3
@test begin
    u_1 = 1.7279825730298626e-5

    rs1 = FullChem().sys
    @unpack O3,O1D = rs1
    o1 = solve(ODEProblem(rs1, [O3 => 20, O1D => 0], tspan, [], combinatoric_ratelaws=false),AutoTsit5(Rosenbrock23()))
    rs2 = FullChem().sys
    @unpack O3,O1D = rs2
    o2 = solve(ODEProblem(rs2, [O3 => 20, O1D => 10], tspan, [], combinatoric_ratelaws=false),AutoTsit5(Rosenbrock23()))
    test1 = o1[O3][end]-o2[O3][end]

    isapprox(test1, u_1, atol=0.001)
end

# Unit Test 2: OH sensitivity to O3
@test begin
    u_2 = 4.688821775289398e-9

    rs1 = FullChem().sys
    @unpack O3,OH = rs1
    o1 = solve(ODEProblem(rs1, [O3 => 20, OH => 0], tspan, [], combinatoric_ratelaws=false),AutoTsit5(Rosenbrock23()))
    rs2 = FullChem().sys
    @unpack O3,OH = rs2
    o2 = solve(ODEProblem(rs2, [O3 => 20, OH => 10], tspan, [], combinatoric_ratelaws=false),AutoTsit5(Rosenbrock23()))
    test2 = o1[O3][end]-o2[O3][end]

    isapprox(test2, u_2, atol=0.001)
end

# Unit Test 3: NO2 sensitivity to O3
@test begin
    u_3 = 3.61666252501891e-12

    rs1 = FullChem().sys
    @unpack O3,NO2 = rs1
    o1 = solve(ODEProblem(rs1, [O3 => 20, NO2 => 20], tspan, [], combinatoric_ratelaws=false),AutoTsit5(Rosenbrock23()))
    rs2 = FullChem().sys
    @unpack O3,NO2 = rs2
    o2 = solve(ODEProblem(rs2, [O3 => 20, NO2 => 40], tspan, [], combinatoric_ratelaws=false),AutoTsit5(Rosenbrock23()))
    test3 = o1[O3][end]-o2[O3][end]

    isapprox(test3, u_3, atol=0.001)
end

# Unit Test 4: HO2 sensitivity to O3
@test begin
    u_4 = 2.617817074224149e-10

    rs1 = FullChem().sys
    @unpack O3,HO2 = rs1
    o1 = solve(ODEProblem(rs1, [O3 => 20, HO2 => 0], tspan, [], combinatoric_ratelaws=false),AutoTsit5(Rosenbrock23()))
    rs2 = FullChem().sys
    @unpack O3,HO2 = rs2
    o2 = solve(ODEProblem(rs2, [O3 => 20, HO2 => 20], tspan, [], combinatoric_ratelaws=false),AutoTsit5(Rosenbrock23()))
    test4 = o1[O3][end]-o2[O3][end]

    isapprox(test4, u_4, atol=0.001)
end