using GasChemMTK
using Test, OrdinaryDiffEq, ModelingToolkit

#=

# ╔═╡ 282a2582-59a2-45f5-abdc-dfd944006290


# ╔═╡ 17a038fa-b4f5-46fb-8107-5861b8fd2c17
md"""
### Example 1
Relation between ozone concentration and laltitude/temperature
"""

# ╔═╡ 7554ae3d-b04a-4596-923f-63b0f3b7d4cf
ozone(220,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30),ozone(250,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30),ozone(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30)  
# Higher temperature, less ozone. Reasonable for this model, because higher temperature means higher chemical rates, so more ozone reacts with other chemicals. Besides, higher altitude means lower temperature in stratrosphere, the result of lowering ozone concentration corresponds with the figure shown in the paper. But it could be a coincidence, because we didn't add relation between temperature and photolysis rates and ozone is produced by NO2 photolysis reaction in this model. 

# ╔═╡ 995d8bd8-0f40-4c85-8026-989a86eadcaa
ozone_(220,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)[1],ozone_(250,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)[1],ozone_(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)[1] 
# balanced concentration of ozone at 220K, 250K and 280K

# ╔═╡ 3a87bc48-c2ab-4633-8114-0aaf09123c11
ozone(220,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30),ozone(220,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0137, 0.00014, 8.9573*10^-6,30) 
# From the J values provided by CEOS-CHEM, the photolysis reaction rates are higher in January than in July, so if we assume that higher temperature causes that and the same happens in troposphere. Then less ozone produced under higher temperature, which also corresponds with the figure shown in the paper.
# But the range of increased concentration of ozone is not as large as shown in the paper, one possible explanation is that we didn't include O3 + hv --> 2OH in the model because we didn't have that reaction rate.

# ╔═╡ c634cc4b-2834-4053-9ba1-52102acc157c
ozone_(220,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)[1],ozone_(220,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0137, 0.00014, 8.9573*10^-6)[1]
# balanced concentration of ozone in January and July with different photolysis rates

# ╔═╡ d83e6078-4ea0-4bd5-98a2-79701fde26ad
md"""
### Example 2
Sensitivity Analysis of DMS concentration
"""

# ╔═╡ ed7014d5-8875-4aee-802d-eb03309ccfbf
# DMS: 120-200 nanogram/m3 --> 46-76ppb
ozone(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 46, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30), ozone(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 76, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30) 
# more DMS (30 ppb, 65%), more SO2 (0.9ppb, 45%) and DMS(29ppb,63%)a 

# ╔═╡ 3484f819-1d4f-4fbf-aca0-82d86ea1792b
ozone_(250,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 76, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)-ozone_(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 46, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)

# ╔═╡ 6f28714c-a7cd-4ebe-bb33-db7be4e32bd4
md"""
### Example 3
Sensitivity Analysis of ISOP concentration
"""

# ╔═╡ 1a1acaa7-49f1-4aa2-b07a-2a44e401e07e
# ISOP: 0.13-0.17ppb up to 0.54 with a lot of people
ozone(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.13, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30), ozone(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.54, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30) 
# More ISOP(0.41,315%), more CH3OOH(0.12 ppb, 80%), CH3O2(0.017 ppb, 70%), CH3OOH(0.47 ppb, 30%), less ozone (-0.16ppb, 1.6%)
# ISOP doesn't produce CH3OOH, while CH3O2 produces CH3OOH

# ╔═╡ 22968cc1-68cf-42f0-89da-ff9ecf9f06de
ozone_(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.54, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)-ozone_(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.13, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)

# ╔═╡ 47de614d-47a7-4ba8-ade0-6fe03f4d8823
md"""
### Example 4
Sensitivity Analysis of NO2 concentration
"""

# ╔═╡ e3bf5bbf-1bfb-4de0-86a5-13e3735728d4
#NO2： 10-100ppb (Normally the NO2 concentration in troposphere is up to 10 ppb, but near the earth, the concentration can be up to 100 ppb.)
ozone(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 100.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,100),ozone(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 50.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,60),ozone(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30)   
#NO2 is critical to the chemical system. More NO2(90ppb, 900%), more ozone(46ppb, 460%), NO（31ppb）， NO2(59ppb, 590%), CH2O(0.11ppb,73%), CH3OOH(-0.67ppb, 42%),H2O2(-2.2ppb, 94%)

# ╔═╡ 7fb99ec2-bdad-4971-9918-9c5a18a99c2a
ozone_(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 100.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 0.1, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)-ozone_(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 0.1, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)

# ╔═╡ 8d957275-6464-49c2-9af3-48b182526363
md"""
### Example 5
Sensitivity Analysis of CO concentration
"""

# ╔═╡ 3a52abeb-3d41-4200-80f5-8567431a7764
#CO:50-500ppb
ozone(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 50.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30),ozone(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 500.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30)  
# The system is not sensitive to CO

# ╔═╡ 238b7903-656a-4b47-b237-7a84a894ddba
ozone_(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 50.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)-ozone_(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 500.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)  

# ╔═╡ 96c84e41-38b3-4b8e-9630-817bf66f61ca
md"""
### Example 6
Sensitivity Analysis of CH4 concentration
"""

# ╔═╡ 3b839b5c-f34c-40b3-9b81-939035424bbb
#CH4: 1600-1900ppb
ozone(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1900.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30),ozone(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1600.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30)
#The system is not sensitive to CH4

# ╔═╡ 66295a42-8303-4953-bd73-af849ee8bd0b
ozone_(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1900.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)-ozone_(280,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1600.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6)   
=#
# ╔═╡ 3e08354e-30ed-428f-b57f-86efab111cbe

tspan = (0.0, 360.0)

# Unit Test 0: Base case
@test begin
    u_0 = [
        19.593406325470127, 4.0773001127828405e-5, 0.00011140076324841333, 2.100000266716814e8, 6.850932471077471, 3.14906752892254,
        1699.9057245133567, 1.8721327694499654e-5, 450.37140500459066, 0.22504141243149042, 274.51136291287236, 1.76706012081502,
        0.04813088503493186, 47.095364149128066, 4.904635850872059, 0.03716992293694597, 6.10272019987819]

    rs = superfast()
    test0 = solve(ODEProblem(rs, [], tspan, []), Tsit5(), saveat=10.0)[37]

    test0 ≈ u_0
end


# Unit Test 1: DMS sensitivity
@test begin
    u_dms = [-0.119651
        -1.91224e-5
        -3.0225e-5
        -0.952411
        0.0153949
        -0.0153949
        0.0182646
        -7.30047e-6
        -0.0788379
        -0.0257533
        0.0961116
        0.00144204
        -0.00974384
        29.033
        0.966983
        0.0109903
        0.00929881]


    rs1 = superfast()
    @unpack DMS = rs1
    o1 = solve(ODEProblem(rs1, [DMS => 76], tspan, []), Tsit5(), saveat=10.0)[37]
    rs2 = superfast()
    @unpack DMS = rs2
    o2 = solve(ODEProblem(rs2, [DMS => 46], tspan, []), Tsit5(), saveat=10.0)[37]
    test1 = o1 - o2

    isapprox(test1, u_dms, atol=0.001)
end

# Unit Test 2: ISOP sensitivity
@test begin
    u_isop = [0.162031
        9.9153e-6
        6.39096e-5
        0.0112869
        -0.0231624
        0.0231624
        0.00147313
        2.79975e-5
        0.0170129
        0.122394
        0.0217178
        0.468289
        0.00558679
        0.0440626
        -0.0440626
        0.104439
        -0.196991]

    rs1 = superfast()
    @unpack ISOP = rs1
    o1 = solve(ODEProblem(rs1, [ISOP => 0.54], tspan, []), Tsit5(), saveat=10.0)[37]
    rs2 = superfast()
    @unpack DMS = rs2
    o2 = solve(ODEProblem(rs2, [ISOP => 0.13], tspan, []), Tsit5(), saveat=10.0)[37]
    test2 = o1 - o2

    isapprox(test2, u_isop, atol=0.001)
end

#  Unit Test 3: NO2 sensitivity
@test begin
    u_no2 = [46.0443
        0.0131283
        -0.00821943
        274.906
        31.0163
        58.9837
        -1.15577
        -0.0016614
        4.02589
        0.114852
        -4.81128
        -0.671062
        0.178781
        -0.0322685
        0.0322685
        -2.22193e-6
        -2.24779]

    rs1 = superfast()
    @unpack NO2, DMS = rs1
    o1 = solve(ODEProblem(rs1, [NO2 => 100.0, DMS => 0.1], tspan, []), Tsit5(), saveat=10.0)[37]
    rs2 = superfast()
    @unpack DMS = rs2
    o2 = solve(ODEProblem(rs2, [DMS => 0.1], tspan, []), Tsit5(), saveat=10.0)[37]
    test3 = o1 - o2

    isapprox(test3, u_no2, atol=0.001)
end

# Unit Test 4: CO sensitivity
@test begin

    u_co = [-0.270119
        -1.07774e-6
        -4.94023e-5
        0.0691338
        0.0326809
        -0.0326809
        -0.00479081
        -8.59785e-7
        0.0154504
        0.0015244
        -449.16
        0.00507079
        0.00245152
        -0.143228
        0.143228
        -0.00262908
        -0.269762]

    rs1 = superfast()
    @unpack CO = rs1
    o1 = solve(ODEProblem(rs1, [CO => 50.0], tspan, []), Tsit5(), saveat=10.0)[37]
    rs2 = superfast()
    @unpack CO = rs2
    o2 = solve(ODEProblem(rs2, [CO => 500.0], tspan, []), Tsit5(), saveat=10.0)[37]
    test4 = o1 - o2

    isapprox(test4, u_co, atol=0.001)
end

# Unit Test 5: CH4 sensitivity
@test begin
    u_ch4 = [0.00678343
        -6.14978e-7
        2.55066e-6
        -0.00222483
        -0.000951862
        0.000951862
        299.983
        1.22716e-6
        0.0166712
        0.00478461
        0.00115292
        0.0110191
        8.07883e-5
        0.00349039
        -0.00349039
        6.41573e-5
        -0.00423393]


    rs1 = superfast()
    @unpack CH4 = rs1
    o1 = solve(ODEProblem(rs1, [CH4 => 1900.0], tspan, []), Tsit5(), saveat=10.0)[37]
    rs2 = superfast()
    @unpack CH4 = rs2
    o2 = solve(ODEProblem(rs2, [CH4 => 1600.0], tspan, []), Tsit5(), saveat=10.0)[37]
    test5 = o1 - o2

    isapprox(test5, u_ch4, atol=0.001)
end
