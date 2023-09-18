using GasChem
using DifferentialEquations, ModelingToolkit

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
        18.86183158341872,
        4.2703017083622305e-5,
        0.00011583779715697954,
        2.1000002643778878e8,
        6.931895809649448,
        3.0681041903505473,
        1699.9160195725997,
        2.0291023439240294e-5,
        450.32824009985256,
        0.23464366198429554,
        274.55619790526,
        1.7535835570470253,
        0.042305366149692464,
        47.40417458095512,
        4.595825419044846,
        0.04328265843279768,
        2.546043110462823e-19,
        6.50621487186388,
    ]

    @parameters t
    rs = SuperFast(t).sys
    test0 = solve(
        ODEProblem(rs, [], tspan, []),
        Tsit5(),
        saveat = 10.0,
    )[37]

    test0 ≈ u_0
end

# Unit Test 1: DMS sensitivity
@test begin
    u_dms = [
        -0.1046858501730874,
        -1.9807122061648487e-5,
        -3.1727909652319645e-5,
        -0.8721528947353363,
        0.014189257766163976,
        -0.01418925776617952,
        0.01575996272640623,
        -7.991805556476097e-6,
        -0.06828251942368979,
        -0.02611920740062834,
        0.08376957284087894,
        0.001355592268675876,
        -0.008276678705923274,
        29.115790709520503,
        0.8842092904795944,
        0.010901015050459976,
        -1.1417920841140374e-21,
        0.005944387165240705,
    ]

    @parameters t
    rs1 = SuperFast(t).sys
    @unpack DMS = rs1
    o1 = solve(
        ODEProblem(rs1, [DMS => 76], tspan, []),
        Tsit5(),
        saveat = 10.0,
    )[37]
    rs2 = SuperFast(t).sys
    @unpack DMS = rs2
    o2 = solve(
        ODEProblem(rs2, [DMS => 46], tspan, []),
        Tsit5(),
        saveat = 10.0,
    )[37]
    test1 = (o1-o2)[1:18]

    isapprox(test1, u_dms, atol = 0.001)
end

# Unit Test 2: ISOP sensitivity
@test begin
    u_isop = [
        0.1938694842738684,
        1.118426789076167e-5,
        7.297766610801546e-5,
        -0.17319664359092712,
        -0.028103556998495094,
        0.028103556998510637,
        0.002062669243059645,
        3.3527857847824344e-5,
        0.009861447993898764,
        0.13707144364696552,
        0.019742953753620895,
        0.42323438434387284,
        0.0035247918766164593,
        0.062111863953681734,
        -0.06211186395370838,
        0.1230086530380752,
        2.0356787887724163e-21,
        -0.18169817774252905,
    ]

    @parameters t
    rs1 = SuperFast(t).sys
    @unpack ISOP = rs1
    o1 = solve(
        ODEProblem(rs1, [ISOP => 0.54], tspan, []),
        Tsit5(),
        saveat = 10.0,
    )[37]
    rs2 = SuperFast(t).sys
    @unpack DMS = rs2
    o2 = solve(
        ODEProblem(rs2, [ISOP => 0.13], tspan, []),
        Tsit5(),
        saveat = 10.0,
    )[37]
    test2 = (o1-o2)[1:18]

    isapprox(test2, u_isop, atol = 0.001)
end

#  Unit Test 3: NO2 sensitivity
@test begin
    u_no2 = [
        45.85362232850652,
        0.013443620451214309,
        -0.008251829362739181,
        273.2852625846863,
        31.238653656377366,
        58.761346343622655,
        -1.1626679626895111,
        -0.001647534014074231,
        4.278662562353361,
        0.10819944788865754,
        -4.841045205104081,
        -0.6656298071604302,
        0.17730281649833546,
        -0.03273608513133203,
        0.03273608513132453,
        -2.828906838233996e-6,
        5.907578300467002e-19,
        -2.1250201139384832,
    ]

    @parameters t
    rs1 = SuperFast(t).sys
    @unpack NO2, DMS = rs1
    o1 = solve(
        ODEProblem(
            rs1,
            [NO2 => 100.0, DMS => 0.1],
            tspan,
            [],
            combinatoric_ratelaws = false,
        ),
        Tsit5(),
        saveat = 10.0,
    )[37]
    rs2 = SuperFast(t).sys
    @unpack DMS = rs2
    o2 = solve(
        ODEProblem(rs2, [DMS => 0.1], tspan, []),
        Tsit5(),
        saveat = 10.0,
    )[37]
    test3 = (o1-o2)[1:18]

    isapprox(test3, u_no2, atol = 0.001)
end

# Unit Test 4: CO sensitivity
@test begin

    u_co = [
        -0.19386314144522743,
        -5.377970599299107e-7,
        -5.068308230428001e-5,
        0.0732552707195282,
        0.024938716038154674,
        -0.024938716038149344,
        -0.004746265539097294,
        -9.492630095183464e-7,
        0.016005078203818357,
        0.0027396748633061463,
        -449.25447756976206,
        0.0047022617615517515,
        0.0023874161918253714,
        -0.14282534835207628,
        0.14282534835213667,
        -0.0030328269313813147,
        -2.1776848880948197e-21,
        -0.2643181260957359,
    ]

    @parameters t
    rs1 = SuperFast(t).sys
    @unpack CO = rs1
    o1 = solve(
        ODEProblem(rs1, [CO => 50.0], tspan, []),
        Tsit5(),
        saveat = 10.0,
    )[37]
    rs2 = SuperFast(t).sys
    @unpack CO = rs2
    o2 = solve(
        ODEProblem(rs2, [CO => 500.0], tspan, []),
        Tsit5(),
        saveat = 10.0,
    )[37]
    test4 = (o1-o2)[1:18]

    isapprox(test4, u_co, atol = 0.001)
end



# Unit Test 5: CH4 sensitivity
@test begin
    u_ch4 = [
        0.006664853560593542,
        3.730835115360904e-7,
        2.303193302226881e-6,
        -0.0018870234489440918,
        -0.0009582041938998032,
        0.000958204193904244,
        299.9852733692228,
        1.1478698813863283e-6,
        0.014830814901642952,
        0.004579581197668298,
        0.0007497166716916581,
        0.00972420095777693,
        5.339295050454246e-5,
        0.0026560916389044564,
        -0.0026560916390483413,
        5.647267142151746e-5,
        7.053477993188078e-23,
        -0.0038984064496272453,
    ]

    @parameters t
    rs1 = SuperFast(t).sys
    @unpack CH4 = rs1
    o1 = solve(
        ODEProblem(rs1, [CH4 => 1900.0], tspan, []),
        Tsit5(),
        saveat = 10.0,
    )[37]
    rs2 = SuperFast(t).sys
    @unpack CH4 = rs2
    o2 = solve(
        ODEProblem(rs2, [CH4 => 1600.0], tspan, []),
        Tsit5(),
        saveat = 10.0,
    )[37]
    test5 = (o1-o2)[1:18]

    isapprox(test5, u_ch4, atol = 0.001)
end
