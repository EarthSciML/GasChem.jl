var documenterSearchIndex = {"docs":
[{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"CurrentModule = GasChem","category":"page"},{"location":"composing_models/#Composing-models","page":"Composing models","title":"Composing models","text":"","category":"section"},{"location":"composing_models/#Illustrative-Example","page":"Composing models","title":"Illustrative Example","text":"","category":"section"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"Here is the complete example of composing, visualizing and solving the SuperFast model and the Fast-JX model, with explanation to follow:","category":"page"},{"location":"composing_models/#Generating-ODESystems","page":"Composing models","title":"Generating ODESystems","text":"","category":"section"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"First, we need to build the two different components we want to compose together.","category":"page"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"@parameters t \n\nsf = superfast(t) \nfj = fast_jx(t) ","category":"page"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"Compose the SuperFast system and the Fast-JX system together.","category":"page"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"connect = compose_fastjx_superfast(fj,sf)","category":"page"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"This is now a differential-algebraic equation (DAE) of 18 variables. We can then define the resulting ODEProblem and send it over to DifferentialEquations.jl:","category":"page"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"tspan = (0.0, 3600*24*2) #simulating for 2 days\nsol = solve(ODEProblem(connect, [], tspan, [], combinatoric_ratelaws=false),Tsit5(), saveat=10.0)\n\nusing Plots\nplot(sol,ylims=(0,20),xlabel=\"Time (second)\", ylabel=\"concentration (ppb)\",legend=:outertopright)","category":"page"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"(Image: Example1 Graph)","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GasChem","category":"page"},{"location":"#GasChem","page":"Home","title":"GasChem","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This atmospheric chemical system model GasChem is built based on the Super Fast Chemical Mechanism, which is one of the simplest representations of atmospheric chemistry. It can efficiently simulate background tropheric ozone chemistry and perform well for those species included in the mechanism. The chemical equations used is included in the supporting table S2 of the paper, \"Evaluating simplified chemical mechanisms within present-day simulations of the Community Earth System Model version 1.2 with CAM4 (CESM1.2 CAM-chem): MOZART-4 vs. Reduced Hydrocarbon vs. Super-Fast chemistry\" (2018), Benjamin Brown-Steiner, Noelle E. Selin, Ronald G. Prinn, Simone Tilmes, Louisa Emmons, Jean-François Lamarque, and Philip Cameron-Smith.","category":"page"},{"location":"#Illustrative-Example","page":"Home","title":"Illustrative Example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Here is a simple example of generating, visualizing and solving the superfast model. We first define all parameters and variables.","category":"page"},{"location":"","page":"Home","title":"Home","text":"begin\n\t@parameters r1 [unit = u\"ppb/s\"]\n\t@parameters r2 [unit = u\"ppb/s\"]\n\t@parameters r3 [unit = u\"ppb/s\"]\n\t@parameters r6 [unit = u\"ppb/s\"]\n\t@parameters r7 [unit = u\"ppb/s\"]\n\t@parameters r9 [unit = u\"ppb/s\"]\n\t@parameters r11 [unit = u\"ppb/s\"]\n\t@parameters r12 [unit = u\"ppb/s\"]\n\t@parameters r13a [unit = u\"ppb/s\"]\n\t@parameters r13b [unit = u\"ppb/s\"]\n\t@parameters r14 [unit = u\"ppb/s\"]\n\t@parameters r15 [unit = u\"ppb/s\"]\n\t@parameters r17 [unit = u\"ppb/s\"]\n\t@parameters r21a [unit = u\"ppb/s\"]\n\t@parameters r21c [unit = u\"ppb/s\"]\n\t@parameters r22 [unit = u\"ppb/s\"]\n\t@parameters rr1 [unit = u\"ppb/s\"]\n\t@parameters rr2 [unit = u\"ppb/s\"]\n\t@parameters rr3 [unit = u\"ppb/s\"]\n\t@parameters rr4 [unit = u\"ppb/s\"]\n\t@parameters rr5 [unit = u\"ppb/s\"]\n\t@parameters rr6  [unit = u\"ppb/s\"]\n\t@parameters r4 [unit = u\"ppb/s\"]\n\t@parameters r5 [unit = u\"ppb/s\"]\n\t@parameters r10 [unit = u\"ppb/s\"]\n\t@parameters t [unit = u\"s\"]\n\t@variables O3(t) [unit = u\"ppb\"]\n\t@variables OH(t) [unit = u\"ppb\"]\n\t@variables HO2(t) [unit = u\"ppb\"]\n\t@variables O2(t) [unit = u\"ppb\"]\n\t@variables H2O(t) [unit = u\"ppb\"]\n\t@variables NO(t) [unit = u\"ppb\"]\n\t@variables NO2(t) [unit = u\"ppb\"]\n\t@variables CH4(t) [unit = u\"ppb\"]\n\t@variables CH3O2(t) [unit = u\"ppb\"]\n\t@variables CH2O(t) [unit = u\"ppb\"]\n\t@variables CO(t) [unit = u\"ppb\"]\n\t@variables CH3OOH(t) [unit = u\"ppb\"]\n\t@variables CH3O(t) [unit = u\"ppb\"]\n\t@variables DMS(t) [unit = u\"ppb\"]\n\t@variables SO2(t) [unit = u\"ppb\"]\n\t@variables ISOP(t) [unit = u\"ppb\"]\n\t@variables H2O2(t) [unit = u\"ppb\"]\nend","category":"page"},{"location":"","page":"Home","title":"Home","text":"Then we define the chemical reaction model.","category":"page"},{"location":"","page":"Home","title":"Home","text":"rxs = [Reaction(r1, [O3,OH], [HO2,O2], [1,1], [1,1]) #O3 + OH --> HO2 + O2\n       Reaction(r2, [HO2,O3], [O2,OH],[1,1],[2,1])#HO2 + O3 --> 2O2 + OH\n\t   Reaction(r3, [HO2,OH], [HO2,O2],[1,1],[1,1]) #HO2 + OH --> H2O + O2\n\t   Reaction(r6, [NO,O3], [NO2,O2],[1,1],[1,1]) #NO + O3 --> NO2 + O2\n\t   Reaction(r7, [HO2,NO], [NO2,OH],[1,1],[1,1]) #HO2 + NO --> NO2 + OH\n\t   Reaction(r9, [CH4, OH], [CH3O2, H2O],[1,1],[1,1]) #CH4 + OH --> CH3O2 + H2O\n\t   Reaction(r11, [CH2O,OH], [CO,H2O,HO2],[1,1],[1,1,1]) \n\t   #CH2O + OH --> CO + H2O + HO2\n\t   Reaction(r12, [CH3O2,HO2], [CH3OOH,O2],[1,1],[1,1]) \n\t   #CH3O2 + HO2 --> CH3OOH + O2\n\t   Reaction(r13a, [CH3OOH,OH], [CH3O2,H2O],[1,1],[1,1]) \n\t   #CH3OOH + OH --> CH3O2 + H2O\n\t   Reaction(r13b, [CH3OOH,OH], [CH3O,H2O,OH],[1,1],[1,1,1]) \n\t   #CH3OOH + OH --> CH3O + H2O + OH\n\t   Reaction(r14, [CH3O2,NO], [CH2O,HO2,NO2],[1,1],[1,1,1])\n\t   #CH3O2 + NO --> CH2O + HO2 + NO2\n\t   Reaction(r15, [CH3O2,CH3O2], [CH2O,H2O],[10,10],[20,8])\n\t   #10CH3O2 + 10CH3O2 --> 20CH2O + 8HO2\n\t   Reaction(r17, [DMS,OH], [SO2],[1,1],[1]) \n\t   #DMS + OH --> SO2\n\t   Reaction(r21a, [ISOP,OH], [CH3O2],[1,1],[2])\n\t   #ISOP +OH --> 2CH3O2\n\t   Reaction(r21c, [ISOP,OH], [ISOP,OH],[2,2],[2,1])\n\t   #2ISOP + 2OH --> 2ISOP + OH\n\t   Reaction(r22, [ISOP,O3], [CH2O,CH3O2,HO2,CO], [1,1.0], [0.87,1.86,0.06,0.05])\n\t   #ISOP + O3 --> 0.87CH2O + 1.86CH3O2 + 0.06HO2 + 0.05CO\n\t   Reaction(rr2, [H2O2], [OH], [1], [2])\n\t   #H2O2 --> 2OH\n\t   Reaction(rr3, [NO2], [NO,O3], [1], [1,1])\n\t   #NO2 --> NO + O3\n\t   Reaction(rr4, [CH2O], [CO,HO2], [1], [1,2])\n\t   #CH2O --> CO + 2HO2\n\t   Reaction(rr5, [CH2O], [CO], [1], [1])\n\t   #CH2O --> CO\n\t   Reaction(rr6, [CH3OOH], [CH2O,HO2,OH], [1], [1,1,1])\n\t   #CH3OOH --> CH2O + HO2 + OH\n\t   Reaction(r4, [HO2],[H2O2,O2],[2],[1,1])\n\t   #HO2 + HO2 = H2O2 + O2\n\t   Reaction(r5,[OH,H2O2],[H2O,HO2],[1,1],[1,1])\n\t   #OH + H2O2 = H2O + HO2\n\t   Reaction(r10,[OH,CO],[HO2],[1,1],[1])\n\t   #OH + CO = HO2\n\t   ] ","category":"page"},{"location":"","page":"Home","title":"Home","text":"@named rs = ReactionSystem(rxs, t)","category":"page"},{"location":"","page":"Home","title":"Home","text":"which in Jupyter notebooks will give the figure that represents the reation networks.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Chemical Network Graph)","category":"page"},{"location":"","page":"Home","title":"Home","text":"We build a function that can predict the change of the concentration of the chemicals in the superfast mechanism with input of temperature, the initial concentrations and reaction rates of photolysis reactions. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"function ozone(T_,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,z)\n\tu₀map_ = [O3 => a, OH => b, HO2 => c, O2 => d, H2O => e, NO => f, NO2 => g, CH4 => h, CH3O2 => i, CH2O => j, CO => k, CH3OOH => l, CH3O => m, DMS => n, SO2 => o, ISOP => p, H2O2 => q]\n\tparammap_ = [r1 => 1.7*10^(-12)*exp(-940/T_)*2.46*10^10, \n\tr2 => 1.0*10^(-14)*exp(-490/T_)*2.46*10^10,\n\tr3 => 4.8*10^(-11)*exp(250/T_)*2.46*10^10,\n\tr6 => 3.0*10^(-12)*exp(-1500/T_)*2.46*10^10,\n\tr7 => 3.5*10^(-12)*exp(250/T_)*2.46*10^10,\n\tr9 => 2.45*10^(-12)*exp(-1775/T_)*2.46*10^10,\n\tr11 => 5.50*10^(-12)*exp(125/T_)*2.46*10^10,\n\tr12 => 4.10*10^(-13)*exp(750/T_)*2.46*10^10,\n\tr13a => 2.70*10^(-12)*exp(200/T_)*2.46*10^10,\n\tr13b => 1.10*10^(-12)*exp(200/T_)*2.46*10^10,\n\tr14 => 2.80*10^(-12)*exp(300/T_)*2.46*10^10,\n\tr15 => 9.50*10^(-14)*exp(390/T_)/10*2.46*10^10,\n\tr17 => 1.10*10^(-11)*exp(-240/T_)*2.46*10^10,\n\tr21a => 2.70*10^(-11)*exp(390/T_)*2.46*10^10,\n\tr21c => 2.70*10^(-11)*exp(390/T_)/2*2.46*10^10,\n\tr22 => 5.59*10^(-15)*exp(-1814/T_)*2.46*10^10,\n\trr2 => r,\n\trr3 => s,\n\trr4 => t,\n\trr5 => t,\n\trr6 => u,\n\tr4 => 3.0*10^(-13)*exp(460/T_)*2.46*10^10,\n\tr5 => 1.8*10^(-12)*2.46*10^10,\n\tr10 => 1.5*10^(-13)*2.46*10^10\n]\n\toprob_ = ODEProblem(rs, u₀map_, tspan, parammap_)\n\tsol_ = solve(oprob_, Tsit5(), saveat=10.)\n\tplot(sol_,ylim=[0,z], lw=2)\nend","category":"page"},{"location":"","page":"Home","title":"Home","text":"For example, below is the graph at 220K, when the initial concentrations are as followed: O3 => 10.0, OH => 10.0, HO2 => 10.0, O2 => 2.1(10^8), H2O => 450.0, NO => 0.0, NO2 => 10.0, CH4 => 1700, CH3O2 => 0.01, CH2O => 0.15, CO => 275, CH3OOH => 1.6, CH3O => 0.0, DMS => 50, SO2 => 2.0, ISOP => 0.15, H2O2 => 2.34, and the potolysis rates for rr2 to rr6 are 1.0097 10^-5, 0.0149, 0.00014, 0.00014 and 8.9573* 10^-6.","category":"page"},{"location":"","page":"Home","title":"Home","text":"ozone(220,10.0, 10.0, 10.0, 2.1*(10^8), 450.0, 0.0, 10.0, 1700.0, 0.01, 0.15, 275.0, 1.6, 0.0, 50, 2.0, 0.15, 2.34, 1.0097*10^-5, 0.0149, 0.00014, 8.9573*10^-6,30)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Example1 Graph)","category":"page"},{"location":"","page":"Home","title":"Home","text":"superfast","category":"page"},{"location":"#GasChem.superfast","page":"Home","title":"GasChem.superfast","text":"superfast()\n\nThis atmospheric chemical system model is built based on the Super Fast Chemical Mechanism, which is one of the simplest representations of atmospheric chemistry. It can efficiently simulate background tropheric ozone chemistry and perform well for those species included in the mechanism. The chemical equations we used is included in the supporting table S2 of the paper:\n\n\"Evaluating simplified chemical mechanisms within present-day simulations of the Community Earth System Model version 1.2 with CAM4 (CESM1.2 CAM-chem): MOZART-4 vs. Reduced Hydrocarbon vs. Super-Fast chemistry\" (2018), Benjamin Brown-Steiner, Noelle E. Selin, Ronald G. Prinn, Simone Tilmes, Louisa Emmons, Jean-François Lamarque, and Philip Cameron-Smith.\n\nThe input of the function is Temperature, concentrations of all chemicals, and reaction rates of photolysis reactions \n\nExample\n\nusing OrdinaryDiffEq, Plots\nrs = superfast()\nsol = solve(ODEProblem(rs, [], (0,360), [], combinatoric_ratelaws=false), Tsit5())\nplot(sol)\n\nWe set combinatoric_ratelaws=false because we are modeling macroscopic rather than microscopic behavior. See here and here.\n\n\n\n\n\n","category":"function"}]
}
