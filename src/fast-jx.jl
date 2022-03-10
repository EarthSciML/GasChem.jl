table_σ_o31D_jx = [[4.842, 4.922, 5.071, 5.228, 6.040, 6.803, 7.190, 7.549, 9.000, 8.989, 8.929, 9.000, 8.901, 4.130, 8.985*0.1, 6.782*0.1, 0.000, 0.000]*0.1 [4.843, 4.922, 5.072, 5.229, 6.040, 6.802, 7.189, 7.549, 9.000, 8.989, 8.929, 9.000, 8.916, 4.656, 1.417, 6.995*0.1, 0.000, 0.000]*0.1 [4.843, 4.922, 5.072, 5.229, 6.040, 6.805, 7.189, 7.550, 9.000, 8.989, 8.929, 9.000, 8.967, 5.852, 2.919, 7.943, 0.000, 0.000]*0.1]
# first column: 200K
# second column: 260K
# third column: 320K

function calc_J_o31D(actinic_flux,t)
	r = zeros(18)
	for i in 1:18
        if t <= 200
		    r[i] = actinic_flux[i] * table_σ_o31D_jx[i,1] * 1.000
        elseif t > 320
            r[i] = actinic_flux[i] * table_σ_o31D_jx[i,3] * 1.000
        elseif 200 < t <= 260
            k = (t-200)/(260-200)
            σ = (table_σ_o31D_jx[i,2]-table_σ_o31D_jx[i,1])*k + table_σ_o31D_jx[i,1]
            r[i] = actinic_flux[i] * σ * 1.000
        elseif 260 < t <= 320
            k = (t-260)/(320-260)
            σ = (table_σ_o31D_jx[i,3]-table_σ_o31D_jx[i,2])*k + table_σ_o31D_jx[i,2]
            r[i] = actinic_flux[i] * σ * 1.000
        end
	end
	return r
end
# for all 18 bins
# for one grid? -> given actinic flux at one bin
# O3 -> O2 + O(1D)
# superfast: O3 -> 2OH  include O(1D) + H2O -> 2OH


table_σ_H2O2_jx =[[2.325, 4.629, 5.394, 5.429, 4.447, 3.755, 3.457, 3.197, 5.346*0.1, 4.855*0.1, 3.423*0.1, 8.407*0.01, 5.029*0.01, 3.308*0.01, 2.221*0.01, 8.598*0.001, 1.807*0.0001,0]*10^-19.0 [2.325, 4.629, 5.394, 5.429, 4.447, 3.755, 3.457, 3.197, 5.465*0.1, 4.966*0.1, 3.524*0.1, 9.354*0.01, 5.763*0.01, 3.911*0.01, 2.718*0.01, 1.138*0.01, 2.419*0.0001, 0]*10^-19.0]
# first column : 200K
# second column: 300K

function calc_J_H2O2(actinic_flux,t)
	r = zeros(18)
	for i in 1:18
        if t <= 200
		    r[i] = actinic_flux[i] * table_σ_H2O2_jx[i,1] * 1.000
        elseif t >= 300
            r[i] = actinic_flux[i] * table_σ_H2O2_jx[i,2] * 1.000
        elseif 200 < t < 300
            k = (t-200)/(300-200)
            r[i] = (table_σ_H2O2_jx[i,2]-table_σ_H2O2_jx[i,1])*k + table_σ_H2O2_jx[i,1]
        end
	end
	return r
end
# for all bins
# for one grid? -> given actinic flux at one bin
# H2O2 -> OH + OH same with superfast

table_σ_CH2Oa_jx = [[0.0, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.143*0.1, 1.021, 1.269, 2.323, 2.498, 1.133, 2.183, 4.746*0.1, 0.000, 0.000]*10^-20.0 [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.147*0.1, 1.018, 1.266, 2.315, 2.497, 1.131, 2.189, 4.751*0.1, 0.000, 0.000]*10^-20.0]
# first column : 223K
# second column: 298K

function calc_J_CH2Oa(actinic_flux, t)
	r = zeros(18)
	for i in 1:18
        if t <= 223
		    r[i] = actinic_flux[i] * table_σ_CH2Oa_jx[i,1] * 1.000
        elseif t >= 298
            r[i] = actinic_flux[i] * table_σ_CH2Oa_jx[i,2] * 1.000
        elseif 223 < t < 298
            k = (t-223)/(298-223)
            r[i] = (table_σ_CH2Oa_jx[i,2]-table_σ_CH2Oa_jx[i,1])*k + table_σ_CH2Oa_jx[i,1]
        end
	end
	return r
end
# for all bins
# for one grid? -> given actinic flux at one bin
# CH2O -> H + HO2 + CO 
# superfast: CH2O -> CO + 2HO2; CH2O -> CO

table_σ_CH3OOH_jx = [0.000, 0.000, 0.000, 0.000, 0.000, 3.120, 2.882, 2.250, 2.716*0.1, 2.740*0.1, 2.143*0.1, 5.624*0.01, 3.520*0.01, 2.403*0.01, 1.697*0.01, 7.230*0.001, 6.973*0.0001, 0.000]*10^-19
#298K

function calc_J_CH3OOH(actinic_flux,t)
	r = zeros(18)
	for i in 1:18
		r[i,1] = actinic_flux[i] * table_σ_CH3OOH_jx[i] * 1.000
	end
	return r
end
# for all bins
# for one grid? -> given actinic flux at one bin
# CH3OOH -> OH + HO2 + CH2O same with superfast

table_σ_NO2_jx = [[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.835*0.1, 4.693*0.1, 7.705*0.1, 1.078, 1.470, 1.832, 2.181, 3.138, 4.321, 1.386*0.001]*10^-19.0 [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.313*0.1, 4.694*0.1, 7.553*0.1, 1.063, 1.477, 1.869, 2.295, 3.448, 4.643, 4.345*0.001]*10^-19] 
# first column: 200K
# second column:294K

function calc_J_NO2(actinic_flux, t)
    r = zeros(18)
	for i in 1:18
		if t <= 200
		    r[i] = actinic_flux[i] * table_σ_NO2_jx[i,1] * 1.000
        elseif t >= 294
            r[i] = actinic_flux[i] * table_σ_NO2_jx[i,2] * 1.000
        elseif 200 < t < 294
            k = (t-200)/(294-200)
            r[i] = (table_σ_NO2_jx[i,2]-table_σ_NO2_jx[i,1])*k + table_σ_NO2_jx[i,1]
        end
	end
	return r
end
# for all bins
#for one grid? -> given actinic flux at one bin
# NO2 -> NO + O 
# superfast: NO2 -> NO + O3  include O + O2 -> O3

act = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
# sample_actinic_flux, not realistic

#Unit Test 1 -- O3 -> O2 + O(1D)
function u1()
    n = 0
    m = ""
    act = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    function calc_J_o31Du1(actinic_flux)
        r = zeros(18,3)
        for i in 1:18, j in 1:3
            r[i,j] = actinic_flux[i] * table_σ_o31D_jx[i,j] * 1.000
        end
        return r
    end
    u_1 = calc_J_o31Du1(act)
    test1 = [calc_J_o31D(act,200) calc_J_o31D(act,260) calc_J_o31D(act,320)]
	if isapprox(test1, u_1, atol=0.001) == true
	    m = "unit test 1 passed"
	else 
		m = "unit test 1 failed"
	end
return m
end

#Unit Test 2 -- H2O2 -> OH + OH
function u2()
    n = 0
    m = ""
    act = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    function calc_J_H2O2u2(actinic_flux)
        r = zeros(18,2)
        for i in 1:18, j in 1:2
            r[i,j] = actinic_flux[i] * table_σ_H2O2_jx[i,j] * 1.000
        end
        return r
    end
    u_2 = calc_J_H2O2u2(act)
    test2 = [calc_J_H2O2(act,200) calc_J_H2O2(act,300)]
	if isapprox(test2, u_2, atol=0.001) == true
	    m = "unit test 2 passed"
	else 
		m = "unit test 2 failed"
	end
return m
end

#Unit Test 3 -- CH2O -> H + HO2 + CO
function u3()
    n = 0
    m = ""
    act = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    function calc_J_CH2Oau3(actinic_flux)
        r = zeros(18,2)
        for i in 1:18, j in 1:2
            r[i,j] = actinic_flux[i] * table_σ_CH2Oa_jx[i,j] * 1.000
        end
        return r
    end
    u_3 = calc_J_CH2Oau3(act)
    test3 = [calc_J_CH2Oa(act,223) calc_J_CH2Oa(act,298)]
	if isapprox(test3, u_3, atol=0.001) == true
	    m = "unit test 3 passed"
	else 
		m = "unit test 3 failed"
	end
return m
end

##Unit Test 4 -- CH3OOH -> OH + HO2 + CH2O
function u4()
    n = 0
    m = ""
    act = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    function calc_J_CH3OOHu4(actinic_flux)
        r = zeros(18)
        for i in 1:18,
            r[i,1] = actinic_flux[i] * table_σ_CH3OOH_jx[i,1] * 1.000
        end
        return r
    end
    u_4 = calc_J_CH3OOHu4(act)
    test4 = calc_J_CH3OOH(act,350)
	if isapprox(test4, u_4, atol=0.001) == true
	    m = "unit test 4 passed"
	else 
		m = "unit test 4 failed"
	end
return m
end

#Unit Test 5 -- NO2 -> NO + O3
function u5()
    n = 0
    m = ""
    act = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
    function calc_J_NO2u5(actinic_flux)
        r = zeros(18,2)
        for i in 1:18, j in 1:2
            r[i,j] = actinic_flux[i] * table_σ_NO2_jx[i,j] * 1.000
        end
        return r
    end
    u_5 = calc_J_NO2u5(act)
    test5 = [calc_J_NO2(act,200) calc_J_NO2(act,294)]
	if isapprox(test5, u_5, atol=0.001) == true
	    m = "unit test 5 passed"
	else 
		m = "unit test 5 failed"
	end
return m
end


