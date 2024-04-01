export date, Emission

function date(t)
    return Dates.unix2datetime(t)
end
@register_symbolic date(t)

"""
This Emission function serves as a fundamental base model for emissions and can be integrated with various datasets in the extension module.
"""
struct Emission <: EarthSciMLODESystem
    sys::ODESystem
    function Emission(t)

        D = Differential(t)
        @variables NO(t) [unit = u"nmol/mol"]
		@variables NO2(t) [unit = u"nmol/mol"]
        @variables CH2O(t) [unit = u"nmol/mol"]
		@variables CH4(t) [unit = u"nmol/mol"]
		@variables CO(t) [unit = u"nmol/mol"]
		@variables SO2(t) [unit = u"nmol/mol"]
		@variables ISOP(t) [unit = u"nmol/mol"]

        eqs = [
            D(NO) ~  0
            D(NO2) ~ 0
			D(CH2O) ~ 0
			D(CH4) ~ 0
			D(CO) ~ 0
			D(SO2) ~ 0
			D(ISOP) ~ 0
        ]

        new(ODESystem(eqs, t, [NO, NO2, CH2O, CH4, CO, SO2, ISOP], []; name=:Emission))
    end
end 