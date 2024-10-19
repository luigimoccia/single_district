# LOG of cost parameter files

- The main scenarios are based on the following file
'costs_2050_plus_others.csv'

## parametric analysis
	- battery capacity cost, Scenario A (high discount rate), by varying 'short_duration_storage_energy,investment' from the base value of 115.7 to 80.
	The cost file is 
		'costs_2050_plus_others_par_battery_A.csv'
		
	- battery capacity cost, scenario \dot A (low discount rate), by varying 'short_duration_storage_energy,investment' from the base value of 115.7 to 35.
	The cost file is 
	'costs_2050_plus_others_par_battery_A_low_discount.csv'

	- nuclear investment cost, Scenario \dot A with and without bio (low discount rate, biogenic gases), by varying 'nuclear,investment' from the base value of 6195 EUR/kW_e to 2200 EUR/kW_e.
	The cost files are 
	'costs_2050_plus_others_par_nuclear_A_low_discount_bio.csv'
	'costs_2050_plus_others_par_nuclear_A_low_discount.csv'
	
	- DAC, Scenario \dot A_bio (low discount rate, biogenic gases), by varying 'direct air capture,investment' from the base value of 4000000 EUR/(tCO2/h) to 12000000 EUR/(tCO2/h).
	The cost file is 
	'costs_2050_plus_others_DAC_high.csv'
	
	- Stagnant wind and solar, 2023 global avg as reported by IRENA2024: onshore wind 1160 €/kW, +20% wrt base case, and solar utility tracking at 752 €/kW, 2.7 times than base level
	The cost file is
	'costs_2050_plus_others_par_stagnant_wind_and_solar.csv'