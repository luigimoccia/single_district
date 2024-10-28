#!/usr/bin/env python
# coding: utf-8

# # PyPSA model of a 1 GWe decarbonized district
# 
# Features:
# - Single zone and two possible demands of energy carriers:
#     - electricty only
#     - electricity and hydrogen
# - Electricity demand under two variants:
#     - constant, at a given value per hour 
#     - dynamic, variable in time with a hourly profile
# - Renewables' capacity factors as per Renewable Ninja, four files are needed:
#     - solar utility scale, ground based PV
#     - solar rooftop, high capacity factor
#     - solar rooftop, low capacity factor
#     - onshore wind
# - Nuclear
# - Thermoelectric generation:
#       - OCGT
#..........- CCGT (disabled in main scenarios)
# - Storage and primary generation technologies can be activated/deactivated
# - Main index: overgeneration vs LCOE vs CO2 wrt to CCGT/OCGT-fossil
# - Thermoelectric generators as Links from a multi-source gas bus
# - The sources for the gas bus are grouped as primary and secondary sources
#     - primary:
#         - fossil gas
#         - bio-CH4
#     - secondary:
#         - e-CH4
#         - H2 and this link can be activated or deactivated
# - Multiple technologies available for storage
#     - Directly connected to the electricity bus:
#         - PHS
#         - A-CAES (disabled in the main scenarios)
#         - Li-ion
#     - Indirect storage via secondary e-fuel
#         - Electrolysis
#         - Methanation
# - Multiple sources of CO2 for methanation
#     - exogenous biogenic source at a fixed price
#     - DAC
# - CO2 limit
#   (NB: For the first run set a sufficiently high value such that the CO2 constraint is not binding
#   and the CO2 shadow price is near-zero)
# 
# Requirements: 
#     - 1 csv file for the techno-economic parameters
#   - 1 csv file for the time-varying electricity demand (not required in the main scenarios that use a constant el. demand)
#   - 4 csv file with wind and solar availability profiles
#     -  1 yaml file for scenario-specific options
# - Sub-folders with input data or required for output:
#     - data (folder with the techno-economic parameters and electricity profile)
#     - data_ren_ninja (folder with the wind and solar CF profiles)
#     - scenarios_and_results (folder of the scenario yaml files and the output files)
# - Functions from custom python libraries:
#     - my_pypsa_functions.py (some modifications to standard PyPSA procedures for processing techno-economic coefficients)
#
# HOW TO LAUNCH
# from the command line:
# python single_district.py -s Scenario_A


import pypsa
import pandas as pd
import numpy as np
from scipy import stats

import yaml

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

# Here the required custom functions
from my_pypsa_functions import load_tech_parameters

# Hard wired coefficients and parameters, these do not changes between scenarios and variants
# CO2 emission factor of fossil gas
fossilCH4_CO2_emission_factor = 0.19 # t_CO2/MWh_th
# CF max
max_capacity_factor_H2 = 0.9
max_capacity_factor_nuclear = 0.9
# Solver specifics
name_solver = 'gurobi'
epsilon = 0.001



# parser code for the executable version
import argparse
parser = argparse.ArgumentParser(description='PyPSA model, single district, v4')
parser.add_argument('-s','--config', help='Scenario config yaml file name (without .yaml extension)',required=True)
args = parser.parse_args()
name_file_run_results = args.config


# Load config file
# Directory of scenario config yaml files and their results
path_file_run_results = "./scenarios_and_results/"


with open(path_file_run_results + name_file_run_results + '.yaml', 'r') as file:
    my_config = yaml.safe_load(file)


# Switch the generator and storage tech set by these boolean values 

#  Generation=
with_onshorewind     = my_config['with_onshorewind']
with_solar_utility_tracker = my_config['with_solar_utility_tracker']
with_solar_rooftop   = my_config['with_solar_rooftop']
with_solar_roof_best = my_config['with_solar_roof_best']  # if my_config[''] active a rooftop profile with a smaller CF
with_OCGT            = my_config['with_OCGT']
with_CCGT            = my_config['with_CCGT']  # CCGT has a very small system benefit, for simplicity can be deactivated
with_nuclear         = my_config['with_nuclear']
with_nuc_fleet       = my_config['with_nuc_fleet']  # If my_config[''] impose the capacity factor as a global constraint, otherwise devise a CF array with summer refueling & maint.
with_bioCH4          = my_config['with_bioCH4']
with_mft             = my_config['with_mft']   # allow OCGT and CCGT to receive as input H2, i.e. multi-fuel turbine

# Storage=
with_SD              = my_config['with_SD']    # Li-ion battery, SD short duration
with_LD1             = my_config['with_LD1']   # ACAES
with_PHS             = my_config['with_PHS']
with_FC              = my_config['with_FC']    # Fuel Cell
with_SynCH4          = my_config['with_SynCH4']    # Methanation 
with_DAC             = my_config['with_DAC']
with_ex_CO2          = my_config['with_ex_CO2']   # Exogenous source of CO2 at a fixed price  

# Demand=
#Electricity demand profile
with_constand_d      = my_config['with_constand_d']   # if my_config[''] a constant 1GW demand otherwise a variable demand with the 2015 profile scaled down to 8.76 TWh/y
with_demand_H2       = my_config['with_demand_H2']  # if larger than zero add a constant demand of H2 with this value multiplied by the average electric hourly load

# Make a flag True only if with_demand_H2 is larger than zero
with_demand_H2_flag = False
if with_demand_H2 > 0. :
    with_demand_H2_flag = True
    
Hourly_Avg_demand    = my_config['Hourly_Avg_demand'] # MWh
file_demand_path     = my_config['file_demand_path']

# Other=
do_parametric_runs = my_config['do_parametric_runs'] # If a zero CO2 constraint is not active and if this parameter is active then it activates a parametric multi-run
co2_limit = my_config['co2_limit']

# Fuel and feedstock costs of exogenous inputs 
fuel_cost_fossilCH4 = my_config['fuel_cost_fossilCH4']  # €/MWh_th
fuel_cost_bioCH4    = my_config['fuel_cost_bioCH4']  # €/MWh_th
co2_feedstock_cost  = my_config['co2_feedstock_cost']  # €/t_CO2   input cost for methanization, this is exogenous


# Files 

name_file_tech_parameters = my_config['name_file_tech_parameters']

# Four Renewable ninja csv files, one onshore wind and three for solar
# Please verify that the parameter year_ninja expresses the same year of the capacity factor files, otherwise an error will occur
year_ninja                = my_config['year_ninja']

wind_file_name_and_path   = my_config['wind_file_name_and_path']

solar_file_name_and_path  = my_config['solar_file_name_and_path']

solar_file_name_and_path2 = my_config['solar_file_name_and_path2']

solar_file_name_and_path3 = my_config['solar_file_name_and_path3']

# Switch bewteen different H2 tec sets
type_of_H2_storage    = my_config['type_of_H2_storage']
type_of_electrolysis  = my_config['type_of_electrolysis']
type_of_H2_compressor = my_config['type_of_H2_compressor']


# construct the dates for the year. It has to match the year in the ren ninja capacity factor files!
start_date = '%d-01-01 00:00' % year_ninja
end_date   = '%d-12-31 23:00' % year_ninja
print(start_date, end_date)


# Define a year with hourly resolution and constant demand

#network = pypsa.Network(override_component_attrs=override_component_attrs)
network = pypsa.Network()
# Set the hours in the same year used for the CF wind and solar profiles
hours_in_year = pd.date_range(start_date, end_date, freq='H')
network.set_snapshots(hours_in_year)
# Num of hours in a year, to account for leap years
Num_hours = len(hours_in_year)
# Add the electricty bus where the electric load will be assigned
network.add("Bus","AC", carrier = "AC")

# Add a H2 bus where an optional load may be assigned
network.add("Bus", "H2", carrier = "H2")


# Define a constant load
# Default, a constant 1 GW
Yearly_demand = Hourly_Avg_demand * Num_hours /1.e+6 # TWh
print('An average hourly demand (MWh):', Hourly_Avg_demand)
print('For a yearly demand (TWh):', Yearly_demand)

# Add the electric load to the electricity bus
if with_constand_d:
    network.add("Load",
            "load", 
            bus="AC", 
            p_set = Hourly_Avg_demand, # This is a constant profile 
           )
else:
    #Alternatively use a scaled version of the 2015 Italy demand
    df_demand = pd.read_csv(file_demand_path, index_col=0)
    # df_demand.index = pd.to_datetime(df_demand.index) #change index to datatime
    df_demand.index = network.snapshots
    my_demand = df_demand.squeeze()
    network.add("Load",
            "load", 
            bus="AC", 
            p_set = my_demand, # This is a demand profile as the 2015 elec. demand
           )
    
# If with_demand_H2_flag = True add the H2 load
if with_demand_H2_flag:
    network.add("Load",
            "load_H2", 
            bus="H2", 
            p_set = with_demand_H2 * Hourly_Avg_demand, # This is a constant profile 
           )




Nyears = network.snapshot_weightings.objective.sum() / Num_hours

tech_parameters = load_tech_parameters(name_file_tech_parameters, my_config["costs"], Nyears)
#print(tech_parameters)
#print(tech_parameters["land_spacing"]["onwind"])
print(tech_parameters["capital_cost"]["electrolysis_IEA"])


df_onshorewind = pd.read_csv(wind_file_name_and_path, sep=',', comment='#', index_col=0)
df_onshorewind.index = pd.to_datetime(df_onshorewind.index)
#print(df_onshorewind)
CF_wind = df_onshorewind['electricity'][[hour.strftime("%Y-%m-%d %H:%M:%S") for hour in network.snapshots]]
CF_wind_y = CF_wind.mean() * 100.
print("CF Wind %.1f%%" % CF_wind_y)


# Add PV utility (1-axis tracking) from a renewable ninja csv file
df_solar = pd.read_csv(solar_file_name_and_path, sep=',', comment='#', index_col=0)
df_solar.index = pd.to_datetime(df_solar.index)
#print(df_solar)
CF_solar = df_solar['electricity'][[hour.strftime("%Y-%m-%d %H:%M:%S") for hour in network.snapshots]]
CF_solar_y = CF_solar.mean() * 100.
print("CF solar %.1f%%" % CF_solar_y)

# Add PV rooftop from a renewable ninja csv file
df_solar2 = pd.read_csv(solar_file_name_and_path2, sep=',', comment='#', index_col=0)
df_solar2.index = pd.to_datetime(df_solar2.index)
#print(df_solar)
CF_solar2 = df_solar2['electricity'][[hour.strftime("%Y-%m-%d %H:%M:%S") for hour in network.snapshots]]
CF_solar2_y = CF_solar2.mean() * 100.
print("CF solar %.1f%%" % CF_solar2_y)

df_solar3 = pd.read_csv(solar_file_name_and_path3, sep=',', comment='#', index_col=0)
df_solar3.index = pd.to_datetime(df_solar3.index)
#print(df_solar)
CF_solar3 = df_solar3['electricity'][[hour.strftime("%Y-%m-%d %H:%M:%S") for hour in network.snapshots]]
CF_solar3_y = CF_solar3.mean() * 100.
print("CF solar %.1f%%" % CF_solar3_y)





# # A naive PV-Wind mix 
# This may be used after the optimization to highlight the usefulness of techno-economic opt.
naive_mix =   CF_wind +  CF_solar
naive_mix_clipped = np.clip(naive_mix, 0, 1)
naive_over_generation_ratio = (naive_mix.sum() - naive_mix_clipped.sum())/(naive_mix_clipped.sum())
naive_mix_demand_fraction = naive_mix_clipped.mean()
plt.plot(np.sort(naive_mix_clipped)[::-1], color='blue', ls=':', label='PV-Wind naïve mix')
plt.legend(fancybox=True, loc='best')
plt.tight_layout()
print("Fraction of the demand covered by mix %.3f" % naive_mix_demand_fraction)
print("Overgeneration fraction %.3f" % naive_over_generation_ratio)



# add the different carriers, only fossil gas emits CO2
network.add("Carrier", "fossilCH4", co2_emissions=fossilCH4_CO2_emission_factor) # in t_CO2/MWh_th
network.add("Carrier", "bioCH4", co2_emissions=0.0) # in t_CO2/MWh_th
network.add("Carrier", "CH4") # it can be from biogenic, synthesis, or fossil origin
network.add("Carrier", "gas") # general gas accounted in MWh 
network.add("Carrier", "onshorewind")
network.add("Carrier", "solar")
network.add("Carrier", "nuclear")
network.add("Carrier", "heat")
network.add("Carrier", "heat_pump")





# Store in this dict a map between the PyPSA keys and the keys of the parameter dataset
# dict with tech components that require a capital investment
dict_pypsa_par_id = {
    "onshorewind" : 'onwind',
    "solar"       : "solar_utility_tracker",
    "solar_rooftop": "solar-rooftop commercial",
    'nuclear'     : "nuclear",
    "OCGT"        : "OCGT",
    "CCGT"        : "CCGT",
    "SD_store"    : "short_duration_storage_energy",
    "SD_in"       : "short_duration_storage_power",
    "LD1_store"   : "long_duration_1_storage_energy",
    "LD1_in"      : "long_duration_1_storage_power_charging",
    "LD1_out"     : "long_duration_1_storage_power_discharging",
    "PHS_store"   : "long_duration_2_storage_energy",
    "PHS_in"      : "long_duration_2_storage_power",
    "compH2_store": type_of_H2_storage,
    "H2_in"       : type_of_electrolysis,
    "compH2_in"   : type_of_H2_compressor,
    "FC"          : "fuel cell",
    "SynCH4_in"   : "methanation",
    "heat_pump_mid" : "industrial heat pump medium temperature",
    "CO2_storage_store"   : "CO2 storage tank",
    "dac"         : "direct air capture",

}


# Add onshore wind generator
this_id = "onshorewind"
network.add("Generator",
            this_id,
            bus="AC",
            p_nom_extendable = with_onshorewind,
            carrier="onshorewind",
            capital_cost = tech_parameters.at[dict_pypsa_par_id[this_id], 'capital_cost'],
            marginal_cost = tech_parameters.at[dict_pypsa_par_id[this_id], 'marginal_cost_no_fuel'],
            p_max_pu = CF_wind)


# Add PV utility (tracking) 
this_id = "solar"
network.add("Generator",
            this_id,
            bus="AC",
            p_nom_extendable = with_solar_utility_tracker,
            carrier="solar",
            capital_cost = tech_parameters.at[dict_pypsa_par_id[this_id], 'capital_cost'],
            marginal_cost = tech_parameters.at[dict_pypsa_par_id[this_id], 'marginal_cost_no_fuel'], #Small cost to prefer curtailment to destroying energy in storage, this is smaller than wind VOM, so solar curtails before wind
            p_max_pu = CF_solar
           )

# Add PV rooftop
this_id = "solar_rooftop"
if with_solar_roof_best :
    CF_solar_roof = CF_solar2
else:
    CF_solar_roof = CF_solar3
    
network.add("Generator",
            this_id,
            bus="AC",
            p_nom_extendable = with_solar_rooftop,
            carrier="solar",
#            capital_cost = tech_parameters.at['solar-utility', 'capital_cost'],
#            marginal_cost = tech_parameters.at['solar-utility', 'marginal_cost_no_fuel'],
            capital_cost = tech_parameters.at[dict_pypsa_par_id[this_id], 'capital_cost'],
            marginal_cost = tech_parameters.at[dict_pypsa_par_id[this_id], 'marginal_cost_no_fuel'], #Small cost to prefer curtailment to destroying energy in storage, this is smaller than wind VOM, so solar curtails before wind
            p_max_pu = CF_solar_roof,
           # p_nom_min = 1000.  # to test mix of utility and rooftop
           )



# Add nuclear generators

# Define CF for a nuclear reactor with average availibity in the year equal to 90%, i.e. 7884 hours
# Assume a continuos block of hours as not available for the remaining hours 

CF_nuclear=np.ones(Num_hours)
# Assume the nuclear reactor is 10% of the time in manintenance/refuelling during the mid of the year, best case for nuclear and PV mix
CF_nuclear[3942:4818]=0.
this_id ="nuclear"

if with_nuc_fleet:
    network.add("Generator",
            this_id,
            bus="AC",
            p_nom_extendable=with_nuclear,
            carrier="nuclear",
            capital_cost = tech_parameters.at[dict_pypsa_par_id[this_id], 'capital_cost'],
            marginal_cost = tech_parameters.at[dict_pypsa_par_id[this_id], 'marginal_cost_no_fuel'],#here no_fel means no *fossil* fuel
    )
else: 
    network.add("Generator",
            this_id,
            bus="AC",
            p_nom_extendable=with_nuclear,
            carrier="nuclear",
            capital_cost = tech_parameters.at[dict_pypsa_par_id[this_id], 'capital_cost'],
            marginal_cost = tech_parameters.at[dict_pypsa_par_id[this_id], 'marginal_cost_no_fuel'],#here no_fuel means no *fossil* fuel
            p_max_pu = CF_nuclear
    )
# Eventually add a second generator with unavalibility in a different period than the previous reactor
# This only under the no-storage no-transmission constraints—some pretend that for renewables, so why not for nuclear? They will not like the results!

#CF_nuclear_B = np.ones(Num_hours)
#CF_nuclear_B[:876] = 0.

#network.add("Generator",
#            "nuclear_B",
#            bus="AC",
#            p_nom_extendable=with_nuclear,
#            carrier="nuclear",
#            capital_cost = tech_parameters.at['nuclear', 'capital_cost'],
#            marginal_cost = tech_parameters.at['nuclear', 'marginal_cost_no_fuel'],#here no_fel means no *fossil* fuel
#            p_max_pu = CF_nuclear_B
#           )


# Create fossilCH4 as a bus, and a store

network.add("Bus", "fossilCH4", carrier = "fossilCH4")

# Add a depletable store of fossil methane, in MWh_th
network.add("Store",
      "fossilCH4_store",
      bus = "fossilCH4",
      e_nom_extendable = False,
      e_cyclic = False,
      e_nom = 250e6, #MWh_th # A sufficient large value so that this is not binding
      e_initial = 250e6, #MWh_th
    )

# Add links from fossilCH4 to the general CH4 bus
network.add("Bus", "CH4", carrier = "CH4")

network.add("Link",
      "from_fossilCH4", 
      bus0 = "fossilCH4",
      bus1 = "CH4",     
      p_nom_extendable = True,
      capital_cost = 0. ,
      marginal_cost = fuel_cost_fossilCH4,
      #p_nom_max = 0,
      efficiency =  1.
      )

# Create bioCH4 as a bus, and a store

network.add("Bus", "bioCH4", carrier = "bioCH4")

# Add a depletable store of biomethane, in MWh_th
network.add("Store",
      "bioCH4_store",
      bus = "bioCH4",
      e_nom_extendable = False,
      e_cyclic = False,
      e_nom = 10e6, #MWh_th 
      e_initial = 10e6, #MWh_th
    )

network.add("Link",
      "from_bioCH4", 
      bus0 = "bioCH4",
      bus1 = "CH4",     
      p_nom_extendable = with_bioCH4,
      capital_cost = 0. ,
      marginal_cost = fuel_cost_bioCH4,
      #p_nom_max = 0,
      efficiency =  1.
      )

network.add("Bus", "gas", carrier = "gas")


network.add("Link",
      "from_CH4_to_gas", 
      bus0 = "CH4",
      bus1 = "gas",     
      p_nom_extendable = True,
      capital_cost = 0. ,
      efficiency =  1.
      )

# Add OCGT, open cycle gas turbine

# NB in the following link the origin bus0 is in MWh_th, 
# hence the capital cost of the generator has to be scaled to the thermal input
# similarly for the marginal cost that was defined for the electric input
# the fuel cost at its calorific value is accounted in the link between the store and the general CH4 bus in €/MWh_th

this_id = "OCGT"
efficiency_OCGT = tech_parameters.at[dict_pypsa_par_id[this_id], 'efficiency'] 
capital_cost_OCGT_th = tech_parameters.at[dict_pypsa_par_id[this_id], 'capital_cost'] * efficiency_OCGT # in €/MW_th
marginal_cost_OCGT_th = tech_parameters.at[dict_pypsa_par_id[this_id], 'marginal_cost_no_fuel'] * efficiency_OCGT # in €/MW_th
network.add("Link",
      this_id, 
      bus0 = "gas",
      bus1 = "AC",     
      p_nom_extendable = with_OCGT,
      capital_cost = capital_cost_OCGT_th ,
      marginal_cost = marginal_cost_OCGT_th,
      efficiency =  efficiency_OCGT
      )

# Add CCGT, combined cycle gas turbine

# NB in the following link the origin bus0 is in MWh_th, 
# hence the capital cost of the generator has to be scaled to the thermal input
# similarly for the marginal cost that was defined for the electric input
# the fuel cost at its calorific value is accounted in the link between the store and the general CH4 bus in €/MWh_th

this_id ="CCGT"
efficiency_CCGT = tech_parameters.at[dict_pypsa_par_id[this_id], 'efficiency'] 
capital_cost_CCGT_th = tech_parameters.at[dict_pypsa_par_id[this_id], 'capital_cost'] * efficiency_CCGT # in €/MW_th
marginal_cost_CCGT_th = tech_parameters.at[dict_pypsa_par_id[this_id], 'marginal_cost_no_fuel'] * efficiency_CCGT # in €/MW_th

network.add("Link",
      this_id, 
      bus0 = "gas",
      bus1 = "AC",     
      p_nom_extendable = with_CCGT,
      capital_cost = capital_cost_CCGT_th,
      marginal_cost = marginal_cost_CCGT_th,
      efficiency =  efficiency_CCGT
      )          
          
    

# CO2 constraint

network.add("GlobalConstraint",
            "co2_limit",
            type="primary_energy",
            carrier_attribute="co2_emissions",
            sense="<=",
            constant=co2_limit)


# Add Storage technologies as stores and links

# add carriers, used by storage tec
network.add("Carrier", "H2")
network.add("Carrier", "water")
network.add("Carrier", "compH2") # Compressed H2, for storage
network.add("Carrier", "SynCH4")
network.add("Carrier", "CO2") # CO2 either biogenic or from direct air capure
network.add("Carrier", "CO2_stored") # CO2 stored



# Add Short duration (SD) storage, for example Li-ion battery
network.add("Bus", "SD_bus", carrier = "AC")
this_id_store = "SD_store"
this_id_in    = "SD_in"
SD_power_capital_cost =  tech_parameters.at[dict_pypsa_par_id[this_id_in], 'capital_cost'] 
SD_energy_capital_cost =  tech_parameters.at[dict_pypsa_par_id[this_id_store], 'capital_cost'] 
SD_VOM = tech_parameters.at[dict_pypsa_par_id[this_id_in], 'VOM'] 
SD_efficiency = tech_parameters.at[dict_pypsa_par_id[this_id_in], 'efficiency']
SD_half_cycle_efficiency = np.sqrt(SD_efficiency)

network.add("Store",
      this_id_store,
      bus = "SD_bus",
      e_nom_extendable = with_SD,
     # e_nom_max = 100e5,
      e_cyclic = True,
      capital_cost = SD_energy_capital_cost
    )


network.add("Link",
      "SD_out", 
      bus0 = "SD_bus",
      bus1 = "AC",
      p_nom_extendable = with_SD,
      efficiency = SD_half_cycle_efficiency,
      marginal_cost = SD_VOM * SD_half_cycle_efficiency # has to be scaled to the input
    )

network.add("Link",
      this_id_in,
      bus0 = "AC",
      bus1 = "SD_bus",
      p_nom_extendable = with_SD,
     # p_nom_max = 80e3,
      efficiency = SD_half_cycle_efficiency,
      capital_cost = SD_power_capital_cost
    )


# Add Long duration (LD1) storage, as in Guerra et al 2021 where LD1 is modeled as CAES
network.add("Bus", "LD1_bus", carrier = "AC")
this_id_store ="LD1_store"
this_id_out   = "LD1_out"
this_id_in    = "LD1_in"

LD1_power_discharging_capital_cost =  tech_parameters.at[dict_pypsa_par_id[this_id_out], 'capital_cost'] 
LD1_power_charging_capital_cost =  tech_parameters.at[dict_pypsa_par_id[this_id_in], 'capital_cost'] 
LD1_energy_capital_cost =  tech_parameters.at[dict_pypsa_par_id[this_id_store], 'capital_cost'] 
LD1_VOM = tech_parameters.at[dict_pypsa_par_id[this_id_out], 'VOM'] 
LD1_efficiency = tech_parameters.at['long_duration_1_storage_power', 'efficiency']
LD1_half_cycle_efficiency = np.sqrt(LD1_efficiency)

network.add("Store",
      this_id_store,
      bus = "LD1_bus",
      e_nom_extendable = with_LD1,
      e_nom_max = 100e5,
      e_cyclic = True,
      capital_cost = LD1_energy_capital_cost
    )


network.add("Link",
      this_id_out, 
      bus0 = "LD1_bus",
      bus1 = "AC",
      p_nom_extendable = with_LD1,
      efficiency = LD1_half_cycle_efficiency,
      capital_cost = LD1_power_discharging_capital_cost * LD1_half_cycle_efficiency,
      marginal_cost = LD1_VOM * LD1_half_cycle_efficiency # has to be scaled to the input
    )

network.add("Link",
      this_id_in,
      bus0 = "AC",
      bus1 = "LD1_bus",
      p_nom_extendable = with_LD1,
      p_nom_max = 80e3,
      efficiency = LD1_half_cycle_efficiency,
      capital_cost = LD1_power_charging_capital_cost 
    )



# Add Pumped Hydro storage modeled as LD2 in Guerra et al 2021
network.add("Bus", "water_PHS", carrier = "water_PHS")
this_id_store = "PHS_store"
this_id_in    = "PHS_in"

PHS_power_capital_cost =  tech_parameters.at[dict_pypsa_par_id[this_id_in], 'capital_cost'] 
PHS_energy_capital_cost =  tech_parameters.at[dict_pypsa_par_id[this_id_store], 'capital_cost'] 
PHS_efficiency = tech_parameters.at[dict_pypsa_par_id[this_id_in], 'efficiency']
PHS_VOM = tech_parameters.at[dict_pypsa_par_id[this_id_in], 'VOM'] 

hydro_turbine_efficiency = np.sqrt(PHS_efficiency)

network.add("Store",
      this_id_store,
      bus = "water_PHS",
      e_nom_extendable = with_PHS,
      e_nom_max = 100e5,
      e_cyclic = True,
      capital_cost = PHS_energy_capital_cost
    )


network.add("Link",
      "PHS_out", 
      bus0 = "water_PHS",
      bus1 = "AC",
      p_nom_extendable = with_PHS,
    #  p_nom_max = 80e3 / hydro_turbine_efficiency,
      efficiency = hydro_turbine_efficiency,
      marginal_cost = PHS_VOM * hydro_turbine_efficiency
    )

network.add("Link",
      this_id_in,
      bus0 = "AC",
      bus1 = "water_PHS",
      p_nom_extendable = with_PHS,
    #  p_nom_max = 80e3,
      efficiency = hydro_turbine_efficiency,
      capital_cost = PHS_power_capital_cost
    )



#Create  H2 bus, store, and links

#network.add("Bus", "H2_store", carrier = "H2") # Uncompressed H2
network.add("Bus", "compH2_store", carrier = "H2")  #Compressed H2

#network.add("Bus", "SynCH4", carrier = "SynCH4")
network.add("Bus", "SynCH4_store", carrier = "SynCH4")




this_id_in = "H2_in"

#Add the link "H2_in", i.e. Electrolysis that transport energy from the electricity bus (bus0) to the H2 bus (bus1)
network.add("Link",
      this_id_in, 
      bus0 = "AC",
      bus1 = "H2",     
      p_nom_extendable = True,
      efficiency = tech_parameters.at[dict_pypsa_par_id[this_id_in], 'efficiency'],
      capital_cost = tech_parameters.at[dict_pypsa_par_id[this_id_in], 'capital_cost']
    )



# Storage of H2 requires intermediate step of compression
# Compressed H2  three  componentes

this_id_in = "compH2_in"
# primary bus of the compressor is electricity
# Hence the efficiency and the fixed cost must be scaled
# comp_el_use =  "compression-electricity-input"  is in MWh_el/MWh_H2
comp_fixed = tech_parameters.at[dict_pypsa_par_id[this_id_in], 'capital_cost'] 
comp_el_use = tech_parameters.at[dict_pypsa_par_id[this_id_in], 'compression-electricity-input']

network.add("Link",
      this_id_in, 
      bus0 = "H2",
      bus1 = "compH2_store",
      bus2 = "AC",
      p_nom_extendable = True,
      efficiency = 1,
      efficiency2 = - comp_el_use,
      capital_cost = comp_fixed,
    )

this_id_store = "compH2_store"

#Connect the store to the bus
network.add("Store",
      this_id_store,
      bus = "compH2_store",
      e_nom_extendable = True,
      e_cyclic = True,
      capital_cost = tech_parameters.at[dict_pypsa_par_id[this_id_store], 'capital_cost'],
      #capital_cost = 0., # assume underground existing storage
      marginal_cost = tech_parameters.at[dict_pypsa_par_id[this_id_store], 'marginal_cost_no_fuel']
      )


network.add("Link",
      "compH2_out", 
      bus0 = "compH2_store",
      bus1 = "H2",     
      p_nom_extendable = True,
      efficiency = 1.,
      capital_cost = 0.
    )

#Add the generator link "FC" that transports energy from the H2 bus (bus0) to the electricity bus (bus1) via Fuel Cell
# NB, the fuel cell capital cost must be scaled to the primary H2 input
this_id = "FC"
fuel_cell_efficiency = tech_parameters.at[dict_pypsa_par_id[this_id], 'efficiency']
fuel_cell_capital_cost = tech_parameters.at[dict_pypsa_par_id[this_id], 'capital_cost'] * fuel_cell_efficiency

network.add("Link",
      this_id, 
      bus0 = "H2",
      bus1 = "AC",     
      p_nom_extendable = with_FC,
      efficiency = fuel_cell_efficiency,
      capital_cost = fuel_cell_capital_cost
    ) 

# Add a link from the H2 to the general gas bus, this allows the OCGT and CCGT to be multifuel
network.add("Link",
      'from_H2_to_gas', 
      bus0 = "H2",
      bus1 = "gas",     
      p_nom_extendable = with_mft,
      efficiency = 1,
    )

# Methanization process

# Add two buses and a depletable store of CO2

network.add("Bus", "exogenous_CO2", carrier = "CO2")
network.add("Bus", "CO2", carrier = "CO2")



network.add("Store",
      "exogenous_CO2_store",
      bus = "exogenous_CO2",
      e_nom_extendable = False,
      e_cyclic = False,
      e_nom = 30e6, # ton
      e_initial = 30e6, # ton
    )

# Add a link from the CO2 store to the bus where the CO2 can be used
network.add("Link",
      "from_CO2", 
      bus0 = "exogenous_CO2",
      bus1 = "CO2",     
      p_nom_extendable = with_ex_CO2,
      capital_cost = 0.,
      efficiency =  1.,
      marginal_cost = co2_feedstock_cost,
      )
# Add a national CH4 store, in this model this is useful only if there is some endogenous elCH4 production
# The max energy stored is equal to the current storage capacity of fossil gas in MWh_th

network.add("Store",
      "SynCH4_store",
      bus = "SynCH4_store",
      e_nom_extendable = True,
      e_cyclic = True,
      e_nom_max = 190e6, #MWh_th
      #e_initial = 0. #MWh_th
      capital_cost = 0.00000001 #Negligible capital cost to avoid selecting the max capacity as optimal capacity
    )

## Add a link that uses H2 and CO2 to produce CH4
# NB, the methanation capital cost must be scaled to the primary H2 input
this_id = "SynCH4_in"

SynCH4_efficiency = tech_parameters.at[dict_pypsa_par_id[this_id], 'efficiency']
SynCH4_capital_cost = tech_parameters.at[dict_pypsa_par_id[this_id], 'capital_cost'] * SynCH4_efficiency
SynCH4_min_relative_power = tech_parameters.at[dict_pypsa_par_id[this_id], 'min_part_load'] / 100.

#cost_CO2_for_CH4 = tech_parameters.at["methanation", 'CO2'] * SynCH4_efficiency * co2_feedstock_cost
# input CO2 from bus usedCO2 is in ton
eta_eq_tCO2_per_MWh_H2 = tech_parameters.at["methanation", 'CO2'] * SynCH4_efficiency 

network.add(
    "Link",
    this_id,
    bus0="H2",
    bus1="SynCH4_store",
    bus2="CO2",
    efficiency=SynCH4_efficiency,
    efficiency2= - eta_eq_tCO2_per_MWh_H2,
    p_min_pu = SynCH4_min_relative_power,
    p_nom_extendable = with_SynCH4,
    capital_cost = SynCH4_capital_cost,
 #   marginal_cost = cost_CO2_for_CH4,
    )

network.add(
    "Link",
    "SynCH4_out",
    bus0="SynCH4_store",
    bus1="CH4",
    efficiency=1.0,
    p_nom_extendable = with_SynCH4,
    marginal_cost = 0.  #SynCH4 costs are accounted elsewhere
    )



# In[25]:


# Add endogenous production of CO2, such as DAC ....

# DAC requires heat, and storage of CO2

network.add("Bus", "heat", carrier="heat")



this_id = "heat_pump_mid"
efficiency_HP120 = tech_parameters.at[dict_pypsa_par_id[this_id], 'efficiency']  # This COP/100
capital_cost_HP120_e = tech_parameters.at[dict_pypsa_par_id[this_id], 'capital_cost'] * efficiency_HP120 # in €/MW_e
marginal_cost_HP120_e = tech_parameters.at[dict_pypsa_par_id[this_id], 'marginal_cost_no_fuel'] * efficiency_HP120 # in €/MW_e

network.add("Link",
             this_id,
             bus0="AC",
             bus1="heat",
             carrier="heat_pump",
             p_nom_extendable=True,
             capital_cost=capital_cost_HP120_e,
             marginal_cost = marginal_cost_HP120_e,
             efficiency=efficiency_HP120,
           )

this_id = "dac"
efficiency_dac_e  = 1. / tech_parameters.at[dict_pypsa_par_id[this_id], 'electricity-input']  # tCO2/MWh_e
efficiency_dac_th = tech_parameters.at[dict_pypsa_par_id[this_id], 'heat-input']  # MWh_th/tCO2
efficiency_th_e   = efficiency_dac_e * efficiency_dac_th  # MWh_th/MWh_e
capital_cost_dac  = tech_parameters.at[dict_pypsa_par_id[this_id], 'capital_cost'] * efficiency_dac_e # in €/MW_e

network.add("Link",
             this_id,
             bus0="AC",
             bus1="CO2",
             bus2="heat",
             carrier=this_id,
             p_nom_extendable = with_DAC,
             capital_cost = capital_cost_dac,
             efficiency   = efficiency_dac_e,
             efficiency2  = - efficiency_th_e,
           )

# Add storage of CO2 as two Links and one Store (it's redundant, but it may be useful for adding marginal costs)
network.add("Bus", "CO2_stored", carrier = "CO2")

this_id = "CO2_storage_store"
capital_cost_CO2_store = tech_parameters.at[dict_pypsa_par_id[this_id], 'capital_cost']
network.add("Store",
            this_id,
            bus="CO2_stored",
            carrier="CO2_stored",
            e_nom_extendable=True,
            e_cyclic=True,
            capital_cost=capital_cost_CO2_store,
           )

network.add("Link",
      "CO2_storage_in", 
      bus0 = "CO2",
      bus1 = "CO2_stored",     
      p_nom_extendable = True,
      efficiency = 1.,
      capital_cost = 0.
    )

network.add("Link",
      "CO2_storage_out", 
      bus0 = "CO2_stored",
      bus1 = "CO2",     
      p_nom_extendable = True,
      efficiency = 1.,
      capital_cost = 0.
    )



# In[26]:


# Refine the model with storage and nuclear CF specific constraints and call the solver to optimize the model

model = network.optimize.create_model(snapshots=network.snapshots)

# Add constraint that power disharge = power charge for PHS and SD
capacity = model.variables["Link-p_nom"]

if with_SD:
    eff_SD = network.links.at["SD_out", "efficiency"]
    lhs_SD = capacity.loc["SD_in"] - eff_SD * capacity.loc["SD_out"]
    model.add_constraints(lhs_SD == 0, name="Link-SD_fix_ratio")

if with_PHS:
    eff_PHS = network.links.at["PHS_out", "efficiency"]
    lhs_PHS = capacity.loc["PHS_in"] - eff_PHS * capacity.loc["PHS_out"]
    model.add_constraints(lhs_PHS == 0, name="Link-PHS_fix_ratio")

#Add a constraint on the yearly average capacity factor if nuclear is modeled as a fleet of reactors
if with_nuc_fleet and with_nuclear:
    P_nuc = model.variables["Generator-p_nom"].loc["nuclear"]
    E_nuc = model.variables["Generator-p"].loc[:, "nuclear"].sum()
    lhs_nuc = E_nuc - P_nuc * Num_hours * max_capacity_factor_nuclear
    model.add_constraints(lhs_nuc <= 0, name="Nuc_CF_fleet")

# Add a constraint for the maximum capacity factor of hydrolisis
if with_demand_H2_flag:
        P_H2 = model.variables['Link-p_nom'].loc['H2_in']
        E_H2 = model.variables['Link-p'].loc[:,'H2_in'].sum()
        lhs_H2 = E_H2 - P_H2 * Num_hours * max_capacity_factor_H2 
        model.add_constraints(lhs_H2 <= 0, name="H2_CF_Max")
    
# Solve
solver_options={'NumericFocus':3}
#network.optimize.solve_model(snapshots=network.snapshots, pyomo=False, solver_name=name_solver, solver_options={'NumericFocus':2}, keep_references = True, keep_shadowprices = True)

network.optimize.solve_model(snapshots=network.snapshots, pyomo=False, solver_name=name_solver, solver_options=solver_options,  keep_references = True, keep_shadowprices = True)



#model.variables["Generator-p_nom"]
#model.variables["Generator-p"].loc[:, "nuclear"].sum()
#model.variables['Link-p'].loc[:,'H2_in']
#model.variables['Link-p_nom'].loc['H2_in']

#network.loads_t.p['load_H2'].sum()/1.e6.loc


def colors_distinguishable_by_color_blind_person(num_of_categories):
    my_min = 0.
    my_max = 1.
    if num_of_categories < 4:
        # then restrict interval to avoid too much contrast
        my_min = 0.15
        my_max = 0.85
        
    category_colors = plt.colormaps['inferno'](np.linspace(my_min, my_max, num_of_categories))
    
    return category_colors[::-1]

def my_stacked_horizontal_bar_plot(results, category_labels, my_title=None, my_x_label=None):
    """
    Parameters
    ----------
    results : dict
        A key entry is a label that describe the results 
        The list of values associated to a key are the results specifics to a category.
        The category labels are in the input category_labels
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_labels*.
    category_labels : list of str
        The category labels.
    """
    result_labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    if my_title is None:
        my_title = ""
    if my_x_label is None:
        my_x_label = ""
    # first check how may categories are active with a tollerance epsilon
    #epsilon = 0.00001
    #cat_displayed = 0
    #map_cat_color = np.zeros(len(category_labels))
    #for i, colname in enumerate(category_labels):
        #if data[:, i].sum() > epsilon:
            #map_cat_color[i] = cat_displayed
            #cat_displayed +=1
    
    # Otherwise, to offer always the same colors under different subsets of techs make this active
    cat_displayed = len(category_labels)
    map_cat_color = np.arange(cat_displayed)
    
    category_colors = colors_distinguishable_by_color_blind_person(cat_displayed)
    fig, ax = plt.subplots(figsize=(12, 3.5))
    ax.invert_yaxis()
    #ax.xaxis.set_visible(False)
    ax.set_xlabel(my_x_label)
    max_value =  np.sum(data, axis=1).max()
    ax.set_xlim(0,max_value)
    pc_min_value_text_displayed = 1.5
    min_value_displayed = max_value * pc_min_value_text_displayed / 100.
    for i, colname in enumerate(category_labels):
        widths = data[:, i]
        # Avoid displaying zero categories, approx to a very small values
        if widths.sum() > 0.00001:
            color = category_colors[int(map_cat_color[i])]
            r, g, b, _ = color
            text_color = 'white' if r * g * b < 0.055 else 'black'
            
            used_col_name = False
            for j in range(len(result_labels)):
                if widths[j] > 0.:
                    start = data_cum[j, i] - widths[j]
                    if not used_col_name:
                        # this avoid duplicate names in the category legend
                        used_col_name = True
                        rect = ax.barh(result_labels[j], widths[j], left=start, height=0.8,
                                    color=color, label=colname)
                    else:
                        rect = ax.barh(result_labels[j], widths[j], left=start, height=0.8,
                                    color=color)
                    if widths[j] > min_value_displayed:
                        ax.bar_label(rect, label_type='center', color=text_color, fmt="%.1f")
            
   
    ax.legend(ncol=8, bbox_to_anchor=(0, 1), loc='lower left', fontsize='small')
    #ax.legend(bbox_to_anchor=(0, 1), loc='lower left', fontsize='small')
    
    fig.suptitle(t=my_title, fontweight='bold', y=1.1)
    plt.tight_layout()


    return fig




def basic_network_statistics(network):
    """
    From a pypsa network return some basic statistics in a dataframe
    Assumptions:
    - the network must contain data in EUR,  MW, and MWh for cost, power, and energy, respectively
    
    """
    Tot_elec_delivered = network.loads_t.p['load'].sum()/1.e6 # TWh
    LCOEwS = network.loads_t.p['load'].multiply(network.snapshot_weightings.generators, axis=0).multiply(network.buses_t.marginal_price["AC"], axis=0).sum()/1.e+6 / Tot_elec_delivered # €/MWh

    Tot_H2_delivered = 0.
    LCOH = 0.
    if with_demand_H2_flag:
        Tot_H2_delivered = network.loads_t.p['load_H2'].sum()/1.e6 # TWh
        LCOH =  network.loads_t.p['load_H2'].multiply(network.snapshot_weightings.generators, axis=0).multiply(network.buses_t.marginal_price["H2"], axis=0).sum()/1.e+6 / Tot_H2_delivered # €/MWh

    Tot_en_delivered = Tot_elec_delivered +   Tot_H2_delivered
    Ann_system_cost = network.objective/1.e9 #in 10^9 €
    LCOS = Ann_system_cost / Tot_en_delivered * 1.e+3 # €/MWh
    mu_CO2 = -float(network.model.dual["GlobalConstraint-co2_limit"])
    
    dict_system_indices = {
                        "Tot_elec_delivered" : ["TWh", float(Tot_elec_delivered)],
                        "Tot_H2_delivered"   : ["TWh", float(Tot_H2_delivered)],
                        "Tot_en_delivered"   : ["TWh", float(Tot_en_delivered)],
                        "Ann_system_cost"  : ["G€", float(Ann_system_cost)],
                        "LCOEwS"           : ["€/MWh", float(LCOEwS)],        
                        "LCOH"             : ["€/MWh", float(LCOH)],
                        "LCOS"             : ["€/MWh", float(LCOS)],
                        "mu_CO2"           : ["€/tCO2", mu_CO2],
                          }

    df_sys_indices = pd.DataFrame(dict_system_indices, index=["unit", "value"])
    return df_sys_indices
    
    
    


def advanced_network_statistics_V2(network, 
                                   non_fossil_gen_labels=None, non_fossil_gen_CF=None, 
                                   fossil_gen_labels=None, links_gen_labels=None, 
                                   primary_fuel_label=None,
                                  indirect_fuel_label = None,
                                   df_sys_indices=None,
                                   link_to_general_gas_labels = None,
                                   bus_general_gas_label = None,
                                  ):
    """
    From a pyspsa network return some basic and advanced statistics on the primary generators
    The generators are differentiated between fossil or non fossil
    Some generators may be modeled as links that derive their input from a store
    
    The advanced statistics require the following input about the underling model:
    - non_fossil_gen_labels: a list of str, the labels of the non fossil generators, if any
    - non_fossil_gen_CF: a list of arrays with the non fossil generators' CF, if any
    - 
    """
    
    
    # derived indices
    tot_non_fossil_output       = 0. # TWh
    tot_pot_non_fossil_output   = 0. # TWh
    tot_fossil_output           = 0. # TWh
    primary_fossil_input        = 0. # TWh
    tot_fossil_CO2              = 0. # Mt
    CO2_intensity_g             = 0.  # g/kWh
    
    # data structures for all generators
    tot_generator_out_t = np.zeros(len(network.snapshots))

    all_gen_plus_link_label = []
    all_gen_label = []
    num_gen_non_fossil = 0
    if non_fossil_gen_labels is not None:
        num_gen_non_fossil = len(non_fossil_gen_labels)
        all_gen_label = non_fossil_gen_labels[:]
    
    num_gen_fossil = 0
    if fossil_gen_labels is not None:
        num_gen_fossil = len(fossil_gen_labels)
        all_gen_label += fossil_gen_labels[:]
    
    num_gen_links = 0
    if links_gen_labels is not None:
        num_gen_links = len(links_gen_labels)
        all_gen_plus_link_label = all_gen_label[:] + links_gen_labels[:]
    
    num_all_gen =  num_gen_non_fossil + num_gen_fossil
    num_all_gen_plus_link = len(all_gen_plus_link_label)
    generator_output                  = np.zeros(num_all_gen_plus_link) # This is only the el output from primary sources not fuel based
    generator_any_output              = np.zeros(num_all_gen_plus_link) # This is the el output not matter the sources
    generator_installed_capacity      = np.zeros(num_all_gen_plus_link)
    generator_capex                   = np.zeros(num_all_gen_plus_link)
    generator_marginal_cost           = np.zeros(num_all_gen_plus_link)
    generator_marg_cost_inc_prim_fuel = np.zeros(num_all_gen_plus_link)
    generator_el_out_primary          = np.zeros(num_all_gen_plus_link)
    generator_capex_opex_primary_fuel = np.zeros(num_all_gen_plus_link)
    generator_LCOE                    = np.zeros(num_all_gen_plus_link)
    # This is similar as network.statistics.capacity_factor() but with the provision of accounting only for primary sources
    generator_primary_CF              = np.zeros(num_all_gen_plus_link)
    generator_output_shares           = np.zeros(num_all_gen_plus_link)
    generator_revenue                 = np.zeros(num_all_gen_plus_link)
    generator_input_cost              = np.zeros(num_all_gen_plus_link)
    investment                        = np.zeros(num_all_gen_plus_link)
    generator_any_CF                  = np.zeros(num_all_gen_plus_link)
    
    all_fuel_labels = []
    
    num_primary_fuel_label = 0
    if primary_fuel_label is not None:
        num_primary_fuel_label = len(primary_fuel_label)
        all_fuel_labels = primary_fuel_label[:]
    
    primary_fuel_costs        = np.zeros(num_primary_fuel_label)
    primary_fuel_out_shares   = np.zeros(num_primary_fuel_label)
    primary_fuel_electric_out = np.zeros(num_primary_fuel_label)
    
    num_indirect_fuel_label = 0
    if indirect_fuel_label is not None:
        num_indirect_fuel_label = len(indirect_fuel_label)
        all_fuel_labels += indirect_fuel_label[:]
        
    print("all_fuel_labels", all_fuel_labels)
    print("primary_fuel_label", primary_fuel_label)
    bus_dict = {} # create an empty dictionary
    bus_dict_price = {}
    bus_dict_fossil_share = {}
    bus_dict_primary_share = {}
    bus_primary_fuel_dict = {}
    bus_primary_fuel_shares_dict = {}
    
    # Process fuels to determine primary shares for generator links
    if len(all_fuel_labels) > 0:
        for i, p_s in enumerate(all_fuel_labels):
            
            this_fuel_output = -network.links_t.p1[p_s].sum()/1e6 # TWh
            # avoid that a near minus zero value from fossil CH4 induces confusion
            if this_fuel_output < 0:
                this_fuel_output = 0.
            
            this_fuel_output_cost = this_fuel_output * 1.e+6 * network.links.marginal_cost[p_s] # This cost takes care only of the last link step, if it's a syn fuel its cost is accounteded elsewhere
            this_fuel_primary_output = 0
            if p_s in primary_fuel_label:
                this_fuel_primary_output = this_fuel_output
                primary_fuel_costs[i] = this_fuel_output_cost
            #print(p_s + " this_fuel_primary_output ", this_fuel_primary_output)   
            this_link_bus_out = network.links.bus1[p_s]
            this_link_bus_in = network.links.bus0[p_s]
            this_link_carrier_in =network.buses.loc[this_link_bus_in, "carrier"]
            this_carrier_CO2_intensity =  network.carriers.loc[this_link_carrier_in, "co2_emissions"]
            this_link_fossil_q = 0.
            if this_carrier_CO2_intensity > 0:
                # It's a fossil carrier 
                tot_fossil_CO2 += this_fuel_output * this_carrier_CO2_intensity
                primary_fossil_input += this_fuel_output
                this_link_fossil_q = this_fuel_output
                
            if this_link_bus_out in bus_dict.keys():
                # update balance
                bus_dict[this_link_bus_out][0] += this_fuel_output
                bus_dict[this_link_bus_out][1] += this_fuel_output_cost 
                bus_dict[this_link_bus_out][2] += this_link_fossil_q
                bus_dict[this_link_bus_out][3] += this_fuel_primary_output
                if p_s in primary_fuel_label:
                    bus_primary_fuel_dict[this_link_bus_out][i] += this_fuel_primary_output
            else:
                # add bus to the dictionary
                bus_dict[this_link_bus_out] = [this_fuel_output, 
                                               this_fuel_output_cost, 
                                               this_link_fossil_q, 
                                               this_fuel_primary_output]
                print(p_s, i, this_fuel_primary_output)
                bus_primary_fuel_dict[this_link_bus_out] = [0] * len(primary_fuel_label)
                if p_s in primary_fuel_label:
                    bus_primary_fuel_dict[this_link_bus_out][i] = this_fuel_primary_output
        
        print("bus_primary_fuel_dict", bus_primary_fuel_dict)   
        
        # Now process dict and determine equivalent energy price per bus
        for k in bus_dict.keys():
            if  bus_dict[k][0] > 0:
                bus_dict_price[k]         = bus_dict[k][1]/ bus_dict[k][0] /1.e+6 # €/;Wh
                bus_dict_fossil_share[k]  = bus_dict[k][2]/ bus_dict[k][0]
                bus_dict_primary_share[k] = bus_dict[k][3]/ bus_dict[k][0]
                bus_primary_fuel_shares_dict[k] = [0] * len(primary_fuel_label)
                for i, p_s in enumerate(primary_fuel_label):
                    bus_primary_fuel_shares_dict[k][i] = bus_primary_fuel_dict[k][i]/ bus_dict[k][0]
            else:
                bus_dict_price[k]         = 0
                bus_dict_fossil_share[k]  = 0
                bus_dict_primary_share[k] = 0
                bus_primary_fuel_shares_dict[k] = [0] * len(primary_fuel_label)
        
   # Now process secondary links to final bus gas
    if bus_general_gas_label is not None and link_to_general_gas_labels is not None:
        tot_e_to_gas = 0
        bus_dict_to_gas = {}
        for k in link_to_general_gas_labels:
            this_link_bus_in = network.links.bus0[k] 
            this_link_e_out = - network.links_t.p1[k].sum()/1e6 # TWh
            tot_e_to_gas += this_link_e_out
            if this_link_bus_in in bus_dict_to_gas.keys():
                bus_dict_to_gas[this_link_bus_in] += this_link_e_out
            else:
                bus_dict_to_gas[this_link_bus_in] = this_link_e_out
        # we are interested only in the relative shares
        print("tot_e_to_gas", tot_e_to_gas)
        for k in bus_dict_to_gas.keys():
            bus_dict_to_gas[k] = bus_dict_to_gas[k]/tot_e_to_gas
        
        bus_dict_price[bus_general_gas_label]         = 0
        bus_dict_primary_share[bus_general_gas_label] = 0
        bus_dict_fossil_share[bus_general_gas_label]  = 0
        bus_primary_fuel_shares_dict[bus_general_gas_label] = np.zeros(len(primary_fuel_label))
        for k in bus_dict_to_gas.keys():
            bus_dict_price[bus_general_gas_label]         += bus_dict_price[k] * bus_dict_to_gas[k]
            bus_dict_primary_share[bus_general_gas_label] += bus_dict_primary_share[k] * bus_dict_to_gas[k]
            bus_dict_fossil_share[bus_general_gas_label]  += bus_dict_fossil_share[k] * bus_dict_to_gas[k]
            bus_primary_fuel_shares_dict[bus_general_gas_label] += np.multiply(np.array(bus_primary_fuel_shares_dict[k][:]), bus_dict_to_gas[k])
        
        
    print("bus_dict (tot_output, cost, fossil_out, primary_out)", bus_dict)
    print("bus_dict_price", bus_dict_price)
    print("bus_dict_fossil_share", bus_dict_fossil_share)
    print("bus_dict_primary_share", bus_dict_primary_share)
    print("bus_primary_fuel_shares_dict", bus_primary_fuel_shares_dict)
    
        
    if non_fossil_gen_labels is not None:
        for i, r_g in enumerate(non_fossil_gen_labels):
            tot_non_fossil_output += float(network.generators_t.p[r_g].sum())/1e6
    
            if non_fossil_gen_CF is not None:
                tot_pot_non_fossil_output += float((network.generators.p_nom_opt[r_g] * non_fossil_gen_CF[i]).sum()) /1.e+6
    
    if fossil_gen_labels is not None:
        for i, f_g in enumerate(fossil_gen_labels):
            this_gen_output = float(network.generators_t.p[f_g].sum()/1e6)
            tot_fossil_output += this_gen_output
            this_gen_primary_fossil = this_gen_output / network.generators.loc[f_g, "efficiency"]
            primary_fossil_input += this_gen_primary_fossil
            this_gen_carrier = network.generators.loc[f_g, "carrier"]
            this_gen_carrier_CO2_factor = network.carriers.loc[this_gen_carrier, "co2_emissions"]
            tot_fossil_CO2 += this_gen_primary_fossil * this_gen_carrier_CO2_factor
    
# Process generators, fossil and fossil-free, and generation links

    for i, my_l in enumerate(all_gen_plus_link_label):
        
        if i < num_all_gen:
            # Then data are stored in the generators component
            this_gen_el_output = network.generators_t.p[my_l].sum()/1e6 # TWh
            generator_output[i] = this_gen_el_output
            generator_any_output[i] = this_gen_el_output
            generator_el_out_primary[i] = this_gen_el_output
            tot_generator_out_t += network.generators_t.p[my_l]
            generator_installed_capacity[i] = network.generators.p_nom_opt[my_l]
            generator_capex[i] = generator_installed_capacity[i] * network.generators.loc[my_l, "capital_cost"] 
            generator_marginal_cost[i] = network.generators_t.p[my_l].multiply(network.snapshot_weightings.generators, axis=0).sum() * network.generators.marginal_cost[my_l]
            generator_marg_cost_inc_prim_fuel[i] = generator_marginal_cost[i]
            generator_revenue[i] = network.generators_t.p[my_l].multiply(network.snapshot_weightings.generators, axis=0).multiply(network.buses_t.marginal_price["AC"], axis=0).sum()/1.e+9 #G€
        else:
            # Then data are in the links component
            # This link whole output may not all be primary
            # Account as primary generation only the primary fraction of its origin bus carrier
            this_gen_el_output =  - network.links_t.p1[my_l].sum()/1e6 # TWh
            if this_gen_el_output > 0:
                generator_any_output[i] = this_gen_el_output
                this_link_eff = network.links.efficiency[my_l]
                generator_installed_capacity[i] = network.links.p_nom_opt[my_l] * this_link_eff
                generator_capex[i] = network.links.capital_cost[my_l] * network.links.p_nom_opt[my_l]
                generator_marginal_cost[i] = network.links.marginal_cost[my_l] * network.links_t.p0[my_l].sum() # NB this may be without fuel
                this_link_bus_in =  network.links.bus0[my_l]
                this_link_energy_in = this_gen_el_output / this_link_eff
                this_link_fuel_cost = this_link_energy_in * bus_dict_price[this_link_bus_in] * 1.e+6 # in € # This cost takes care only of the last link step, if it's a syn fuel its cost is accounteded elsewhere
                # if accounted as fuel do not attribute it at the generator level
                generator_marg_cost_inc_prim_fuel[i] += this_link_fuel_cost 
                tot_fossil_output += bus_dict_fossil_share[this_link_bus_in] * this_gen_el_output
                tot_generator_out_t += - network.links_t.p1[my_l] * bus_dict_primary_share[this_link_bus_in]
                # if accounted at the fuel level do not assign output at the generator level
                generator_el_out_primary[i] = this_gen_el_output *  bus_dict_primary_share[this_link_bus_in]
                for j in range(len(primary_fuel_label)):  
                    primary_fuel_electric_out[j] += this_gen_el_output * bus_primary_fuel_shares_dict[this_link_bus_in][j]
                    #print("Verify on bus, output, j, share", this_link_bus_in, this_gen_el_output, j, bus_primary_fuel_shares_dict[this_link_bus_in][j])
                #print("verify bus-link-gen (fuel cost and primary share) ", my_l, this_link_bus_in, this_link_fuel_cost, bus_dict_primary_share[this_link_bus_in])
                generator_revenue[i] = -network.links_t.p1[my_l].multiply(network.snapshot_weightings.generators, axis=0).multiply(network.buses_t.marginal_price["AC"], axis=0).sum()/1.e+9 # G€
                generator_input_cost[i] = network.links_t.p0[my_l].multiply(network.snapshot_weightings.generators, axis=0).multiply(network.buses_t.marginal_price[this_link_bus_in], axis=0).sum()/1.e+9 #G€
        print("my_l, cap, label, unit_investment", my_l, generator_installed_capacity[i], dict_pypsa_par_id[my_l],tech_parameters["investment"][dict_pypsa_par_id[my_l]])
        investment[i] = generator_installed_capacity[i] * tech_parameters["investment"][dict_pypsa_par_id[my_l]]/1.e+9 # G€
            
    # Bring everything to G€ and TWh
    generator_capex = generator_capex /1.e+9 # G€
    generator_marginal_cost = generator_marginal_cost /1.e+9 # G€
    generator_marg_cost_inc_prim_fuel = generator_marg_cost_inc_prim_fuel /1.e+9 # G€
    generator_capex_plus_opex = generator_capex + generator_marginal_cost #G€
    primary_fuel_costs = primary_fuel_costs /1.e+9 # G€
    Ann_system_cost_g = generator_capex_plus_opex.sum() + primary_fuel_costs.sum() # G€
    
    generator_capex_opex_primary_fuel = generator_capex + generator_marg_cost_inc_prim_fuel # G€
    

    Tot_el_generated = generator_output.sum() + primary_fuel_electric_out.sum() # TWh

    
    for i, my_l in enumerate(all_gen_plus_link_label):
        if generator_el_out_primary[i] > 0:
            generator_LCOE[i] = generator_capex_opex_primary_fuel[i] / generator_el_out_primary[i] * 1.e3 # €/MWh
            generator_primary_CF[i] = generator_el_out_primary[i] * 1.e+6 / generator_installed_capacity[i] / Num_hours 
            generator_any_CF[i] = generator_any_output[i] * 1.e+6 / generator_installed_capacity[i] / Num_hours 
        else:
            generator_LCOE[i] = 0.
            generator_primary_CF[i] = 0.
            if generator_installed_capacity[i] > 0. and generator_any_output[i] > 0. :
                generator_any_CF[i] = generator_any_output[i] * 1.e+6 / generator_installed_capacity[i] / Num_hours 
            else:
                generator_any_CF[i] = 0.
            
    
    
    # these shares could also be accessed with (NB as base 1000) 
    # network.statistics.supply()
    generator_output_shares =  generator_output/float(Tot_el_generated)
    primary_fuel_out_shares = primary_fuel_electric_out /float(Tot_el_generated) 
    print("This in output primary_fuel_out_shares", primary_fuel_out_shares)
    print("Comes from this el out TWh", primary_fuel_electric_out)

    LCOS_g = generator_capex_opex_primary_fuel.sum() / float(Tot_el_generated) * 1.e3 # €/MWh

    non_fossil_fraction = tot_non_fossil_output / Tot_el_generated
    Curtailment = tot_pot_non_fossil_output - tot_non_fossil_output # TWh
    
    curtailment_pc = 0.
    if tot_non_fossil_output > epsilon:
        Overgeneration_factor = tot_pot_non_fossil_output / tot_non_fossil_output
        curtailment_pc = (Overgeneration_factor - 1.) * 100
    else:
        Overgeneration_factor = 0.
    if tot_fossil_CO2 > 0:
        CO2_intensity_g = tot_fossil_CO2 / Tot_el_generated * 1.e+3 # gCO2/kWh
    

    # add new system indices
    df_sys_indices["Tot_el_generated"] = ["TWh", Tot_el_generated]
    df_sys_indices["tot_non_fossil_output"] = ["TWh", tot_non_fossil_output] 
    df_sys_indices["tot_fossil_output"] = ["TWh", tot_fossil_output] 
    df_sys_indices["tot_fossil_CO2"] = ["Mt", tot_fossil_CO2] 
    df_sys_indices["CO2_intensity_g"] = ["gCO2/kWh", CO2_intensity_g]
    df_sys_indices["tot_pot_non_fossil_output"] = ["TWh", tot_pot_non_fossil_output]
    df_sys_indices["Ann_system_cost_g"] = ["G€", Ann_system_cost_g]
    df_sys_indices["LCOS_g"] = ["€/MWh", LCOS_g]
    df_sys_indices["Curtailment"] = ["TWh", Curtailment]
    df_sys_indices["curtailment_pc"] = ["%", curtailment_pc]
    df_sys_indices["tot_generator_out_t"] = ["MW", tot_generator_out_t]
    
    # Create a dataframe with the results that may need further processing
    df_gen_results = pd.DataFrame(
                    {
                        "generator_installed_capacity" : pd.Series(generator_installed_capacity, index=all_gen_plus_link_label, copy=True),
                        "investment" : pd.Series(investment, index=all_gen_plus_link_label, copy=True),
                        "generator_output" : pd.Series(generator_output, index=all_gen_plus_link_label, copy=True),
                        "generator_any_output" : pd.Series(generator_any_output, index=all_gen_plus_link_label, copy=True),                   
                        "generator_capex" : pd.Series(generator_capex, index=all_gen_plus_link_label, copy=True),
                        "generator_marginal_cost" : pd.Series(generator_marginal_cost, index=all_gen_plus_link_label, copy=True),
                        "generator_capex_plus_opex" : pd.Series(generator_capex_plus_opex, index=all_gen_plus_link_label, copy=True),
                        "generator_marg_cost_inc_prim_fuel" : pd.Series(generator_marg_cost_inc_prim_fuel, index=all_gen_plus_link_label, copy=True),
                        "generator_el_out_primary" : pd.Series(generator_el_out_primary, index=all_gen_plus_link_label, copy=True),
                        "generator_capex_opex_primary_fuel" : pd.Series(generator_capex_opex_primary_fuel, index=all_gen_plus_link_label, copy=True),
                        "generator_output_shares" : pd.Series(generator_output_shares, index=all_gen_plus_link_label, copy=True),
                        "generator_LCOE" : pd.Series(generator_LCOE, index=all_gen_plus_link_label, copy=True),
                        "generator_primary_CF" : pd.Series(generator_primary_CF, index=all_gen_plus_link_label, copy=True),
                        "generator_any_CF" : pd.Series(generator_any_CF, index=all_gen_plus_link_label, copy=True),
                        "generator_revenue" : pd.Series(generator_revenue, index=all_gen_plus_link_label, copy=True),
                        "generator_input_cost" : pd.Series(generator_input_cost, index=all_gen_plus_link_label, copy=True),
                
                    }
                )
    df_gen_results = df_gen_results.assign(price = df_gen_results.generator_revenue / df_gen_results.generator_any_output * 1.e+3) # €/MWh
    df_gen_results = df_gen_results.assign(generator_market_value = (df_gen_results.generator_capex_plus_opex +df_gen_results.generator_input_cost)/ df_gen_results.generator_any_output * 1.e+3) # €/MWh
    df_gen_results = df_gen_results.assign(capex_unit = df_gen_results.generator_capex / df_gen_results.generator_any_output * 1.e+3) # €/MWh
    df_gen_results = df_gen_results.assign(opex_unit = df_gen_results.generator_marginal_cost / df_gen_results.generator_any_output * 1.e+3) # €/MWh
    df_gen_results = df_gen_results.assign(input_cost_unit = df_gen_results.generator_input_cost / df_gen_results.generator_any_output * 1.e+3) # €/MWh

    # dataframe for the primary fuels
    df_fuel_primary_results = pd.DataFrame({
                            "primary_fuel_costs" : pd.Series(primary_fuel_costs, index=primary_fuel_label, copy=True),
                            "primary_fuel_out_shares" : pd.Series(primary_fuel_out_shares, index=primary_fuel_label, copy=True),
       
    })
    
    #Dataframe for the bus statistics
    df_bus_results = pd.DataFrame({
                    "primary_share" : pd.Series(bus_dict_primary_share, index=bus_dict_primary_share.keys())
    })
    
    return df_bus_results, df_fuel_primary_results, df_gen_results, df_sys_indices

   
    
    
    
def storage_el_balance(network, 
                       all_storage_tech_list = None,
                       list_syn_fuels = None,
                       df_sys_indices = None,
                      storage_single_link_list = None,
                       storage_with_bus2 = None,
                      ):
    """"
    Parameters:
     - all_storage_tech_list: list of strings of the identifier of storage techs connected 
         both as input and output via links to a store
         the convention is to provide the storage id only, and the to derive as id_store, id_in, and id_out
         the three components of the storage tech in the PyPSA network
         Es. the id is "PHS" and then we have "PHS_store", "PHS_in", "PHS_out" as keys
         
     - list_syn_fuels list of strings of the identifier of synthetic fuels which output is *not* accounted
       in the balance of the elctric bus. The production of these syn fuels may require electric input,
       and in that case this input is accounted in the electric input total profile
    
    """
    time_steps = len(network.snapshots)
    tot_storage_out_t = np.zeros(time_steps)
    tot_storage_in_t = np.zeros(time_steps)
    fraction_gen_to_storage_t = np.zeros(time_steps)
    fraction_gen_to_dispatch_t = np.ones(time_steps)
    storage_to_storage_t = np.zeros(time_steps)
    fraction_storage_to_dispatch_t = np.zeros(time_steps)
    tot_storage_dispatched_t = np.zeros(time_steps)
    
    map_indices_syn_fuel = {}
    if list_syn_fuels is not None:
        num_syn_fuels = len(list_syn_fuels)
        for i, l in enumerate(list_syn_fuels):
            map_indices_syn_fuel[l] = i
    
    syn_fuel_out                        = np.zeros(num_syn_fuels)
    syn_fuel_capex_plus_opex_plus_input = np.zeros(num_syn_fuels)
    syn_fuel_market_value               = np.zeros(num_syn_fuels)
    
    if all_storage_tech_list is not None:
        num_all_storage_tech_list = len(all_storage_tech_list)
        storage_in                         = np.zeros(num_all_storage_tech_list) # Only electricty
        storage_any_in                     = np.zeros(num_all_storage_tech_list) # Any energy vector
        storage_out                        = np.zeros(num_all_storage_tech_list) # Only electricty
        storage_dispatched                 = np.zeros(num_all_storage_tech_list) # Only electricty, net of storage-to-storage flow
        storage_any_out                    = np.zeros(num_all_storage_tech_list) # Any energy vector
        storage_any_out                    = np.zeros(num_all_storage_tech_list)
        storage_cap                        = np.zeros(num_all_storage_tech_list)
        storage_in_pnom                    = np.zeros(num_all_storage_tech_list)
        storage_in_peff                    = np.zeros(num_all_storage_tech_list)
        storage_out_pnom                   = np.zeros(num_all_storage_tech_list)
        storage_out_peff                   = np.zeros(num_all_storage_tech_list)
        storage_max_discharge_duration     = np.zeros(num_all_storage_tech_list)
        storage_full_cycles_eq_per_period  = np.zeros(num_all_storage_tech_list)
        storage_capex                      = np.zeros(num_all_storage_tech_list)
        storage_opex                       = np.zeros(num_all_storage_tech_list)
        storage_input_cost                 = np.zeros(num_all_storage_tech_list)
        storage_out_revenue                = np.zeros(num_all_storage_tech_list)
        storage_capex_plus_opex            = np.zeros(num_all_storage_tech_list)
        storage_capex_plus_opex_plus_input = np.zeros(num_all_storage_tech_list)
        storage_market_value               = np.zeros(num_all_storage_tech_list)
        balance                            = np.zeros(num_all_storage_tech_list)
        storage_CF_in                      = np.zeros(num_all_storage_tech_list)
        storage_CF_out                     = np.zeros(num_all_storage_tech_list)
        investment                         = np.zeros(num_all_storage_tech_list)
        investment_power                   = np.zeros(num_all_storage_tech_list)
        investment_energy                  = np.zeros(num_all_storage_tech_list)

        
        for i, s_l in enumerate(all_storage_tech_list):
            if s_l in storage_single_link_list:
                # then this tech is represented as a single link
                s_l_store = None
                s_l_in = s_l
                s_l_out = s_l
            else:
                s_l_store = s_l + "_store"
                s_l_in = s_l + "_in"
                s_l_out = s_l + "_out"
            
            b_l_in   = network.links.bus0[s_l_in] # label of the Bus in 
            b_l_in1   = network.links.bus1[s_l_in] # label of the Bus in1
            b_l_out0   = network.links.bus0[s_l_out] # label of the Bus out0
            b_l_out  = network.links.bus1[s_l_out] # label of the Bus out
            b2_id = ""
            
            storage_out[i] = -network.links_t.p1[s_l_out].sum()/1e6     # TWh
            storage_in[i] = network.links_t.p0[s_l_in].sum()/1e6        # TWh
            this_tech_aux_in = 0.
            this_tech_aux_cost = 0.
            if s_l in storage_with_bus2:
                # then there is an auxiliary input
                # it could be an electric input e.g.for H2 compression 
                # or a feedstock input
                b2_id = network.links.bus2[s_l_in]
                if b2_id == "AC":
                    # account the electric input only as TWh
                    this_tech_aux_in = network.links_t.p2[s_l_in].sum()/1e6        # TWh
                this_tech_aux_cost  = network.links_t.p2[s_l_in].multiply(network.snapshot_weightings.generators, axis=0).multiply(network.buses_t.marginal_price[b2_id], axis=0).sum()/1.e+9 # in G€

            storage_any_in[i] = storage_in[i] + this_tech_aux_in # TWh
            storage_any_out[i] = storage_out[i] # TWh
            storage_in_pnom[i] = network.links.p_nom_opt[s_l_in]/1.e+3  # GW
            storage_in_peff[i] = storage_in_pnom[i] * network.links.efficiency[s_l_in] # GW
            storage_out_pnom[i] = network.links.p_nom_opt[s_l_out]/1.e+3  # GW 
            storage_discharge_eff = network.links.efficiency[s_l_out]
            storage_out_peff[i] = storage_out_pnom[i] * storage_discharge_eff  # GW
            # check if power in has a capital cost
            if s_l_in in dict_pypsa_par_id.keys():
                if b_l_in == 'AC' and b_l_out != "heat":
                    # then the investment is proportional to the P nominal in
                    investment_power[i] = (storage_in_pnom[i] * 1.e+3 *  tech_parameters["investment"][dict_pypsa_par_id[s_l_in]]) /1.e+9 #G€
                else:
                    # otherwise you have to scale it by eff
                    investment_power[i] = (storage_in_peff[i] * 1.e+3 *  tech_parameters["investment"][dict_pypsa_par_id[s_l_in]]) /1.e+9 #G€
            if s_l_out in dict_pypsa_par_id.keys() and s_l not in storage_single_link_list:
                investment_power[i] += (storage_out_peff[i] * 1.e+3 *  tech_parameters["investment"][dict_pypsa_par_id[s_l_out]]) /1.e+9 #G€
            storage_CF_in[i]   = storage_any_in[i] * 1.e+3 / storage_in_pnom[i] / time_steps
            storage_CF_out[i]   = storage_any_out[i] * 1.e+3 / storage_out_peff[i] / time_steps
            if s_l_store is not None:
                Cap_available     = network.stores.e_nom_opt[s_l_store]/1e+3 #GWh
                Cap_max_used      = network.stores_t.e[s_l_store].max()/1e+3 #GWh
                # Check if storage is significantly used (it can be a near zero value)
                if Cap_max_used < 1.e-9:
                    Cap_max_used = 0.
                storage_cap[i] = Cap_max_used  # GWh
                if s_l_store in dict_pypsa_par_id.keys():
                    investment_energy[i] = (storage_cap[i] *1.e+3 * tech_parameters["investment"][dict_pypsa_par_id[s_l_store]]) /1.e+9 #G€
                    print("s_l_store", s_l_store, storage_cap[i],  tech_parameters["investment"][dict_pypsa_par_id[s_l_store]])
                investment[i] = investment_power[i] + investment_energy[i]
                storage_max_discharge_duration[i] = storage_cap[i] / storage_out_peff[i] #h
            
                if Cap_max_used > 0. :
                    storage_full_cycles_eq_per_period[i] = storage_out[i] *1.e+3 / storage_cap[i]
            
                storage_capex[i] = (storage_cap[i] * network.stores.capital_cost[s_l_store] * 1.e+3 
                                + storage_in_pnom[i] * network.links.capital_cost[s_l_in] * 1.e+3 
                                + storage_out_pnom[i] * network.links.capital_cost[s_l_out] * 1.e+3) /1.e+9 # in G€
                
                storage_opex[i] = (storage_cap[i] * network.stores.marginal_cost[s_l_store] * 1.e+3 
                               + storage_in[i] * network.links.marginal_cost[s_l_in] * 1.e+6 
                               + network.links_t.p0[s_l_out].sum() * network.links.marginal_cost[s_l_out]) /1.e+9 # in G€
            else:
                # This is a single link tech
                storage_capex[i] = (storage_in_pnom[i] * network.links.capital_cost[s_l_in] * 1.e+3) /1.e+9 # in G€
                storage_opex[i]  = (storage_in[i]  * 1.e+6  * network.links.marginal_cost[s_l_in] )/1.e+9 # in G€
                investment[i]    = investment_power[i] 
            
            # Avoid that a near zero or negative value because of precision propagates an error
            if storage_capex[i] < 0.0000000001:
                storage_capex[i] = 0.     
            
            storage_input_cost[i]  = this_tech_aux_cost + network.links_t.p0[s_l_in].multiply(network.snapshot_weightings.generators, axis=0).multiply(network.buses_t.marginal_price[b_l_in], axis=0).sum()/1.e+9 # in G€
            storage_out_revenue[i] =  -network.links_t.p1[s_l_out].multiply(network.snapshot_weightings.generators, axis=0).multiply(network.buses_t.marginal_price[b_l_out], axis=0).sum()/1.e+9 # in G€          
            storage_capex_plus_opex[i] = storage_capex[i] + storage_opex[i] # G€          
            storage_capex_plus_opex_plus_input[i] = storage_capex_plus_opex[i] + storage_input_cost[i]
            balance[i] = storage_capex_plus_opex_plus_input[i] - storage_out_revenue[i]
            
            if storage_any_out[i] > 0:
                storage_market_value[i] = storage_capex_plus_opex_plus_input[i] / storage_any_out[i] * 1.e+3 # €/MWh
            
            if storage_capex_plus_opex_plus_input[i] > 0:
                if np.abs(balance[i]) > 0.00001:
                    print("verify balance at storage tech " + s_l, balance[i], storage_capex_plus_opex_plus_input[i],  storage_out_revenue[i])
                else:
                    print("Ok revenue-cost balance for tech " + s_l + " market value (€/MWh) ", storage_market_value[i])
    
            
            if b_l_in not in ["AC"]:
                # then it's indirect storage, not related to the electricity bus
                # therefore set to zero the data for the electric balance
                storage_in[i] = this_tech_aux_in
                # eventually add aux power
                if s_l in storage_with_bus2 and b2_id == 'AC':
                     tot_storage_in_t +=  network.links_t.p2[s_l_in] # MWh
            else: 
                # only for storage directly connected to the electric bus account for in out profiles
                tot_storage_in_t +=  network.links_t.p0[s_l_in] # MWh
            
            # store output data in other array to distinquish between electric and general energy vector
            storage_any_out[i] = storage_out[i]
            if b_l_out not in ["AC"]:
                # then it's indirect storage, not related to the electricity bus
                # therefore set to zero the data for the electric balance
                storage_out[i] = 0
            else: 
                # only for storage directly connected to the electric bus account for in out profiles
                tot_storage_out_t +=  -network.links_t.p1[s_l_out] # MWh
            
    # Compute fractions from generation to storage and dispatch
    tot_generator_out_t = df_sys_indices.loc["value"]["tot_generator_out_t"]
    for i in range(time_steps):
        if tot_generator_out_t[i] >= tot_storage_in_t[i] and tot_storage_in_t[i] > 0:
            fraction_gen_to_storage_t[i] = tot_storage_in_t[i]/tot_generator_out_t[i]
        elif tot_generator_out_t[i] < tot_storage_in_t[i] and tot_storage_in_t[i] > 0:
            fraction_gen_to_storage_t[i] = 1.
            from_storage = tot_storage_in_t[i] - tot_generator_out_t[i]
            if from_storage <= tot_storage_out_t[i]:
                storage_to_storage_t[i] = from_storage
                fraction_storage_to_dispatch_t[i] = from_storage/tot_storage_out_t[i]
                tot_storage_out_t[i] -= from_storage
            else:
                print("missing", from_storage)
        else:
            fraction_gen_to_storage_t[i] = 0.
    if storage_to_storage_t.sum()>0:
        print("Check storage to storage (MWh)", storage_to_storage_t.sum())
    fraction_gen_to_dispatch_t -= fraction_gen_to_storage_t
    
    for i, s_l in enumerate(all_storage_tech_list):
        if s_l in storage_single_link_list:
                # then this tech is represented as a single link
            s_l_out = s_l
        else:
            s_l_out = s_l + "_out"
        b_l_out  = network.links.bus1[s_l_out] # label of the Bus out
        if b_l_out in ["AC"]:       
            this_storage_out_t        = -network.links_t.p1[s_l_out]
            this_storage_dispatched_t = np.subtract(this_storage_out_t, np.multiply(this_storage_out_t, fraction_storage_to_dispatch_t))
            storage_dispatched[i]     = this_storage_dispatched_t.sum()/1e6     # TWh
            tot_storage_dispatched_t += this_storage_dispatched_t
            
    Tot_el_from_storage = storage_out.sum() #TWh
    Tot_el_to_storage = storage_in.sum() #TWh, Tot_el_to_storage
    print("Total electricity from and to storage (TWh)", Tot_el_from_storage, Tot_el_to_storage)
    print("Avg efficiency %.2f%%" % round(Tot_el_from_storage/Tot_el_to_storage * 100, 3))
    print("storage_capex_plus_opex (G€)", storage_capex_plus_opex)
    

    df_storage_results = pd.DataFrame({
                        "storage_el_in" : pd.Series(storage_in, index=all_storage_tech_list, copy=True),
                        "storage_in_pnom" : pd.Series(storage_in_pnom, index=all_storage_tech_list, copy=True),
                        "storage_out_peff" : pd.Series(storage_out_peff, index=all_storage_tech_list, copy=True),
                        "storage_full_cycles_eq_per_period" : pd.Series(storage_full_cycles_eq_per_period, index=all_storage_tech_list, copy=True),
                        "storage_cap" : pd.Series(storage_cap, index=all_storage_tech_list, copy=True),
                        "storage_max_discharge_duration" : pd.Series(storage_max_discharge_duration, index=all_storage_tech_list, copy=True),
                        "storage_el_out" : pd.Series(storage_out, index=all_storage_tech_list, copy=True),
                        "storage_dispatched" : pd.Series(storage_dispatched, index=all_storage_tech_list, copy=True),
                        "storage_capex" : pd.Series(storage_capex, index=all_storage_tech_list, copy=True),
                        "storage_opex" : pd.Series(storage_opex, index=all_storage_tech_list, copy=True),
                        "storage_capex_plus_opex" : pd.Series(storage_capex_plus_opex, index=all_storage_tech_list, copy=True),
                        "storage_input_cost" : pd.Series(storage_input_cost, index=all_storage_tech_list, copy=True),
                        "storage_out_revenue" : pd.Series(storage_out_revenue, index=all_storage_tech_list, copy=True),
                        "storage_any_out" : pd.Series(storage_any_out, index=all_storage_tech_list, copy=True),
                        "storage_market_value" : pd.Series(storage_market_value, index=all_storage_tech_list, copy=True),
                        "storage_CF_in" : pd.Series(storage_CF_in, index=all_storage_tech_list, copy=True),
                        "storage_CF_out" : pd.Series(storage_CF_out, index=all_storage_tech_list, copy=True),
                        "investment"        : pd.Series(investment, index=all_storage_tech_list, copy=True),
                        "investment_power"  : pd.Series(investment_power, index=all_storage_tech_list, copy=True),
                        "investment_energy" : pd.Series(investment_energy, index=all_storage_tech_list, copy=True),


    })
    
    df_storage_results = df_storage_results.assign(price = df_storage_results.storage_out_revenue/df_storage_results.storage_any_out *1.e+3) # €/MWh
    df_storage_results = df_storage_results.assign(capex_unit = df_storage_results.storage_capex/df_storage_results.storage_any_out *1.e+3) # €/MWh
    df_storage_results = df_storage_results.assign(opex_unit = df_storage_results.storage_opex/df_storage_results.storage_any_out *1.e+3) # €/MWh
    df_storage_results = df_storage_results.assign(input_cost_unit = df_storage_results.storage_input_cost/df_storage_results.storage_any_out *1.e+3) # €/MWh

    df_sys_indices["tot_storage_in_t"] = ["MW", tot_storage_in_t]
    df_sys_indices["tot_storage_out_t"] = ["MW", tot_storage_out_t]
    df_sys_indices["fraction_gen_to_storage_t"] = ["/", fraction_gen_to_storage_t]
    df_sys_indices["fraction_gen_to_dispatch_t"] = ["/", fraction_gen_to_dispatch_t]
    df_sys_indices["storage_to_storage_t"] = ["/", storage_to_storage_t]
    df_sys_indices["tot_storage_dispatched_t"] = ["/", tot_storage_dispatched_t]
    
    return df_sys_indices, df_storage_results  



    


# Write down the solved network for further processing in other notebooks
#network.export_to_netcdf("./networks/Calabria_i08_noPHS_NoCAES_noSynCH4_1GWcost.nc")
#myc = colors_distinguishable_by_color_blind_person(5)
#print(myc)






def determine_gen_directly_dispatched(network, 
                                      non_fossil_gen_labels=None,  
                                      fossil_gen_labels=None, 
                                      links_gen_labels=None, 
                                      primary_fuel_label=None, 
                                      df_sys_indices=None,
                                      df_gen_results=None,
                                      df_bus_results=None):
    """
    -
    - 
    """
    
    fraction_gen_to_dispatch_t = df_sys_indices.loc["value"]["fraction_gen_to_dispatch_t"]
    # data structures for all generators
    tot_gen_directly_out_t          = np.zeros(len(network.snapshots))
    tot_gen_to_storage_t            = np.zeros(len(network.snapshots))
    tot_gen_from_indirect_storage_t = np.zeros(len(network.snapshots))

    all_gen_plus_link_label = []
    all_gen_label = []
    num_gen_non_fossil = 0
    if non_fossil_gen_labels is not None:
        num_gen_non_fossil = len(non_fossil_gen_labels)
        all_gen_label = non_fossil_gen_labels[:]
    
    num_gen_fossil = 0
    if fossil_gen_labels is not None:
        num_gen_fossil = len(fossil_gen_labels)
        all_gen_label += fossil_gen_labels[:]
    
    num_gen_links = 0
    if links_gen_labels is not None:
        num_gen_links = len(links_gen_labels)
        all_gen_plus_link_label = all_gen_label + links_gen_labels
    
    num_all_gen =  num_gen_non_fossil + num_gen_fossil
    num_all_gen_plus_link = len(all_gen_plus_link_label)
    
    generator_directly_dispatched   = np.zeros(num_all_gen_plus_link)
    generator_to_storage            = np.zeros(num_all_gen_plus_link)
    generator_from_indirect_storage = np.zeros(num_all_gen_plus_link)
    
    
    
    for i, my_l in enumerate(all_gen_plus_link_label):
        
        this_gen_all_out_t                 = np.zeros(len(network.snapshots))
        this_gen_primary_out_t             = np.zeros(len(network.snapshots))
        this_gen_directly_out_t            = np.zeros(len(network.snapshots))
        this_gen_to_storage_t              = np.zeros(len(network.snapshots))
        this_gen_from_indirect_storage_t   = np.zeros(len(network.snapshots))
        
        if i < num_all_gen:
            # Then data are stored in the generators component
            this_gen_all_out_t      = network.generators_t.p[my_l]
            this_gen_primary_out_t  = this_gen_all_out_t
            this_gen_directly_out_t = np.multiply(this_gen_all_out_t, fraction_gen_to_dispatch_t)
            this_gen_to_storage_t   = np.subtract(this_gen_all_out_t, this_gen_directly_out_t)
        else:
            # Then data are in the links component and we have to discern between primary and indirect cshares
            this_link_bus_in        =  network.links.bus0[my_l]
            this_gen_all_out_t      = - network.links_t.p1[my_l]
            this_gen_primary_out_t  =   this_gen_all_out_t * df_bus_results.loc[this_link_bus_in]["primary_share"]
            this_gen_directly_out_t = np.multiply(this_gen_primary_out_t, fraction_gen_to_dispatch_t)
            this_gen_to_storage_t   = np.subtract(this_gen_primary_out_t, this_gen_directly_out_t)
            # Then add the non primary share, if any
            this_gen_from_indirect_storage_t = this_gen_all_out_t * (1. - df_bus_results.loc[this_link_bus_in]["primary_share"])
            this_gen_directly_out_t += this_gen_from_indirect_storage_t
            
        tot_gen_directly_out_t        += this_gen_directly_out_t
        tot_gen_to_storage_t          += this_gen_to_storage_t
        tot_gen_from_indirect_storage_t += this_gen_from_indirect_storage_t
        
        generator_directly_dispatched[i]   =  this_gen_directly_out_t.sum()/1.e+6 # TWh
        generator_to_storage[i]            =  this_gen_to_storage_t.sum()/1.e+6 # TWh
        generator_from_indirect_storage[i] =  this_gen_from_indirect_storage_t.sum()/1.e+6 # TWh
        
        
        df_sys_indices["tot_gen_directly_out_t"] = ["MW", tot_gen_directly_out_t]
        df_sys_indices["tot_gen_to_storage_t"]   = ["MW", tot_gen_to_storage_t]
        df_sys_indices["tot_gen_from_indirect_storage_t"] = ["MW", tot_gen_from_indirect_storage_t]
        df_sys_indices["tot_gen_directly_out_t"] = ["MW", tot_gen_directly_out_t]
        
        df_gen_results["generator_directly_dispatched"] = generator_directly_dispatched
        df_gen_results["generator_to_storage"] = generator_to_storage
        df_gen_results["generator_from_indirect_storage"] = generator_from_indirect_storage
        

    return df_gen_results, df_sys_indices



def plot_generators_and_storage(df_storage_results, df_gen_results):

    
    cat_displayed = 3
    #category_colors = ["y", "g", "b"]
    category_colors = colors_distinguishable_by_color_blind_person(cat_displayed)

    categories = ["Power input", "Energy or mass", "Power or mass flow output"]
    used_col_name = [False] * (len(categories) + 1)

    temp_df  = df_gen_results.loc[:]['generator_installed_capacity']
    temp_df  = temp_df[temp_df > 0]
    temp_df  = temp_df.sort_values()
    # count how many non-zero gen you have
    nz_gen = temp_df.count()
    
    temp_df2  = pd.concat([df_storage_results.loc[:]['storage_in_pnom'],df_storage_results.loc[:]['storage_out_peff'],df_storage_results.loc[:]['storage_cap'],df_storage_results.loc[:]['storage_max_discharge_duration']],axis=1)
    temp_df2  = temp_df2[temp_df2['storage_in_pnom'] > 1.e-6] #in_pnom is in GW hence threshold is 1 kW
    temp_df2  = temp_df2.sort_values(by=["storage_cap"])
    # count how many non-zero storage you have
    nz_sto = temp_df2['storage_in_pnom'].count()
    multiplier = 2
    bar_widht = 0.5

    y_range = np.arange(0, nz_gen)
    y_range2 = np.arange(nz_gen, nz_gen + nz_sto * multiplier, multiplier)
    y_range2_mid = y_range2 + bar_widht
    y_range_all = [*y_range[:], *y_range2_mid[:]]
    
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.invert_yaxis()
    ax2 = ax.twiny()
    
    ax.set_xscale("log")
    
    ax.set_xlabel(r"MWh or tCO$_2$")
    ax2.set_xlabel(r"MW or tCO$_2$ per hour")


    nice_tec_labels = []
    
    for i, rowname in enumerate(temp_df.index):
        p_installed = temp_df.loc[rowname]
        nice_name = dict_all_tech_nice_labels[rowname]
        if not used_col_name[2]:
            used_col_name[2] = True
            rect = ax2.barh(y_range[i], p_installed, left=0, height=bar_widht, 
                                    color=category_colors[2], label=categories[2])
        else:
            rect = ax2.barh(y_range[i], p_installed, left=0, height=bar_widht, color=category_colors[2])
        nice_tec_labels += [nice_name]
    
    # Now do also the storage techs
    for i, rowname in enumerate(temp_df2.index):
        p_in = temp_df2.loc[rowname]['storage_in_pnom'] * 1.e+3
        e_cap = temp_df2.loc[rowname]['storage_cap'] * 1.e+3
        p_out = temp_df2.loc[rowname]['storage_out_peff'] * 1.e+3 
        max_discharge =  np.round(temp_df2.loc[rowname]['storage_max_discharge_duration'],0)
        max_discharge_label = '(%.d h)' % max_discharge
        nice_name = dict_all_tech_nice_labels[rowname]
        if not used_col_name[0]:
            used_col_name[0] = True
            rect = ax2.barh(y_range2[i], p_in, left=0, height=bar_widht, color=category_colors[0], label=categories[0])
        else:
            rect = ax2.barh(y_range2[i], p_in, left=0, height=bar_widht, color=category_colors[0])
        if not used_col_name[1]:
            used_col_name[1] = True
            rect = ax.barh(y_range2[i] + bar_widht, e_cap, left=0, height=bar_widht, color=category_colors[1], label=categories[1])
        else:
            rect = ax.barh(y_range2[i] + bar_widht, e_cap, left=0, height=bar_widht, color=category_colors[1])
        ax.bar_label(rect, labels=[max_discharge_label], label_type='edge', padding = 2)
        
        
        
        rect = ax2.barh(y_range2[i] + 2 * bar_widht, p_out, left=0, height=bar_widht, color=category_colors[2])
        nice_tec_labels += [nice_name]

    ax.legend(ncol=1, bbox_to_anchor=(0, -0.05), loc='upper left', fontsize='small')
    ax2.legend(ncol=2, bbox_to_anchor=(0, 1.1), loc='lower left', fontsize='small')

    ax.set_yticks(y_range_all, nice_tec_labels)
    my_title = "Installed capacities: \n power or mass flow (top axis, linear scale) \n and energy or mass (bottom axis, log scale, \n and in parenthesis the maximum discharge duration in hours, rounded to the nearest integer)"
    fig.suptitle(t=my_title, fontweight='bold', y=1.1)
    plt.tight_layout()


    return fig


def cost_vs_revenue_chart(df_storage_results, df_gen_results
                         ):
    
    """
    Parameters
    ----------
    df_storage_results : 
    
    """
    
    cat_displayed = 4
    map_cat_color = np.arange(cat_displayed)
    
    category_colors = colors_distinguishable_by_color_blind_person(cat_displayed)
    
    fig, ax = plt.subplots(figsize=(12, 3.5))
    ax.invert_yaxis()


    #ax.xaxis.set_visible(False)
    ax.set_xlabel("M€/y")
    
    temp_df = pd.concat([df_storage_results.loc[:]["storage_capex"], df_storage_results.loc[:]["storage_opex"], df_storage_results.loc[:]["storage_input_cost"]],axis=1)
    temp_df = pd.concat([temp_df, df_storage_results.loc[:]["storage_out_revenue"]],axis=1)
    temp_df.columns = ["capex","opex","input","revenue"]
    
    temp_df2 = pd.concat([df_gen_results.loc[:]["generator_capex"], df_gen_results.loc[:]["generator_marginal_cost"]], axis=1)
    temp_df2 = pd.concat([temp_df2, df_gen_results.loc[:]["generator_input_cost"], df_gen_results.loc[:]["generator_revenue"]], axis=1)
    temp_df2.columns = ["capex","opex","input","revenue"]
    temp_df = pd.concat([temp_df, temp_df2], axis=0)
    temp_df = temp_df.sort_values(by=["revenue"])
    
    categories = ["CAPEX + FO&M", "VO&M", "input"]
    used_col_name = [False] * (len(categories) + 1)
    epsilon = 0.0001
    for i, rowname in enumerate(temp_df.index):
        # Avoid displaying zero categories, approx to a very small values
        if temp_df.capex[rowname] > epsilon:
            start = 0.
            nice_name = dict_all_tech_nice_labels[rowname]
            for j in range(len(categories)):
                width = temp_df.loc[rowname][j] * 1.e+3 #M€/y
                if not used_col_name[j]:
                    used_col_name[j] = True
                    rect = ax.barh(nice_name, width, left=start, height=0.8,
                                    color=category_colors[j], label=categories[j])
                else:
                    rect = ax.barh(nice_name, width, left=start, height=0.8,
                                    color=category_colors[j])
                start+= width
            if not used_col_name[3]:
                used_col_name[3] = True
                ax.scatter(temp_df.loc[rowname]["revenue"] * 1.e+3, nice_name, marker="D", color=category_colors[3], label="revenue")
            else:
                ax.scatter(temp_df.loc[rowname]["revenue"] * 1.e+3, nice_name, marker="D", color=category_colors[3])

    
    ax.legend(ncol=cat_displayed, bbox_to_anchor=(0, 1), loc='lower left', fontsize='small')
    my_title = "Yearly revenue and cost components"
    fig.suptitle(t=my_title, fontweight='bold', y=1.1)
    plt.tight_layout()


    return fig



def cost_vs_price_chart(df_storage_results, df_gen_results
                         ):
    
    """
    Parameters
    ----------
    df_storage_results : 
    
    """
    
    cat_displayed = 4
    map_cat_color = np.arange(cat_displayed)
    
    category_colors = colors_distinguishable_by_color_blind_person(cat_displayed)
    
    fig, ax = plt.subplots(figsize=(12, 3.5))
    ax.invert_yaxis()
    #ax.xaxis.set_visible(False)
    ax.set_xlabel(r"€/MWh or €/tCO$_2$")
    
    temp_df = pd.concat([df_storage_results.loc[:]["capex_unit"], df_storage_results.loc[:]["opex_unit"], df_storage_results.loc[:]["input_cost_unit"]],axis=1)
    temp_df = pd.concat([temp_df, df_storage_results.loc[:]["price"], df_storage_results.loc[:]["storage_CF_in"]],axis=1)
    temp_df.columns = ["capex","opex","input","price","CF"]
    
    temp_df2 = pd.concat([df_gen_results.loc[:]["capex_unit"], df_gen_results.loc[:]["opex_unit"]], axis=1)
    temp_df2 = pd.concat([temp_df2, df_gen_results.loc[:]["input_cost_unit"], df_gen_results.loc[:]["price"], df_gen_results.loc[:]["generator_any_CF"]], axis=1)
    temp_df2.columns = ["capex","opex","input","price","CF"]
    temp_df = pd.concat([temp_df, temp_df2], axis=0)
    temp_df = temp_df.sort_values(by=["price"])
    
    categories = ["CAPEX + FO&M", "VO&M", "input"]
    used_col_name = [False] * (len(categories) + 1)
    epsilon = 0.0001
    for i, rowname in enumerate(temp_df.index):
        # Avoid displaying zero categories, approx to a very small values
        if temp_df.capex[rowname] > epsilon:
            start = 0.
            nice_name = dict_all_tech_nice_labels[rowname]
            for j in range(len(categories)):
                if not used_col_name[j]:
                    used_col_name[j] = True
                    rect = ax.barh(nice_name, temp_df.loc[rowname][j], left=start, height=0.8,
                                    color=category_colors[j], label=categories[j])
                else:
                    rect = ax.barh(nice_name, temp_df.loc[rowname][j], left=start, height=0.8,
                                    color=category_colors[j])
                start+= temp_df.loc[rowname][j]
            this_CF = np.round(temp_df.loc[rowname]["CF"] * 100,0)
            CF_label = '(%.0f %%)' % this_CF
            ax.bar_label(rect, labels=[CF_label], label_type='edge', padding = 5)
            if not used_col_name[3]:
                used_col_name[3] = True
                ax.scatter(temp_df.loc[rowname]["price"], nice_name, marker="D", color=category_colors[3], label="average price")
            else:
                ax.scatter(temp_df.loc[rowname]["price"], nice_name, marker="D", color=category_colors[3])
    
    ax.legend(ncol=cat_displayed, bbox_to_anchor=(0, 1), loc='lower left', fontsize='small')
    #ax.legend(bbox_to_anchor=(0, 1), loc='lower left', fontsize='small')
    my_title = "Average price and cost components per unit of service delivered, \n in parenthesis the net capacity factor rounded to the nearest integer"
    fig.suptitle(t=my_title, fontweight='bold', y=1.1)
    plt.tight_layout()


    return fig


def make_some_nice_charts_statistics(network,
                                     links_gen_labels = None,
                                     primary_fuel_label = None,
                                     indirect_fuel_label = None,
                                     all_storage_tech_list=None,
                                     storage_single_link_list = None,
                                     list_syn_fuels = None,
                                     storage_with_bus2 = None,
                                     link_to_general_gas_labels = None,
                                     bus_general_gas_label = None,
                                    ):
    
    
    df_sys_indices   = basic_network_statistics(network)
    
    Tot_el_delivered = df_sys_indices.loc["value"]["Tot_elec_delivered"]
    Ann_system_cost  = df_sys_indices.loc["value"]["Ann_system_cost"]
    LCOS             = df_sys_indices.loc["value"]["LCOS"]
    mu_CO2           = df_sys_indices.loc["value"]["mu_CO2"]
    
    df_bus_results, df_fuel_primary_results, df_gen_results, df_sys_indices =  advanced_network_statistics_V2(network, 
                            non_fossil_gen_labels=non_fossil_gen_labels, 
                            non_fossil_gen_CF=non_fossil_gen_CF, 
                            fossil_gen_labels=fossil_gen_labels,
                            links_gen_labels = links_gen_labels,
                            primary_fuel_label = primary_fuel_label,
                            indirect_fuel_label = indirect_fuel_label,
                            df_sys_indices = df_sys_indices,
                            link_to_general_gas_labels = link_to_general_gas_labels,
                            bus_general_gas_label = bus_general_gas_label,
                          )
    Tot_el_generated = df_sys_indices.loc["value"]["Tot_el_generated"]
    curtailment_pc   = df_sys_indices.loc["value"]["curtailment_pc"]
    CO2_intensity    = df_sys_indices.loc["value"]["tot_fossil_CO2"] / Tot_el_delivered * 1.e+3 # gCO2/kWh or kg/MWh or Mt/TWh, whatever
    df_sys_indices["CO2_intensity"] = ["gCO2/kWh", CO2_intensity]
    
    
    #print(network.generators.p_nom_opt) # in MW

    # Info from storage techs balances
    df_sys_indices, df_storage_results = storage_el_balance(network, all_storage_tech_list = all_storage_tech_list,
                                                list_syn_fuels = list_syn_fuels,
                                                df_sys_indices = df_sys_indices,
                                                storage_single_link_list = storage_single_link_list,
                                                storage_with_bus2 = storage_with_bus2,
                                        )



 
    df_gen_results, df_sys_indices = determine_gen_directly_dispatched(network, 
                              non_fossil_gen_labels=non_fossil_gen_labels,  
                              fossil_gen_labels=fossil_gen_labels, 
                              links_gen_labels=links_gen_labels, 
                              primary_fuel_label=primary_fuel_label,
                              df_sys_indices=df_sys_indices,
                              df_gen_results=df_gen_results,
                              df_bus_results=df_bus_results)
    
    print(df_sys_indices)

    
    # Update storage generation with gen from syn fuels, if any (this is computed in the generation procedure)
    all_tech_dispatch_shares = [*df_gen_results[:]["generator_directly_dispatched"] /float(Tot_el_delivered) * 100, *df_storage_results.loc[:]["storage_dispatched"]/float(Tot_el_delivered) * 100, *np.zeros(len(primary_fuel_label))]

    Reconstructed_Sys_Cost = df_storage_results.loc[:]["storage_capex_plus_opex"].sum() + df_sys_indices.loc["value"]["Ann_system_cost_g"] # G€
    threshold = np.abs(Reconstructed_Sys_Cost - Ann_system_cost)
    if threshold > 0.001:
        print("Verify Syst Cost, reconstructed %.5f vs objective function %.5f, abs diff %.5f" % (Reconstructed_Sys_Cost, Ann_system_cost, threshold))
    else:
        print("Ok Syst Cost, reconstructed %.5f vs objective function %.5f, abs diff %.5f" % (Reconstructed_Sys_Cost, Ann_system_cost, threshold))

    # Verify dispatching balances
    method_v1 = df_sys_indices.loc["value"]["tot_generator_out_t"].sum() - (df_sys_indices.loc["value"]["tot_gen_directly_out_t"].sum() - df_sys_indices.loc["value"]["tot_gen_from_indirect_storage_t"].sum()) - df_sys_indices.loc["value"]["tot_storage_in_t"].sum() + df_sys_indices.loc["value"]["storage_to_storage_t"].sum()
    method_v2 = df_gen_results[:]["generator_directly_dispatched"].sum() * 1.e+6 + df_sys_indices.loc["value"]["tot_storage_dispatched_t"].sum() - Tot_el_delivered * 1.e+6
    print("tot_generator_out_t", df_sys_indices.loc["value"]["tot_generator_out_t"].sum())
    print("tot_gen_directly_out_t", - df_sys_indices.loc["value"]["tot_gen_directly_out_t"].sum())
    print("tot_gen_from_indirect_storage_t", df_sys_indices.loc["value"]["tot_gen_from_indirect_storage_t"].sum()) 
    print("tot_storage_in_t", - df_sys_indices.loc["value"]["tot_storage_in_t"].sum())
    print("storage_to_storage_t", df_sys_indices.loc["value"]["storage_to_storage_t"].sum())
    print("method_v1 ", method_v1)
    print("generator_directly_dispatched", df_gen_results[:]["generator_directly_dispatched"].sum() * 1.e+6 )
    print("tot_storage_dispatched_t", df_sys_indices.loc["value"]["tot_storage_dispatched_t"].sum())
    print("Tot_el_delivered ", - Tot_el_delivered * 1.e+6)
    print("method_v2 ", method_v2 )
    if abs(method_v1) > 1:
        print("Verify dispatched V1", method_v1)
    else:
        print("Ok dispatch balance v1")
    if abs(method_v2) > 1:
        print("Verify dispatched V2", method_v2)
    else:
        print("Ok dispatch balance v2") 
    generator_LCOS_shares = np.array(df_gen_results.loc[:]["generator_capex_plus_opex"])/Ann_system_cost * 100.
    storage_LCOS_shares = df_storage_results.loc[:]["storage_capex_plus_opex"]/Ann_system_cost * 100.  
    primary_fuel_LCOS_shares = df_fuel_primary_results.loc[:]["primary_fuel_costs"] / Ann_system_cost * 100.
    all_tech_LCOS_shares = [*generator_LCOS_shares, *storage_LCOS_shares, *primary_fuel_LCOS_shares]
    #print(all_tech_LCOS_shares)  

    all_tech_primary_output_shares = [*np.array(df_gen_results.loc[:]["generator_output_shares"]) * 100., *np.zeros(len(all_storage_tech_list)), *df_fuel_primary_results.loc[:]["primary_fuel_out_shares"] * 100.]
    print("all_tech_primary_output_shares", all_tech_primary_output_shares) 

    results = {
        'Primary generation \n (%.2f TWh/y)' %float(Tot_el_generated) : np.round(all_tech_primary_output_shares,2),
        'Delivered electricity  \n (%.2f TWh/y)' %   float(Tot_el_delivered)     : np.round(all_tech_dispatch_shares,2),
        'Levelized Cost of System \n LCOS  (%.f €/MWh) ' % np.round(LCOS,0) : np.round(all_tech_LCOS_shares,2)
    }

    
    
    #Verify results up to a threshold
    threshold =  np.abs(np.array(all_tech_primary_output_shares).sum() - 100)
    if threshold > 0.01:
        print("Check Gen ", all_tech_primary_output_shares, threshold)
    else:
        print("Ok shares Gen ", all_tech_primary_output_shares, threshold)
    threshold =  np.abs(np.array(all_tech_dispatch_shares).sum() - 100)
    if threshold > 0.01:
        print("Check Disp ", all_tech_dispatch_shares, threshold)
    else:
        print("Ok shares Disp ", all_tech_dispatch_shares, threshold)
    threshold =  np.abs(np.array(all_tech_LCOS_shares).sum() - 100)
    if threshold > 0.01:
        print("Check LCOS ", all_tech_LCOS_shares, threshold)
    else:
        print("Ok shares LCOS ", all_tech_LCOS_shares, threshold)
    
    #print("results", results)
    
    my_title = r"Primary generation, delivered electricity, and LCOS values and shares at %.f gCO$_2$/kWh, curtailment %.f%%, CO$_2$ shadow price %.f €/t " % (np.uintc(np.round(CO2_intensity,0)), np.round(curtailment_pc,0), np.round(mu_CO2,0))
    if with_demand_H2_flag:
        # expand the title to include info on H2
        LCOEwS           = df_sys_indices.loc["value"]["LCOEwS"]
        LCOH             = df_sys_indices.loc["value"]["LCOH"]
        Tot_H2_delivered = df_sys_indices.loc["value"]["Tot_H2_delivered"]
        my_title += "\n"+ r"Delivered %.2f TWh/y of electricity at %.1f €/MWh" %  (np.round(Tot_el_delivered,2), np.round(LCOEwS,1))
        my_title +=  r" and %.2f TWh/y of H$_2$ at %.1f €/MWh = %.2f €/kg" %(np.round(Tot_H2_delivered,2), np.round(LCOH,1), np.round(LCOH/33.333,2)) 
    
    my_x_label = "(%)"
    fig = my_stacked_horizontal_bar_plot(results, all_tech_nice_labels, my_title, my_x_label)
    #plt.show()
    
    investment_results = {
        
        "Investment" : [*np.array(df_gen_results.investment[:]), *np.array(df_storage_results.investment[:])]
        
    }
    investment_labels = [*df_gen_results.index, *df_storage_results.index]
    nice_inv_labels = list(map(lambda x : dict_all_tech_nice_labels[x], investment_labels))
    i_g = np.round(df_gen_results.investment.sum(),1)
    i_s = np.round(df_storage_results.investment.sum(),1)
    i_t = i_g + i_s
    my_title = r"Total investment %.1f mld€, of which %.1f for generation and %.1f for storage" % (i_t, i_g, i_s)
    my_x_label = "G€"
    fig2 = my_stacked_horizontal_bar_plot(investment_results, nice_inv_labels, my_title, my_x_label)
    #plt.show()
     
    
    
    return fig, fig2, df_gen_results, df_sys_indices, df_storage_results



# List with the "nice" id are used only for charts
# Ensure that the order matchs with the list of keys used to access network components

generator_nice_labels = ['solar utility', 
                         'solar rooftop',
                         'onshore wind', 
                          'nuclear',
                        #  'nuclear B',
                          'OCGT',
                          'CCGT',
                          'fuel cell',
                         ]
            
non_fossil_gen_labels = ["solar", 'solar_rooftop', "onshorewind", 'nuclear']
non_fossil_gen_CF = [CF_solar, CF_solar2, CF_wind, CF_nuclear]
fossil_gen_labels = None
links_gen_labels = ["OCGT", "CCGT", 'FC']
primary_fuel_label = ["from_fossilCH4", "from_bioCH4", "from_CO2"]  # These are links id related to primary  fuels, not necessarily fossil, e.g. biogas
link_to_general_gas_labels = ["from_H2_to_gas", "from_CH4_to_gas"]
bus_general_gas_label = "gas"

# To be done
feedstock_labels = ["from_CO2"]

primary_fuel_nice_labels = [r"fossil CH$_4$", r"bio-CH$_4$", r"CO$_2$ feedstock"]

indirect_fuel_label = ["SynCH4_out", 'compH2_out'] # These are links that provide syn fuels

all_storage_tech_list = ["SD", "LD1", "PHS", "H2_in", "compH2", "SynCH4", "dac", "heat_pump_mid", "CO2_storage"] # these are all sub.systems that are necessary for storage
# Storage techs are either represented as two links and store or as a single link
# The subset of storage tech representedas single link are:
storage_single_link_list = ["H2_in", "dac", "heat_pump_mid"]
# The subset of storage techs with a bus2 input
storage_with_bus2 = ["compH2", "SynCH4", "dac"]


all_nice_storage_tech_list = ["Li-ion", "ACAES", "PHS", r"H$_2$, electrolysis", r"H$_2$, comp. & storage", r"e-CH$_4$","DAC","heat pump",r"CO$_2$ storage"]
list_syn_fuels = ["H2", "SynCH4"]

all_tech_nice_labels = generator_nice_labels[:] + all_nice_storage_tech_list[:] + primary_fuel_nice_labels[:]

dict_all_tech_nice_labels = {
    "solar"          :"solar utility",
    "solar_rooftop"  :"solar rooftop",
    "onshorewind"    :"onshore wind", 
    'nuclear'        :"nuclear",
    "OCGT"           :"OCGT", 
    "CCGT"           :"CCGT", 
    'FC'             :"fuel cell",
    "SD"             :"Li-ion", 
    "LD1"            :"CAES", 
    "PHS"            :"PHS", 
    "H2_in"          :r"H$_2$, electrolysis", 
    "compH2"         :r"H$_2$, comp. & storage",
    "SynCH4"         :r"e-CH$_4$",
    "from_fossilCH4" :r"fossil CH$_4$", 
    "from_bioCH4"    :r"bio-CH$_4$", 
    "from_CO2"       :r"CO$_2$ feedstock",
    "dac"            : "DAC",
    "CO2_storage"    :r"CO$_2$ storage",
    "heat_pump_mid"  : "heat pump",
}

fig, fig2, df_gen_results, df_sys_indices, df_storage_results = make_some_nice_charts_statistics(network=network,
                                                                                    links_gen_labels = links_gen_labels,
                                                                                    primary_fuel_label = primary_fuel_label,
                                                                                    indirect_fuel_label = indirect_fuel_label,
                                                                                    all_storage_tech_list=all_storage_tech_list,
                                                                                    storage_single_link_list = storage_single_link_list,
                                                                                    list_syn_fuels = list_syn_fuels,
                                                                                    storage_with_bus2 = storage_with_bus2,
                                                                                    link_to_general_gas_labels = link_to_general_gas_labels,
                                                                                    bus_general_gas_label = bus_general_gas_label,
                                                                                    )

fig.savefig(path_file_run_results + name_file_run_results +'_run00.png', dpi=300, bbox_inches = "tight")
fig2.savefig(path_file_run_results + name_file_run_results +'_investment.png', dpi=300, bbox_inches = "tight")



# In[35]:


fig = plot_generators_and_storage(df_storage_results, df_gen_results)
this_name_file = path_file_run_results + name_file_run_results +'_gen_capacities.png'
fig.savefig(this_name_file, dpi=300, bbox_inches = "tight")


fig = cost_vs_price_chart(df_storage_results, df_gen_results)
this_name_file = path_file_run_results + name_file_run_results +'_cost_price.png'
fig.savefig(this_name_file, dpi=300, bbox_inches = "tight")

fig = cost_vs_revenue_chart(df_storage_results, df_gen_results)
this_name_file = path_file_run_results + name_file_run_results +'_cost_revenue.png'
fig.savefig(this_name_file, dpi=300, bbox_inches = "tight")


# In[36]:


pd.set_option('display.max_columns', None)
print(df_gen_results)
print(df_storage_results)
#print(df_fuel_primary_results)



df_gen_results.loc[:]["generator_primary_CF"]





# Save results to an Excel file
this_name_file = path_file_run_results + name_file_run_results +'.xlsx'

with pd.ExcelWriter(this_name_file,engine='xlsxwriter') as writer:  
    df_gen_results.to_excel(writer, sheet_name='df_gen_results')
    df_storage_results.to_excel(writer, sheet_name='df_storage_results')



def plot_all_storage_SOC_log(network, all_storage_tech_list, storage_single_link_list, df_storage_results):
    import matplotlib.dates as mdates

    
    fig,ax = plt.subplots(1,1,figsize=(14,5))
    plt.yscale("log")  
    ax.set_xlabel("Hour of the year")
    ax.set_ylabel(r"MWh or tCO$_2$")
    plt.title("State of Charge  \n in the legend the storage rotation rate (number of full cycles per year)")
    for my_l in all_storage_tech_list:
        if my_l not in storage_single_link_list:
            nice_name = dict_all_tech_nice_labels[my_l]
            label_name = nice_name + " (" + str(np.round(df_storage_results.loc[my_l]["storage_full_cycles_eq_per_period"],1)) + ")"
            my_l_store = my_l + "_store"
            if network.stores_t.e[my_l_store].sum() > 0.1:
                plt.plot(network.stores_t.e[my_l_store], label=label_name)
            #
            
    plt.tight_layout()
    # Define the date format
    date_form = DateFormatter("%b")
    ax.xaxis.set_major_locator(plt.MaxNLocator(6))
    ax.xaxis.set_major_formatter(date_form)
    ax.set_xlim(left=network.snapshots.min(), right=network.snapshots.max())
    
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.25,
                 box.width, box.height * 0.65])

    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=False, ncol=5)
    
    #ax.legend()
    
    #plt.show()
    return fig



fig = plot_all_storage_SOC_log(network, all_storage_tech_list, storage_single_link_list, df_storage_results)
this_name_file = path_file_run_results + name_file_run_results +'_storage_SOC.png'
fig.savefig(this_name_file, dpi=300, bbox_inches = "tight")



network.links.loc["SynCH4_in", "p_nom_opt"]



fig,ax = plt.subplots(1,1,figsize=(12,4))
#plt.yscale("log")  
ax2 = ax.twinx()
ax.set_xlabel("Hour of the year")
ax.set_ylabel(r"MW of H$_2$")
ax.set_ylim(0, network.links.loc["SynCH4_in", "p_nom_opt"] * 1.1)
ax2.set_ylim(0, network.links_t.p2["SynCH4_in"].max() * 1.1)

ax2.set_ylabel(r"ton CO$_2$")
CF_SynCH4 = network.links_t.p0["SynCH4_in"].sum() /network.links.loc["SynCH4_in", "p_nom_opt"]/Num_hours
ax.plot(network.links_t.p0["SynCH4_in"], color="g", ls=":")
ax2.plot(network.links_t.p2["SynCH4_in"], color="y", ls="-.")
E_H2_in = network.links_t.p0["SynCH4_in"].sum()/1.e+6
T_CO2_in = network.links_t.p2["SynCH4_in"].sum()/1.e+6
CH4_out = - network.links_t.p1["SynCH4_in"].sum()/1.e+6
# Temp must be redone for a general case
E_compH2_out = - network.links_t.p1["compH2_out"].sum()/1.e+6
E_electrolysis = - network.links_t.p1["H2_in"].sum()/1.e+6
share_comp = E_compH2_out / E_electrolysis * 100.
plt.title(r"Methanation, H$_2$ input (left y-axis), CO$_2$ (right y-axis), CF %.2f, input H$_2$ %.2f (TWh/y) of which comp. %.2f (%%), plus CO$_2$ %.2f (Mt/y), output CH$_4$ %.2f (TWh/y)" % (CF_SynCH4, E_H2_in, share_comp, T_CO2_in, CH4_out))
print(share_comp)
plt.tight_layout()
# Define the date format
date_form = DateFormatter("%b")
ax.xaxis.set_major_locator(plt.MaxNLocator(6))
ax.xaxis.set_major_formatter(date_form)
ax.set_xlim(left=network.snapshots.min(), right=network.snapshots.max())
ax.legend()
#plt.show()
this_name_file = path_file_run_results + name_file_run_results +'_Methanation.png'
fig.savefig(this_name_file, dpi=300, bbox_inches = "tight")


fig,ax = plt.subplots(1,1,figsize=(12,4))
#plt.yscale("log")  
ax.set_xlabel("Day of the year")
ax.set_ylabel("MW of H2")
TWh_compH2 = network.links_t.p0["compH2_in"].sum()/1.e+6
GWh_used = network.links_t.p2["compH2_in"].sum() /1.e+3
this_tech_aux_cost  = network.links_t.p2["compH2_in"].multiply(network.snapshot_weightings.generators, axis=0).multiply(network.buses_t.marginal_price["AC"], axis=0).sum()/1.e+6 # in M€

plt.plot(network.links_t.p2["compH2_in"])

plt.title("H2 input at the compressor %.2f TWh/y, el. used %.2f GWh/y with a cost of %.2f M€/y" % (TWh_compH2, GWh_used, this_tech_aux_cost) )

plt.tight_layout()
# Define the date format
date_form = DateFormatter("%b")
ax.xaxis.set_major_locator(plt.MaxNLocator(6))
ax.xaxis.set_major_formatter(date_form)
ax.set_xlim(left=network.snapshots.min(), right=network.snapshots.max())
ax.legend()
#plt.show()
this_name_file = path_file_run_results + name_file_run_results +'_H2compressor.png'
fig.savefig(this_name_file, dpi=300, bbox_inches = "tight")





T_CO2_source = network.links_t.p0["CO2_storage_out"].sum()/1.e+6
T_CO2_in = network.links_t.p2["SynCH4_in"].sum()/1.e+6
print(T_CO2_source, T_CO2_in)



my_prices = np.array(network.buses_t.marginal_price["AC"].values)
my_prices.sort()

fig, ax = plt.subplots(1, 1)
ax.semilogy(my_prices)
ax.set_ylabel("€/MWh")
ax.set_xlabel("Hours per year")
ax.set_xlim(left=0, right=Num_hours)
plt.tight_layout()
plt.title("Price duration curve, electricity")
#plt.show()


this_name_file = path_file_run_results + name_file_run_results +'_price_duration.png'
fig.savefig(this_name_file, dpi=300, bbox_inches = "tight")



my_bus = "CO2"
network.buses_t.marginal_price[my_bus].plot()
network.buses_t.marginal_price[my_bus].mean()



network.buses_t.marginal_price["H2"].mean()




network.statistics.capacity_factor().loc[:][:]



#network.statistics.dispatch()


network.statistics.curtailment()


network.statistics.withdrawal()


network.statistics.optimal_capacity()


network.statistics.supply()



network.statistics.capex()

network.statistics.energy_balance()


current_CO2 = df_sys_indices.loc["value"]["CO2_intensity"] * df_sys_indices.loc["value"]["Tot_elec_delivered"] * 1.e+3 # ton CO2
iterate_on_CO2 = False
num_of_runs = 1

if current_CO2 > epsilon:
    iterate_on_CO2 = True
    list_CO2_limits = np.flip(np.linspace(0, current_CO2 * 0.6, num=4))
    print("list_CO2_limits", list_CO2_limits)
    num_of_runs = len(list_CO2_limits) + 1

all_runs_LCOS = np.zeros(num_of_runs)
all_runs_CO2 = np.zeros(num_of_runs)
all_runs_mu_CO2 = np.zeros(num_of_runs)

all_runs_LCOS[0]    = df_sys_indices.loc["value"]["LCOS"]
all_runs_CO2[0]     = df_sys_indices.loc["value"]["CO2_intensity"]
all_runs_mu_CO2[0]  = df_sys_indices.loc["value"]["mu_CO2"]    

print("from current_CO2", current_CO2)
if not do_parametric_runs:
    iterate_on_CO2 = False






if iterate_on_CO2:
    for i, co2_limit in enumerate(list_CO2_limits):
    
        network.global_constraints.at["co2_limit", "constant"] = co2_limit
        print(network.global_constraints.at["co2_limit", "constant"])


        network.optimize(snapshots=network.snapshots, pyomo=False, solver_name=name_solver, solver_options=solver_options, keep_references = True, keep_shadowprices = True)
    
        fig, fig2, df_gen_results, df_sys_indices, df_storage_results = make_some_nice_charts_statistics(network=network,
                                                                                              primary_fuel_label = primary_fuel_label,
                                                                                              indirect_fuel_label = indirect_fuel_label,
                                                                                              links_gen_labels = links_gen_labels,
                                                                                              all_storage_tech_list=all_storage_tech_list,
                                                                                              storage_single_link_list = storage_single_link_list,
                                                                                              list_syn_fuels = list_syn_fuels,
                                                                                              storage_with_bus2 = storage_with_bus2,
                                                                                              link_to_general_gas_labels = link_to_general_gas_labels,
                                                                                              bus_general_gas_label = bus_general_gas_label,
                                                                                                 )

        all_runs_LCOS[i + 1] = df_sys_indices.loc["value"]["LCOS"]
        all_runs_CO2[i + 1]  = df_sys_indices.loc["value"]["CO2_intensity"]
        all_runs_mu_CO2[i + 1]  = df_sys_indices.loc["value"]["mu_CO2"]
        this_run_name_file = path_file_run_results + name_file_run_results +'_run{:02}.png'.format(i + 1)
        fig.savefig(this_run_name_file, dpi=300, bbox_inches = "tight")



# Save some results
df = pd.DataFrame({'CO2':all_runs_CO2, 
                   'LCOS':all_runs_LCOS,
                   'mu_CO2': all_runs_mu_CO2
                  })



path_and_name_file = path_file_run_results + name_file_run_results +'.csv'
import os  
os.makedirs(path_file_run_results, exist_ok=True)  
df.to_csv(path_and_name_file, index=False) 


fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax2 = ax.twinx()

x = df['CO2']
y = df['LCOS']
y2 = df['mu_CO2']

l1, = ax.plot(x,y, c="y", label="LCOS")
l2, = ax2.plot(x,y2,  c="b", ls='-.', label=r"$\mu$ CO$_2$")

ax2.set_yscale("log")
ax2.legend(handles=[l1, l2])


ax.set_ylabel("LCOS (€/MWh)")
ax2.set_ylabel(r"$\mu_{CO2}$ (€/t)")
ax.set_xlabel(r"CO$_2$ (kg/MWh)")

plt.tight_layout()

#plt.show()
print(y)



import subprocess, os, platform
from glob import glob

if platform.system() == 'Darwin':       # macOS
    all_files_path = path_file_run_results + name_file_run_results +'*.png'
    list_files = glob(all_files_path)
    subprocess.call(('open', *list_files))



tech_parameters.at['onwind', 'marginal_cost_no_fuel']


tech_parameters.at['solar-utility', 'marginal_cost_no_fuel']


