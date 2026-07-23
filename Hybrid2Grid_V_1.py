'''
-------------------------------------------------------------------------------
Hybrid power system modelling
Copyright (C) 2026 James Manwell
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
See <http://www.gnu.org/licenses/>.
Contact: James Manwell, University of Massachusetts Amherst
Email: manwell@umass.edu
-------------------------------------------------------------------------------
'''

"""
This models a hybrid power system with ideal storage
It includes a minimum power level for the conventional generator and allows the generator
to be required to run all the time (or not).  It also includes a linear fuel use estimator
multiple non-identical generators are possible
Default assumption is for power units of MW and storage in MWh
(outputs are in GW or TWh for convenience) 
This is test version, based on <hybrid_with_storage.py> with 
procedures from <fuel_multi_gen.py> and <dispatch_test.py> incorporated in it
the original version was 6/3/25; it has been considerably upgraded
This version (based on 11) includes solar PV ana is derived from version 10, 
the primary difference being cleaning up extraneous comments and 
improving the formatting
Results match those of xlsm file from OSES 2023
Further testing and improvements to the dispatch model (now imported as dispatch_test_10)
would still be appropriate
This version was used in OSES 2026
"""
# the first block of imports are from python libraries

# these imports are py files that accompany hybrid2grid

import time
import csv
from datetime import datetime
from scipy.optimize import root_scalar
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import yaml
import solar_functions as SolarConvert
import dispatch_test_10 as dispatchFileName
pd.set_option('display.max_columns', None)
# Set expand_frame_repr to False to force all columns onto one line
# (if the output area width allows it, otherwise it might still wrap)
pd.set_option('display.expand_frame_repr', False)


# YNshutoff = 'N'   # Y means that all generators may be shutoff unless overwritten
# code_on = 'xxx'  # this is likely no longer needed, used for testing

# Get the file with the Hybrid2 parameters
this_dir = os.path.dirname(__file__)
parameter_file = os.path.join(this_dir, 'Hybrid2_parameters_OSES_4a.yaml')

# get inputs from yaml parameter file
inputs = yaml.safe_load(open(parameter_file))

parameter_file = os.path.basename(parameter_file)  # shorten for print out

# Read in the PV panel parameters
panel_parameters = inputs['panel_params']
k_boltz = panel_parameters['k']             # Boltzmanns constant, J/K
e_charge = panel_parameters['e']            # Electron energy, Coulombs
GT_ref = panel_parameters['GT_ref']         # Reference insolation, W/m^2
TC_ref = panel_parameters['TC_ref']         # Reference temperature, C
N_cells = panel_parameters['N_cells']       # NUmber of cells, -
Isc_ref = panel_parameters['Isc']           # Short circuit current, A
IL_ref = Isc_ref
V_oc_ref = panel_parameters['V_oc']         # Open circuit voltage, V
I_mp = panel_parameters['I_mp']             # Maximum power point current, A
V_mp = panel_parameters['V_mp']             # Maximum power point voltage, V
m = panel_parameters['m']                   # Ideality factor, -

Vt_ref = k_boltz*(TC_ref+273.15)/e_charge  # referebce terminal voltage
I0_ref = Isc_ref*np.exp(-V_oc_ref/(N_cells*m*Vt_ref))  # reference I0
Rs = (m*Vt_ref*(np.log((IL_ref-I_mp-I0_ref)
                       / I0_ref))-V_mp/N_cells)/I_mp  # reference resistance, ohms

# Short circuit current temperature coefficient, A/C
Isc_coeff = panel_parameters['Isc_coeff']

# Open circuit voltage temperature coefficient, V/C
Voc_coeff = panel_parameters['Voc_coeff']

# Panel heat loss coefficient, W/m^2 C
PV_UL = panel_parameters['PV_UL']
PV_Length = panel_parameters['PV_Length']   # PV panel length, m
PV_Width = panel_parameters['PV_Width']     # PV panel width, m
PV_panel_rated = V_mp*I_mp                  # rated PV panel power

# Read in the wind turbine parameters
turbine_parameters = inputs['turbine_params']

U_in = turbine_parameters['U_in']           # cut-in wind speed, m/s

U_rated = turbine_parameters['U_rated']     # rated wind speed, m/s

U_out = turbine_parameters['U_out']         # cut-out wind speed, m/s

P_rated = turbine_parameters['P_rated']     # rated power, MW
# turbine power curve
turbineFileName = turbine_parameters['powerCurveFileName']

# Read in the wind plant parameters

plant_params = inputs['plant_params']       # header for wind plant parameters
windFileName = plant_params['windFile']     # wind data file name
nTurbines = plant_params['n_turbines']      # number of wind turbines

# Read in solar site parameters

solar_site = inputs['solar_site_params']
# code for whether solar PV is considered
YNSolar = solar_site['YNSolar']
latitude = solar_site['latitude']           # latitude, degrees
slope = solar_site['slope']                 # panel slope, degrees
azimuth_w = solar_site['azimuth_w']         # panel azimuth angle, degrees
reflec = solar_site['reflec']               # surface reflectance, -
day_1 = solar_site['first_day']             # first day of solar data
first_hour = solar_site['first_hour']       # first hour of solar data
PV_system_rated = solar_site['PV_rated']    # rated PV system power, MW
solarFileName = solar_site['solar_data']    # solar data file name
# temperature data file name
temperatureFileName = solar_site['temperature_data']


# Read in load parameters
load_params = inputs['load_params']         # header for load parameters
loadFileName = load_params['loadFile']      # load data file name

# Get generator parameters
gen_params = inputs['gen_params']

# generator_dispatch_code = gen_params['gen_dispatch'] # removed 5/31/26
generator_dispatch_code = 1                 # default  added 5/31/26

if generator_dispatch_code != 1:

    '''# This block is no longer used as of 5/31/26
    '''

else:
    # generator parameters are from separate csv file
    # get generator parameter file
    generator_file = gen_params['gen_file']  # generator parameter file
    # generator_file = gen_params['gen_parameters'] #"generators_test_4.csv"
    file_path = generator_file
    _, ext = os.path.splitext(file_path.lower())

    # data starts on 2nd row
    data = np.loadtxt(generator_file, delimiter=",", skiprows=1)

    # get number of generators from generator paramater file
    nGen = data.shape[0]

    # code for generators that must be on
    GenOn = np.zeros(nGen)
    count = np.zeros(nGen)                    # counter for generator run time

    # code for generators that must be on, genertors only
    GenOn_gen_only = np.zeros(nGen)
    # counter for generator run time, generators only
    count_gen_only = np.zeros(nGen)

    NConfig = 2**nGen  # nGen**2                         # number of configurations
    ncols = nGen+8
    GenConfig = np.zeros((NConfig, ncols))
    genParams = np.array(data)      # generator paramaters from input file

    # initialize generator parameters, names self evident
    GenRatedPower = np.zeros(nGen)
    GenMinPower = np.zeros(nGen)
    GeneratorFullLoadFuel = np.zeros(nGen)
    GeneratorNoLoadFuel = np.zeros(nGen)
    GenMinFuel = np.zeros(nGen)
    GenFuelCost = np.zeros(nGen)
    GenMinRunTime = np.zeros(nGen)

    # code for whether all generators can be off
    YNshutoff = gen_params['YNshutoff']

    for i in range(nGen):
        GenRatedPower[i] = genParams[i, 0]
        GenMinPower[i] = genParams[i, 1]
        GeneratorFullLoadFuel[i] = genParams[i, 2]
        GeneratorNoLoadFuel[i] = genParams[i, 3]
        GenFuelCost[i] = genParams[i, 4]
        GenMinRunTime[i] = genParams[i, 5]
        GenMinFuel[i] = GenMinPower[i] * (GeneratorFullLoadFuel[i]
                                          - GeneratorNoLoadFuel[i]) / GenRatedPower[i] + GeneratorNoLoadFuel[i]

    result = dispatchFileName.initialize(
        genParams, nGen, NConfig, GenFuelCost, GenMinRunTime)

    GenConfig = result[0]

    n_configs = result[1]  # number of configurations

print_x = 10000   # for testing only!!!

output_parameters = inputs['time_series']  # header for output time series
# code for printing time series
YNTimeSeries = output_parameters['YNTimeSeries']

output_file_csv = output_parameters['OutputFile']  # name of output csv file


def InitializeStorage(inputs):
    # Read in the storage parameters from the yaml file
    store_params = inputs['store_params']
    # maximum storage capacity, MWh
    S0 = float(store_params['StorageLevelMax'])
    # initial storage level, MWh
    S1 = store_params['StoreInitRatio']*S0
    # rated charger output power
    S2 = float(store_params['P_chg_rated'])
    # fixed loss ratio for charging (rectifier)
    S3 = float(store_params['fxd_loss_ratio_chg'])
    # variable loss ratio for charging (rectifier)
    S4 = float(store_params['variable_loss_ratio_chg'])
    # rated discharge output power
    S5 = float(store_params['P_dis_rated'])
    # fixed loss ratio for discharge (inverter)
    S6 = float(store_params['fxd_loss_ratio_dis'])
    # variable loss ratio for discharge (inverterer
    S7 = float(store_params['variable_loss_ratio_dis'])
    return S0, S1, S2, S3, S4, S5, S6, S7


def powerCurve(powerCurve):
    # This reads the wind turbine power curve data
    nPtsPC = 0                  # initalize power curve points
    U_Pc = []                   # wind speed vector
    P_Pc = []                   # power vector
    powerCurve.readlines(1)
    for line in powerCurve:
        fields = line.split(',')
        U_Pc.append(float(fields[0]))
        P_Pc.append(float(fields[1]))
        nPtsPC = nPtsPC + 1
    powerCurve.close()
    return nPtsPC, U_Pc, P_Pc


def getNdata(fileName):
    # This gets the number of data points in a file
    file = open(fileName, 'r')
    file.readlines(1)  # skip first line
    nData = 0
    for line in file:
        nData += 1
    file.close()
    return nData


def ReadFile(fileName, n):
    # This is for reading load and wind speed data files
    data = np.zeros(n)
    file = open(fileName, 'r')
    file.readlines(1)  # skip first line
    for i in range(n):
        data[i] = float(file.readline())
    file.close()
    return data


def WindPower(U, U_Pc, P_Pc):
    # this functions calculates the power from multiple wind turbines
    # given the wind speed

    if U < U_in:
        P_wind = 0
    elif U < U_rated:
        P_wind = np.interp(U, U_Pc, P_Pc)
    elif U <= U_out:
        P_wind = P_rated
    else:
        P_wind = 0
    P_wind = P_wind*nTurbines                              # total wind power, MW

    return P_wind


def eff(converter_type, P, P_rated, fxd_loss_ratio, variable_loss_ratio):
    """
    This finds the efficiency and power out of a conversion device
    P: input power for charging
    P_rated: rated output power
    fxd_loss_ratio: ratio of the fixed loss to the rated input power
    variable_loss_ratio: the ratio of the variable loss to the input power
    Note that if the input power is less than the fixed loss, the output is zero
    Also note that converter type = 1 is for charging,  = 2 for discharging
    and = 3 for when storage is near full
    """
    P_rated_in = P_rated/(1-(fxd_loss_ratio + variable_loss_ratio))

    if converter_type == 1:
        # This is for charging
        if P != 0:
            # efficency = 1-variable_loss_ratio - P_rated_in*fxd_loss_ratio/P_in
            loss = P*variable_loss_ratio + fxd_loss_ratio*P_rated_in
            P_out = P - loss
            efficency = P_out/P
        else:
            efficency = 0
            P_out = 0
            loss = P

        if efficency < 0:
            efficency = 0
            P_out = 0
            loss = P
        return P_out, loss

    elif converter_type == 2:
        # this is for discharging unless there is a storage constraint

        if P > 0:
            # efficency = 1-variable_loss_ratio - P_rated_in*fxd_loss_ratio/P_in
            loss = P*variable_loss_ratio + fxd_loss_ratio*P_rated_in
            P_from_store = P + loss
            efficency = P/P_from_store
        else:
            efficency = 0
            P_from_store = 0
            loss = P
        return P_from_store, loss

    else:
        # this is for discharging when there is a storage constraint
        loss = P*variable_loss_ratio + fxd_loss_ratio*P_rated_in
        P_from_store = P - loss

        if P_from_store <= 0:
            P_from_store = 0
            efficency = 0
        else:
            efficency = P_from_store/P

    return P_from_store, loss


def IdealStorage_2(E_required, E_StorageLevelMax, E_StorageLevel, P_chg_rated, fxd_loss_ratio_chg, variable_loss_ratio_chg, P_dis_rated, fxd_loss_ratio_dis, variable_loss_ratio_dis):
    # This is to model ideal storage
    # Charging/discharging are 100% efficient
    # Note that in charging E_required is negative
    E_P_store_in = 0                          # energy into storage, kWh
    E_P_store_out = 0                         # energy into storage, kWh
    E_not_met = 0                             # load not met, kWh
    E_excess = 0                              # excess energy, kWh
    P_loss_chg = 0                            # charging loss
    P_loss_dis = 0                            # discharging loss

    if E_required > 0:                      # try to take energy from storage
        # check loss
        dum = eff(2, E_required, P_dis_rated,
                  fxd_loss_ratio_dis, variable_loss_ratio_dis)
        loss_temp = dum[1]
        E_P_store_out_temp = dum[0]

        if E_P_store_out_temp < P_dis_rated:
            # inverter is large enough for required power
            # enough storage to supply required and cover loss
            if E_P_store_out_temp < E_StorageLevel:
                if i == print_x:
                    print("")
                # dum = eff(E_required,P_dis_rated,fxd_loss_ratio_dis, variable_loss_ratio_dis)
                # E_P_store_out = dum[0]
                P_loss_dis = loss_temp
                E_P_store_out = E_required  # this is useful power actually delivered from storage
                '''
                This was used for testing, can be ignored
                if i == print_x:
                    print("")
                    print("###", i)
                    print("E_required", E_required)
                    print("E_P_store_out, P_loss_dis",
                          E_P_store_out, P_loss_dis)
                '''

            else:
                # NOT enough storage to supply required and cover loss, take as much as possible
                dum = eff(3, E_StorageLevel, P_dis_rated,
                          fxd_loss_ratio_dis, variable_loss_ratio_dis)
                E_P_store_out = dum[0]
                P_loss_dis = dum[1]
                if E_P_store_out == 0:
                    P_loss_dis = 0

                # E_not_met = E_required - E_StorageLevel + P_loss_dis
                # power must come from generator(s)
                E_not_met = E_required - E_P_store_out
                '''
                This was used for testing, can be ignored
                if i == print_x:
                    print('')
                    print("i", i)
                    print("E_StorageLevel", E_StorageLevel)
                    print("@@@")
                    print("E_required", E_required)
                    print("E_P_store_out", E_P_store_out)
                    print('P_loss_dis', P_loss_dis)
                    print('E_not_met', E_not_met)
                '''
        else:
            dum = eff(3, P_dis_rated, P_dis_rated,
                      fxd_loss_ratio_dis, variable_loss_ratio_dis)
            E_P_store_out_temp = dum[0]
            loss_temp = dum[1]
            if E_P_store_out == 0:
                loss_temp = 0
            '''
            This was used for testing, can be ignored            
            if E_P_store_out_temp < .01:
                print("E_P_store_out_temp", E_P_store_out_temp)
                print('loss_temp', loss_temp)'''

            if E_P_store_out_temp < E_StorageLevel:

                # dum = eff(E_required,P_dis_rated,fxd_loss_ratio_dis, variable_loss_ratio_dis)
                # E_P_store_out = dum[0]
                P_loss_dis = loss_temp
                # this is useful power actually delivered from storage
                E_P_store_out = E_P_store_out_temp
                # power must come from generator(s)
                E_not_met = E_required - E_P_store_out
                '''
                This was used for testing, can be ignored                            
                if i == print_x:
                    print("")
                    print("###", i)
                    print("E_required", E_required)
                    print("E_P_store_out, P_loss_dis",
                          E_P_store_out, P_loss_dis)
                '''
            else:
                # NOT enough storage to supply required and cover loss, take as much as possible
                # dum = eff(3,P_dis_rated,P_dis_rated,fxd_loss_ratio_dis, variable_loss_ratio_dis)
                # E_P_store_out = dum[0]
                # P_loss_dis = dum[1]
                # E_not_met = E_required - E_StorageLevel + P_loss_dis
                # power must come from generator(s)
                E_not_met = E_required - E_P_store_out

                '''
                This was used for testing, can be ignored   
                if i == print_x:
                    print('')
                    print("i", i)
                    print("E_StorageLevel", E_StorageLevel)
                    print("@@@")
                    print("E_required", E_required)
                    print("E_P_store_out", E_P_store_out)
                    print('P_loss_dis', P_loss_dis)
                    print('E_not_met', E_not_met)
                    '''

            # NewE_StorageLevel = E_StorageLevel - E_P_store_out  # storage level is now lower
        NewE_StorageLevel = E_StorageLevel - E_P_store_out - \
            P_loss_dis  # storage level is now lower

    else:
        # try to put all available energy into storage
        if E_StorageLevel - E_required < E_StorageLevelMax:

            # 1st test of charging with losses and limit on charge capacity, ok so far
            if P_chg_rated > - E_required:

                # no limit due to charger
                dum = eff(1, -E_required, P_chg_rated,
                          fxd_loss_ratio_chg, variable_loss_ratio_chg)
                E_P_store_in = dum[0]
                # E_P_store_in = - E_required
                P_loss_chg = dum[1]
                NewE_StorageLevel = E_StorageLevel + E_P_store_in  # updated storage level
                E_excess = 0

            else:
                # put power in up to P_chg_rated
                dum = eff(1, P_chg_rated, P_chg_rated,
                          fxd_loss_ratio_chg, variable_loss_ratio_chg)
                E_P_store_in = dum[0]
                # E_P_store_in = - E_required
                P_loss_chg = dum[1]

                NewE_StorageLevel = E_StorageLevel + E_P_store_in  # updated storage level
                E_excess = - E_required - E_P_store_in - P_loss_chg
                """
                This was used for testing
                print('E_P_store_in',E_P_store_in)
                print("E_required",E_required)
                print('E_excess',E_excess)
                print('P_loss_chg',P_loss_chg)
                """

        else:  # storage cannot take all that is available
            E_excess = E_StorageLevel - E_StorageLevelMax - E_required
            E_P_store_in = -E_required - E_excess
            NewE_StorageLevel = E_StorageLevelMax

    return E_P_store_in, E_P_store_out, E_not_met, E_excess, NewE_StorageLevel, P_loss_chg, P_loss_dis

# Define function for fuel use


def Dispatch_multi_gens_only(n, P_min, P_net, GenConfig, n_configs, GenOn, count):
    # This is for fuel use multiple non-indentical generators
    # without wind, solar or storage
    fuel = 0.
    P_extra = 0.
    P_unMet = 0.
    extraMinRun = 0.
    P_gen = P_net
    # print('P_net zzz',P_net)
    # print('GenConfig[n_configs-1,n+1] zzz',GenConfig[n_configs-1,n+1])
    if P_gen > GenConfig[n_configs-1, n+1]:
        P_unMet = P_gen-GenConfig[n_configs-1, n+1]
        P_gen = GenConfig[n_configs-1, n+1]
    # print('P_net zzz',P_net)
    result = dispatchFileName.new_fuel_use(genParams, P_gen, n, GenConfig,
                                           GenOn, GenRatedPower,
                                           GenMinPower, GeneratorFullLoadFuel,
                                           GeneratorNoLoadFuel, GenFuelCost, GenMinRunTime, GenMinFuel)
    fuel = result[1]
    config = result[0]
    extraMinRun = result[4]

    # account for minimum run time
    for j in range(n):
        # print()
        # print('config zzz 1',config)
        # print('j,config,GenConfig[config,j] zzz',j,config,GenConfig[config,j])
        if GenConfig[config, j] == 1:
            # print('GenConfig[config,j] zzz',GenConfig[config,j])
            if GenMinRunTime[j] > count[j]:
                # generator must still be running
                count[j] += 1
                GenOn[j] = 1

            else:
                # reset minimum run counter
                count[j] = 0
                GenOn[j] = 0

    return fuel, P_gen, P_extra, P_unMet, extraMinRun, GenOn, count


def Dispatch_multi(YNshutoff, n, P_min, P_net, E_StorageLevelMax,
                   E_StorageLevel, GenConfig, n_configs, GenOn, count, P_chg_rated, fxd_loss_ratio_chg,
                   variable_loss_ratio_chg, P_dis_rated, fxd_loss_ratio_dis, variable_loss_ratio_dis):
    # This is storage/generator dispatch for multiple non-indentical generators
    fuel = 0.
    P_curtail = 0.
    p_excess = 0.
    P_extra = 0.
    P_loss_chg = 0.
    P_loss_dis = 0.
    P_unMet = 0.
    extraMinRun = 0.

    if YNshutoff == "GGG":                         # Generator can be shut off
        ''' This block is no longer used
        if P_net < 0:
            # net power is negative, store if possible, generator is off
            p_excess = - P_net
            # get values from IdealStorage_2
            storage = IdealStorage_2(-p_excess, E_StorageLevelMax, E_StorageLevel, P_chg_rated,
                                     fxd_loss_ratio_chg, variable_loss_ratio_chg, P_dis_rated, fxd_loss_ratio_dis,
                                     variable_loss_ratio_dis)
            P_store_in = storage[0]
            P_store_out = storage[1]
            P_curtail = storage[3]
            E_StorageLevel = storage[4]
            P_loss_chg = storage[5]
            P_loss_dis = storage[6]
            P_gen = 0                        # generator can sbe hut off
            fuel = 0.
            """
            # P_curtail_av = P_curtail_av + P_curtail

            #if i == print_x:
             #   print("1 Dispatch",i)
            #    print("P_curtail",P_curtail)"""

        else:  # net power > 0, use store if possible, generator for remaining

            # Check on amount of storage available
            if P_net < E_StorageLevel:
                # may be enough storage, so take from storage; generator can be shut off, subject to rated inverter power
                # get values from IdealStorage_2
                storage = IdealStorage_2(P_net, E_StorageLevelMax, E_StorageLevel, P_chg_rated,
                                         fxd_loss_ratio_chg, variable_loss_ratio_chg, P_dis_rated, fxd_loss_ratio_dis,
                                         variable_loss_ratio_dis)
                P_store_in = storage[0]
                P_store_out = storage[1]
                P_not_met = storage[2]
                P_curtail = storage[3]
                E_StorageLevel = storage[4]
                P_loss_chg = storage[5]
                P_loss_dis = storage[6]

                if P_not_met > 0:
                    # P_gen = P_min   # generator must be on and run at minimum allowed power why???
                    P_gen = P_not_met  # why not this???

                     change here
                    # fuel = FuelUse_2(P_gen)

                    # fuel = 0 # temporary!!!
                    if P_gen > 0:
                        temp = dispatchFileName.new_fuel_use(genParams, P_net, n, GenConfig, GenOn, GenRatedPower,
                                                             GenMinPower, GeneratorFullLoadFuel, GeneratorNoLoadFuel, GenFuelCost, GenMinRunTime, GenMinFuel)
                        fuel = temp[1]
                        extraMinRun = temp[4]

                    if P_gen > (P_net - P_store_out):
                        P_extra = P_gen - (P_net - P_store_out)

                else:
                    P_gen = 0                        # shut off generator
                    fuel = 0.
                    P_extra = 0
                # P_curtail_av = P_curtail_av + P_curtail

            else:
                # not enough in storage; at least one generator must be running
                # p_excess = P_min - P_net
                # get values from IdealStorage_2
                storage = IdealStorage_2(P_net,
                                         E_StorageLevelMax, E_StorageLevel, P_chg_rated,
                                         fxd_loss_ratio_chg, variable_loss_ratio_chg, P_dis_rated,
                                         fxd_loss_ratio_dis, variable_loss_ratio_dis)
                P_store_in = storage[0]
                P_store_out = storage[1]
                P_curtail = storage[3]
                E_StorageLevel = storage[4]
                P_loss_chg = storage[5]
                P_loss_dis = storage[6]

                 change here
                # P_gen = P_min   # generator must be on run at minimum allowed power why is this true???

                P_gen = P_net - P_store_out  # isn;t this more reasonable 3/13/26

                 change here
                # fuel = FuelUse_2(P_gen)
                # fuel = 0 # temporary!!!
                if P_net > 0:
                    temp = dispatchFileName.new_fuel_use(genParams, P_gen, n, GenConfig, GenOn, GenRatedPower,
                                                         GenMinPower, GeneratorFullLoadFuel, GeneratorNoLoadFuel, GenFuelCost, GenMinRunTime, GenMinFuel)
                    fuel = temp[1]
                    extraMinRun = temp[4]

                if P_gen > (P_net - P_store_out):
                    # P_extra is from minimum power level
                    P_extra = P_gen - (P_net - P_store_out)
                """print("")
                print(i)
                print('P_extra',P_extra)
                print('P_gen',P_gen)
                print('P_net',P_net)
                print("P_store_out",P_store_out)"""

                # There nay be extra power due to minimum power level; try to store it

                # P_curtail_av = P_curtail_av + P_curtail
                # if i== print_x:
                #    print("P_gen 2",P_gen)

                # fuel = fuel_temp[0]
                # P_extra = fuel_temp[1]
        '''
    else:
        # this block is now used

        if P_net < P_min:

            OK_off = 0  # no shutoff unless explicitly allowed

            if P_net <= 0 and YNshutoff == 'Y':
                # net load < 0, all generators can be shut off
                OK_off = 1
                for i in range(n):
                    if GenOn[i] != 0:
                        OK_off = 0

            if OK_off == 1:
                # shutoff allowed
                fuel = 0  # temporary!!!
                P_gen = 0
                p_excess = P_net  # extra power (negative), try to store it

            else:
                """shutoff is allowed find fuel use due to minimum power of one generator"""
                P_gen = P_min
                result = dispatchFileName.new_fuel_use(genParams, P_min, n, GenConfig,
                                                       GenOn, GenRatedPower, GenMinPower, GeneratorFullLoadFuel,
                                                       GeneratorNoLoadFuel, GenFuelCost, GenMinRunTime, GenMinFuel)

                fuel = result[1]
                config = result[0]
                extraMinRun = result[4]
                # extra power (negative), try to store it
                p_excess = P_net - P_min

            # get values from IdealStorage_2
            storage = IdealStorage_2(p_excess, E_StorageLevelMax, E_StorageLevel,
                                     P_chg_rated, fxd_loss_ratio_chg,
                                     variable_loss_ratio_chg, P_dis_rated, fxd_loss_ratio_dis,
                                     variable_loss_ratio_dis)
            P_store_in = storage[0]
            P_store_out = storage[1]
            P_curtail = storage[3]
            E_StorageLevel = storage[4]
            P_loss_chg = storage[5]
            P_loss_dis = storage[6]

            ''' ignore here
            # P_gen = P_min
            # P_curtail_av = P_curtail_av + P_curtail
            move black below 6/7/25
                fuel_temp = FuelUse_2(P_gen)
                fuel = fuel_temp[0]
                P_extra = fuel_temp[1]
            # print("P_extra 1",P_extra)'''

        else:

            """net load is greater than minimum generator power level, take what is possible"""

            if YNshutoff != 'Y':
                # shut off NOT equal to 1; at least one generator must be on
                P_required = P_net - P_min
            else:
                P_required = P_net

            # get values from IdealStorage_2
            storage = IdealStorage_2(P_required, E_StorageLevelMax, E_StorageLevel,
                                     P_chg_rated, fxd_loss_ratio_chg,
                                     variable_loss_ratio_chg, P_dis_rated,
                                     fxd_loss_ratio_dis, variable_loss_ratio_dis)
            P_store_in = storage[0]
            P_store_out = storage[1]

            if YNshutoff != 'Y':
                P_gen = storage[2] + P_min
            else:
                P_gen = storage[2]

            E_StorageLevel = storage[4]
            P_loss_chg = storage[5]
            P_loss_dis = storage[6]

            ''' change here'''
            if P_gen > 0:
                # print('zzz')
                # print('GenOn zzz',GenOn)
                # print('n_configs',n_configs)
                if P_gen > GenConfig[n_configs-1, n+1]:
                    P_unMet = P_gen-GenConfig[n_configs-1, n+1]
                    P_gen = GenConfig[n_configs-1, n+1]
                result = dispatchFileName.new_fuel_use(genParams, P_gen, n, GenConfig,
                                                       GenOn, GenRatedPower, GenMinPower, GeneratorFullLoadFuel,
                                                       GeneratorNoLoadFuel, GenFuelCost, GenMinRunTime, GenMinFuel)
                fuel = result[1]
                config = result[0]
                extraMinRun = result[4]

            # account for minimum run time
            for j in range(n):
                # print()
                # print('config zzz 1',config)
                # print('j,config,GenConfig[config,j] zzz',j,config,GenConfig[config,j])
                if P_gen > 0:
                    if GenConfig[config, j] == 1:
                        # print('GenConfig[config,j] zzz',GenConfig[config,j])
                        if GenMinRunTime[j] > count[j]:
                            # generator must still be running
                            count[j] += 1
                            GenOn[j] = 1

                        else:
                            # reset minimum run counter
                            count[j] = 0
                            GenOn[j] = 0

    return (fuel, P_gen, E_StorageLevel, P_store_in, P_store_out,
            P_curtail, P_extra, P_loss_chg, P_loss_dis,
            P_unMet, extraMinRun, GenOn, count)


'''
no longer used 5/31/26
def FuelUse_2(load):
    # 'load' here refers to power required from generator(s)
    fuel = 0.  # initialize fuel use
    P_extra = 0.  # initialize extra power, MW
    # number of generators that must be on
    i_gen = int(load/P_gen_rated) + 1

    if load > 0:
        # if i==print_x:
        #   print("fuel",i,load)
        # Find power from load following generator
        if i_gen == load/P_gen_rated + 1:
            # if i== print_x:
            #    print("dum")
            P_load_follow = load - (i_gen - 2)*P_gen_rated
        else:
            P_load_follow = load - (i_gen - 1)*P_gen_rated
        # if i== print_x:
        #    print('P_load_follow',P_load_follow)

        if P_load_follow == P_gen_rated:
            # load following generator will run at rated power
            fuel = (i_gen-1)*(F_no_load + F_slope*P_gen_rated)
          #  if i== print_x:
          #      print("1 FuelUse_2",i)

        elif P_load_follow > P_min:
            # load following generator can follow load, others are at rated power
            fuel = (i_gen-1)*(F_no_load + F_slope*P_gen_rated)\
                + (F_no_load + F_slope*P_load_follow)
            # if i== print_x:
            #    print('fuel',fuel)
            #    print("2 FuelUse_2")

        # elif P_load_follow < P_min and load > P_min:

        elif load >= P_min:
            # load > P_min so power is divided equally among all generatrs
            fuel = i_gen*F_no_load + F_slope*load
            # if i== print_x:
            #    print("3")
           # print('load/i_gen',load/i_gen)
           # print('F_slope*P_min',F_slope*P_min)

        else:
            fuel = i_gen*F_no_load + F_slope*P_min
            P_extra = P_min-load
            print("Not Needed???")

    # block below likely no longer needed
    else:
        print("Needed???")
        code_on = 'N'
        if code_on == "Y":
            # print("1",i)
            # one generator must be on
            fuel = i_gen*F_no_load + F_slope*P_min
            P_extra = P_min
        else:
            fuel = 0
            # print("2",i)

    return fuel
'''


def TCell(GT_, T_ambC):
    # This finds the cell temperature from the ambient temperature, the solar insolation
    # the heat loss coefficient and the assumed cell efficiency
    eff_cell = 0.15  # This assumes an average efficiency of 0.15
    T_cellC = T_ambC + GT_ * (1 - eff_cell) / PV_UL
    return T_cellC


def solve_pv_current(V, T_cellC, GT_, left_start, right_start, tolerance):
    """
    Uses a root-finding method (Brent's method) to find the current I
    such that pv_iv_zero() = 0 for given V, T_cellC, and GT_.

    Inputs:
        V: Voltage across the PV panel
        T_cellC: Cell temperature in Celsius
        GT_ : Solar irradiance (W/m^2)
        left_start: Start of current search range (e.g., 0)
        right_start: End of current search range (e.g., Isc max)
        tolerance: Convergence threshold

    Returns:
        I: Estimated current 
    """

    def f(I):
        # root when residual = 0 (not absolute value)
        delta_temp = T_cellC - TC_ref
        V_thermal = k_boltz * (T_cellC + 273.15) / e_charge

        IL = (GT_ / GT_ref) * (IL_ref + delta_temp * Isc_coeff)

        VOC_adj = V_oc_ref + delta_temp * Voc_coeff
        I0 = (GT_ / GT_ref) * (IL_ref + delta_temp * Isc_coeff) * np.exp(
            -VOC_adj / (N_cells * m * V_thermal)
        )

        exp_term = np.exp((V / N_cells + I * Rs) / (m * V_thermal))
        I_model = IL - I0 * (exp_term - 1)

        return I - I_model

    try:
        sol = root_scalar(
            f,
            bracket=[left_start, right_start],
            method='brentq',
            xtol=tolerance
        )
        return sol.root if sol.converged else 0.0
    except ValueError:
        # No sign change in bracket → fallback
        return 0.0


def pv_iv_zero(I, V, T_cellC, GT_):
    """
    Computes residual between guessed current and modeled current from PV diode equation
    considering temperature and irradiance effects.

    Inputs:
        I: Guessed current
        V: Voltage input
        T_cellC: Cell temperature in Celsius
        GT_: Solar irradiance (W/m^2)

    Returns:
         Absolute error between I and model current
    """
    delta_temp = T_cellC - TC_ref
    V_thermal = k_boltz * (T_cellC + 273.15) / e_charge

    IL = (GT_ / GT_ref) * (IL_ref + delta_temp * Isc_coeff)

    VOC_adj = V_oc_ref + delta_temp * Voc_coeff
    I0 = (GT_ / GT_ref) * (IL_ref + delta_temp * Isc_coeff) * np.exp(
        -VOC_adj / (N_cells * m * V_thermal)
    )

    exp_term = np.exp((V / N_cells + I * Rs) / (m * V_thermal))

    I_model = IL - I0 * (exp_term - 1)

    return abs(I - I_model)


def pv_max_power(T_cellC, GT_, tolerance):
    """
    Finds the maximum power (PMax) from a PV panel using ternary search
    across voltage range (from 0 to VOC-like max voltage).

    Inputs:
        T_cellC (float): Cell temperature in Celsius
        GT_ (float): Solar irradiance
        left_start (float): Start of voltage search range
        right_start (float): End of voltage search range
        tolerance (float): Convergence tolerance
        IL_ref (float): Reference light-generated current (used in TernarySearchPVIV)

    Returns:
        float: Maximum power (PMax)
    """
    left = 0            # left_start - starting voltage
    right = V_oc_ref    # right_start- ending voltage
    PMax = 0.0

    for i in range(100):
        left_third = left + (right - left) / 3
        right_third = right - (right - left) / 3

        I1 = solve_pv_current(left_third, T_cellC, GT_, 0, IL_ref, tolerance)
        P1 = I1 * left_third

        I2 = solve_pv_current(right_third, T_cellC, GT_, 0, IL_ref, tolerance)
        P2 = I2 * right_third

        if P1 < P2:
            left = left_third
        else:
            right = right_third

        if abs(right - left)/abs(right) < tolerance:
            # Use average voltage and last computed current for final power
            V_avg = (left + right) / 2
            I_avg = solve_pv_current(V_avg, T_cellC, GT_, 0, IL_ref, tolerance)
            PMax = V_avg * I_avg
            break
    else:
        # If loop didn't converge
        PMax = 0.0
        V_avg = 0
    return PMax, V_avg


'''
no longer used 5/31/26
def Dispatch(YNshutoff, P_min, P_net, E_StorageLevelMax, E_StorageLevel, P_chg_rated,
             fxd_loss_ratio_chg, variable_loss_ratio_chg, P_dis_rated, fxd_loss_ratio_dis, variable_loss_ratio_dis):
    fuel = 0.
    P_curtail = 0.
    p_excess = 0.
    P_extra = 0.
    P_loss_chg = 0.
    P_loss_dis = 0.
    if i == print_x:
        print('')

    if YNshutoff == "Y":                         # Generator can be shut off

        if P_net < 0:
            # net power is negative, store if possible, generator is off
            p_excess = - P_net
            # get values from IdealStorage_2
            storage = IdealStorage_2(-p_excess, E_StorageLevelMax, E_StorageLevel, P_chg_rated,
                                     fxd_loss_ratio_chg, variable_loss_ratio_chg, P_dis_rated,
                                     fxd_loss_ratio_dis, variable_loss_ratio_dis)
            P_store_in = storage[0]
            P_store_out = storage[1]
            P_curtail = storage[3]
            E_StorageLevel = storage[4]
            P_loss_chg = storage[5]
            P_loss_dis = storage[6]
            P_gen = 0                        # generator can sbe hut off
            fuel = 0.
            """
            # P_curtail_av = P_curtail_av + P_curtail
            
            #if i == print_x:
             #   print("1 Dispatch",i)
            #    print("P_curtail",P_curtail)"""

        else:  # net power > 0, use store if possible, generator for remaining

            # Check on amount of storage available
            if P_net < E_StorageLevel:
                # may be enough storage, so take from storage; generator can be shut off, subject to rated inverter power
                # get values from IdealStorage_2
                storage = IdealStorage_2(P_net,
                                         E_StorageLevelMax, E_StorageLevel, P_chg_rated,
                                         fxd_loss_ratio_chg, variable_loss_ratio_chg, P_dis_rated,
                                         fxd_loss_ratio_dis, variable_loss_ratio_dis)
                P_store_in = storage[0]
                P_store_out = storage[1]
                P_not_met = storage[2]
                P_curtail = storage[3]
                E_StorageLevel = storage[4]
                P_loss_chg = storage[5]
                P_loss_dis = storage[6]

                if P_not_met > 0:
                    P_gen = P_min   # generator must be on and run at minimum allowed power
                    fuel = FuelUse_2(P_gen)
                    if P_gen > (P_net - P_store_out):
                        P_extra = P_gen - (P_net - P_store_out)

                else:
                    P_gen = 0                        # shut off generator
                    fuel = 0.
                    P_extra = 0
                # P_curtail_av = P_curtail_av + P_curtail

            else:
                # not enough in storage; at least one generator must be running
                # p_excess = P_min - P_net
                # get values from IdealStorage_2
                storage = IdealStorage_2(P_net,
                                         E_StorageLevelMax, E_StorageLevel, P_chg_rated,
                                         fxd_loss_ratio_chg, variable_loss_ratio_chg, P_dis_rated,
                                         fxd_loss_ratio_dis, variable_loss_ratio_dis)
                P_store_in = storage[0]
                P_store_out = storage[1]
                P_curtail = storage[3]
                E_StorageLevel = storage[4]
                P_loss_chg = storage[5]
                P_loss_dis = storage[6]
                P_gen = P_min   # generator must be on run at minimum allowed power
                fuel = FuelUse_2(P_gen)
                if P_gen > (P_net - P_store_out):
                    # P_extra is from minimum power level
                    P_extra = P_gen - (P_net - P_store_out)
                """print("")
                print(i)
                print('P_extra',P_extra)
                print('P_gen',P_gen)
                print('P_net',P_net)
                print("P_store_out",P_store_out)"""

                # There nay be extra power due to minimum power level; try to store it

                # P_curtail_av = P_curtail_av + P_curtail
                # if i== print_x:
                #    print("P_gen 2",P_gen)

                # fuel = fuel_temp[0]
                # P_extra = fuel_temp[1]

    else:
        """no shutoff allowed, at least one generator must be on"""

        if P_net < P_min:
            # if i== print_x:
            #   print("4",i)
            """find fuel use due to minimum power of one generator"""
            fuel = FuelUse_2(P_min)
            # fuel = fuel_temp[0]

            p_excess = P_net - P_min  # extra power (negative), try to store it
            # get values from IdealStorage_2
            storage = IdealStorage_2(p_excess, E_StorageLevelMax, E_StorageLevel, P_chg_rated,
                                     fxd_loss_ratio_chg, variable_loss_ratio_chg, P_dis_rated, fxd_loss_ratio_dis,
                                     variable_loss_ratio_dis)
            P_store_in = storage[0]
            P_store_out = storage[1]
            P_curtail = storage[3]
            E_StorageLevel = storage[4]
            P_loss_chg = storage[5]
            P_loss_dis = storage[6]
            P_gen = P_min
            # P_curtail_av = P_curtail_av + P_curtail
            """ move black below 6/7/25
            fuel_temp = FuelUse_2(P_gen)
            fuel = fuel_temp[0]
            P_extra = fuel_temp[1]
            """
            # print("P_extra 1",P_extra)
        else:
            # if i== print_x:
            #    print("5",i)
            """net load is greater than minimum generator power level, take what is possible"""
            P_required = P_net - P_min
            # get values from IdealStorage_2
            storage = IdealStorage_2(P_required, E_StorageLevelMax, E_StorageLevel,
                                     P_chg_rated, fxd_loss_ratio_chg, variable_loss_ratio_chg,
                                     P_dis_rated, fxd_loss_ratio_dis, variable_loss_ratio_dis)
            P_store_in = storage[0]
            P_store_out = storage[1]
            P_gen = storage[2] + P_min
            E_StorageLevel = storage[4]
            P_loss_chg = storage[5]
            P_loss_dis = storage[6]
            fuel = FuelUse_2(P_gen)
            """fuel = fuel_temp[0]
            P_extra = fuel_temp[1]
            print("P_extra 2",P_extra)"""
    print('fuel', fuel)
    return fuel, P_gen, E_StorageLevel, P_store_in, \
        P_store_out, P_curtail, P_extra, P_loss_chg, P_loss_dis
'''


def resolve_output_filename(filename, overwrite=False, auto_rename=True):
    """
    Resolves output filename safely.

    overwrite=False   → raise error if file exists
    auto_rename=True  → append _1, _2, ... if file exists
    """

    if not os.path.exists(filename):
        return filename

    if overwrite:
        return filename

    if not auto_rename:
        raise FileExistsError(
            f"Output file '{filename}' already exists.\n"
            f"Set overwrite=True or auto_rename=True."
        )

    base, ext = os.path.splitext(filename)
    i = 1
    while True:
        new_name = f"{base}_{i}{ext}"
        if not os.path.exists(new_name):
            return new_name
        i += 1


def safe_write_time_series_csv(filename, data, header=None,
                               overwrite=False, auto_rename=True,
                               retries=5, delay=0.5):

    filename = resolve_output_filename(filename, overwrite, auto_rename)
    # print("data",data)
    for attempt in range(retries):
        try:
            with open(filename, "w", newline="") as f:
                if header:
                    for line in header.split("\n"):
                        f.write(f"# {line}\n")

                writer = csv.writer(f)

                # --- Handle column-wise data ---
                if isinstance(data, (list, tuple)) and len(data) > 0:
                    # Check if this is column-wise data
                    if hasattr(data[0], "__len__"):
                        columns = data
                        ncols = len(columns)
                        nrows = len(columns[0])

                        # Optional: check all columns same length
                        for col in columns:
                            if len(col) != nrows:
                                raise ValueError(
                                    "All columns must have same length")

                        # Header
                        col_names = [f"Value_{j+1}" for j in range(ncols)]
                        col_names = ['wind', 'load', 'net',
                                     'generator', 'store', 'fuel']
                        writer.writerow(["Index"] + col_names)

                        # Write rows by zipping columns
                        for i, row_vals in enumerate(zip(*columns), start=1):
                            writer.writerow([i] + list(row_vals))

                    else:
                        # Single column fallback
                        writer.writerow(["Index", "Value"])
                        for i, val in enumerate(data, start=1):
                            writer.writerow([i, val])

                else:
                    raise ValueError("Data format not recognized")

            print("")
            print(f"CSV time series written to: {filename}")
            return filename

        except PermissionError:
            if attempt == retries - 1:
                raise PermissionError(
                    f"Cannot write '{filename}'.\n"
                    f"File is open in another program (Excel, etc.)."
                )
            time.sleep(delay)


def display_results(Hour, P_wind, P_load, P_PV, P_gen, P_extraMinRun, P_curtail, store, P_net, fuel_wtg,
                    P_load_av, P_wtg_av, P_net_av, P_gen_av, P_curtail_av, P_extra_av,
                    P_loss_chg_av, P_loss_dis_av, P_extraRunTime_av, P_PV_av, P_unMet_av, Fuel_total,
                    Fuel_wtg_total, wtg_ratio, Fuel_save_ratio, balance_total, hours_storage,
                    E_StorageLevel, E_StorageLevelMax):
    # print output to screen
    print()
    print('Input files:')
    print('yaml file:', parameter_file)
    print('wind file:', windFileName)
    print('load file:', loadFileName)
    if YNSolar == 'Y':
        print('solar file:', solarFileName)
        print('temperature file:', temperatureFileName)
    print('turbine file:', turbineFileName)
    if generator_dispatch_code == 1:
        print('generator file:', generator_file)
    else:
        print('generator parameters from yaml file')
    print("")
    print(f"Average load = {P_load_av/1000:,.4} GW")
    print(f"Average wind power = {P_wtg_av/1000:,.4} GW")
    if YNSolar == 'Y':
        print(f"Average PV power = {P_PV_av/1000:,.3} GW")
    # P_wtg_av, P_net_av,P_gen_av,P_curtail_av,P_extra_av,P_loss_chg_av,P_loss_dis_av,P_extraRunTime_av,Fuel_total
    print(f"Average net load = {P_net_av/1000:,.3} GW")
    print("")

    if YNshutoff == 'Y':
        print('All generators can be shut off')
    else:
        print("At least one generator must be on")

    print('')
    print(f"Average conventional power = {P_gen_av/1000:.4} GW")
    if YNSolar == 'Y':
        print(f"Average curtailed wind/PV power = {P_curtail_av/1000:,.3} GW")
    else:
        print(f"Average curtailed wind power = {P_curtail_av/1000:,.3} GW")
    print(f"Average extra power = {P_extra_av/1000:,.3} GW")
    print(f"Average charging loss = {P_loss_chg_av/1000:,.3} GW")
    print(f"Average discharging loss = {P_loss_dis_av/1000:,.3} GW")
    print(f"Average unmet power = {P_unMet_av/1000:,.3} GW")
    print(f"Average extra min run time power = {
          P_extraRunTime_av/1000:,.3} GW")
    if YNSolar == 'Y':
        print(f"Fuel use without wind/PV = {Fuel_total:,.3} units/yr")
    else:
        print(f"Fuel use without wind = {Fuel_total:,.3} units/yr")
    if generator_dispatch_code == 1:
        if YNSolar == 'Y':
            print(f"Fuel use with wind/PV = {Fuel_wtg_total:,.3} units/yr")
            print(f"Fuel savings ratio = {Fuel_save_ratio:,.3}")
        else:
            print(f"Fuel use with wind = {Fuel_wtg_total:,.3} units/yr")
            print('')
            print(f"Wind power to load ratio = {wtg_ratio:,.3}")
            print(f"Fuel savings ratio = {Fuel_save_ratio:,.3}")

    else:
        if YNSolar == 'Y':
            print(f"Fuel use with wind /PV= {Fuel_wtg_total} units/yr")
            print(f"Fuel savings ratio = {Fuel_save_ratio:,.3}")
        else:
            print(f"Fuel use with wind = {Fuel_wtg_total} units/yr")
            print('')
            print(f"Wind power to load ratio = {wtg_ratio:,.3}")
            print(f"Fuel savings ratio = {Fuel_save_ratio:,.3}")
            print('')

    print('')
    print(f"balance overall = {balance_total:,.4}")
    print("")
    print(f"Maximum storage capacity = {E_StorageLevelMax/1e6:,.4} TWh")
    print(f"(hours of storage at average P_load = {hours_storage:,.4})")
    print(f"Final storage Level = {E_StorageLevel/1e6:,.3} TWh")

    # plot time series in GW
    plt.plot(Hour, P_wind/1000)
    plt.plot(Hour, P_load/1000)
    plt.legend(["Wind power", "Load"])
    plt.title("Load and wind power")
    plt.xlabel("Time, hrs")
    plt.ylabel("Power, GW")
    # plt.ylim(0,30)

    plt.figure(2)
    plt.plot(Hour, (P_gen+P_extraMinRun)/1000)
    plt.title("Conventional power")
    plt.xlabel("Time, hrs")
    plt.ylabel("Power, GW")
    # plt.ylim(0,30)

    plt.figure(3)
    plt.plot(Hour, P_curtail/1000)
    plt.title("Curtailed power")
    plt.xlabel("Time, hrs")
    plt.ylabel("Power, GW")
    # plt.ylim(0,30)

    plt.figure(4)
    plt.plot(Hour, store/1e6)
    plt.title("Stored energy")
    plt.xlabel("Time, hrs")
    plt.ylabel("Energy, TWh")
    # plt.ylim(0,30)

    if YNSolar == "Y":
        plt.figure(5)
        plt.plot(Hour, P_PV/1000)
        plt.plot(Hour, P_load/1000)
        plt.legend(["PV power", "Load"])
        plt.title("Load and PV power")
        plt.xlabel("Time, hrs")
        plt.ylabel("Power, GW")

    plt.show()

    if YNTimeSeries == 'Y':
        data = P_wind, P_load, P_net, P_gen, store, fuel_wtg

        summary = 'Simulation run: '+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'\n'
        summary += 'parameter file = ' + parameter_file + \
            '; generator file = ' + generator_file + '\n'+'load file = ' + loadFileName +\
            '; wind file = ' + windFileName
        if YNSolar == "Y":
            summary += '; solar file = '+solarFileName
        if YNshutoff == 'Y':
            summary += '; all generators may be shut off'
        else:
            summary += '; at least one generator must be on'
        safe_write_time_series_csv(
            output_file_csv,
            data,
            header=summary,
            overwrite=False,      # protect existing file
            auto_rename=True      # timeseries_1.csv, timeseries_2.csv, ...
        )

# Main loop


def main(YNshutoff, nGen, GenMinPower,
         GenConfig, n_configs, GenOn, GenOn_gen_only, count, count_gen_only):
    # Read in load data, get the number of data points
    nData1 = getNdata(loadFileName)

    # Read in wind data, get the number of data points
    nData2 = getNdata(windFileName)

    # use lower number of data points if different
    if nData1 < nData2:
        nData = nData1
    else:
        nData = nData2

    # nData = 100  # override for testing only!!!
    print()
    print(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print('number of simulation steps:', nData)

    # Initialize arrays
    # initialize load array
    P_load = np.zeros(nData)
    # initialize wind array
    U = np.zeros(nData)
    # initialize net load array
    P_net = np.zeros(nData)
    # initialize wind power array
    P_wind = np.zeros(nData)
    # initialize time (hours) array
    Hour = np.zeros(nData)
    # initialize curtailed generation array
    P_curtail = np.zeros(nData)
    # initialize conventional generators array
    P_gen = np.zeros(nData)
    # initialize storage level array (for plot)
    store = np.zeros(nData)
    P_extraMinRun = np.zeros(nData)
    fuel_wtg = np.zeros(nData)
    Hour = np.zeros(nData)
    I = np.zeros(nData)
    T = np.zeros(nData)
    P_PV = np.zeros(nData)

    # Initialize averages
    # initialize average wind turbine power
    P_wtg_av = 0
    # initialize average net load
    P_net_av = 0
    # initialize average conventional power, MW
    P_gen_av = 0
    # initialize average curtail power, MW
    P_curtail_av = 0
    # initialize fuel use w/out turbines
    Fuel_total = 0
    # initialize fuel use with turbines
    Fuel_wtg_total = 0
    P_extra_av = 0
    P_loss_chg_av = 0
    P_loss_dis_av = 0
    P_extraRunTime_av = 0
    P_PV_av = 0
    P_unMet_av = 0
    # total balance, should sum to 0
    balance_total = 0
    # dum = 0                                                     # dummy variable
    # initializing, for converting solar data
    hour_count = -0.5
    day = day_1

    # Read in load data
    P_load = ReadFile(loadFileName, nData)
    P_load_av = np.average(P_load)               # calculate average load,

    # Read in wind data
    U = ReadFile(windFileName, nData)

    # Get power curve
    turbineFile = open(turbineFileName, 'r')
    temp = powerCurve(turbineFile)
    nPoints = temp[0]
    U_Pc = np.zeros(nPoints)
    P_Pc = np.zeros(nPoints)
    U_Pc = temp[1]
    P_Pc = temp[2]

    # Read in solar site data
    if YNSolar == 'Y':
        I = ReadFile(solarFileName, nData)  # solar insolation data
        T = ReadFile(temperatureFileName, nData)  # temperature data

    # Initialize Storage
    temp = InitializeStorage(inputs)
    E_StorageLevelMax = temp[0]
    E_StorageLevel = temp[1]
    P_chg_rated = temp[2]
    fxd_loss_ratio_chg = temp[3]
    variable_loss_ratio_chg = temp[4]
    P_dis_rated = temp[5]
    fxd_loss_ratio_dis = temp[6]
    variable_loss_ratio_dis = temp[7]

    '''azimuth_w = 0 # test!
    slope = 40 # test!
    latitude = 40 #test!'''

    for i in range(nData):
        # print('i', i)  # !!!$$$
        hour_count += 1
        if hour_count > 23:
            hour_count += -24
            day += 1

        # hour (for plotting)
        Hour[i] = i
        # calculate power at U[i], MW
        P_wind[i] = WindPower(U[i], U_Pc, P_Pc)

        if YNSolar == 'Y':
            if I[i] > 0:
                I_slope = SolarConvert.SunOnSlope(
                    I[i], latitude, day, hour_count, slope, azimuth_w, reflec)
                T_panel = TCell(I_slope[0], T[i])  # fixed 5/31/26
                # print('*************')
                # print('I_slope,T_panel',I_slope,T_panel)
                # I_slope = 726 # test!
                # T_panel = 56 # test!

                temp = pv_max_power(T_panel, I_slope[0], .01)  # fixed 5/31/26

                P_PV[i] = temp[0]*PV_system_rated/PV_panel_rated
                '''if i<25:
                    print()
                    print('hour_count',hour_count)
                    print('T[i],T_panel',T[i],T_panel)
                    print('I[i]',I[i])
                    print('I_slope,P_PV[i]',I_slope,P_PV[i])'''
            else:
                P_PV[i] = 0

        # calculate net load, MW
        P_net[i] = P_load[i] - P_wind[i] - P_PV[i]
        P_gen[i] = 0.                                        # initialize generator power, MW
        P_curtail[i] = 0.                                    # initialize curtailed power, MW
        # P_required = 0.                                      # initialize required power, MW

        # for testing
        # if P_net[i]>0:
        #      fuel_wtg =0

        if generator_dispatch_code == 1:

            # if i == 20:
            #    dum = dum  # zzz for testing only
            outs = Dispatch_multi(YNshutoff, nGen, GenMinPower[1], P_net[i], E_StorageLevelMax,
                                  E_StorageLevel, GenConfig, n_configs, GenOn, count, P_chg_rated, fxd_loss_ratio_chg,
                                  variable_loss_ratio_chg, P_dis_rated, fxd_loss_ratio_dis, variable_loss_ratio_dis)
            # outs = Dispatch(YNshutoff, P_min, P_net[i], E_StorageLevelMax, E_StorageLevel)
            # print('P_net[i]', P_net[i])  # !!!
            count = outs[12]
            GenOn = outs[11]

            # if i < 5000:
            #    dum = dum  # zzz for testing only
            # print('count',count)
            outs_gen_only = Dispatch_multi_gens_only(nGen, GenMinPower[1], P_load[i], GenConfig,
                                                     n_configs, GenOn_gen_only, count_gen_only)
            fuel_gen_only = outs_gen_only[0]
            count_gen_only = outs_gen_only[6]
            GenOn_gen_only = outs_gen_only[5]

        else:
            '''# no longer used 5/31/26
            #outs = Dispatch(YNshutoff, P_min,
            #                P_net[i], E_StorageLevelMax, E_StorageLevel)'''

        # print('i',i)
        fuel_wtg[i] = outs[0]
        P_gen[i] = outs[1]
        E_StorageLevel = outs[2]
        P_store_in = outs[3]
        P_store_out = outs[4]
        P_curtail[i] = outs[5]
        P_extra = outs[6]
        P_loss_chg = outs[7]
        P_loss_dis = outs[8]

        if generator_dispatch_code == 1:
            P_unMet = outs[9]
        else:
            P_unMet = 0
        P_extraMinRun[i] = outs[10]

        if P_loss_chg > 0 and P_loss_dis > 0:
            print(i)

        balance = P_wind[i] + P_PV[i] - P_load[i] - P_store_in - P_curtail[i] \
            + P_store_out + P_gen[i] - P_loss_chg - P_extra + P_unMet

        balance_total = balance_total + balance
        # sum for average net load
        P_net_av = P_net_av + P_net[i]
        # sum for average wind turbine power
        P_wtg_av = P_wtg_av + P_wind[i]
        P_extra_av = P_extra_av + P_extra
        P_loss_chg_av = P_loss_chg_av + P_loss_chg
        P_loss_dis_av = P_loss_dis_av + P_loss_dis
        P_curtail_av = P_curtail_av + P_curtail[i]
        P_extraRunTime_av += P_extraMinRun[i]
        P_gen_av = P_gen_av + P_gen[i]
        Fuel_wtg_total += fuel_wtg[i]
        Fuel_total += fuel_gen_only
        P_extra_av = P_extra_av + P_extra
        P_PV_av += P_PV[i]
        P_unMet_av += P_unMet
        # fill storage array for plotting
        store[i] = E_StorageLevel

    # total fuel savings, units
    Fuel_save = Fuel_total - Fuel_wtg_total
    # average power from wind turbines, MW
    P_wtg_av = P_wtg_av/nData
    # average net P_load, MW
    P_net_av = P_net_av/nData
    # average generator power, MW
    P_gen_av = P_gen_av/nData
    # average curtailed power, MW
    P_curtail_av = P_curtail_av/nData
    P_extra_av = P_extra_av/nData
    P_loss_chg_av = P_loss_chg_av/nData
    P_loss_dis_av = P_loss_dis_av/nData
    P_extraRunTime_av = P_extraRunTime_av/nData
    P_PV_av = P_PV_av/nData
    P_unMet_av = P_unMet_av/nData
    '''if P_load_av>0:
        Net_load_ratio = P_net_av/np.average(P_load)        # net load to load ratio
    else:
        P_load_av'''  # not used now

    if Fuel_total > 0:
        Fuel_save_ratio = Fuel_save/Fuel_total              # fuel savings ratio
    else:
        Fuel_save_ratio = 0.

    if P_load_av > 0:
        wtg_ratio = P_wtg_av/P_load_av                     # wind power to load ratio
        hours_storage = E_StorageLevelMax/P_load_av
    else:
        hours_storage = 0.
        wtg_ratio = 0.

    display_results(Hour, P_wind, P_load, P_PV, P_gen, P_extraMinRun, P_curtail, store, P_net, fuel_wtg,
                    P_load_av, P_wtg_av, P_net_av, P_gen_av, P_curtail_av, P_extra_av,
                    P_loss_chg_av, P_loss_dis_av, P_extraRunTime_av, P_PV_av, P_unMet_av, Fuel_total,
                    Fuel_wtg_total, wtg_ratio, Fuel_save_ratio, balance_total, hours_storage, E_StorageLevel, E_StorageLevelMax)


if __name__ == "__main__":
    main(YNshutoff, nGen, GenMinPower,
         GenConfig, n_configs, GenOn, GenOn_gen_only, count, count_gen_only)
