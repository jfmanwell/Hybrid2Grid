'''
-------------------------------------------------------------------------------
Generator dispatch testing code, used as input to Hybrid2Grid.py
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
Created on Mon Mar  2 14:09:27 2026
This is for testing multiple generators
The previous version _3 seemed to work pretty well
This version (_4) is to add minimum run time (not there yet)
much seems to work, next task is to add fuel use and verify, being
careful to distinguish single generator case from multi-generator
current version (_5) is to clean up and add fuel
Version _5 worked pretty well, now clean up and put more features in functions
to make merging with hybrid2_test_2.py easier
Version _6 is much better organized than _5, now continue clean up and get rid of old code
that is probably useless; this is saved in version _6 in case it is useful
The current version is #10 and is subject to update as necessary
"""
'''
pd.set_option('display.max_columns', None)
# Set expand_frame_repr to False to force all columns onto one line
# (if the output area width allows it, otherwise it might still wrap)
pd.set_option('display.expand_frame_repr', False)

input_file = "generators_test_3.csv"

file_path = input_file
_, ext = os.path.splitext(file_path.lower())

data = np.loadtxt(input_file, delimiter=",",skiprows=1) # data starts on 2nd row

nGen = data.shape[0]       # get number of generators from generator paramater file
NConfig=nGen**2            # number of configurations
genParams = np.array(data)      # generator paramaters from input file
LocationVector = np.zeros(nGen)'''

# genParams[i,0]        # rated power
# genParams[i,1]        # minimum power
# genParams[i,2]        # GeneratorFullLoadFuel
# genParams[i,3]        # full load fuel
# genParams[i,4]        # no load fuel
# genParams[i,5]        # minimum run time

'''
for i in range(0,nGen):
   LocationVector[i] = i  # location code for generator
GenOn = np.zeros(nGen)  # test!

for i in range(0,nGen):
   LocationVector[i] = i  # location code for generator
GenOn = np.zeros(nGen)  # test!
'''




import os
import numpy as np
import pandas as pd
import sys
def Swap(x, y):
    # This swaps x and y
    tmp = x
    x = y
    y = tmp
    return x, y


def Sort2D3(A_matrix, LocationVector, LocationIndex, N_rows, N_columns, OrderCode):
    '''This is a bubble sorting routine
    A_matrix() is returned reordered, from smallest to largest or largest to smallest
    LocationIndex = index specify location in A_matrix on which sort is based
    LocationVector() gives the original position
    N_rows = number of rows in matrix
    N_columns = number of columns in matrix
    OrderCode = 0 for smallest to largest
    OrderCode = 1 for largest to smallest
    This version looks at each COLUMN in the matrix
    Input: A_matrix(),LocationIndex,N_rows, N_columns
    Output: A_matrix(),LocationVector()
    adapted from version in Hybrid2 2/11/17 JM'''

    flip = 1
    while flip:  # Flip is "true"
        flip = 0
        for i in range(0, N_rows - 1):  # cycle through all values of vector A(_vector)
            # print('A_matrix[i, LocationIndex]',A_matrix[i, LocationIndex])
            # print('A_matrix[i + 1, LocationIndex]',A_matrix[i + 1, LocationIndex])
            if OrderCode == 0:
                if A_matrix[i, LocationIndex] > A_matrix[i + 1, LocationIndex]:
                    # print('A_matrix[i, LocationIndex]',A_matrix[i, LocationIndex])
                    for j in range(0, N_columns):

                        A_temp = Swap(A_matrix[i, j], A_matrix[i + 1, j])
                        A_matrix[i, j] = A_temp[0]
                        A_matrix[i+1, j] = A_temp[1]

                    L_temp = Swap(LocationVector[i], LocationVector[i + 1])
                    LocationVector[i] = L_temp[0]
                    LocationVector[i+1] = L_temp[1]
                    flip = 1   # rows have been swaoped; allow next row to be checked and possibly swapped
            elif OrderCode == 1:
                if A_matrix[i, LocationIndex] < A_matrix[i + 1, LocationIndex]:
                    # print('A_matrix[i, LocationIndex]',A_matrix[i, LocationIndex])
                    for j in range(0, N_columns):

                        A_temp = Swap(A_matrix[i, j], A_matrix[i + 1, j])
                        A_matrix[i, j] = A_temp[0]
                        A_matrix[i+1, j] = A_temp[1]

                    L_temp = Swap(LocationVector[i], LocationVector[i + 1])
                    LocationVector[i] = L_temp[0]
                    LocationVector[i+1] = L_temp[1]
                    flip = 1   # rows have been swaoped; allow next row to be checked and possibly swapped
            else:
                print()
                print("incorrect OrderCode")
                print('Run terminated')
                print()
                sys.exit()

    return A_matrix, LocationVector


'''temp = Sort2D3(x,LocationVector, 4, n, 6,1)
x = temp[0]  # replace original parameter matrix with one reordered by maximum fuel usage

    #for testing!!! 
print('A_matrix')
print(temp[0])
print()
print('LocationVector',temp[1])'''
'''for i in range(n):
    k = LocationVector[i]
    print('LocationVector ')
    print('k',k)'''


# below generator parameters that will go in configuration file are initialized
# GeneratorNumber = np.zeros(n)
'''
GenRatedPower = np.zeros(nGen)
GenMinPower = np.zeros(nGen)
GeneratorFullLoadFuel = np.zeros(nGen)
GeneratorNoLoadFuel = np.zeros(nGen)
GenMinFuel = np.zeros(nGen)
GenFuelCost = np.zeros(nGen)
GenMinRunTime = np.zeros(nGen)

GenOrder = np.zeros((nGen,6))  # this may be unnecessary

# generator parameters from input file are renamed and put in vectors
for i in range(nGen):
    #GeneratorNumber[i] = int(i+1)
    GenRatedPower[i] = genParams[i,0]
    GenMinPower[i] = genParams[i,1]
    GeneratorFullLoadFuel[i] = genParams[i,2]
    GeneratorNoLoadFuel[i] = genParams[i,3]
    GenFuelCost[i] = genParams[i,4]
    GenMinRunTime[i] = genParams[i,5]
                
    GenMinFuel[i] = GenMinPower[i] * (GeneratorFullLoadFuel[i] \
            - GeneratorNoLoadFuel[i]) / GenRatedPower[i] + GeneratorNoLoadFuel[i]
    '''
# else:
#    GenMinFuel[i] = 0.

# Initialize generator according to order in x

'''
this routine is not used now
def FindGenOrder(n):  
    for i in range(n):
        GenOrder[i, 0] = GeneratorNumber[i] #Generator number
        GenOrder[i, 1] = GenRatedPower[i] #rated power
        GenOrder[i, 2] = GenMinPower[i] #minimum power
        GenOrder[i, 3] = GeneratorFullLoadFuel[i] #full load fuel
        GenOrder[i, 4] = GeneratorNoLoadFuel[i] #no load fuel
        GenOrder[i, 5] = GenFuelCost[i] #fuel cost
        GenOrder[i,7] = GenMinRunTime[i]
        
    #Let#s try ordering generators according cheapest fuel use
    #per kWh at rated
        #GenOrder[i, 6] = GenFuelCost[i] * GeneratorFullLoadFuel[i] / GenRatedPower[i] #fuel cost
        #print('GenRatedPower[i]',GenOrder[i, 6])
    
    
    #Now find the order, most expensive first
    #Call Sort2D3[GenOrder[], 6, nGen, 6, dumvector![], 0]
    TotalGeneratorRated = 0
    for i in range(n):
        #GenOrder[i, 6] = GenOrder[i, 6]
        #GenOrder[i, 1] = GenOrder[i, 1]
        TotalGeneratorRated = TotalGeneratorRated + GenRatedPower[i]
    print('TotalGeneratorRated',TotalGeneratorRated)
    return GenOrder
'''
'''NConfig=n**2  # number of configurations
ncols = n + 8  # number of columns in configuration matrix
GenConfig = np.zeros((NConfig, ncols))
new_GenConfig = np.zeros((NConfig, ncols))'''


def OnOffMatrix(n, NConfig, GenConfig):
    '''This produces a matrix of generators that might be on at any given time
    based on the total number of generators available'''
    # ==========================================================
    # THIS PART GENERATES MATRIX CORRESPONDING TO CONFIGURATION
    # OF Generators THAT MAY BE ON
    # ==========================================================
    nsq = NConfig  # number of configurations
    for i in range(nsq):
        # print()
        # print('i',i)
        x = i
        # print('x',x%2)
        for j in range(n):
            GenConfig[i, j] = int(x % 2)
            if GenConfig[i, j] > 0:
                GenConfig[i, n] = GenConfig[i, n] + 1
            if x >= 2:
                x = x/2
            else:
                # print("sssss")
                break
            # print('i,j',i,j,GenConfig[i, j])
        # loop below indicates reasonable result!
        # but note that first row and column will not be used
        '''for i in range(1,nsq+1):
            print()
            print('i',i)
            for j in range(1,n+2):
                print('j',GenConfig[i, j])'''

    return GenConfig


'''temp = OnOffMatrix(n,GenConfig)
GenConfig = temp   # initial version of configuration'''

# df = np.array(temp)
# print(df)  #useful for illustration


def SortedConfig(n, NConfig, GenConfig, genParams, GenFuelCost):
    # ==========================================================
    # This next part finds rated power and total minimum fuel use
    # or lowest cost of energy for each configuration
    # ==========================================================

    nsq = NConfig
    for L in range(nsq):
        for k in range(n):

            # TOTAL RATED POWER OF LTH CONFIGURATION
            GenConfig[L, n + 1] += (
                # rated power of configuration
                GenConfig[L, k] * genParams[k, 0]
            )
            # if L == 2:
            #   print('L,k,GenRatedPower[k-1]',k,GenConfig[L, k]*GenRatedPower[k-1])

            # TOTAL MINIMUM FUEL USE OF LTH CONFIGURATION
            GenConfig[L, n + 2] += (
                GenConfig[L, k] * genParams[k, 1]*(
                    # GenMinFuel[k]
                    genParams[k, 2]-genParams[k, 3])/genParams[k, 0] + genParams[k, 3]
            )

            # TOTAL MINIMUM POWER OF LTH CONFIGURATION
            GenConfig[L, n + 3] += (
                # GenMinPower[k] of configuration
                GenConfig[L, k] * genParams[k, 1]
            )

            '''GenMinPower[i] * (GeneratorFullLoadFuel[i] \
                    - GeneratorNoLoadFuel[i]) / GenRatedPower[i] + GeneratorNoLoadFuel[i]'''

            # Cost of energy at rated power
            GenConfig[L, n + 4] += (
                # GenRatedPower[k]* GenFuelCost[k])
                GenConfig[L, k] * genParams[k, 0] * genParams[k, 4])

            # Cost of energy at minimum power
            GenConfig[L, n + 5] += (
                # * GenMinPower[k]* GenFuelCost[k]
                GenConfig[L, k]*genParams[k, 1]*genParams[k, 4])

            # print('k,genParams[k,4] zzz',k,genParams[k,4])
            # print()
            # minimu run time - is this correct? zzz
            GenConfig[L, n + 7] = genParams[k, 5]

            # 11/11/22: code for whether generator is in configuration
            # GenConfig[L, n + 6 + k] = GenConfig[L, k]

        # Average cost of energy from configuration
        '''if GenConfig[L, n + 2] > 0:
            GenConfig[L, n + 5] /= GenConfig[L, n + 2]
        else:
            GenConfig[L, n + 2] = 0'''
    # ==========================================================
    # Find most expensive generator in each configuration
    # ==========================================================
    for L in range(NConfig):
        # print()
        # print('L zzz',L)
        i_expensive = 0
        temp_coe = 0
        # !!! note that GenFuelCost etc. still keep their originally indexed positio!!!
        for k in range(n):
            # print('k zzz',k,genParams[k,4])
            # if genParams[k,4] > temp_coe:
            if GenConfig[L, k] * GenFuelCost[k] > temp_coe:
                # print()
                # print('k,GenFuelCost[k] zzz',k,genParams[k,4])
                temp_coe = genParams[k, 4]
                i_expensive = k
                # print('k,i_expensive zzz', k, i_expensive)
        GenConfig[L, n + 6] = i_expensive
        '''print()
        print('L', L)
        print('i_expensive zzz', i_expensive)
        print('genParams[k, 4]', genParams[k, 4])
        print('GenConfig[L, n + 6]', GenConfig[L, n + 6])'''

   # ==========================================================
   # Reorder configurations by:
   # total minimum fuel use (OrderCode = 3)
   # or minimum cost of energy (OrderCode = 5)
   # ==========================================================

    # print('OrderCode',OrderCode)
    # sort_column = n + OrderCode   # VBA index shift
    # dum = np.zeros(nsq)

    # Below sorts GenConfig according to the sort_column
    # REPLACE Sort2D3(GenConfig,dum,sort_column,nsq,ncols)
    # GenConfig = temp
    # print(temp[0])

    '''chkk = 0
    tempNConfig = NConfig
    for L in range(1,NConfig):   #= 1 To NConfig - 1    #check each configuration
        #print('L',L)
        #See if maximum power of configuration = that of another
        
        #print('GenConfig[L, n + 2]',GenConfig[L, n + 2])
        
         if GenConfig[L + 1, n + 2] <= GenConfig[L, n + 2]:
         #Redundant configuration, eliminate it
             chkk = 1
             for ll in range[L,NConfig-1]: #= L To NConfig - 1
         
                 if ll < NConfig - 1:
                     GenConfig[ll + 1, k] == GenConfig[ll + 2, k]
                 else:
                     GenConfig[ll + 1, k] == 0

         tempNConfig = tempNConfig - 1'''

    # if chkk = 1 Then GoTo check
    return GenConfig


'''print()
print('@@@@',GenConfig)'''
# OrderCode = 5
# temp = SortedConfig(n,OrderCode,GenConfig)
# ttt=pd.DataFrame(temp)
# print('7777@@@@',ttt)


def ReduceRows(rows, n_row, n_col, col):
    # This is used to eliminate reduncznt rows when values in col are the same
    checkAgain = True
    count = 0
    # col = column to check
    while checkAgain:
        count += 1
        for i in range(n_row-1):  # look at each row
            if rows[i, col] >= rows[i+1, col]:  # check if value is same in next row
                # print('!!! x[i,col]',i,x[i,col],x[i+1,col])
                index_to_delete = i+1
                reduced_matrix = np.delete(rows, index_to_delete, axis=0)
                rows = reduced_matrix
                checkAgain = True
                n_row = rows.shape[0]
                break
            else:
                checkAgain = False

        if count > 100:
            checkAgain = False
            print('safety stop')
            break
    return rows


'''print()
print('filled and sorted GenConfig')
zz = pd.DataFrame(qq)
#zz=np.array(qq)
print(zz)  # this is sorted configuration matric=x
print()'''
'''n_row=GenConfig.shape[0]  # initial number of rows
tempConfig = SortedConfig(n,GenConfig)
new_GenConfig = ReduceRows(tempConfig,n_row,ncols,5)
# n_row = new_GenConfig.shape[0]  # final number of rows not needed

#new_GenConfig = GenConfig  #TEST ONLY!
rr = pd.DataFrame(new_GenConfig)
print()
print('final configurations')
print(rr)'''


def fuel_Use(load, P_rated, fuel_rated, fuel_no_load):
    # This finds the fuel for a given power, assuming a linear relation
    fuel = load*(fuel_rated - fuel_no_load)/P_rated + fuel_no_load
    return fuel


def allocate(n, config, GenConfig, load, GenRatedPower, GenMinPower):
    # This is to allocate load among generators, with ratios
    # find sum of power differences between minimum and rated
    # for all active generators except the one skipped
    # n = number of generators
    loadMax = 0  # maximum load that can be supplied by generators
    loadMaxInc = 0
    loadMin = 0
    tempLoad = np.zeros(n+1)
    extra = 0
    skip = GenConfig[config, n+6]
    # print()
    # print('config',config)
    # print('skip',skip)
    for i in range(n):
        # print('i',i)
        # print('GenRatedPower[i]',GenRatedPower[i])
        if i != skip:
            loadMax += GenConfig[config, i]*GenRatedPower[i]
            loadMin += GenConfig[config, i]*GenMinPower[i]
            loadMaxInc += GenConfig[config, i] * \
                (GenRatedPower[i] - GenMinPower[i])

    if load > loadMax:
        # print("too much power asked for")
        loadMax = loadMax  # placeholder for print statement
    else:

        for i in range(n):
            # find proportional allocation of generators before considering minimum power level
            if i != skip:
                # print('i xxx',i,loadMax)
                if loadMax > 0:
                    tempLoad[i] = load*GenConfig[config, i] * \
                        GenRatedPower[i]/loadMax
                else:
                    loadMax = 0
                # print()
                # print('GenRatedPower[i]',GenRatedPower[i])
                # print('i,tempLoad[i]',i,tempLoad[i])
            else:
                tempLoad[i] = 0
        # for i in range(n):
        #    print()
        #    print('i',i)

        for j in range(0, n):
            # print('j',j)
            if j != skip and GenConfig[config, j] != 0:

                if tempLoad[j] < GenMinPower[j]:
                    # print('tempLoad[j] 2',tempLoad[j])
                    extra += GenMinPower[j] - tempLoad[j]
                    # print('i,extra',i,extra)
                    tempLoad[j] = GenMinPower[j]

                    # print('j,tempLoad[j]',j,tempLoad[j])

                if j < n:
                    # print('j,i',j)
                    # print('j,tempLoad[j+1] ',j,tempLoad[j+1] )
                    if j+1 != skip and j+1 < n:
                        if tempLoad[j+1] > extra:
                            tempLoad[j+1] = tempLoad[j+1]-extra
                            extra = 0
                    # print('j',j)
                    # print('!!!! tempLoad[j]',tempLoad[j])
                    # print('!!!!!tempLoad[j+1] ',tempLoad[j+1] )
                    #
                    # extra = 0

                    # print('GenRatedPower[j+1]',GenRatedPower[j+1])
            elif GenConfig[config, j] != 0:
                # print('j',j)
                # print('j+1,tempLoad[j+1]',j+1,tempLoad[j+1])
                tempLoad[j+1] = tempLoad[j+1]-extra
                extra = 0
                # print('tempLoad[j+1]',tempLoad[j+1])
                # print()
                # print('***')
                # print('tempLoad[j+1]',tempLoad[j+1])

# update tempLoad according to whther the generator is on or not
    for i in range(n):
        tempLoad[i] = tempLoad[i]*GenConfig[config, i]

    # print('temp',temp)
    # print('loadMaxInc',loadMaxInc)

    # print('fuel',fuel)
    # print(tempLoad)
    return loadMax, tempLoad, extra


def new_fuel_use(genParams, netLoad, n, GenConfig, GenOn, GenRatedPower, GenMinPower,
                 GeneratorFullLoadFuel, GeneratorNoLoadFuel, GenFuelCost, GenMinRunTime,
                 GenMinFuel):
    nGen = n
    extraMinRun = 0
    N_config = GenConfig.shape[0]
    for i in range(n):
        GenRatedPower[i] = genParams[i, 0]
        GenMinPower[i] = genParams[i, 1]
        GeneratorFullLoadFuel[i] = genParams[i, 2]
        GeneratorNoLoadFuel[i] = genParams[i, 3]
        GenFuelCost[i] = genParams[i, 4]
        GenMinRunTime[i] = genParams[i, 5]

        GenMinFuel[i] = GenMinPower[i] * (GeneratorFullLoadFuel[i]
                                          - GeneratorNoLoadFuel[i]) / GenRatedPower[i] + GeneratorNoLoadFuel[i]

    # print('N_config',N_config)
    config = N_config   # default for load greater than maximum rated of configuration
    unmet = 0  # this many not be used!
    # print('N_config',N_config)
    '''key configuration parameters 
        n+1 = column for rated power
        n+2 = column for cost fuel use/unit power at minimum power
        n+3 = column for minimum power
        n+4 = column for cost for fuel at rated power
        n+5 = column for cost for fuel at minimum power
        n+6 = column for index of most expensive generator
        
        key for generators
        GenRatedPower[i] 
        GenMinPower[i] 
        GeneratorFullLoadFuel[i] 
        GeneratorNoLoadFuel[i]
        GenFuelCost[i]
        MinRunTime[i]
    '''

    # OK = 1 #test
    for i in range(N_config):
        # OK = 1 #test
        if netLoad > GenConfig[i, n+1]:
            config = i
        else:
            config += 1
            # print("config",config)
            break

            # print('test min run')
            '''for j in range(n):
                
                if GenConfig[i,j]==0: #and GenOn[j]==1:
                    #print()
                    #print('j',j)
                    OK = 0
                    #print(GenConfig[i,j])
                    #print('GenOn[j]',GenOn[j])
                    #print('OK',OK)
                    break
            #print('OK 2',OK)
            if OK == 1:
                config = i
                break'''

            # fix this here!
    tempArray = np.array(GenConfig)  # numpy array of configurations
    # columns with generators on in configuration
    columns = np.arange(0, nGen)
    # array of generators all configurations
    GensInConfig = np.array(tempArray[:, columns])

    GensMinRun = np.array(GenOn)  # GensMinRun
    GensMinRun = GensMinRun.astype(int)
    if netLoad > 0:
        GensDueToLoad = GensInConfig[config, :]  # GensDueToLoad
        GensDueToLoad = GensDueToLoad.astype(int)
    # print('GensDueToLoad',GensDueToLoad)
    # print('GensMinRun',GensMinRun)
    if netLoad > 0:
        GensAll = GensMinRun | GensDueToLoad  # GensAll
    else:
        GensAll = GensMinRun

    newConfig = np.flatnonzero(np.equal(GensInConfig, GensAll).all(1))
    '''print('newConfig',newConfig)
    #config = newConfig[0]  # update
    print('GensAll')
    print(GensAll)
    print(GensInConfig)
    print('newConfig')
    print(newConfig)'''
    if netLoad > 0:
        # print('newConfig[0]', newConfig[0])  # !!!
        if config > 1:  # added 5/30/26 to try to eliminate bug !!!
            if newConfig[0] != config:
                '''Note!! at the moment extraMinRun only affects the fuel; it is not available for storage'''
                # print('newConfig,config zzz',newConfig,config)
                # find extra power due to minimum run time
                extraMinRun = GenConfig[newConfig[0], nGen+3]\
                    - GenConfig[config, nGen+3]
                # print('extraMinRun zzz',extraMinRun)
                netLoad += extraMinRun
    '''print()
    print('GenConfig[newConfig[0],nGen+3] zzz',GenConfig[newConfig[0],nGen+3])
    print('GenConfig[config,nGen+3] zzz',GenConfig[config,nGen+3])
    print('extraMinRun zzz',extraMinRun)
    print('netLoad zzz',netLoad)'''
    if config > 1:  # added 5/30/26 to try to eliminate bug !!!
        config = newConfig[0]
    # print(config)

    # index=[k for k,row in enumerate(GensInConfig)if row==D]
    # print('index',index)
    # g=index[0]
    # print(g)

    '''    for j in range(i,N_config):
                #check all other configurations to see if they
                # include generators that must be on
                for k in range(nGen):
                    if GenConfig[j,k]==GenOn[k]:
                        config=j'''

    '''for j in range(n):                
                if GenConfig[i,j]==0 and GenOn[j]==1:
                    # this applies when generator is required to be on
                    OK = 0
                    #break
            if OK == 1:
                config = i
                break'''
    # print()
    # print('config fuel use zzz',config)
    # print('GenOn zzz',GenOn)

    ''' Try this alternative below
    for i in range(N_config):
        OK = 1 #test
        for j in range(n):
            #print('zzz GenOn[j]',GenOn[j])                
            if GenOn[j]==1:
                # this applies when generator is required to be on
                OK = 0
                break
        if OK == 1:
            config = i
            break    
        print('i, config zzz',i,config)
    Try this alternative above'''

    if config == N_config:
        config = config - 1
        unmet = netLoad - GenConfig[i, n+1]
        # met = GenConfig[i,n+1]

    # method below assumes all generators except most expensive
    # are run at full power
    fuel = 0
    nonPeak = 0
    extra = 0
    # print('config',config)

    if GenConfig[config, n] == 1:   # check to see if there is only one generator on
        # only one generator is on
        for i in range(n-1):   # note change 5/29/26 from range(n)
            # print('i,GenConfig[config,i]',i,GenConfig[config,i])
            if GenConfig[config, i] == 1:
                k = i
                # print()
                # print('k',k)
                # print('GenMinPower[k]',GenMinPower[k])

        peak = netLoad
        if peak <= GenMinPower[k]:
            extra = GenMinPower[k]-peak
            peak = GenMinPower[k]
            # print('peak 2',peak)
            # print('2 extra',extra)
        fuel = GenConfig[config, k]*(peak*(GeneratorFullLoadFuel[k]
                                           - GeneratorNoLoadFuel[k])/GenRatedPower[k] + GeneratorNoLoadFuel[k])
    else:
        # more than one generator
        # look at first columns to see which generator is on # note change 5/29/26 now n-1
        for i in range(n-1):

            # ensures most expensive runs at lowest power possible
            if i != GenConfig[config, n+6]:
                '''
                print('i',i)
                print('GenConfig[config,n+6]',GenConfig[config,n+6])
                print('GenRatedPower[i]',GenRatedPower[i])
                print('GeneratorFullLoadFuel[i]',i,GeneratorFullLoadFuel[i])
                print('GenConfig[config,i]',GenConfig[config,i])'''
                nonPeakMax = GenConfig[config, i] * \
                    GenRatedPower[i]  # test change removed 5/30/26
                nonPeakMax = GenRatedPower[i]  # replaced with this 5/30/26
                # fuel += GenConfig[config,i]*GeneratorFullLoadFuel[i]
            else:
                k = i  # code for peak generator
                # print('peak,k',k)

        # compare load to nonPeak
        diff = netLoad-nonPeakMax   # note change 5/29/26  now reversed
        if diff < GenMinPower[k]:
            peak = GenMinPower[k]
            nonPeak = netLoad-peak

        else:
            peak = diff
            nonPeak = netLoad-peak

        #!!! more tests!!! zzz
        '''if netLoad<GenMinPower[k]:
            peak=GenMinPower[k]
            nonPeak = netLoad-peak'''
        # else:

        # fuel use by peak generator
        fuel = GenConfig[config, k]*(peak*(GeneratorFullLoadFuel[k]
                                           - GeneratorNoLoadFuel[k])/GenRatedPower[k] + GeneratorNoLoadFuel[k])

        temp = allocate(n, config, GenConfig, nonPeak,
                        GenRatedPower, GenMinPower)
        extra = temp[2]
        genLoad = temp[1]

        for i in range(n):
            if i != k:
                fuelTemp = GenConfig[config, i]*fuel_Use(genLoad[i],
                                                         GenRatedPower[i], GeneratorFullLoadFuel[i], GeneratorNoLoadFuel[i])
                fuel += fuelTemp
            # print('fuelTemp',fuelTemp)

        '''peak = netLoad - nonPeak - unmet
        if peak<GenMinPower[k]:
                
                nonPeak = nonPeak +(peak -GenMinPower[k])
                peak = GenMinPower[k]'''

        # from allocate: loadMax,tempLoad,extra
        # print('loadMax',temp[0])
        # print('tempLoad',temp[1])
        # print('extra',temp[2])

    # print('peak',peak)
    # print('nonPeak',nonPeak)
    '''dum = peak*(GeneratorFullLoadFuel[k]\
                -GeneratorNoLoadFuel[k])/GenRatedPower[k]+GeneratorNoLoadFuel[k]
    print('dum',dum)'''
    '''Below assumes non-peak generators operated at rated
That will not be true when certain generators are forced to stay on'''

    # peak or only generator fuel
    # fuel = GenConfig[config,k]*(peak*(GeneratorFullLoadFuel[k]\
    #        -GeneratorNoLoadFuel[k])/GenRatedPower[k] + GeneratorNoLoadFuel[k])
    # print('peak fuel',fuel)

    # find additional fuel when there is more than one generator
    '''if GenConfig[config,n] != 1:
        genLoad =temp[1]
        
        for i in range(1,n):
            fuelTemp = fuel_Use(genLoad[i],GenRatedPower[i],GeneratorFullLoadFuel[i],GeneratorNoLoadFuel[i])
            fuel+=fuelTemp
            #print('fuelTemp',fuelTemp)'''
    # print('config zzz new_fuel_use',config)
    return config, fuel, unmet, extra, extraMinRun


'''GenOn[0]=1  # test!
GenOn[1]=1  # test!
GenOn[2]=1  # test!'''

# netLoad = 150   #for testing!!!
# temp = new_fuel_use(netLoad, n, new_GenConfig,GenOn)


def initialize(genParams, n, NConfig, GenFuelCost, GenMinRunTime):
    # NConfig=n**2  # number of configurations
    ncols = n + 8  # number of columns in configuration matrix
    GenConfig = np.zeros((NConfig, ncols))
    new_GenConfig = np.zeros((NConfig, ncols))

    temp = OnOffMatrix(n, NConfig, GenConfig)
    GenConfig = temp   # initial version of configuration
    dum = np.zeros(n)
    temp = Sort2D3(genParams, dum, 4, n, 6, 1)
    # replace original parameter matrix with one reordered by maximum fuel usage
    genParams = temp[0]

    genParams = np.array(genParams)

    n_row = GenConfig.shape[0]  # initial number of rows
    GenConfig_print = pd.DataFrame(GenConfig)

    tempConfig = SortedConfig(n, NConfig, GenConfig, genParams, GenFuelCost)
    # check min run times to avoid possibility of problem from reducing redundant rows
    # !!! eliminating reduce causes a problem; remove for time being 5/30/26
    noReduce = 0
    '''for i in range(n):
        if GenMinRunTime[i] > 1:
            noReduce = 1
            break'''

    if noReduce == 0:
        new_GenConfig = ReduceRows(tempConfig, n_row, ncols, n+2)
    else:
        new_GenConfig = tempConfig

    # n_row = new_GenConfig.shape[0]  # final number of rows not needed

    # new_GenConfig = GenConfig  #TEST ONLY!
    '''GenConfig_print = pd.DataFrame(new_GenConfig)
    print()
    print('final configurations')
    print(GenConfig_print)'''
    n_row_final = new_GenConfig.shape[0]
    return new_GenConfig, n_row_final


def main():

    pd.set_option('display.max_columns', None)
    # Set expand_frame_repr to False to force all columns onto one line
    # (if the output area width allows it, otherwise it might still wrap)
    pd.set_option('display.expand_frame_repr', False)

    input_file = "generators_test_3.csv"

    file_path = input_file
    _, ext = os.path.splitext(file_path.lower())

    data = np.loadtxt(input_file, delimiter=",",
                      skiprows=1)  # data starts on 2nd row

    # get number of generators from generator paramater file
    nGen = data.shape[0]
    NConfig = nGen**2            # number of configurations
    genParams = np.array(data)      # generator paramaters from input file

    GenOn = np.zeros(nGen)

    n = nGen

    GenRatedPower = np.zeros(nGen)
    GenMinPower = np.zeros(nGen)
    GeneratorFullLoadFuel = np.zeros(nGen)
    GeneratorNoLoadFuel = np.zeros(nGen)
    GenMinFuel = np.zeros(nGen)
    GenFuelCost = np.zeros(nGen)
    GenMinRunTime = np.zeros(nGen)

    # GenOrder = np.zeros((nGen,6))  # this may be unnecessary

    # generator parameters from input file are renamed and put in vectors
    for i in range(nGen):
        # GeneratorNumber[i] = int(i+1)
        GenRatedPower[i] = genParams[i, 0]
        GenMinPower[i] = genParams[i, 1]
        GeneratorFullLoadFuel[i] = genParams[i, 2]
        GeneratorNoLoadFuel[i] = genParams[i, 3]
        GenFuelCost[i] = genParams[i, 4]
        GenMinRunTime[i] = genParams[i, 5]

        GenMinFuel[i] = GenMinPower[i] * (GeneratorFullLoadFuel[i]
                                          - GeneratorNoLoadFuel[i]) / GenRatedPower[i] + GeneratorNoLoadFuel[i]

    temp = initialize(genParams, n, NConfig, GenFuelCost, GenMinRunTime)
    new_GenConfig = temp[0]

    pd.set_option('display.max_columns', None)
    # Set expand_frame_repr to False to force all columns onto one line
    # (if the output area width allows it, otherwise it might still wrap)
    pd.set_option('display.expand_frame_repr', False)

    input_file = "generators_test_3.csv"

    file_path = input_file
    _, ext = os.path.splitext(file_path.lower())

    data = np.loadtxt(input_file, delimiter=",",
                      skiprows=1)  # data starts on 2nd row

    # get number of generators from generator paramater file
    nGen = data.shape[0]
    NConfig = 2**nGen   # nGen**2            # number of configurations
    genParams = np.array(data)      # generator paramaters from input file

    count = np.zeros(n)

    netLoad = [50, 50, 50, 150, 150, 50, 50, 150, 150, 50, 50, 50]
    cfg = 0
    for i in range(12):
        '''print()
        print('i',i)
        print('netLoad',netLoad[i])'''
        for j in range(n):
            '''print()
            print('j',j)
            print('new_GenConfig',new_GenConfig[cfg,j])
            print('count[j]',count[j])
            print('GenMinRunTime[j]',GenMinRunTime[j])'''
            if new_GenConfig[cfg, j] == 1:
                if GenMinRunTime[j] > count[j]:
                    # generator must still be running

                    count[j] += 1
                    GenOn[j] = 1
                    # print('zzz',j,count[j])
                else:
                    # reset minimum run counter
                    count[j] = 0
                    GenOn[j] = 0
                # print('j,GenOn',j,GenOn[j])
        temp = new_fuel_use(genParams, netLoad[i], n, new_GenConfig, GenOn, GenRatedPower,
                            GenMinPower, GeneratorFullLoadFuel, GeneratorNoLoadFuel, GenFuelCost, GenMinRunTime, GenMinFuel)
        cfg = temp[0]

        # print('fuel',fuel)  # xxx

        # print('configuration',temp[0])
        # print('fuel',temp[1])
        # print('unmet',temp[2])
        # print('extra',temp[3])

        return genParams, n, GenRatedPower, GenMinPower, GeneratorFullLoadFuel, \
            GeneratorNoLoadFuel, GenFuelCost, GenMinRunTime, GenMinFuel


if __name__ == "__main__":
    temp = main()
    genParams = temp[0]
    n = temp[1]
    GenRatedPower = temp[2]
    GenMinPower = temp[3]
    GeneratorFullLoadFuel = temp[4]
    GeneratorNoLoadFuel = temp[5]
    GenFuelCost = temp[6]
    GenMinRunTime = temp[7]
    GenMinFuel = temp[8]
