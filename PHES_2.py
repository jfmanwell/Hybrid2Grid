'''
-------------------------------------------------------------------------------
Script for pumped hydroelectric energy storage
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
Created on Wed Jun  3 14:42:52 2026
This script is to characterize pumped hydroelectric energy storage
Most recent update 6/30/26
"""


# common inputs
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
g = 9.81                    # gravitational constant, m/s^2
aRou = 0.0005               # Absolute roughness, m
rho_w = 1000                # density of water, kg/m^3
mu = 1.002e-3               # viscosity of water, kg m^-1, s^-1
pi = np.pi                  # pi


def Reynolds(rho, mu, V, D):
    # This finds the Reynolds number; uses SI units
    # Inputs:
    #    rho_w: water density, kg/m^3
    #    V: velocity, m/s
    #    D:  diameter, m
    Re = rho*V*D/mu
    return Re


def f_friction(D, Re, aRou):
    # This finds the Swamee Jain friction factor
    # Inputs:
    #   D       : Inner diameter of the pipe , m
    #   Re      : Reynolds Number , -
    #   aRou    : Absolute roughness of pipe, m

    f = (0.25/(np.log10((aRou/D)/3.7+5.74/Re**0.9))**2)

    return f


def reservoir_height(D_res, Vol, head_0, height):
    '''
    Inputs
    ----------
    D_res : diameterof reservoir, m
    Vol : volume of water in reservoir, m^3
    head_0 : height of bottom of reservoir, m
    height : height of water surface, m

    Returns
    height : final height, m
    '''
    A = pi*D_res**2/4               # surface area, assuming circular, m^2
    delta_h = Vol/A                 # change in height, m
    height = delta_h + head_0       # new height, m
    if height < head_0:
        height = head_0
    return height


def solve_charge(P, V, D, head, L, V_guess=None):

    if V_guess is None:
        V_guess = 2.0  # m/s reasonable velocity guess

    def flowCharge(V):
        Q = V*np.pi*D**2/4  # volume flow rate
        Re = Reynolds(rho_w, mu, V, D)  # Reynolds number
        f = f_friction(D, Re, aRou)  # friction factor
        head_loss = f*L*V**2/(D*2*g)  # frictional head loss in pipe, m
        P_temp = Q*rho_w*g*head  # useful power, W
        P_loss = Q*rho_w*g*head_loss  # power loss, W
        return P - (P_temp + P_loss)  # this should be close to zero

    # final value for velocity, m/s
    V_final = fsolve(flowCharge, V_guess).item()
    Q_final = V_final*np.pi*D**2/4   # final value for volume flow rate, m^3/s
    P_useful = Q_final*rho_w*g*head     # useful power, W

    return {'Q_final': Q_final, 'P_useful': P_useful}


def solve_discharge(P, V, D, head, L, eff_turbine, Vol_stored, V_guess=None):
    # This will be used when power out is to be specified rather than flow
    if V_guess is None:
        V_guess = 2.0  # reasonable velocity first guess, m/s

    def flowDischarge(V):
        Q = V*np.pi*D**2/4                  # volume flow rate
        Re = Reynolds(rho_w, mu, V, D)      # Reynolds number
        f = f_friction(D, Re, aRou)         # friction factor
        head_loss = f*L*V**2/(D*2*g)        # frictional head loss in pipe, m
        P_temp = Q*rho_w*g*(head)             # useful power, W
        P_loss = Q*rho_w*g*head_loss        # power loss, W
        # this should be close to zero
        return P+P_loss - (P_temp + P_loss)

    # final value for velocity

    V_final = fsolve(flowDischarge, V_guess).item()
    Q_final = V_final*np.pi*D**2/4          # final value for volume flow rate, m^3/s
    if Q_final*3600 < Vol_stored:
        P_useful = Q_final*rho_w*g*head*eff_turbine      # useful power, W
        Re = Reynolds(rho_w, mu, V_final, D)             # Reynolds number
        f = f_friction(D, Re, aRou)                      # friction factor
        # frictional head loss in pipe, m
        head_loss = f*L*V_final**2/(D*2*g)
        P_loss = Q_final*rho_w*g*head_loss               # power loss, W
    else:
        Q_final = 0
        P_useful = 0
        P_loss = 0

    # print()
    # print('P_useful dis', P_useful)
    return {'Q_final': Q_final, 'P_useful': P_useful, 'P_loss': P_loss}


def flow_Loss(Q, D, L):
    # This is to find loss for a given flow
    V = Q/(np.pi*D**2/4)  # volume flow rate
    if V > 0:
        Re = Reynolds(rho_w, mu, V, D)  # Reynolds number
        f = f_friction(D, Re, aRou)  # friction factor
        head_loss = f*L*V**2/(D*2*g)  # frictional head loss in pipe, m
        P_loss = Q*rho_w*g*head_loss
    else:
        P_loss = 0
    return P_loss


def IdealReservoir(Vol_required, Vol_stored, Vol_storedMax):
    # This is to model a storage reserovoir
    # Filling/discharging are 100% efficient
    # Note that in filling Vol_required is negative
    '''
    Inputs
        Vol_required: volume required from reservoir, m^3
        Vol_stored: volume in reservoir, m^3
        Vol_storedMax: maxmimum capacity of reservoir, m^3

    Returns
        Vol_in: volume in, m^3
        Vol_out: volume out, m^3
        Vol_not_met: volume not met, m^3
        Vol_excess: excess volume,  m^3
        NewVol_stored: updated stored volume, m^3

    '''
    Vol_in = 0                              # volume into storage, m^3
    Vol_out = 0                             # volume into storage, m^3
    Vol_not_met = 0                         # load not met, m^3
    Vol_excess = 0                          # excess volume, m^3

    if Vol_required > 0:                    # this is for discharge
        if Vol_required < Vol_stored:       # enough storage to supply required volume
            Vol_out = Vol_required
        else:
            Vol_out = Vol_stored      # not enough volume; take what there is
            Vol_not_met = Vol_required - Vol_stored
        NewVol_stored = Vol_stored - Vol_out  # volume level is now lower

    else:
        # this is for filling the reservoir
        if Vol_stored - Vol_required < Vol_storedMax:  # try to put all available volume into reservoir
            Vol_in = -Vol_required
            NewVol_stored = Vol_stored + Vol_in  # updated reservoir volume
            Vol_excess = 0

        else:  # reservoir cannot take all that is available
            Vol_excess = Vol_stored - Vol_storedMax - Vol_required
            Vol_in = -Vol_required - Vol_excess
            NewVol_stored = Vol_storedMax

    return {'Vol_in': Vol_in,
            'Vol_out': Vol_out,
            'Vol_not_met': Vol_not_met,
            'Vol_excess': Vol_excess,
            'Vol_stored': NewVol_stored}


def main():
    # sample inputs provided below
    D = 2                               # penstock diameter, m
    D_res = 2500                         # diameter of reservoir, m
    head_0 = 30                         # elevation of reservoir w.r.t. intake, m
    L = 70                              # penstock length, m
    eff_pump = .9                        # pump efficiency, -
    eff_turbine = .9                     # turbine efficiency, -

    E_in = 0
    E_out = 0
    E_excess = 0
    E_not_met = 0
    Vol_stored = 0
    E_loss_chg = 0
    E_loss_dis = 0
    Vol_storedMax = 20000  # 27500
    nPts = 10  # numper of points, this is for testing
    Vol = np.zeros(nPts)            # volume for plotting, m^3
    time = np.zeros(nPts)           # time for plotting, hrs
    z = np.zeros(nPts)              # head for plotting, m

    # test power, kW
    P = [-500, -500, -500, -500, -500, 500, 500, 500, 500, 500]

    head = head_0  # initial head (empty reservoir), m
    for i in range(nPts):
        time[i] = i

        if P[i] < 0:
            # fill reservoir when P <0
            V_guess = 1                     # velocity guess, m/s
            charge = solve_charge(-P[i]*eff_pump*1e3,
                                  V_guess, D, head, L, V_guess)
            # volume flow into reservoir, m^3/s
            Q_in = charge['Q_final']
            reservoir = IdealReservoir(-Q_in*3600, Vol_stored, Vol_storedMax)
            NewVol = reservoir['Vol_stored']  # new volume in reservoir, m^3

            if reservoir['Vol_excess'] < .001:             # Check if there is any excess
                # no excess, reservoir can take all thw water available
                E_in += -P[i]*1e3                           # input energy, Wh
                # flow rate in, m^3/s
                Q_in = reservoir['Vol_in']/3600
                # charging flow loss, W
                flow_loss_chg = flow_Loss(Q_in, D, L)
                P_loss_chg = -P[i]*(1-eff_pump)*1000+flow_loss_chg
            else:
                # reservoir full or near full
                E_in += -P[i]*1e3  # input is same, even is some is excess
                Q_in_red = reservoir['Vol_in']/3600  # reduced flow rate, m^3/s
                # flow power loss in filling, W
                flow_loss_chg = flow_Loss(Q_in_red, D, L)

                P_in_reduced = Q_in_red*rho_w*g*head  # reduced power in, W

                P_loss_reduced = (P_in_reduced*(1/eff_pump-1) +
                                  flow_loss_chg)  # loss, W

                P_excess = -P[i]*1e3-P_in_reduced-P_loss_reduced
                # net1 = P_in_reduced+P_loss_reduced+P_excess
                E_excess += P_excess
                P_loss_chg = P_loss_reduced
            Vol_stored = NewVol

            # P_loss_chg = -P[i]*(1-eff_pump)*1000+flow_loss_chg
            E_loss_chg += P_loss_chg

        elif P[i] > 0:
            # This is for discharging

            V_guess = 1
            # solve_discharge(P, V, D, head, L, V_guess=None)
            discharge = solve_discharge(
                P[i]*1e3/eff_turbine, V_guess, D, head, L, eff_turbine, Vol_stored, V_guess)
            Q_out = discharge['Q_final']

            '''Note!  If there is not enough water to fullfill the requirement
               the load is not met and no water is taken is taken out'''
            if Q_out > 0:
                reservoir = IdealReservoir(
                    Q_out*3600, Vol_stored, Vol_storedMax)
                NewVol = reservoir['Vol_stored']
                E_out += P[i]*1e3
                P_loss_dis = P[i]*1e3*(1/eff_turbine-1)+discharge['P_loss']
                E_loss_dis += P_loss_dis

            else:
                NewVol = Vol_stored
                E_not_met += P[i]*1e3
            Vol_stored = NewVol

        head = reservoir_height(D_res, Vol_stored, head_0, head)
        Vol[i] = Vol_stored         # for plotting
        z[i] = head                 # for plotting

    E_store = rho_w*NewVol*.5*(head+head_0)*g/3.6e3  # convert to Wh
    # energy balance below, should be close to 0
    balance = E_in-E_out-E_loss_chg-E_loss_dis-E_excess-E_store

    print()
    print(f"final volume: {NewVol:,.1f} m^3")
    print(f"stored energy: {E_store/1000:,.1f} kWh")
    print(f"E_in: {E_in/1000:,.1f} kWh")
    print(f"E_out: {E_out/1000:,.1f} kWh")
    print(f"E_not_met: {E_not_met/1000:,.1f} kWh")
    print(f"E_excess: {E_excess/1000:,.1f} kWh")
    print(f"charging loss: {E_loss_chg/1000:,.1f} kWh")
    print(f"discharging loss: {E_loss_dis/1000:,.1f} kWh")
    print(f"balance: {balance/1e3:,.1f} kWh")
    plt.figure(0)
    plt.plot(time, Vol)
    plt.title("Volume stored")
    plt.xlabel("Time, hrs")
    plt.ylabel("Volume, m^3")

    plt.figure(1)
    plt.plot(time, z)
    plt.title("Head")
    plt.xlabel("Time, hrs")
    plt.ylabel("Head, m")

    plt.show()

    return


if __name__ == "__main__":
    main()
