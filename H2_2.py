'''
-------------------------------------------------------------------------------
Script for investigating hydrogen synthesis and compression
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
Created on Wed Jan 21 16:51:15 2026
@author: jfman
This is for the production and compression of hydrogen
most recent update: 6/30/26
"""
import numpy as np
from scipy.optimize import fsolve
from datetime import datetime
pi = np.pi
Ru = 8.31446                                # ideal gas consant, m^3 kPa/kmol K
M = 2.016                                   # molecular weight of hydrogen, kg/kmol
k = 1.407                                   # ratio of Cp to Cv for hydrogen
# energy content of H2 (HHV), kWh/kg
e_H2 = 39.4


def H2_mass_work(P1, P2, T_K, v_spec, eff, E_wind, flag):
    """ specific compressor work"""
    # isentropic specific compressor work, kWh/kg

    w_c_isen = P1*v_spec*(k/(k-1))*((P2/P1)**((k-1) / k)-1) / \
        3600

    # isothermal specific compressor work kWh/kg
    w_c_iso = P1*v_spec*np.log(P2/P1)/3600

    # specific work required for electrolysis, kWh/kg
    a = e_H2/eff

    """ use search to find mass that can be produced (via electolysis) 
   and compressed with energy available"""

    # isentropic case

    def f2(m_H2, a=a, b=w_c_isen):
        return abs(E_wind - m_H2*(a + b))

    # Isothermal case

    def f1(m_H2, a=a, b=w_c_iso):
        return abs(E_wind - m_H2*(a + b))

    if flag == 1:
        # this is for isothermal compression
        result = fsolve(f1, 1)                           # use fsolve here
        # mass hydrogen produced (isothermal)
        m_H2 = float(result[0])
        V1 = v_spec*m_H2
        # volume after isothermal compression
        V2 = V1*(P1/P2)
        T2 = T_K                                # temperature after isothermal compression
        # electrical energy used for electrolysis (isothermal)
        E_elec = m_H2*a
        # electrical energy used for compression (isothermal)
        W_c = m_H2*w_c_iso

    else:
        # this is for isentropic compression
        result = fsolve(f2, 1)                           # use fsolve here
        # mass hydrogen produced (isentropic)
        m_H2 = float(result[0])
        V1 = v_spec*m_H2                       # volume of H2 at initial temperature
        # volume after isentropic compression
        V2 = V1*(P2/P1)**(-1/k)
        # temperature after isentropic compression
        T2 = T_K*(P2/P1)**(1-1/k)
        E_elec = m_H2*a                   # isentropic electrolysis work, kWh
        W_c = m_H2*w_c_isen               # isentropic compressor work, kWh

    return m_H2, V2, T2, E_elec, W_c


def main(flag):

    P1 = 100                                     # initial pressure, kPa
    P2 = 40000.                                 # final pressure, kPa

    # maximum practical H2 storage pressure is ~ 70,000 kPa

    T_C = 15                                      # initialtemperature, C
    T1 = T_C + 273.15                            # initial temperature, K
    # specific volume of H2 @ P1, T
    v_spec = (Ru/M)*T1/P1
    Pr = 5                                      # rated power of wind turbine, MW
    cf = .5                                     # assumed capacity factor of turbine
    t = 1.                                      # time period of interest, h
    # energy available to produce and compress H2, kWh
    E_wind = cf*Pr*t*1000
    # assumed electrolysis conversion efficiency, -
    eff = 0.65

    # Get results from H2 production
    m_H2, V2, T2,  E_elec, W_c = H2_mass_work(
        P1, P2, T1, v_spec, eff, E_wind, flag)

    print()
    print(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("initial specific volume H2: %5.2f" % v_spec, "m^3/kg")
    print("Energy from wind: %5.2f" % E_wind, "kWh")
    print()

    if flag == 1:
        print("Isothermal compression")
        # print(f"initial volume H2: = {v_spec*m_H2_iso:,.4} m^3")
        print(f"final volume H2 isothermal = {V2:,.2} m^3")
        print(f"mass of hydrogen isothermal = {m_H2:,.3} kg")
        # print("specific compressor work: % 5.2f" % w_c_iso, "kWh/kg")
        print("Electrolysis work: %5.2f" % E_elec, "kWh")
        print("Compressor work: %5.2f" % W_c, "kWh")
        print("Final temperature: %5.2f" % T2, "K")

    else:
        print("Isentropic compression")
        print(f"initial volume H2: = {v_spec*m_H2:,.4} m^3")
        print(f"final volume H2 isentropic = {V2:,.2} m^3")
        print(f"mass of hydrogen isenstropic = {m_H2:,.3} kg")
        # print("specific compressor work: % 5.2f" % w_c_isen, "kWh/kg")
        print("Electrolysis work: %5.2f" % E_elec, "kWh")
        print("Compressor work: %5.2f" % W_c, "kWh")
        print("Final temperature: %5.2f" % T2, "K")

    # find dimension assuming spherical tank, V = (4/3)*pi*R**3

    if flag == 1:
        # radus of spherical tank, isothermal
        R2_sphere_iso = np.cbrt(V2/(4*pi/3))
    else:
        # radus of spherical tank, isentropic
        R2_sphere_isen = np.cbrt(V2/(pi*4/3))

    print()
    # Tank dimensions at T and V
    print("Dimensions at initial conditions")
    print()
    print("Dimensions at final pressure")

    if flag == 1:
        # Tank dimensions with isothermal compression
        print("R2_sphere_isothermal %5.2f" % R2_sphere_iso, "m")
        m_check_iso = V2*P2/(Ru*T2/M)
        print('check on mass, isothermal:% 5.2f' % m_check_iso, 'kg')

    else:
        # Tank dimensions with isentropic compression
        print("R2_sphere_isen %5.2f" % R2_sphere_isen, "m")
        m_check_isen = V2*P2/(Ru*T2/M)
        print('check on mass, isentropic:% 5.2f' % m_check_isen, 'kg')


# flag for type of compression, = 1 for isothermal
compressor_flag = 2

if __name__ == "__main__":
    main(compressor_flag)

"""from original problem about hindenburg
m_H2_isen = P1*V*M/(Ru*(T + 273.1))           # mass of hydrogen, kg
E = m_H2*e_H2/eff                       # energy required, kWh
print(f"Required energy = {E/1e6:,.4} GWh")
print(f"mass of hydrogen = {m_H2/1000:,.4} tonnes")
"""
