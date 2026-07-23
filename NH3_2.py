'''
-------------------------------------------------------------------------------
Script for ammmonia production
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
This is for the production of ammonia from hydrogen
It uses functions in H2_2.py, which is imported
most recent update: 6/30/26
"""

import numpy as np
from datetime import datetime
import H2_2
pi = np.pi
Ru = 8.31446                            # ideal gas consant, m^3 kPa/kmol K
M_H2 = 2.016                            # molecular weight of hydrogen, kg/kmol
k = 1.407                               # ratio of Cp to Cv for hydrogen
e_H2 = 39.4                             # energy content of H2, kWh/kg
M_NH3 = 17.03                           # molecular weight of NH3, kg/kmol
rho_NH3 = 681                           # density liquid NH3 at -33.3 C, kg/m^3

# set compressor_flag = 1 for isothermal compression of H2, 0 otherwise
compressor_flag = 1


def NH3_H2(m_H2):
    '''This routine finds mass, energy and losses
    when creating NH3 from a given mass of hydrogen
    input: mass of H2, kg
    output: mass of NH3, kg; heat released, kWh; volume, m^3'''
    # recall: 3 moles of H2 produce 2 of NH3
    deltaH_NH3 = 46.1               # enthalpy of reaction, kJ/mol NH3
    m_NH3 = m_H2*2*M_NH3/(3*M_H2)   # mass of NH3, kg
    mol_kg_H2 = 1000/M_H2           # mol H2 per kg H2
    kJ_mol_NH3 = 2*deltaH_NH3/3     # heat released, kJ/mol
    Q_loss_kg_J = mol_kg_H2*kJ_mol_NH3  # heat released, kJ/kg
    Q_loss_kWh = m_H2*Q_loss_kg_J/3600  # heat released, kWh
    # volume of liquid NH3 at -33.3 C and 100 kPa, m^3
    Vol_NH3 = m_NH3/rho_NH3
    return m_NH3, Q_loss_kWh, Vol_NH3


def main(flag):
    # inputs here are for testing
    P1 = 100                                    # initial pressure, kPa
    P2 = 4000.                                 # final pressure, kPa

    # maximum practical H2 storage pressure is ~ 70,000 kPa

    T_C = 15                                     # initialtemperature, C
    T1 = T_C + 273.15                            # initial temperature, K
    # specific volume of H2 @ P1, T
    v_spec = (Ru/M_H2)*T1/P1
    Pr = 5                                      # rated power of wind turbine, MW
    cf = .5                                     # assumed capacity factor of turbine
    t = 1.                                      # time period of interest, h
    # energy available to produce and compress H2, kWh
    E_wind = cf*Pr*t*1000

    eff = 0.65  # assumed electrolysis conversion efficiency, -

    # Get results from H2 production
    m_H2, V2, T2,  E_elec, W_c = H2_2.H2_mass_work(
        P1, P2, T1, v_spec, eff, E_wind, flag)
    E_H2 = e_H2*m_H2
    # find energy in H2
    kWh_kg_H2 = 39.4
    E_H2 = m_H2*kWh_kg_H2
    # eff_H2 = E_H2/E_wind

    # find mass of NH3 and energy loss

    m_NH3, Q_loss, Vol_NH3 = NH3_H2(m_H2)
    # find energy in NH3

    kWh_kg_NH3 = 6.24
    E_NH3 = m_NH3*kWh_kg_NH3
    eff_NH3 = E_NH3/E_wind
    print()
    print(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("initial specific volume of H2: %5.2f" % v_spec, "m^3/kg")
    print("Energy from wind: %5.2f" % E_wind, "kWh")
    print(f"mass of hydrogen = {m_H2:,.3} kg")
    print(f"energy in hydrogen = {E_H2:1.1f} kWh")
    print()
    print(f"mass of ammonia = {m_NH3:1.1f} kg")
    print(f"energy in ammonia = {E_NH3:1.1f} kWh")
    print(f"energy loss = {Q_loss:1.1f} kWh")
    print(f"volume NH3 = {Vol_NH3:1.3f} m^3")
    print(f"efficiency (from wind) = {eff_NH3:,.3} -")

    print()

    ''' hydrogen values may be printed out here if desired
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
    '''


if __name__ == "__main__":
    main(compressor_flag)
