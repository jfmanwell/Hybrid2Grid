'''
-------------------------------------------------------------------------------
Compressed air energy storage modelling
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
This is for idealized compressed air energy storage.  It assumes that air goes into 
a submerged bladder at constant pressure
It allows adiabatic, polytropic or polytropic compression
Future versions should include heat transfer and thermal storage model
This version based on CAES_final_6,py; last changes 7/1/26
"""

import numpy as np
R_air = 287.05                              # ideal constant for air, J/kg K
J_kWh = 3600000                             # conversion from J to kWh
# constant pressure heat capacity, J/kg C
Cp_air = 1005.0
# constant volume heat capacity, J/kg C
Cv_air = 718.0


def polytropic_compression(
    W_in_kWh,
    P1,
    T1,
    P2,
    n
):
    """
    Polytropic compression model.

    Works correctly for:
        n -> 1   (isothermal)
        1 < n < 1.4
        n = 1.4  (isentropic)

    PARAMETERS
    ----------
    W_in_kWh : compressor work input, kWh
    P1 : initial pressure, Pa
    T1 : initial temperature, K
    P2 : final pressure, Pa
    n : polytropic exponent
    """

    Win_J = W_in_kWh * 3.6e6

    rp = P2 / P1

    # ------------------------------------------------------
    # Near-isothermal special case
    # ------------------------------------------------------
    eps = 1e-4

    if abs(n - 1.0) < eps:

        # Isothermal compression
        T2 = T1

        w = (
            R_air * T1 * np.log(rp)
        )

        # all compressor work rejected as heat
        Qreject_J = Win_J

    else:

        # --------------------------------------------------
        # General polytropic compression
        # --------------------------------------------------
        T2 = (
            T1
            * rp**((n - 1)/n)
        )

        w = (
            (n/(n - 1))
            * R_air
            * T1
            * (
                rp**((n - 1)/n)
                - 1
            )
        )

        # compressed mass
        # (temporarily needed below)
        m_temp = Win_J / w

        # enthalpy rise
        deltaH_J = (
            m_temp
            * Cp_air
            * (T2 - T1)
        )

        # rejected heat
        Qreject_J = (
            Win_J - deltaH_J
        )

    # ------------------------------------------------------
    # Final compressed mass
    # ------------------------------------------------------
    m = Win_J / w

    return {
        "pressure_ratio": rp,
        "final_temperature_K": T2,
        "final_temperature_C": T2 - 273.15,
        "specific_work_Jpkg": w,
        "air_mass_kg": m,
        "heat_rejected_kWh": (
            Qreject_J / 3.6e6
        )
    }

    # ------------------------------------------------------
    # Bladder cooling
    # ------------------------------------------------------


def bladder_cooling(T_hot, T_ambient, cooling_fraction):
    """
    cooling_fraction = fraction of original
    temperature rise remaining after cooling.

    Example:
        cooling_fraction = 0.5

    means:

        T_new - Tamb =
            0.5*(T_hot - Tamb)
    """

    return (
        T_ambient
        + cooling_fraction * (T_hot - T_ambient)
    )


def storage_system_analysis(
    compressor_input_kWh,
    m,
    P_storage,
    T_storage,
    thermal_storage_kWh,
    thermal_storage_recovery_fraction,
    P0,
    T0,
    gamma=1.4,
    eta_turbine=1  # 0.90
):
    """
    CAES analysis with non-isentropic turbine.
    """

    # ------------------------------------------------------
    # Ideal isentropic outlet temperature
    # ------------------------------------------------------
    T_iso = (
        T_storage
        * (P0 / P_storage)**((gamma - 1)/gamma)
    )

    # ------------------------------------------------------
    # Ideal isentropic turbine work
    # ------------------------------------------------------
    w_iso_Jpkg = (
        Cp_air
        * (T_storage - T_iso)
    )

    # ------------------------------------------------------
    # Actual turbine work
    # ------------------------------------------------------
    w_actual_Jpkg = (
        eta_turbine
        * w_iso_Jpkg
    )

    # ------------------------------------------------------
    # Actual outlet temperature
    # ------------------------------------------------------
    T_actual = (
        T_storage
        - w_actual_Jpkg / Cp_air
    )

    # ------------------------------------------------------
    # Total turbine work
    # ------------------------------------------------------
    W_turbine_kWh = (
        m * w_actual_Jpkg / 3.6e6
    )

    # ------------------------------------------------------
    # Reheat required to reach ambient
    # ------------------------------------------------------
    q_reheat_Jpkg = (
        Cp_air * (T0 - T_actual)
    )

    q_reheat_Jpkg = max(0.0, q_reheat_Jpkg)

    Q_required_kWh = (
        m * q_reheat_Jpkg / 3.6e6
    )

    # ------------------------------------------------------
    # Recovered thermal storage
    # ------------------------------------------------------
    Q_recovered_kWh = (
        thermal_storage_kWh
        * thermal_storage_recovery_fraction
    )

    # ------------------------------------------------------
    # External heating required
    # ------------------------------------------------------
    Q_external_kWh = max(
        0.0,
        Q_required_kWh - Q_recovered_kWh
    )

    # ------------------------------------------------------
    # Overall efficiency
    # ------------------------------------------------------
    eta_overall = (
        W_turbine_kWh
        / (
            compressor_input_kWh
            + Q_external_kWh
        )
    )

    return {
        "isentropic_exit_temperature_K": T_iso,
        "actual_exit_temperature_K": T_actual,
        "required_reheat_kWh": Q_required_kWh,
        "external_heat_required_kWh": Q_external_kWh,
        "turbine_work_kWh": W_turbine_kWh,
        "overall_efficiency": eta_overall
    }


# ------------------------------------------------------
# COMPLETE EXAMPLE
# ------------------------------------------------------
if __name__ == "__main__":

    # ------------------------------------------------------
    # Compression conditions
    # ------------------------------------------------------
    W_in_kWh = 10000.0  # # work in, kWh
    P1 = 100e3  # initial pressure, Pa
    T1 = 293.15  # initial temperature, K
    P2 = 7000e3  # final pressure, Pa
    # superseded T_thermal_storage = 500  # K
    n = 1.   # polytropic exponent, -
    thermal_storage_recovery_fraction = 1
    # 0.5 # fraction of original temperature rise remaining after cooling.
    # 1: no cooling, 0: maximal cooling
    cooling_fraction = 1
    # ------------------------------------------------------
    # Compression
    # ------------------------------------------------------
    comp = polytropic_compression(
        W_in_kWh=W_in_kWh,
        P1=P1,
        T1=T1,
        P2=P2,
        n=n
    )

    print("\n===== Compression Results =====\n")

    for key, value in comp.items():
        print(f"{key}: {value:,.3f}")

    # ------------------------------------------------------
    # Bladder cooling
    # ------------------------------------------------------

    Tcooled = bladder_cooling(
        T_hot=comp["final_temperature_K"],
        T_ambient=T1,
        cooling_fraction=cooling_fraction
    )

    print("\n===== After Bladder Cooling =====\n")

    print(f"Cooled temperature: {Tcooled:.2f} K")

    # ------------------------------------------------------
    # Storage analysis
    # ------------------------------------------------------
    results = storage_system_analysis(
        compressor_input_kWh=W_in_kWh,
        m=comp["air_mass_kg"],
        P_storage=P2,
        T_storage=Tcooled,
        thermal_storage_kWh=(
            comp["heat_rejected_kWh"]
        ),
        thermal_storage_recovery_fraction=thermal_storage_recovery_fraction,
        P0=P1,
        T0=T1
    )

    print("\n===== Storage System Results =====\n")

    for key, value in results.items():

        if "efficiency" in key:
            print(f"{key}: {100*value:.2f} %")
        else:
            print(f"{key}: {value:,.2f}")
