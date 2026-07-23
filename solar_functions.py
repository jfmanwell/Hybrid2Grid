'''
-------------------------------------------------------------------------------
Solar functions
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
Solar angles
updated 4/7/26
subject to additional updates as needed
"""

import numpy as np
pi = np.pi
deg2rad = pi/180.
Gsc = 1367.    # solar constant, W/m^2


def declination(day):
    # This finds the solar declination angle in degrees
    # declination angle, deg
    decl_angle = 23.45 * np.sin(2 * pi * (284 + day) / 365)
    return decl_angle


def extraTerrestrial(latitude, day, hour):
    # this finds the extraterrestrial radiation, W/m^2
    G_ext = Gsc*(1 + 0.033*np.cos(day*2*pi/365))\
        * np.cos(zenith(latitude, day, hour)*deg2rad)
    return G_ext


def clearness(I, latitude, day, hour):
    """This calculates the clearness index from the measured data 
    and the extraterrestrial radiation for the corresponding
    day and time, both on horizontal"""
    G_ext = extraTerrestrial(latitude, day, hour)
    if G_ext > 0:
        kT = I/G_ext
        if kT > 1:
            kT = 1
    else:
        kT = 0

    return kT                                                   # clearness index


def diffuse(kT, I):
    """This is for diffuse radiation, given total radiation on the horizontal
    The main inputs are the clearness index kT and radiaion I (W/m^2)
    Id_I is diffuse fraction.  The method is from Erbs"""

    if kT < 0.22:
        Id_I = 1 - 0.09*kT

    elif kT <= 0.8:
        Id_I = 0.9511 - 0.1604*kT + 4.39 * kT**2 - 16.64 * kT**3\
            + 12.336 * kT**4
    else:
        Id_I = 0.165

    # diffuse radiation, W/m^2
    Id = Id_I*I

    return Id, Id_I


def hourAngle(hour):
    # This returns the hour angle in degrees, given the hour
    hour_angle = float(15 * (hour - 12))
    return hour_angle             # hour angle, deg


def zenith(latitude, day, hour):
    # This finds solar zenith angle in degrees
    # It must be equal to or less than 90 deg
    delta = declination(day)
    hour_angle = hourAngle(hour)
    cos_zenith = np.cos(deg2rad*delta) * np.cos(deg2rad*hour_angle)\
        * np.cos(deg2rad*latitude)+np.sin(deg2rad*delta)\
        * np.sin(deg2rad*latitude)                          # cos of zenith angle

    if cos_zenith >= 0:
        zenith_Angle = np.arccos(cos_zenith)/deg2rad
    else:
        zenith_Angle = 90.
    return zenith_Angle                                         # zenith angle, deg


def SolarAzimuth(latitude, day, hour):
    """This finds the solar azimuth in degrees
    The method is from Duffie and Beckman, 2nd. ed, 1991
    Angles initially converted to radians for convenience"""

    # declination in radians
    decl_r = declination(day)*deg2rad
    # latitude in radians
    phi_r = latitude*deg2rad
    # hour angle in radians
    omega_r = hourAngle(hour)*deg2rad
    # zenith angle in radians
    theta_z_r = zenith(latitude, day, hour)*deg2rad
    # hour sun is directly E or W
    omega_ew_r = np.arccos(np.tan(decl_r)/np.tan(phi_r))

    if abs(omega_r) <= omega_ew_r:
        C1 = 1
    else:
        C1 = -1
    if phi_r - decl_r >= 0:
        C2 = 1
    else:
        C2 = -1
    if omega_r >= 0:
        C3 = 1
    else:
        C3 = -1
    gamma_p_r = np.arcsin(np.sin(omega_r)
                          * np.cos(decl_r)/np.sin(theta_z_r))

    azimuth = (C1*C2*gamma_p_r + C3*(1-C1*C2)*pi/2) / \
        deg2rad     # solar azimuth angle, deg

    return azimuth


def sunSet(latitude, day):
    """This finds the sunset hour angle
    Angles initially converted to radians for convenience"""

    # declination in radians
    decl_r = declination(day)*deg2rad
    # latitude in radians
    phi_r = latitude*deg2rad

    if abs(np.tan(decl_r)*np.tan(phi_r)) <= 1:
        hour_angle_sunset = np.arccos(-np.tan(decl_r) * np.tan(phi_r))/deg2rad
        # sunset hour angle, deg
        hour_sunset = hour_angle_sunset/15 + 12
    else:
        hour_sunset = 0.                                        # sunset hpur

    return hour_sunset


def incidence(latitude, day, hour, slope, azimuth_w):
    """This finds solar incidence angle in degrees
    as a function of latitude,day,hour, slope and wall surface azimuth azimuth_w"""

    """# An alternate method, which needs to be double checked but can be ignored:
    zenith_ = zenith(latitude,day,hour)                         # zenith angle, deg
    azimuth_s = SolarAzimuth(latitude,day,hour)                 # solar azimuth, deg
    cos_theta = np.cos(zenith_*deg2rad)*np.cos(slope*deg2rad)\
        + np.sin(zenith_*deg2rad)*np.sin(slope*deg2rad)\
            *np.sin((azimuth_s-azimuth_w)*deg2rad) # cosine of incidence angle
    """

    delta = declination(day)  # declination, deg
    omega = hourAngle(hour)  # hour angle, deg

    # format of cos_theta equation below is from Duffie  Beckman, 1991
    cos_theta = np.sin(delta*deg2rad) * np.sin(latitude*deg2rad) * np.cos(slope*deg2rad)\
        - np.sin(delta*deg2rad) * np.cos(latitude*deg2rad) * np.sin(slope*deg2rad) * np.cos(azimuth_w*deg2rad)\
        + np.cos(delta*deg2rad) * np.cos(latitude*deg2rad) * np.cos(slope*deg2rad) * np.cos(omega*deg2rad)\
        + np.cos(delta*deg2rad)*np.sin(latitude*deg2rad) * np.sin(slope*deg2rad) * np.cos(azimuth_w*deg2rad) * np.cos(omega*deg2rad)\
        + np.cos(delta*deg2rad) * np.sin(slope*deg2rad) * \
        np.sin(azimuth_w*deg2rad) * np.sin(omega*deg2rad)

    theta = np.arccos(cos_theta)/deg2rad  # incidence angle, deg

    return theta


def SunOnSlope(I, latitude, day, hour, slope, azimuth_w, reflec):
    """This finds the solar radiation on a sloped surface, which need not
    be south facing, given the inputs"""

    if I > 0:
        # radation must be > 0
        cos_theta_1 = np.cos(
            incidence(latitude, day, hour, slope, azimuth_w)*deg2rad)
        cos_theta_z = np.cos(zenith(latitude, day, hour)*deg2rad)

        # clearness index
        kT = clearness(I, latitude, day, hour)

        # diffuse radiation
        Id = diffuse(kT, I)[0]
        I_B = I - Id                                        # beam radiation
        # the zenith angle and incidence angle must both be > 0 to avoid difficulties
        if cos_theta_1 > 0 and cos_theta_z > 0:
            Rb = cos_theta_1/cos_theta_z                    # beam radiation ratio

            if Rb > 10:
                # ensure that Rb is not unrealistically large
                Rb = 10.

            I_slope = I_B*Rb + (1 + np.cos(slope*deg2rad))*Id/2\
                + reflec*(1 - np.cos(slope*deg2rad))*I / \
                2               # radiation on the sloped surface, W/m^2
        else:
            Rb = 0.
            I_slope = Id
    else:
        I_slope = 0.

    return I_slope, I_B, Id, Rb


"""------------------------------------------
# Main program for testing
"""


def main():
    latitude = 40.
    day = 180.
    hour = 17.
    slope = 40.
    azimuth_w = 0
    reflec = 0.2
    I_test = 500
    # make sure I_test is not impossibly big
    G0 = extraTerrestrial(latitude, day, hour)
    if I_test > G0:
        print()
        print(f"I too high, set to G0 = {G0:.1f} W/m^2")
        I_test = G0

    theta = incidence(latitude, day, hour, slope, azimuth_w)
    Zenith = zenith(latitude, day, hour)
    sun = SunOnSlope(I_test, latitude, day,
                     hour, slope, azimuth_w, reflec)
    I_slope = sun[0]
    I_B = sun[1]
    Id = sun[2]
    Rb = sun[3]
    Sunset_hour = sunSet(latitude, day)
    print("")
    print(f"theta = {theta:.1f} deg")
    if theta > 90:
        print('sunlight is from behind panel')
    print(f"zenith = {Zenith:.3} deg")
    print(f"I beam = {I_B:.1f} W/m^2")
    print(f"I diffuse = {Id:.1f} W/m^2")
    print(f"Rb = {Rb:.3f} ")
    print(f"I_slope = {I_slope:.1f} W/m^2")
    print(f"sunset hour = {Sunset_hour:.3} deg")


if __name__ == "__main__":
    main()
