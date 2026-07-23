'''
-------------------------------------------------------------------------------
Script for synthesizing time series data
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
Created on Wed Jan 14 16:30:52 2026
This script provides multiple options for time series data synthesis
One of them uses McNerney's method
This code was adapted and expanded from Hybrid2's data synthesizer
Solar data synthesizer now improved over original in Hybrid2
most recent update: 6/30/26
"""


import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import os
import time
import random
from scipy.stats import norm
from datetime import datetime
pi = np.pi                          # pi
Gsc = 1367                          # solar constant, W/m^2
deg2rad = pi / 180                  # conversion from degrees to radians
bin_i, pdf_2, x = [], [], []

output_file_csv = "synth_data.csv"

# flag =1 for wind, 2 for solar, 3 for shifted Weibull, 4 for normal, 5 for log normal
flag = 2
target_mean = 10                    # overwritten for solar
target_std = 5                      # overwritten for solar
diurnal_ratio = 1.                  # not used for solar
diurnal_hour_max = 12               # not used for solar
long_period = 8760                  # not used for solar
long_ratio = 1                      # not used for solar
long_hour_max = 4300                # not used for solar
max_val = 50                        # overwritten for solar
min_val = 0                         # overwritten for solar
target_max = 50                     # overwritten for solar
target_min = 0                      # overwritten for solar

# used primarily for load offset, not used for solar
offset = target_min

# bin width for pdf; >= 1 for wind, load; < 1 for solar, overwritten for solar
bin_width = 1
first_day = 0                       # first Julian day for data; most important for solar
last_day = 364                      # last Julian day for data; most important for solar
# number of points, initially assume hourly
N_points = (last_day - first_day+1) * 24
autocor = .95                       # target autocorrelation at lag number lagN
lagN = 1                            # target autocorrelation lag number lagN

kt_av = .5                          # average clearness index, only used for solar
latitude = 45                       # latitude (degrees) only used for solar

if flag == 2:                       # this is for clearness index which is less than 1.0
    target_mean = 1                 # overwritten for solar
    target_std = 0                  # overwritten for solar
    bin_width = .05                 # overwritten for solar
    max_val = 1                     # overwritten for solar
    min_val = 0                     # overwritten for solar
    target_max = 1                  # overwritten for solar
    target_min = 0                  # overwritten for solar


def control_para(pdf_flag, autocor):
    '''This returns initial values for starting synthesis
    adapted from Gannong Deng's Basic program'''
    if pdf_flag in (1, 2, 3):
        b1, b2 = 1.0, 3.0
        dR = autocor
        pass1 = 0.02
    return b1, b2, dR, pass1


def correl2(data, mean, std, lag_start, lag_end):
    '''This performs an autocorrelation up to a specified lag number

    'Inputs
    data      : Array of data
    mean  : mean of synthesized data
    std   : standard deviation of synthesized data
    lag_start : start calculate autocorrelation with this lag
    lag_end   : stop calculate autocorrelation with this lag

    Output
    aCorr    : Array of autocorrelation value; index starts with LagStart to LagEnd'''

    n = len(data)  # number of data points to calculate autocorrelation
    aCorr = np.zeros(lag_end + 1)

    for lag in range(lag_start, lag_end + 1):
        s = 0.0
        for i in range(n - lag):
            s += (data[i] - mean) * (data[i + lag] - mean)
        aCorr[lag] = s / ((n - lag) * std**2)

    return aCorr


def decay_mat(B, pdfN):
    i = np.arange(pdfN + 1)
    return B ** (-np.abs(i[:, None] - i[None, :]))


def get_norm_r(decay, p0):
    pdfN = len(p0) - 1
    r = np.zeros(pdfN + 1)

    for i in range(pdfN + 1):
        for j in range(pdfN + 1):
            r[i] += decay[i, j] * p0[i] * p0[j]

    r /= np.sum(r)
    return r


def get_cumul(decay, p0):
    '''This generates the predicted cumulative pdf ("wheel of fortune")
    from the starting vector and decay matrix
    Adapted from Gannong Deng's Basic code

    Inputs
    decay = decay matrix
    p0  = starting vector

    Outputs
    cumul =  cumulative TPM'''

    pdfN = len(p0) - 1  # dimension of pdf
    # print("pdfN",pdfN)
    den = decay @ p0
    tpm = np.zeros((pdfN + 1, pdfN + 1))
    cumul = np.zeros_like(tpm)

    for i in range(pdfN + 1):
        for j in range(pdfN + 1):
            if den[i] > 0:
                tpm[i, j] = decay[i, j] * p0[j] / den[i]
            cumul[i, j] = tpm[i, j] if j == 0 else cumul[i, j-1] + tpm[i, j]
            # if i ==10:
            # print("cumul[i, j]",i/(pdfN+1),j/(pdfN+1),cumul[i, j])
    return cumul


def initial_p(decay, pdfN, target_pdf):
    relax = 0.5
    diff = 0.001
    itermax = 2 * pdfN

    # Initialize p0 with target pdf
    p0 = target_pdf.copy()

    for _ in range(itermax + 1):

        r = get_norm_r(decay, p0)
        grad = target_pdf - r
        p0 = p0 + relax * grad
        outNrm = np.sqrt(np.sum(grad**2))
        if outNrm <= diff:
            break

    return p0


def get_pdf(flag, offset, pdfN, kt_av, bin_width):
    pdf = np.zeros(pdfN + 1)
    bin_j = np.zeros(pdfN + 1)
    total = 0.0

    if flag == 1:   # Wind (Weibull/Rayleigh)
        if target_std == -1:
            # Rayleigh
            k = 2.0
            c = 2 * target_mean / np.sqrt(pi)
        else:
            # Weibull
            k = (target_std / (target_mean - offset)) ** -1.086
            c = (target_mean-offset)/math.gamma(1+1/k)

        for j in range(pdfN + 1):
            wd = (j * bin_width - offset) / c

            if wd >= 0:
                pdf[j] = (k / c) * wd**(k - 1) * np.exp(-wd**k) * bin_width
            total += pdf[j]
            bin_j[j] = j

        if total < 1:
            pdf[pdfN] += (1 - total)

    elif flag == 2:
        # this is for solar data via the clearness index
        dum = clearness_frequency(kt_av, pdfN)
        pdf = dum
        # for i in range(pdfN):
        #   print(i/pdfN,pdf[i])

    elif flag == 3:  # Load (shifted Rayleigh)
        sm = target_mean - target_min
        k = 0.5 * pi / sm**2

        for j in range(pdfN + 1):
            ld = j * bin_width
            pdf[j] = k * ld * np.exp(-0.5 * k * ld**2) * bin_width
            total += pdf[j]
            bin_j[j] = j

        if total < 1:
            pdf[0] += (1 - total)

    return pdf, bin_j


def mean(x):
    return np.mean(x)


def stdev(x):
    return np.std(x, ddof=1)


def mtimesv(mat, v):
    return mat @ v


def search_table(cumul, i, r):
    jl, ju = -1, cumul.shape[1]

    while ju - jl > 1:
        j = (ju + jl) // 2
        if r > cumul[i, j]:
            jl = j
        else:
            ju = j
    return jl


def rwalk(bin_width, min_val, target_mean, cumul, num_points):
    pdfN = cumul.shape[0] - 1
    tally = np.zeros(pdfN + 1)
    data = np.zeros(num_points)

    i1 = 1
    for i in range(num_points):
        r = np.random.rand()

        j2 = search_table(cumul, i1, r)

        i1 = j2 + 1
        tally[i1] += 1
        data[i] = (i1 + (np.random.rand() - 0.5)) * \
            bin_width  # + target_mean#min_val

    return data


def sine_scale(period, data, mag, start_point, dt):
    n = len(data)
    for i in range(n):
        hr = (i + 1) * dt
        data[i] *= 1 + (mag - 1) * np.cos(2 * pi * (hr - start_point) / period)
        if data[i] < 0:
            data[i] = 0


def hr_ang(hr, dt):
    return (hr - 12 - dt / 2) * 15


def get_time_and_decl(counter, first_day, time_step, state):
    if state["first"]:
        state["hour"] = time_step
        state["day"] = first_day
        state["first"] = False

    decl = 23.45 * np.sin(0.01721 * (284 + state["day"])) * deg2rad

    state["hour"] += time_step
    if state["hour"] > 24:
        state["hour"] -= 24
        state["day"] += 1

    hour_angle = hr_ang(state["hour"], time_step) * deg2rad
    return state["day"], decl, hour_angle


'''def write_time_series(filename, data, header=None):
    not used now, save for later
    with open(filename, "w") as f:
        if header:
            f.write(header + "\n")
        for i, val in enumerate(data, start=1):
            f.write(f"{i:6d}  {val:14.6f}\n")'''


def clearness_frequency(kt_av, n_bins):
    """
    This provides a relation between hourly
    expected frequency of solar radiation as a function
    of monthly average clearness index.  It is based on
    Bendt's correlations of the Liu and Jordan method
    as described in Beckman and Duffie, 1991, p. 79
    Originally written JM 12/24/98

    Inputs
     kt_av = average clearness index, 0-1
     n_bins     = number of bins

    Output
     kt_freq()  = pdf of clearness index
    """

    # Locals
    # gamma, eta = terms used by Huget in estimating
    # kTMax = maximum expected kT
    # ktMin = minimum expected kT
    # kt_av = bin midpoint clearness index
    # binWidth = pdf bin width
    # lowBin  = lower value of bin
    # hiBin   = upper value of bin
    # tot = cumulative probability, for verifying
    # dumktav = calculated value of kTaverage, based on estimate of gamma
    # i = bin number

    kt_freq = np.zeros(n_bins)

    kt_min = 0.05  # proposed by Bendt et al.

    # ktMax proposed by Hollands and Huget, 1983
    kt_max = 0.6313 + 0.267 * kt_av - 11.9 * (kt_av - 0.75) ** 8

    bin_width = 1.0 / n_bins

    # gamma proposed by Huget, 1985
    eta = (kt_max - kt_min) / (kt_max - kt_av)
    gamma = (
        -1.498
        + (1.184 * eta - 27.182 * np.exp(-1.5 * eta))
        / (kt_max - kt_min)
    )

    # Search below improves estimate of gamma
    dumktav = 0.0
    while abs(dumktav - kt_av) > 0.01:
        dumktav = (
            ((kt_min - 1 / gamma) * np.exp(gamma * kt_min)
             - (kt_max - 1 / gamma) * np.exp(gamma * kt_max))
            / (np.exp(gamma * kt_min) - np.exp(gamma * kt_max))
        )

        if dumktav < kt_av:
            gamma += 0.1
        else:
            gamma -= 0.1

    tot = 0.0

    for i in range(n_bins):
        # kt_av = bin midpoint
        kt_av = (i + 1) / n_bins - bin_width / 2

        if kt_av - bin_width / 2 < kt_min:
            kt_freq[i] = 0.0
        else:
            low_bin = kt_av - bin_width / 2
            hi_bin = kt_av + bin_width / 2

            # cumulative for testing
            kt_freq[i] = (
                kt_cumulative(kt_min, kt_max, gamma, hi_bin)
                - kt_cumulative(kt_min, kt_max, gamma, low_bin)
            )

        tot += kt_freq[i]

    # Make sure frequencies sum to 1

    if tot < 1.0:
        kt_freq[0] += 1.0 - tot
    # print("tot",tot)

    return kt_freq


def kt_cumulative(kt_min, kt_max, gamma, kt_av):
    """
    This provides the cumulative distribution of the
    clearness index between ktMin and kTMax

    Inputs:
     kt_min =
     kt_max =
     gamma  =
     kt_av  =
    """

    dumktC = (
        np.exp(gamma * kt_min) - np.exp(gamma * kt_av)
    ) / (
        np.exp(gamma * kt_min) - np.exp(gamma * kt_max)
    )

    if dumktC > 1.0:
        dumktC = 1.0

    return dumktC


def declination(day):
    # This finds the solar declination angle in degrees
    # declination angle, deg
    decl_angle = 23.45 * np.sin(2 * pi * (284 + day) / 365)
    return decl_angle


def hourAngle(hour):
    # This returns hour angle in degrees, given the hour
    hour_angle = float(15 * (hour - 12))
    return hour_angle                                           # hour angle, deg


def incidence(latitude, day, hour, slope, azimuth_w):
    """This finds solar incidence angle in degrees
    as a function of latitude,day,hour, slope and wall surface azimuth"""

    """#alternate method, which needs to be double checked but can be ignored
    zenith_ = zenith(latitude,day,hour)                         # zenith angle, deg
    azimuth_s = SolarAzimuth(latitude,day,hour)                 # solar azimuth, deg
    cos_theta = np.cos(zenith_*deg2rad)*np.cos(slope*deg2rad)\
        + np.sin(zenith_*deg2rad)*np.sin(slope*deg2rad)\
            *np.sin((azimuth_s-azimuth_w)*deg2rad) # cosine of incidence angle
    """

    # declination, deg
    delta = declination(day)
    # hour angle, deg
    omega = hourAngle(hour)

    # format of cos_theta equation below is from Duffie  Beckman, 1991
    cos_theta = np.sin(delta*deg2rad) * np.sin(latitude*deg2rad) * np.cos(slope*deg2rad)\
        - np.sin(delta*deg2rad) * np.cos(latitude*deg2rad) * np.sin(slope*deg2rad) * np.cos(azimuth_w*deg2rad)\
        + np.cos(delta*deg2rad) * np.cos(latitude*deg2rad) * np.cos(slope*deg2rad) * np.cos(omega*deg2rad)\
        + np.cos(delta*deg2rad)*np.sin(latitude*deg2rad) * np.sin(slope*deg2rad) * np.cos(azimuth_w*deg2rad) * np.cos(omega*deg2rad)\
        + np.cos(delta*deg2rad) * np.sin(slope*deg2rad) * \
        np.sin(azimuth_w*deg2rad) * np.sin(omega*deg2rad)

    theta = np.arccos(cos_theta)/deg2rad        # incidence angle, deg

    return theta


def ExtraTer(latitude, day, hour):
    '''
    This finds the extraterrestrial radiation on a horizontal surface
    Input:
        latitude = latitude, degrees
        day = day of year
        hour = hour of day
     '''

    # horizontal surface assumed
    azimuth_w = 0
    # horizontal surface assumed
    slope = 0

    # incidence angle, degrees
    theta = incidence(latitude, day, hour, slope, azimuth_w)

    # incidence angle cannot be greater than zero (below the surface)
    if theta >= 90:
        G0 = 0
    else:
        G0 = Gsc * (1 + 0.033 * np.cos(pi ** 2 * day / 365)) * \
            np.cos(theta*deg2rad)
    if G0 < 0:
        G0 = 0
    return G0


'''def write_time_series_csv(filename, data, header=None):
   !!! replaced by safe_write_time_series_csv !!!
    with open(filename, mode="w", newline="") as f:
        if header:
            for line in header.split("\n"):
                f.write(f"{line}\n")

        writer = csv.writer(f)
        writer.writerow(["Index", "Value"])

        for i, val in enumerate(data, start=1):
            writer.writerow([i, val])'''


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

    for attempt in range(retries):
        try:
            with open(filename, "w", newline="") as f:
                if header:
                    for line in header.split("\n"):
                        f.write(f"# {line}\n")

                writer = csv.writer(f)
                writer.writerow(["Index", "Value"])

                # for i, val in enumerate(data, start=1):
                # writer.writerow([i, val])
                for i, val in enumerate(data, start=1):
                    writer.writerow([i, f"{val:.2f}"])

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


def time_series_main(flag, day1, offset, target_pdf, num_points, bin_width, time_step, output_file_csv):
    t = np.zeros(num_points)  # for tracking time
    data = np.zeros(num_points)

    if flag < 4:  # wind, solar or load
        # when solar data is to be created, first make time series of extraterrestrial radiation
        # assume day starts at midnight
        # initialize extraterrestrial radiation
        G0 = np.zeros(num_points)
        # I = np.zeros(num_points)                # initialize solar radiation (no longer needed)
        i = -1                                  # counter for hour in time period
        # total number of days in time period
        days = int(num_points/24)
        for j in range(days):
            day = j+day1
            for hour in range(24):
                i += 1
                G0[i] = ExtraTer(latitude, day, hour)

        pdfN = len(target_pdf) - 1

        # ---- Control parameters (from ControlPara) ----
        b1, b2, dR, pass1 = control_para(flag, autocor)

        itermax = 20

        # ---- Iterative Markov synthesis ----
        for _ in range(itermax):
            B = 0.5 * (b1 + b2)

            decay = decay_mat(B, pdfN)
            p0 = initial_p(decay, pdfN, target_pdf)
            cumul = get_cumul(decay, p0)

            data = rwalk(bin_width, min_val, target_mean, cumul, num_points)

            # initial mean and st dev of synthetic data
            synth_mean = mean(data)
            synth_std = stdev(data)

            aCorr = correl2(data, synth_mean, synth_std, 1, lagN)
            sDr = aCorr[lagN]

            if abs(sDr - dR) <= pass1:
                break

            if sDr < dR:
                b1 = B
            else:
                b2 = B
        for i in range(num_points):
            t[i] = i

        # ---- Adjust mean & standard deviation ----
        # print("synth_mean",synth_mean)
        if flag == 1:
            zero = data - synth_mean
            if target_std > 0:
                zero *= target_std / synth_std
            data = target_mean + zero  # + offset
        else:
            data = data - synth_mean + kt_av

        data[data < 0] = 0

        # ---- Diurnal & long-term modulation ----
        if flag == 1 or flag == 3:
            if diurnal_ratio != 1:
                sine_scale(24, data, diurnal_ratio, diurnal_hour_max, 1)

            if long_ratio != 1:
                sine_scale(long_period, data, long_ratio, long_hour_max, 1)

        if flag == 2:
            # solar synthesis here, data so far is clearness index
            # first find extraterrestrial radiation
            # replace original data output (kt) with estimate of I
            data = data*G0

    elif flag == 4:
        # this is for uncorrelated normal time series
        for i in range(num_points):
            y = random.random()
            x_temp = target_std*norm.ppf(y) + target_mean
            t[i] = i
            data[i] = x_temp
    else:
        # this is for uncorrelated log normal time series
        sigma2 = np.log(1 + (target_std/target_mean)**2)
        av = np.log(target_mean) - 0.5*sigma2
        std = np.sqrt(sigma2)

        for i in range(num_points):
            y = random.random()
            x_temp = std*norm.ppf(y) + av
            t[i] = i
            data[i] = np.exp(x_temp)
        # data = np.exp(x)

    # final mean and st dev of synthesized data
    synth_mean = mean(data)
    synth_std = stdev(data)

   # ---- Optional file output ----
    if output_file_csv:
        summary = f"Synthesized Time Series | Mean = {
            mean(data):.2f}; Std = {stdev(data):.2f}"
        safe_write_time_series_csv(
            output_file_csv,
            data,
            header=summary,
            overwrite=False,      # protect existing file
            auto_rename=True      # timeseries_1.csv, timeseries_2.csv, ...
        )

    '''if output_file:
        not used now, save for later
        header = f"# Synthesized Time Series | Mean = {mean(data):.2f}, Std = {stdev(data):.2f}"
        write_time_series(output_file, data, header)'''

    return t, data, synth_mean, synth_std


def main(flag):
    print("")
    if flag == 1:
        print("Synthetic wind data")
    elif flag == 2:
        print("Synthetic solar data")
    elif flag == 3:
        print("Synthetic load data")
    elif flag == 4:
        print("Synthetic normal data (uncorrelated)")
    else:
        print("Synthetic log normal data (uncorrelated)")
    print(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    N_bins = int((max_val - min_val)/bin_width)
    # print("N_bins",N_bins)

    # Get the target pdf
    results_1 = get_pdf(flag=flag, offset=offset, pdfN=N_bins,
                        kt_av=kt_av, bin_width=bin_width)
    target_pdf = results_1[0]

    # Get the time series
    results = time_series_main(flag=flag, day1=first_day, offset=offset, target_pdf=target_pdf,
                               num_points=N_points, bin_width=bin_width, time_step=1,
                               output_file_csv=output_file_csv)

    t = results[0]                                      # this is time
    data = results[1]                                  # data time series

    synth_mean = results[2]
    synth_std = results[3]

    print("")
    print("data average: %5.2f" % synth_mean)
    print("data st dev: %5.2f" % synth_std)

    plt.figure
    plt.plot(t, data, label="Synthetic time series)")
    plt.xlabel("Time")
    plt.ylabel("Value")
    plt.legend()
    plt.grid()
    # print("N_bins",N_bins)

    plt.figure()
    sns.histplot(data, bins=N_bins, kde=True,
                 color='lightgreen', edgecolor='red')
    # sns.distplot(data, kde=False, norm_hist=True)
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.title('Histogram')

    plt.show()


if __name__ == "__main__":
    main(flag)
