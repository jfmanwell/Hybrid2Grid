'''
-------------------------------------------------------------------------------
Script for filling gaps in time series data in a plausible method
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
Created by JM on Sun Jan 18 12:15:38 2026
Gap filler, derived from <GapFiller V2 2023.xlsm> 
Good data is any data which is greater than the "missing_value"
The input file name must be specified in code
The default filename is <test_gap_fill_2.csv>
The next update should be to make the file choosable
most recent update: 6/30/26
"""


# -------------------- Configuration -------------------------

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
from datetime import datetime
import random
input_file = "test_gap_fill_2.csv"

root, ext = os.path.splitext(input_file)
if ext.lower() != ".csv":
    raise ValueError("input_file must be a .csv file")

output_file = root + "_fxd.csv"

missing_value = -10.0
Nbins = 20
# RANDOM_SEED = 1234 superseded
# np.random.seed(RANDOM_SEED) superseded

# -------------------- Input/output Utilities -------------------------


def load_data(filename):
    """Load single-column CSV"""
    data = np.loadtxt(filename, delimiter=",")
    return data.astype(float)


def save_data(filename, data):
    """Save single-column CSV with timestamp header"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with open(filename, "w") as f:
        f.write(f"# Generated on {timestamp}\n")
        np.savetxt(f, data, delimiter=",", fmt="%.2f")

# -------------------- Markov Utilities ----------------------


def build_markov_matrix(data, nbins):
    """
    Build Markov transition matrix from valid data
    """

    valid = data[data > missing_value]
    vmin, vmax = valid.min(), valid.max()
    vrange = vmax - vmin

    # Normalize and bin
    bins = np.clip(
        ((valid - vmin) / vrange * nbins).astype(int),
        0, nbins - 1
    )

    T = np.zeros((nbins, nbins))

    for i in range(len(bins) - 1):
        T[bins[i], bins[i + 1]] += 1

    # Normalize rows
    row_sums = T.sum(axis=1)
    for i in range(nbins):
        if row_sums[i] > 0:
            T[i] /= row_sums[i]
        else:
            T[i] = 1.0 / nbins

    # Cumulative wheel
    wheel = np.cumsum(T, axis=1)

    return wheel, vmin, vmax


def generate_markov_series(wheel, length, start_bin):
    """
    Generate a Markov time series of given length
    """
    # nbins = wheel.shape[0] (no longer needed)
    series_bins = np.zeros(length, dtype=int)
    series_bins[0] = start_bin

    for i in range(1, length):
        # r = np.random.rand()
        r = random.random()
        series_bins[i] = np.searchsorted(wheel[series_bins[i - 1]], r)

    return series_bins

# -------------------- Gap Filling ---------------------------


def fill_gaps(data, nbins):
    """
    Fill gaps using Markov chain interpolation
    """

    filled = data.copy()
    wheel, vmin, vmax = build_markov_matrix(data, nbins)
    vrange = vmax - vmin

    # Identify gaps
    # is_gap = data == missing_value superseded
    is_gap = data <= missing_value
    gap_indices = np.where(is_gap)[0]

    if len(gap_indices) == 0:
        return filled

    # Group contiguous gaps
    gaps = np.split(
        gap_indices,
        np.where(np.diff(gap_indices) != 1)[0] + 1
    )

    for gap in gaps:
        g0, g1 = gap[0], gap[-1]

        # Boundary values
        if g0 == 0 or g1 == len(data) - 1:
            continue

        left = filled[g0 - 1]
        right = filled[g1 + 1]

        # Determine starting bin
        start_bin = int((left - vmin) / vrange * nbins)
        start_bin = np.clip(start_bin, 0, nbins - 1)

        # Generate normalized fill
        bins = generate_markov_series(wheel, len(gap), start_bin)
        norm_fill = (bins + 0.5) / nbins

        # Scale to match endpoints
        lin = np.linspace(left, right, len(gap) + 2)[1:-1]
        scaled_fill = lin * norm_fill / norm_fill.mean()

        filled[gap] = scaled_fill

    return filled

# -------------------- Plotting ------------------------------


def plot_results(original, filled):
    """
    Plot original and filled data vs time
    """

    t = np.arange(len(original))

    plt.figure(figsize=(12, 5))
    plt.plot(t, filled, "r-", linewidth=1.5, label="Filled")
    plt.plot(t, original, "k-", linewidth=0.5, label="Original (with gaps)")

    plt.legend()
    plt.xlabel("Index")
    plt.ylabel("Value")
    plt.title("Markov Chain Gap Filling")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# -------------------- Main function ---------------------------


def main():
    if not Path(input_file).exists():
        raise FileNotFoundError(input_file)

    data = load_data(input_file)
    filled = fill_gaps(data, Nbins)

    save_data(output_file, filled)
    plot_results(data, filled)

    print()
    print("Gap filling complete")
    print(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Input file: {input_file}")
    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    main()
