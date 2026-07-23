'''
-------------------------------------------------------------------------------
Time series block averaging
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
Created by JM on Sun Jan 18 22:52:44 2026
Block average, based on old xlsm version
"""

# ============================================================
# -------------------- MAIN INPUTS-----------------------
# ============================================================

import numpy as np
from pathlib import Path
import csv
input_file = "test_block.csv"

# Equivalent to Excel named ranges
step_x_in = 1.0      # Input time step
step_x_out = 24.0     # Output (block) time step

yes_no_stdev = "No"   # "No" or "Yes" for whether the file includes the standard deviation

# -------------------- Input/output  -------------------------


def make_output_filename(input_file):
    p = Path(input_file)
    return p.with_name(f"{p.stem}_fxd{p.suffix}")


def load_input_data(filename):
    """
    Load single-column CSV data
    """
    data = np.loadtxt(filename, delimiter=",")
    return data.astype(float)


def save_output_data(filename, data):
    """
    Save block-averaged data to CSV
    """
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["index", "block_average"])
        for i, val in enumerate(data, start=1):
            writer.writerow([i, f"{val:.2f}"])


def block_average(x_in, step_x_in, step_x_out, yes_no_stddev="No"):
    """
    This procedure is used to block average a set of data
    on trial basis is also a calculation for standard deviation.
    Adapted from VBA version (JM 9/1/09), converted to Python 1/22/26
    """

    n_points_xin = len(x_in)

    # Number of points per block
    block_points = int(step_x_out / step_x_in)

    # Number of output blocks
    n_points_xout = int(n_points_xin / block_points)

    if yes_no_stddev == "No":
        x_out = np.zeros(n_points_xout)
    else:
        x_out = np.zeros((n_points_xout, 2))

    # Perform block averaging
    for i in range(n_points_xout):
        start = i * block_points
        end = start + block_points
        block = x_in[start:end]

        if yes_no_stddev == "No":
            x_out[i] = block.mean()
        else:
            mean_val = block[:, 0].mean()
            std_val = np.sqrt(
                (np.sum(block[:, 1] ** 2)
                 + np.sum(block[:, 0] ** 2)
                 - (np.sum(block[:, 0]) ** 2) / block_points)
                / block_points
            )
            x_out[i, 0] = mean_val
            x_out[i, 1] = std_val

    return x_out

# -------------------- Main ---------------------------------


def main():
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    output_path = make_output_filename(input_path)

    # Load input data
    x_in = load_input_data(input_path)

    # Run block average
    x_out = block_average(
        x_in,
        step_x_in,
        step_x_out,
        yes_no_stdev
    )

    # Save results; note that overwriting previous file is allowed here
    save_output_data(output_path, x_out)
    print()
    print("Block averaging complete.")
    print(f"Input file : {input_path}")
    print(f"Output file: {output_path}")


if __name__ == "__main__":
    main()
