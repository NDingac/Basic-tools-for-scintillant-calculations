import glob
import os
import subprocess
import numpy as np
import argparse

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('-name', '-n', dest='user_name', help='User name for job name', type=str)
args, _ = parser.parse_known_args()
user_name = args.user_name
if not user_name or not user_name.strip():
    try:
        user_name = input("user name please: ").strip()
        while not user_name:
            user_name = input("user name can not be empty, input again: ").strip()
    except EOFError:
        raise SystemExit("user name not be provided, program exited.")
user_name = user_name.strip().replace(' ', '_')

def parse_cif_cell(filename):
    def clean_value(value):
        return float(value.split('(')[0])

    cell = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip().startswith("_cell_length_a"):
                a = clean_value(line.split()[-1])
            elif line.strip().startswith("_cell_length_b"):
                b = clean_value(line.split()[-1])
            elif line.strip().startswith("_cell_length_c"):
                c = clean_value(line.split()[-1])
            elif line.strip().startswith("_cell_angle_alpha"):
                alpha = clean_value(line.split()[-1])
            elif line.strip().startswith("_cell_angle_beta"):
                beta = clean_value(line.split()[-1])
            elif line.strip().startswith("_cell_angle_gamma"):
                gamma = clean_value(line.split()[-1])
        alpha, beta, gamma = np.radians([alpha, beta, gamma])
        cell = [
            [a, 0, 0],
            [b * np.cos(gamma), b * np.sin(gamma), 0],
            [
                c * np.cos(beta),
                c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
                c * np.sqrt(1 - np.cos(beta)**2 - ((np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma))**2),
            ],
        ]
    return np.array(cell)

cif_files = glob.glob("*.cif")
if not cif_files:
    raise FileNotFoundError("cif file not found in the current directory.")
title = os.path.splitext(cif_files[0])[0]

cell = parse_cif_cell(cif_files[0])
cell_lengths = np.linalg.norm(cell, axis=1)
k_lengths = [19, 19, 19]
k_points = np.ceil(k_lengths / cell_lengths).astype(int)
k_A, k_B, k_C = k_points


input_data_opt = f"""cp2k
{title}_opt.inp
-1
4
10
1
3
2
4
8
{k_A},{k_B},{k_C}
0
q
"""

with open(f"{title}_inpfile.log", "w") as output_file:
    process = subprocess.run(
        ["Multiwfn", f"{title}.cif"],
        input=input_data_opt,
        text=True,
        stdout=output_file,
        stderr=subprocess.PIPE
    )

OPT_inp_file = f"{title}_opt.inp"


with open(OPT_inp_file, "r") as file:
    lines = file.readlines()
updated_lines = []

for line in lines:
    if "#     EXTRAPOLATION USE_PREV_P #Use converged density matrix of last geometry as initial guess" in line:
        line = "      EXTRAPOLATION USE_PREV_P #Use converged density matrix of last geometry as initial guess\n"
    updated_lines.append(line)

with open(OPT_inp_file, "w") as file:
    file.writelines(updated_lines)


sbatch_script = f"""#!/bin/bash
#SBATCH -J {user_name}-{title}_cellopt
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH -o {title}_opt.log
#SBATCH -e {title}_opt.err

export cp2kroot=/home/think/app/cp2k-2024.1
source $cp2kroot/tools/toolchain/install/setup
export PATH=$cp2kroot/exe/local:$PATH

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=512m
export OMP_NUM_THREADS=2

echo The job start time at: `date`

mpirun -np 64 -map-by ppr:4:L3cache:pe=2 cp2k.psmp {title}_opt.inp > {title}_opt.out

echo The job end time at: `date`
"""

script_opt_path = f"{title}_opt.sh"
with open(script_opt_path, "w") as script_file:
    script_file.write(sbatch_script)
subprocess.run(["sbatch", script_opt_path])