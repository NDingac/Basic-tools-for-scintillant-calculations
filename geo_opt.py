import glob
import os
import subprocess
import numpy as np
import argparse
import re

def get_user_name():
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
    return user_name.strip().replace(' ', '_')

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

def get_title_and_cif():
    cif_files = glob.glob("*.cif")
    if not cif_files:
        raise FileNotFoundError("cif file not found in the current directory.")
    title = os.path.splitext(cif_files[0])[0]
    return title, cif_files[0]

def get_k_points(cell):
    cell_lengths = np.linalg.norm(cell, axis=1)
    k_lengths = [19, 19, 19]
    k_points = np.ceil(np.array(k_lengths) / cell_lengths).astype(int)
    return k_points

def run_multiwfn(title, cif_file, k_A, k_B, k_C):
    input_data_opt = f"""
cp2k
{title}_opt.inp
-1
3
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

def modify_inp_file(file_path):

    with open(file_path, "r") as file:
        lines = file.readlines()

    has_mn_kind = False
    for L in lines:
        if L.strip().startswith("&KIND"):
            parts = L.split()
            if len(parts) > 1 and parts[1].strip() == "Mn":
                has_mn_kind = True
                break

    updated_lines = []
    current_kind = None
    magnetic_element = False

    in_coord = False
    mn_toggle = False
    coord_elem_pattern = re.compile(r'^(\s*)Mn(\b)(\s+)')

    for line in lines:
        if has_mn_kind:
            if line.strip().startswith("&COORD"):
                in_coord = True
                mn_toggle = False
                updated_lines.append(line)
                continue
            if in_coord:
                if line.strip().startswith("&END COORD"):
                    in_coord = False
                    updated_lines.append(line)
                    continue
                m = coord_elem_pattern.match(line)
                if m:
                    if mn_toggle:
                        new_line = f"{m.group(1)}Mn_1{m.group(3)}" + line[m.end():]
                        updated_lines.append(new_line)
                    else:
                        updated_lines.append(line)
                    mn_toggle = not mn_toggle
                    continue

        if "&KIND" in line:
            parts = line.split()
            current_kind = parts[1].strip() if len(parts) > 1 else None
            updated_lines.append(line)
            continue

        if current_kind and "POTENTIAL GTH-PBE" in line:
            if current_kind == "Mn":
                magnetic_element = True
                updated_lines.append(line)
                updated_lines.extend([
                    "      &BS\n",
                    "        &ALPHA\n",
                    "           N  3 4\n",
                    "           L  2 0\n",
                    "           NEL +5 -2\n",
                    "        &END ALPHA\n",
                    "        &BETA\n",
                    "           N  3 4\n",
                    "           L  2 0\n",
                    "           NEL -5 -2\n",
                    "        &END BETA\n",
                    "      &END BS\n",
                    "    &END KIND\n",
                    "    &KIND Mn_1\n",
                    "      ELEMENT Mn\n",
                    "      BASIS_SET DZVP-MOLOPT-SR-GTH-q15\n",
                    "      POTENTIAL GTH-PBE\n",
                    "      &BS\n",
                    "        &ALPHA\n",
                    "           N  3 4\n",
                    "           L  2 0\n",
                    "           NEL -5 -2\n",
                    "        &END ALPHA\n",
                    "        &BETA\n",
                    "           N 3 4\n",
                    "           L  2 0\n",
                    "           NEL +5 -2\n",
                    "        &END BETA\n",
                    "      &END BS\n",
                ])
                continue

        if "&END KIND" in line:
            current_kind = None

        if magnetic_element and "MULTIPLICITY" in line:
            updated_lines.append(line)
            updated_lines.append("    UKS\n")
            continue

        if "#     EXTRAPOLATION USE_PREV_P #Use converged density matrix of last geometry as initial guess" in line:
            line = "      EXTRAPOLATION USE_PREV_P #Use converged density matrix of last geometry as initial guess\n"

        updated_lines.append(line)

    with open(file_path, "w") as file:
        file.writelines(updated_lines)

def write_sbatch_script(user_name, title):
    sbatch_script = f"""#!/bin/bash
#SBATCH -J {user_name}-{title}_geoopt
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 128
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
    return script_opt_path

def submit_job(script_path):
    # If SINGLE_JOB is set, run the generated script inside the current job;
    # otherwise submit it as a separate sbatch job (original behavior).
    if os.getenv('SINGLE_JOB') == '1':
        subprocess.run(["bash", script_path])
    else:
        subprocess.run(["sbatch", script_path])

def main():
    user_name = get_user_name()
    title, cif_file = get_title_and_cif()
    cell = parse_cif_cell(cif_file)
    k_A, k_B, k_C = get_k_points(cell)
    run_multiwfn(title, cif_file, k_A, k_B, k_C)
    OPT_inp_file = f"{title}_opt.inp"
    modify_inp_file(OPT_inp_file)
    script_opt_path = write_sbatch_script(user_name, title)
    submit_job(script_opt_path)

if __name__ == "__main__":
    main()