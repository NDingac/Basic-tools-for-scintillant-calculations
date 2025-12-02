import subprocess
import glob
import os
from ase.io.cp2k import read_cp2k_restart
import numpy as np
import argparse
import re
def parse_args():
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

def run_multiwfn_generate_inputs(title, k_A, k_B, k_C):
	input_data = f"""300
7
19
{k_A}
{k_B}
{k_C}
-10
0
cp2k
{title}_OT.inp
2
3
3
2
4
0
cp2k
{title}_wfn.inp
4
-9
12
450
0
0
cp2k
{title}_HSE06.inp
1
-8
-2
4
-9
10
2
0
0
q
"""

	with open(f"{title}_inpfile.log", "w") as output_file:
		process = subprocess.run(
			["Multiwfn", f"{title}_bs.inp"],
			input=input_data,
			text=True,
			stdout=output_file,
			stderr=subprocess.PIPE
		)

def submit_sbatch(script_path):
	if os.getenv('SINGLE_JOB') == '1':
		subprocess.run(["bash", script_path])
	else:
		subprocess.run(["sbatch", script_path])

def modify_inp_file(title, file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()

    has_mn_kind = any(
        line.strip().startswith("&KIND") and len(line.split()) > 1 and line.split()[1].strip() == "Mn"
        for line in lines
    )

    updated_lines = []
    current_kind = None
    magnetic_element = False

    in_coord = False
    mn_toggle = False
    coord_elem_pattern = re.compile(r'^(\s*)Mn(\s+)')

    for line in lines:
        if has_mn_kind and line.strip().startswith("&COORD"):
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
                    new_line = f"{m.group(1)}Mn_1{m.group(2)}" + line[m.end():]
                    updated_lines.append(new_line)
                else:
                    updated_lines.append(line)
                mn_toggle = not mn_toggle
                continue
            updated_lines.append(line)
            continue

        if line.strip().startswith("&KIND"):
            parts = line.split()
            current_kind = parts[1].strip() if len(parts) > 1 else None
            updated_lines.append(line)
            continue

        if current_kind and "BASIS_SET" in line:
            if current_kind == "I":
                line = line.replace("TZVP-MOLOPT-GTH-q7", "TZVP-MOLOPT-SR-GTH-q7")
            elif current_kind == "Cu":
                line = line.replace("TZVP-MOLOPT-GTH-q11", "TZVP-MOLOPT-SR-GTH-q11")
            elif current_kind == "Mn":
                line = line.replace("TZVP-MOLOPT-GTH-q15", "TZVP-MOLOPT-SR-GTH-q15")
            updated_lines.append(line)
            continue

        if current_kind and "POTENTIAL GTH-PBE" in line:
            if current_kind == "Mn":
                magnetic_element = True
                updated_lines.append(line)
                extra = []
                if file_path == f"{title}_OT.inp" or file_path == f"{title}_wfn.inp":
                    extra = [
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
                        "      BASIS_SET TZVP-MOLOPT-SR-GTH-q15\n",
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
                    ]
                elif file_path == f"{title}_HSE06.inp":
                    extra = [
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
                        "      BASIS_SET TZVP-MOLOPT-SR-GTH-q15\n",
                        "      BASIS_SET AUX_FIT admm-dzp-q15\n",
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
                    ]
                if extra:
                    updated_lines.extend(extra)
                continue

        if "&END KIND" in line:
            current_kind = None

        if "    BASIS_SET_FILE_NAME  BASIS_MOLOPT" in line:
            updated_lines.append(line)
            updated_lines.append("    BASIS_SET_FILE_NAME  BASIS_MOLOPT_UCL\n")
            continue
        
        if f"#   WFN_RESTART_FILE_NAME {title}_wfn-RESTART.wfn" in line:
            updated_lines.append(f"    WFN_RESTART_FILE_NAME {title}_OT-RESTART.wfn\n")
            continue
                        
        if f"WFN_RESTART_FILE_NAME {title}_HSE06-RESTART.wfn" in line:
            updated_lines.append(f"    WFN_RESTART_FILE_NAME {title}_wfn-RESTART.wfn\n")
            continue
                    
        if magnetic_element and "MULTIPLICITY" in line:
            updated_lines.append(line)
            updated_lines.append("    UKS\n")
            continue
        
        if  file_path == f"{title}_wfn.inp" and "#     SCF_GUESS RESTART #Use wavefunction from WFN_RESTART_FILE_NAME file as initial guess" in line:
            updated_lines.append("      SCF_GUESS RESTART #Use wavefunction from WFN_RESTART_FILE_NAME file as initial guess\n")
            continue
                
        if "MAX_MEMORY" in line and "3000" in line:
            line = line.replace("3000", "30000")
            
        updated_lines.append(line)

    with open(file_path, "w") as file:
        file.writelines(updated_lines)

def main():
    user_name = parse_args()

    bs_files = glob.glob("*_bs.inp")
    if not bs_files:
        raise FileNotFoundError("*_bs.inp file not found in the current directory.")
    title = os.path.splitext(bs_files[0])[0].split("_bs")[0]

    with open(f"{title}_bs.inp", "r") as file:
        atoms = read_cp2k_restart(file)
    cell = atoms.get_cell()
    cell_lengths = np.linalg.norm(cell, axis=1)
    k_lengths = np.array([19, 19, 19], dtype=float)
    k_points = np.ceil(k_lengths / cell_lengths).astype(int)
    k_A, k_B, k_C = k_points

    run_multiwfn_generate_inputs(title, k_A, k_B, k_C)

    modify_inp_file(title, f"{title}_OT.inp")
    modify_inp_file(title, f"{title}_wfn.inp")
    modify_inp_file(title, f"{title}_HSE06.inp")

    sbatch_script = f"""#!/bin/bash
#SBATCH -J {user_name}-{title}_wfn
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -o {title}_wfn.log
#SBATCH -e {title}_wfn.err

export cp2kroot=/home/think/app/cp2k-2024.1
source $cp2kroot/tools/toolchain/install/setup
export PATH=$cp2kroot/exe/local:$PATH
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=512m
export OMP_NUM_THREADS=2

echo The job start time at: `date`

mpirun -np 64 -map-by ppr:4:L3cache:pe=2 cp2k.psmp {title}_OT.inp > {title}_OT.out
mpirun -np 64 -map-by ppr:4:L3cache:pe=2 cp2k.psmp {title}_wfn.inp > {title}_wfn.out

export OMP_NUM_THREADS=8
mpirun -np 16 -map-by ppr:1:L3cache:pe=8 cp2k6A.psmp {title}_HSE06.inp > {title}_HSE06.out
"""

    script_path = f"{title}_wfn.sh"
    with open(script_path, "w") as script_file:
        script_file.write(sbatch_script)

    submit_sbatch(script_path)

if __name__ == "__main__":
	main()


