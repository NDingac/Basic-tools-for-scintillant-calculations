import glob
import os
from seekpath import get_path
import numpy as np
from ase.io import read
from util import atoms_num_dict
import subprocess
import argparse
import re

def read_structure(file_path):
    atoms = read(file_path, format="cp2k-restart")
    lattice = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    species = atoms.get_chemical_symbols()
    return np.array(lattice), np.array(positions), species

def write_cp2k_input(title,file_path, lattice, positions, species, band_path):
    num_to_symbol = {v: k for k, v in atoms_num_dict.items()}
    species_symbols = [num_to_symbol[num] for num in species]

    cartesian_positions = np.dot(positions, lattice)

    lattice = adjust_lattice_vectors(lattice)

    with open(file_path, 'w') as f:
        f.write("&GLOBAL\n  PROJECT {}\n  RUN_TYPE BAND\n&END GLOBAL\n\n".format(file_path.split('_seekpath.inp')[0]))
        f.write("&CELL\n")
        for label, vec in zip(["A", "B", "C"], lattice):
            f.write("  {} {:16.10f} {:16.10f} {:16.10f}\n".format(label, *vec))
        f.write("&END CELL\n\n")
        f.write("&COORD\n")
        for specie, pos in zip(species_symbols, cartesian_positions):
            f.write("  {:2} {:16.10f} {:16.10f} {:16.10f}\n".format(specie, *pos))
        f.write("&END COORD\n\n")
        f.write("&PRINT\n")
        f.write("   &BAND_STRUCTURE\n")
        f.write(f"      FILE_NAME {title}.bs\n")
        for start, end in band_path['path']:
            f.write("      &KPOINT_SET\n")
            f.write("      UNITS B_VECTOR\n")
            f.write("         SPECIAL_POINT {:2} {:16.10f} {:16.10f} {:16.10f}\n".format(start, *band_path['kpoints'][start]))
            f.write("         SPECIAL_POINT {:2} {:16.10f} {:16.10f} {:16.10f}\n".format(end, *band_path['kpoints'][end]))
            f.write("         NPOINTS 20\n")
            f.write("      &END KPOINT_SET\n")
        f.write("   &END BAND_STRUCTURE\n")
        f.write("&END PRINT\n")

def adjust_lattice_vectors(lattice):
    for i in range(len(lattice)):
        max_component_index = np.argmax(np.abs(lattice[i]))
        if lattice[i][max_component_index] < 0:
            lattice[i] = -lattice[i]
    return lattice

def modify_inp_file(title, file_path):
    with open(f"{title}_seekpath.inp", "r") as src:
        src_lines = src.readlines()
    band_path_lines = []
    band_path = False
    for line in src_lines:
        if "&PRINT" in line:
            band_path = True
        if band_path:
            band_path_lines.append(line)
        if "&END PRINT" in line:
            band_path = False
            break

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
                    new_line = f"{m.group(1)}Mn_1{m.group(3)}" + line[m.end():]
                    updated_lines.append(new_line)
                else:
                    updated_lines.append(line)
                mn_toggle = not mn_toggle
                continue
            updated_lines.append(line)
            continue

        if "&KIND" in line:
            parts = line.split()
            current_kind = parts[1].strip() if len(parts) > 1 else None
            updated_lines.append(line)
            continue

        elif current_kind and "BASIS_SET" in line:
            if current_kind == "I":
                line = line.replace("TZVP-MOLOPT-GTH-q7", "TZVP-MOLOPT-SR-GTH-q7")
            elif current_kind == "Cu":
                line = line.replace("TZVP-MOLOPT-GTH-q11", "TZVP-MOLOPT-SR-GTH-q11")
            elif current_kind == "Mn":
                line = line.replace("TZVP-MOLOPT-GTH-q15", "TZVP-MOLOPT-SR-GTH-q15")
            updated_lines.append(line)
            continue

        elif current_kind and "POTENTIAL GTH-PBE" in line:
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
                ])
                continue
        if "&END KIND" in line:
            current_kind = None

        if "    BASIS_SET_FILE_NAME  BASIS_MOLOPT" in line:
            updated_lines.append(line)
            updated_lines.append("    BASIS_SET_FILE_NAME  BASIS_MOLOPT_UCL\n")
            continue

        if magnetic_element and "MULTIPLICITY" in line:
            updated_lines.append(line)
            updated_lines.append("    UKS\n")
            continue

        if "    &END POISSON" in line:
            updated_lines.append(line)
            updated_lines.extend(band_path_lines)
            continue

        updated_lines.append(line)

    with open(file_path, "w") as file:
        file.writelines(updated_lines)


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

def run_multiwfn_seekpath(title, k_A, k_B, k_C):
	input_data_bs = f"""cp2k
{title}_bs.inp
2
3
3
2
8
{k_A},{k_B},{k_C}
-9
12
450
0
0
q
"""
	with open(f"{title}_Multiwfn.log", "w") as log_file:
		process = subprocess.run(
			["Multiwfn", f"{title}_seekpath.inp"],
			input=input_data_bs,
			text=True,
			stdout=log_file,
			stderr=subprocess.PIPE
		)

def submit_sbatch(script_path):
	if os.getenv('SINGLE_JOB') == '1':
		subprocess.run(["bash", script_path])
	else:
		subprocess.run(["sbatch", script_path])

def main():
	user_name = parse_args()

	restart_files = glob.glob("*_opt-1.restart")
	if not restart_files:
		raise FileNotFoundError("*_opt-1.restart file not found in the current directory.")
	title = os.path.splitext(restart_files[0])[0].split("_opt-1")[0]

	input_file = f"{title}_opt-1.restart"
	lattice, positions, species = read_structure(input_file)

	atomic_numbers = [atoms_num_dict[atom] for atom in species]

	result = get_path((lattice, positions, atomic_numbers))
	primitive_lattice = result["primitive_lattice"]
	primitive_positions = result["primitive_positions"]
	primitive_species = result["primitive_types"]
	band_path = {
		'path': result["path"],
		'kpoints': result["point_coords"]
	}

	output_file = f"{title}_seekpath.inp"
	write_cp2k_input(title, output_file, primitive_lattice, primitive_positions, primitive_species, band_path)

	cell_lengths = np.linalg.norm(primitive_lattice, axis=1)
	k_lengths = [19, 19, 19]
	k_points = np.ceil(k_lengths / cell_lengths).astype(int)
	k_A, k_B, k_C = k_points

	run_multiwfn_seekpath(title, k_A, k_B, k_C)

	modify_inp_file(title, f"{title}_bs.inp")

	sbatch_script = f"""#!/bin/bash
#SBATCH -J {user_name}-{title}_bs
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -o {title}_bs.log
#SBATCH -e {title}_bs.err

export cp2kroot=/home/think/app/cp2k-2024.1
source $cp2kroot/tools/toolchain/install/setup
export PATH=$cp2kroot/exe/local:$PATH
export OMP_NUM_THREADS=1

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=512m
export OMP_NUM_THREADS=2

echo The job start time at: `date`

mpirun -np 64 -map-by ppr:4:L3cache:pe=2 cp2k.psmp {title}_bs.inp > {title}_bs.out

echo The job end time at: `date`
"""

	script_path = f"{title}_bs.sh"
	with open(script_path, "w") as script_file:
		script_file.write(sbatch_script)

	submit_sbatch(script_path)

if __name__ == "__main__":
	main()