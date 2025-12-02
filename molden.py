import subprocess
import os
import glob
import shutil
import argparse
def title_recognition():

    global title
    file_path = glob.glob("*_bs.inp")
    if not file_path:
        raise FileNotFoundError("*_bs.inp file not found in the current directory")
    title = os.path.splitext(file_path[0])[0].split("_bs")[0]

def primitive_cell_recognition(file_path):

    global metal, halogen, cell_A, cell_B, cell_C, cell_AA, cell_BB, cell_CC, atom_numbers
    current_kind = None
    in_cell = False

    if file_path == f"{title}_bs.inp":
        metal = None
        halogen = None
        cell_A = cell_B = cell_C = None
        with open(file_path, "r") as file:
            lines = file.readlines()

        for line in lines:
            if "&CELL" in line:
                in_cell = True
                continue
            if in_cell:
                if line.strip().startswith("A"):
                    parts = line.split()
                    if len(parts) == 4:
                        cell_A = f"{parts[1]},{parts[2]},{parts[3]}"
                elif line.strip().startswith("B"):
                    parts = line.split()
                    if len(parts) == 4:
                        cell_B = f"{parts[1]},{parts[2]},{parts[3]}"
                elif line.strip().startswith("C"):
                    parts = line.split()
                    if len(parts) == 4:
                        cell_C = f"{parts[1]},{parts[2]},{parts[3]}"
                elif line.strip().startswith("&END CELL"):
                    in_cell = False
            elif "&KIND" in line:
                parts = line.split()
                if len(parts) > 1:
                    current_kind = parts[1].strip()
            elif current_kind:
                metals = {"Cu", "Sb", "In", "Mn"}
                halogens = {"Cl", "Br", "I"}
                            
                if current_kind in metals:
                    metal = current_kind
                elif current_kind in halogens:
                    halogen = current_kind

            if "&END KIND" in line:
                current_kind = None
                
    elif file_path == f"{title}_bs.out":
        atom_numbers = None

        with open(file_path, "r") as file:
            for line in file:
                if "- Atoms:" in line:
                    atom_numbers = line.strip().split(":")[1].strip()
                    break
    
    elif file_path == f"{title}_HSE06.inp":
        cell_AA = cell_BB = cell_CC = None
        with open(file_path, "r") as file:
            lines = file.readlines()

        for line in lines:
            if "&CELL" in line:
                in_cell = True
                continue
            if in_cell:
                if line.strip().startswith("A"):
                    parts = line.split()
                    if len(parts) == 4:
                        cell_AA = f"{parts[1]},{parts[2]},{parts[3]}"
                elif line.strip().startswith("B"):
                    parts = line.split()
                    if len(parts) == 4:
                        cell_BB = f"{parts[1]},{parts[2]},{parts[3]}"
                elif line.strip().startswith("C"):
                    parts = line.split()
                    if len(parts) == 4:
                        cell_CC = f"{parts[1]},{parts[2]},{parts[3]}"
                elif line.strip().startswith("&END CELL"):
                    in_cell = False

def magnet_recognition(file_path):

    global magnet
    magnet = False

    with open(file_path, "r") as file:
        lines = file.readlines()

    for line in lines:
        if "    UKS" in line:
                magnet = True

def modify_molden_file(file_path):

    add_data = f"""[Cell]
{cell_AA.replace(',', ' ')}
{cell_BB.replace(',', ' ')}
{cell_CC.replace(',', ' ')}
[Nval]
O 6
N 5
C 4
H 1
Cu 11
Sb 5
In 13
Mn 15
P 5
I 7
Br 7
Cl 7
"""
    temp_path = f"{file_path}.tmp"
    
    with open(temp_path, 'w') as new_file, open(file_path, 'r') as orig_file:
        new_file.write(add_data)
        shutil.copyfileobj(orig_file, new_file)

    shutil.move(temp_path, file_path)


def molden_run_and_bs_merge_and_correct(file_path, user_name):

    magnet_input_data_molden_1 = f"""10
-1
1
cond
{metal}
a
S
q
2
cond
{metal}
a
P
q
3
cond
{metal}
a
D
q
4
cond
{halogen}
a
S
q
5
cond
{halogen}
a
P
q
6
e {metal}
e {halogen}
q
7
e C
e H
e N
e O
e P
q
0
3 //FWHM
0.25
-6
H
-0
16 //name
1
{metal}_S
2
{metal}_P
3
{metal}_D
4
{halogen}_S
5
{halogen}_p
6
{metal}{halogen}
7
Organic ligand
0
2 //pdf
3 //data
0
"""
    magnet_input_data_molden_2 = f"""10
-1
1
cond
{metal}
a
S
q
2
cond
{metal}
a
P
q
3
cond
{metal}
a
D
q
4
cond
{halogen}
a
S
q
5
cond
{halogen}
a
P
q
6
e {metal}
e {halogen}
q
7
e C
e H
e N
e O
e P
q
0
3 //FWHM
0.25
-6
H
6
2
-0
16 //name
1
{metal}_S
2
{metal}_P
3
{metal}_D
4
{halogen}_S
5
{halogen}_p
6
{metal}{halogen}
7
Organic ligand
0
2 //pdf
3 //data
0
-10
6
-3 //Keep selected atomic wavefunctions
1-{atom_numbers}
-1
300
7
17 //Keep selected atoms
1-{atom_numbers}
26 //Set cell information
7 //Vector 1
{cell_A}
8
{cell_B}
9
{cell_C}
0
-10
0
200
3
ha
9
0,0,0
0,0,0
0.25
3
la
9
0,0,0
0,0,0
0.25
3
hb
9
0,0,0
0,0,0
0.25
3
lb
9
0,0,0
0,0,0
0.25
"""
    input_data_molden = f"""10
-1
1
cond    
{metal}
a
S
q
2
cond
{metal}
a
P
q
3
cond
{metal}
a
D
q
4
cond
{halogen}
a
S
q
5
cond
{halogen}
a
P
q
6
e {metal}
e {halogen}
q
7
e C
e H
e N
e O
e P
q
0
3 //FWHM
0.25
-6
H
-0
16 //name
1
{metal}_S
2
{metal}_P
3
{metal}_D
4
{halogen}_S
5
{halogen}_p
6
{metal}{halogen}
7
Organic ligand
0
2 //pdf
3 //data
0
-10
6
-3 //Keep selected atomic wavefunctions
1-{atom_numbers}
-1
300
7
17 //Keep selected atoms
1-{atom_numbers}
26 //Set cell information
7 //Vector 1
{cell_A}
8
{cell_B}
9
{cell_C}
0
-10
0
200
3
h
9
0,0,0
0,0,0
0.25
3
l
9
0,0,0
0,0,0
0.25
3
h-1
9
0,0,0
0,0,0
0.25
3
l+1
9
0,0,0
0,0,0
0.25
"""
    
    if magnet:
        sbatch_script_molden = f"""#!/bin/bash
#SBATCH -J {user_name}-{title}_molden
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -o {title}_molden.log
#SBATCH -e {title}_molden.err

Multiwfn {file_path} << EOF > {title}_Multiwfn_molden.out
{magnet_input_data_molden_1}
EOF
mv dislin.pdf {title}_dos_Alpha.pdf
mv DOS_line.txt DOS_line_Alpha.txt
mv DOS_curve.txt DOS_curve_Alpha.txt

Multiwfn {file_path} << EOF > {title}_Multiwfn_molden.out
{magnet_input_data_molden_2}
EOF
mv dislin.pdf {title}_dos_Beta.pdf
mv DOS_line.txt DOS_line_Beta.txt
mv DOS_curve.txt DOS_curve_Beta.txt

python /home/think/Desktop/ningding/script/bs_correct.py
mv dislin.pdf {title}_bs.pdf
"""
    else:
        sbatch_script_molden = f"""#!/bin/bash
#SBATCH -J {user_name}-{title}_molden
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -o {title}_molden.log
#SBATCH -e {title}_molden.err

Multiwfn {file_path} << EOF > {title}_Multiwfn_molden.out
{input_data_molden}
EOF
mv dislin.pdf {title}_dos.pdf

python /home/think/Desktop/ningding/script/bs_correct.py
mv dislin.pdf {title}_bs.pdf
"""

    script_path = f"{title}_molden.sh"
    with open(script_path, "w") as script_file:
        script_file.write(sbatch_script_molden)

    if os.getenv('SINGLE_JOB') == '1':
        subprocess.run(["bash", script_path])
    else:
        subprocess.run(["sbatch", script_path])



def main():
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
    title_recognition()
    primitive_cell_recognition(f"{title}_bs.inp")
    primitive_cell_recognition(f"{title}_bs.out")
    primitive_cell_recognition(f"{title}_HSE06.inp")
    modify_molden_file(f"{title}_HSE06-MOS-1_0.molden")
    magnet_recognition(f"{title}_bs.inp")
    molden_run_and_bs_merge_and_correct(f"{title}_HSE06-MOS-1_0.molden", user_name)

if __name__ == "__main__":

    main()