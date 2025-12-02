import argparse
import os
import subprocess

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

def write_and_submit(user_name):
    sbatch_script = f"""#!/bin/bash
#SBATCH -J {user_name}-all
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -o {user_name}_all.log
#SBATCH -e {user_name}_all.err

# run all steps inside one Slurm job
export SINGLE_JOB=1

echo The combined job start time at: `date`

# Ensure current working directory contains the scripts and input files.
python3 $dnpy/geo_opt.py -n {user_name}
python3 $dnpy/bs.py -n {user_name}
python3 $dnpy/wfn.py -n {user_name}
python3 $dnpy/molden.py -n {user_name}

echo The combined job end time at: `date`
"""
    script_path = "submit_all.sh"
    with open(script_path, "w") as f:
        f.write(sbatch_script)
    subprocess.run(["sbatch", script_path])

def main():
    user = get_user_name()
    write_and_submit(user)

if __name__ == "__main__":
    main()
