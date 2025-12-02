import pandas as pd
import sys
import glob
import os
import subprocess
import re

def merge_and_correct(file_pattern, output_file, adjust_value):
    file_list = sorted(glob.glob(file_pattern))
    if not file_list:
        return False
    dfs = []
    for file in file_list:
        df = pd.read_csv(file, sep='\s+', header=None)
        dfs.append(df)
    merged_df = pd.concat(dfs, ignore_index=True)
    y_columns = merged_df.columns[1:]
    for col in y_columns:
        merged_df[col] = merged_df[col].apply(lambda y: y if y <= 0 else y + adjust_value)
    merged_df.to_csv(output_file, index=False, header=False)
    print(f"Correction completed! Result saved to: {output_file}")
    return True

def magnet_recognition(file_path):

    global magnet
    magnet = False

    with open(file_path, "r") as file:
        lines = file.readlines()

    for line in lines:
        if "    UKS" in line:
                magnet = True


def title_recognition():

    global title
    file_path = glob.glob("*.bs")
    if not file_path:
        raise FileNotFoundError("No *.bs file found in the current directory")
    title = os.path.splitext(file_path[0])[0]
def PBE_bandgap(file_path):
    global Eg_PBE
    Eg_PBE = None
    Eg_PBE_Alpha = None
    Eg_PBE_Beta = None

    magnet_input_data = f"""cp
2
1
10
H
1
1
-3
-1
q
"""

    input_data = f"""cp
2
10
H
1
-3
-1
q
"""

    with open("tmp_Multiwfn_input.txt", "w") as f:
        if magnet:
            f.write(magnet_input_data)
        else:
            f.write(input_data)
    subprocess.run(f"Multiwfn {file_path} < tmp_Multiwfn_input.txt > {title}_Multiwfn_bs.out", shell=True)
    os.remove("tmp_Multiwfn_input.txt")

    float_re = re.compile(r"(\d+\.\d+)")

    with open(f"{title}_Multiwfn_bs.out", "r") as f:
        lines = f.readlines()

    if magnet:
        for line in lines:
            if "Alpha band gap" in line:
                m = float_re.search(line)
                if m:
                    Eg_PBE_Alpha = float(m.group(1))
            elif "Beta band gap" in line:
                m = float_re.search(line)
                if m:
                    Eg_PBE_Beta = float(m.group(1))
        Eg_PBE = min(Eg_PBE_Alpha, Eg_PBE_Beta)
    else:
        for line in lines:
            if " Band gap" in line:
                m = float_re.search(line)
                if m:
                    Eg_PBE = float(m.group(1))
                    break
            
def HSE06_bandgap():
    global Eg_HSE06
    Eg_HSE06 = None

    if magnet:
        with open("DOS_line_Alpha.txt", "r") as f:
            lines = f.readlines()
            for line in lines:
                if not line.strip():
                    continue
                parts = line.split()
                try:
                    value = float(parts[0])
                    if value > 0:
                        Eg_HSE06 = value
                        break
                except (ValueError, IndexError):
                    continue
    else:
        with open("DOS_line.txt", "r") as f:
            lines = f.readlines()
            for line in lines:
                if not line.strip():
                    continue
                parts = line.split()
                try:
                    value = float(parts[0])
                    if value > 0:
                        Eg_HSE06 = value
                        break
                except (ValueError, IndexError):
                    continue
def main():

    adjust_value = Eg_correct

    patterns_outputs = [
        ("path[0-9].txt", "path_merged_corrected.csv"),
        ("path[0-9]_Alpha.txt", "path_merged_corrected_Alpha.csv"),
        ("path[0-9]_Beta.txt", "path_merged_corrected_Beta.csv"),
    ]
    found = False
    for pattern, output in patterns_outputs:
        if merge_and_correct(pattern, output, adjust_value):
            found = True
    if not found:
        print("No path files found")
        sys.exit(1)

if __name__ == "__main__":
    title_recognition()
    magnet_recognition(f"{title}_bs.inp")
    PBE_bandgap(f"{title}.bs")
    HSE06_bandgap()
    Eg_correct = float(Eg_HSE06) - float(Eg_PBE)

    main()