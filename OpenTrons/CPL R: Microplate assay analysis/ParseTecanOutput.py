!pip install --ignore-installed jedi pandas numpy==1.23 numba openpyxl 'setuptools<=64.0.2' cvxpy

import pandas as pd
import re

##########################################
USER DEFINED VARIABLES
filename = "13.03.23-tecan.xlsx" # Enter the path to the file produced after running the microplate reader
##########################################

# Define functions
def parse_chunk(chunk):
    d = {'well': chunk[0],
         'time_hours': chunk[1],
         'mean_od': chunk[2]}
    return pd.DataFrame(d)

def parse_tecan(file):
    df = pd.read_excel(file)
    df.to_csv('tmp.csv', index=False)
    chunks = []
    chunk = []
    with open('tmp.csv', 'rt') as handle:
        for line in handle.readlines():
            line = line.strip()
            if line.startswith('Cycles / Well'):
                if chunk:
                    chunks.append(chunk)
                chunk = []
                continue
            rgx = re.compile('(^[A-H][0-9]+)')
            m = rgx.search(line)
            if m:
                chunk.append(m.group(1))
            elif line.startswith('Time [s]'):
                chunk.append([round(float(x)/60/60, 2) for x in line.split(',')[1:]])
            elif line.startswith('Mean'):
                chunk.append([round(float(x), 4) for x in line.split(',')[1:]])
        chunks.append(chunk)
    pd.concat([parse_chunk(x) for x in chunks]).to_csv('reformatted.csv', index=False)

# Call function
df = parse_tecan(filename)
