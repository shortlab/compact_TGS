import pandas as pd
pname = 'RSIData/MainSetup/Tungsten_Calibration-2024-04-04-07.6 um-spot_postprocessing.txt'
with open(pname) as f:
    lines = f.readlines()
arr = []
for line in lines:
    if "Users" in line:
        arr.append(line.split(' ')[3:10])

pd.DataFrame(arr).to_csv('MainSetupData.csv')
