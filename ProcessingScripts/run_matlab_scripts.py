import matlab.engine
import numpy as np
import pandas as pd
import os 

eng = matlab.engine.start_matlab()
eng.cd('/Users/jonasrajagopal/Desktop/ShortUROP/TGS_Data_Analysis/TGS-Processing-Scripts/MATLAB')
COLUMNS = ["run_name", "frequency","frequency error", "SAW speed", "thermal diffusivity", "thermal diffusivity error", "acoustic damping", "acoustic damping error", "A", "A error", "beta", "beta error", "B", "B error", "theta", "theta error", "C", "C error"]

def process_matlab_output(raw_output, sample_name):
    raw_output = list(raw_output)
    raw_output.insert(0, sample_name)
    raw_output[6] = np.asarray(raw_output[6])[0][2]
    return pd.DataFrame(np.array([raw_output]), columns=COLUMNS)

def data_process_mini_tgs(pos_file, pos_baseline, verbose=0):
    grating = 7.6
    start_phase = 2
    baseline_subtract = 1
    two_SAW_frequencies = 0
    mono_heterodyne = 1
    delimiter = ','
    plot_outputs = 0
    header_length = 5
    sample_name = pos_file.split('/')[-1]
    raw_output = eng.TGSPhaseAnalysis(pos_file,0,grating,start_phase,two_SAW_frequencies,baseline_subtract, pos_baseline, 0, delimiter, verbose, header_length, plot_outputs, mono_heterodyne, nargout=17)
    return process_matlab_output(raw_output, sample_name)

def data_process_main_tgs(pos_file, neg_file, grating, verbose=0):
    start_phase = 2
    baseline_subtract = 0
    two_SAW_frequencies = 0
    mono_heterodyne = 0
    delimiter = ''
    plot_outputs = 0
    header_length = 16
    sample_name = pos_file.split('/')[-1]
    raw_output = eng.TGSPhaseAnalysis(pos_file,neg_file,grating,start_phase,two_SAW_frequencies,baseline_subtract, 0, 0, delimiter, verbose, header_length, plot_outputs, mono_heterodyne, nargout=17)
    return process_matlab_output(raw_output, sample_name)

data_directory = os.path.join(os.getcwd(), 'Data')
run_name = 'Mar21'
output_name = 'C2_wdata.csv'
data = pd.DataFrame([], columns=COLUMNS)

for i in range(25):
    posFile = os.path.join(data_directory, run_name, f'C2w000{f"{i:02}"}.txt')
    posBaseline = os.path.join(data_directory, run_name, f'C2wbase000{f"{i:02}"}.txt')
    df = data_process_mini_tgs(posFile, posBaseline)
    data = pd.concat([data, df])

# for i in range(25):
#     posFile = os.path.join(data_directory, run_name, f'Tungsten_Calibration-2024-04-04-07.6 um-spot{i}-POS-1.txt')
#     negFile = os.path.join(data_directory, run_name, f'Tungsten_Calibration-2024-04-04-07.6 um-spot{i}-NEG-1.txt')
#     df = data_process_main_tgs(posFile, negFile, grating=7.6)
#     data = pd.concat([data, df])

data.to_csv(os.path.join(data_directory, run_name, output_name))

