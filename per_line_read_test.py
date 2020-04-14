import numpy as np
delta_H_arr = np.zeros(10000)
with open("per_line_equi_start_delta_h_energy_tracking.dat") as infile:
    for i,line in enumerate(infile):
        try:
            delta_H_arr[i] = float(line)
        except ValueError as e:
            print(str(e))
