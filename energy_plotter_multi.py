import matplotlib.pyplot as plt
import numpy as np
import os

cur_path = os.getcwd()
print(cur_path)
#new_path = os.path.join(cur_path,"ordered_starts//m_eff_0.1_dt_0.001_lf_steps_10000_av_track_fourier_acc_s_1d_scalar_field_mod_sampling.txt")

"N_steps+1 values of S for each run of the trajectory. I want to average all 100 runs. "
num_steps = 10001
sweeps = 100
#infile = open("var_t_s_traj_outfile_energy_tracking.dat", "r")
#infile = open("ordered_var_t_s_traj_outfile_energy_tracking.dat", "r")
#infile = open("av_s_track_fourier_acc_s_1d_scalar_field_mod_sampling.txt", "r")
#infile = open("av_FFT_accelerated_t_h_unaccel_trajenergy_tracking.dat", "r")
#infile = open(new_path, "r")
file_name_list = [
"1.0_lambda_ordered_start_m_eff_1.2_dt_0.001_lf_steps_10000_av_s_track_FA_AHO.txt",
"1.0_lambda_ordered_start_m_eff_2.0_dt_0.001_lf_steps_10000_av_s_track_FA_AHO.txt",
"1.0_lambda_ordered_start_m_eff_3.0_dt_0.001_lf_steps_10000_av_s_track_FA_AHO.txt",
"1.0_lambda_ordered_start_m_eff_1.0_dt_0.001_lf_steps_10000_av_s_track_FA_AHO.txt", 
"1.0_lambda_ordered_start_m_eff_0.9_dt_0.001_lf_steps_10000_av_s_track_FA_AHO.txt"
]

handl = []
for count,fname in enumerate(file_name_list):
    infile = open(fname, "r")
    s = [float(S) for S in infile.read().split()]
    print(len(s))
    traj_arr = np.zeros((sweeps, num_steps))
    print(len(s))
    for i in range(sweeps):
        for j in range(num_steps):
        #for j in range(num_steps-1):    
            #ind = num_steps*i + j
            #traj_arr[i][j] = s[(num_steps-1)*i + j]
            traj_arr[i][j] = s[(num_steps)*i + j]
            
            #try:
            #    s[ind] = 1.0
            #except IndexError:
            #    print(ind)
        #try:
                #traj_arr[i][j] = s[num_steps*i + j]
            #except:
            #   print(num_steps*i + j)
            #  exit()
    avg_traj_arr = np.mean(traj_arr,0) 
    
    print(len(avg_traj_arr))
    
    plt_param, = plt.plot((1/10.0)*avg_traj_arr + count*1.5*np.ones(len(avg_traj_arr)), label=fname.split("m_eff_")[1].split("_")[0])
    handl.append(plt_param)

plt.title(r"Intra-trajectory evolution of the action for $\phi^4$ theory")
plt.ylabel(r"$\langle S \rangle$")
plt.xlabel(r"Hamiltonian dynamics time")
plt.legend(handles=handl, title=r"$M$", loc=4)    
plt.show()