import matplotlib.pyplot as plt
import numpy as np
"""
I take one equilibrated config, average over all different momentum initializations, store this averaged trajectory in an array. I do this for each of the 100 equilibrated configs 
and then average over these.
N_steps+1 values of S for each run of the trajectory. I want to average over . 
The input file lists the data as 1 value of S per line.
    100 x s_traj_data (10001 steps) for equilib state 1
    100 x s_traj_data for equilib state 2     
    100 x s_tr aj_data for equilib state 3
so I need to do if i, line in enumerate(line) %10001==0: mom_samp_arr
s_arr[i%10001] += line_val and at the end divide by 100x100=10000 as averaging over 
the momenta inits first and then averaging over the equilibrated field initializations
is equivalent to averaging over both.

The file is to large to do what I have been up until now with reading it entirely into one array.
In this case I will read 100x10001 lines into an array and average over the 100 momentum initializations to 
"""
num_steps = 5001
sweeps = 100 # poor choice of name but 1 sweep is 1 momentum initialization
equi_configs = 100
#infile = open("var_t_s_traj_outfile_energy_tracking.dat", "r")
#infile = open("ordered_var_t_s_traj_outfile_energy_tracking.dat", "r")
#infile = open("av_s_track_fourier_acc_s_1d_scalar_field_mod_sampling.txt", "r")
#infile = open("av_FFT_accelerated_t_h_unaccel_trajenergy_tracking.dat", "r")

# infile = open("equi_start_ordered_var_t_s_traj_outfile_energy_tracking.dat", "r")
outfile = open("50_sampls_outfile_equi_start_FA_HO_s_traj.dat", "w")
traj_arr = np.zeros(num_steps)
mom_av_traj_arr = np.zeros(num_steps)
#with open("per_line_equi_start_ordered_var_t_s_traj_outfile_energy_tracking.dat") as infile:
count = 0
# sum over 100 randomized momentum trajectories 
#with open("equi_start_m_eff_1.0_dt_0.001_lf_steps_10000_av_track_fourier_acc_s_1d_scalar_field_mod_sampling.txt") as infile:    
with open("equi_start_500equis_m_eff_1.0_dt_0.001_lf_steps_10000_av_track_fourier_acc_s_1d_scalar_field_mod_sampling.txt") as infile:    
    for i,line in enumerate(infile):
        try:
            # average over batches of 100 momenta inits 
            if count >=0 and count <5000*num_steps :
                traj_arr[i%num_steps] += (1/5000.0)*float(line)
                count += 1
            else:
                for s in traj_arr:
                    outfile.write("{}\n".format(s))
                outfile.close()
                exit()
                #count = 0
                #mom_av_traj_arr += (1/100.0)*traj_arr
                #traj_arr = np.zeros(num_steps)
        except ValueError as e:
            print(i)
            #exit()
            print(str(e))
#tokenize line. count up to 



# print(len(s))
# for n in range(equi_configs):
#     for i in range(sweeps):
#         for j in range(num_steps):
#         #for j in range(num_steps-1):    
#             #ind = num_steps*i + j
#             #traj_arr[i][j] = s[(num_steps-1)*i + j]
#             traj_arr[n][i][j] = s[(num_steps)*i*n + j]
            
#             #try:
#             #    s[ind] = 1.0
#             #except IndexError:
#             #    print(ind)
#         #try:
#                 #traj_arr[i][j] = s[num_steps*i + j]
#             #except:
#              #   print(num_steps*i + j)
#               #  exit()
# avg_traj_arr = np.mean(traj_arr,1) # average over momenta inits
# avg_traj_arr = np.mean(avg_traj_arr,0) # average over equi starts 
for s in mom_av_traj_arr:
    outfile.write("{}\n".format(s))
#print(len(avg_traj_arr))
t_arr = np.arange(10001)
#y = 15.5*np.sin(0.001*t_arr)**2
#z = np.sum((np.sin(w_p*t_arr)**2))
#z = w_p*t_arr
#plt.plot(t_arr, traj_arr[9])
plt.plot(traj_arr)
#plt.plot(t_arr, y)
plt.show()
outfile.close()
infile.close()
