import matplotlib.pyplot as plt
import numpy as np

#infile = open("unaccelerated_s_outfile.txt", "r")
#infile = open("var_t_s_traj_outfile.dat", "r")

#infile = open("FFT_accelerated_t_h_unaccel_trajenergy_tracking.dat", "r")
infile = open("h_track_fourier_acc_s_1d_scalar_field_mod_sampling.txt", "r")
#infile=open("ordered_var_t_h_unaccel_trajenergy_tracking.dat","r")
#infile = open("ordered_var_t_s_traj_outfile_energy_tracking.txt", "r")
#infile = open("var_t_h_unaccel_traj.dat", "r")
s = infile.read().split()
s = [float(s_i) for s_i in s]
#s = s[0:200]
t_arr = np.arange(len(s))
#s_arr = []
#
#for i in range(10): 
#    s_arr.append(float(s[i]))
#plt.ylim(0, 100)    
#print(s)
infile.close()
plt.plot(t_arr[:-1], s[:-1])
#plt.ticklabel_format(useOffset=False)
plt.show()