import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""
Post-processing for 2d scalar field

The data in my files is in the format 

[[1,2,3][4,5,6]][[7,8,9][10,11,12]]


A 2d array like this shall have elements indexed like [0][1] = 2
"""

N_traj = 1000000
N_x = 10
N_tau = 10

infile = open("2d_scalar_field_S_avg_", "r")

phi = np.zeros((N_traj, N_tau, N_x))
s = infile.read()
s = re.sub("\n", " ", s)
phi_list = s.split()


for k in range(N_traj):
    for i in range(N_tau):
        for j in range(N_x):
            phi[k][i][j] = phi_list[k + (i*N_tau) + (j*N_x)]
    

phi = phi[100:]            # Discard first 100 trajectories (thermalization)

"""
Plot of the average of phi**2 (averaged over both the img time and x axes) for each trajectory
"""
reshaped_phi = phi.reshape((phi.shape[0],-1)) # reshape phi array into an array of N_traj 1 dimensional arrays. 
phi_sqr_avg = np.mean(reshaped_phi**2, axis=1) # now average over the flattened 1d arrays for each trajectory 
print("<phi^2> avg = {}".format(np.mean(phi_sqr_avg))) # print average value of phi^2 for simulation
traj_lst = range(100,N_traj)

plt.scatter(traj_lst, phi_sqr_avg)
plt.show()



"""
The commented code plots  the average (over img time axis) of the square of the field as a function of x position for each trajectory 
"""

################################################################################################################################################################
# phi_sqr_avg = np.flatten(np.mean(phi**2,1)) # Obtain the average (over img time axis) of the square of the field as a function of x position for each trajectory #
#                                                                                                                                                              #
# xs = [i for i in range(N_x) for j in range(len(phi_sqr_avg))]                                                                                                  #
# ys = np.tile(np.arange(N_x), len(phi_sqr_avg))                                                                                                                 #
# z = phi_sqr_avg                                                                                                                                                #
#                                                                                                                                                              #
# fig = plt.figure()                                                                                                                                           #
# ax = fig.add_subplot(111,projection="3d")                                                                                                                    #
#                                                                                                                                                              #
# ax.scatter(xs, ys, z)                                                                                                                                        #
# plt.show()                                                                                                                                                   #
################################################################################################################################################################
