import numpy as np
import random as rand
import math
import matplotlib
import matplotlib.pyplot as plt
import sys 
matplotlib.rcParams['text.usetex'] = True
np.set_printoptions(threshold=sys.maxsize)
"""
Steps:

propose a number u which is distributed according to a uniform distro between -h and h. The value of h is adjusted to meet a predefined acceptance ratio.

"""


# initialize (thermalize) the system (cold first). The system is an array of N_tau position values corresponding to N_tau values of time.

N_tau = 1000  # number of timeslices
N_sep = 10 # Every 1000th sweep is used to calculate observable quantities
m = 0.1 # mass of particle
omega = 1.0
h = 1.0 # this parameter is adjusted to meet the predefined acceptance ratio idrate
num_sweeps = 100000          # 100000 sweeps
idrate = 0.4 # this is the ideal acceptance rate. It is used to adjust h
accrate = 0.0

def init(mode, N_tau, h):
    # either hot, cold or approx equilibrated
    if mode=="cold":
        return np.zeros(N_tau)
    elif mode=="hot":
        np.random.uniform(-0.5, 0.5, N_tau) # init lattice sites to random numbers between -1/2 and 1/2

        
def MCMC_init(N_tau,x,h,idrate,m,omega):
    """
    Performs one sweep through the lattice
    """
    accrate = 0.0 # this is the acceptance rate.
    for i,x_i in enumerate(x):
        
        x_new = x_i + np.random.uniform(-h,h)     # we now propose an update to position x_i. x_i -> x_i + u
        if accept(x_new,x,i):
            x[i] = x_new
            # we stop modifying the acceptance rate after equilibration 
            accrate += 1.0/N_tau
    h = h*(accrate/idrate)
    return h,x
    
def MCMC(N_tau,x,h,idrate,m,omega):
    """
    Performs one sweep through the lattice
    """
    #accrate = 0.0 # this is the acceptance rate.
    global accrate
    for i,x_i in enumerate(x):
        
        x_new = x_i + np.random.uniform(-h,h)     # we now propose an update to position x_i. x_i -> x_i + u
        if accept(x_new,x,i):
            x[i] = x_new
            # we stop modifying the acceptance rate after equilibration 
            accrate += 1.0
    return x   

def accept(x_new, x, i):
    tau_plus = (i+1)%N_tau
    tau_minus = (i-1)%N_tau
    x_plus = x[tau_plus]
    x_minus = x[tau_minus]
    S_old = 0.5*m*(x_plus - x[i])**2 + 0.5*m*(x[i]-x_minus)**2 + 0.5*m*(omega**2)*(x[i]**2)     # we need only consider positions x_(i-1), x_i, x_i+1
    S_new = 0.5*m*(x_plus - x_new)**2 + 0.5*m*(x_new-x_minus)**2 + 0.5*m*(omega**2)*(x_new**2)   # as all others cancel in the difference delta_S
    delta_S = S_new - S_old

    p_accept =  math.exp(-delta_S)
    # accept update with probability p_accept. If the update lowers the action, then exp(-delta_S)>1, and the update is definitely made.
    # if the update doesn't lower the action, acceptance occurs with probability exp(-delta_S) as ensured by comparison with a unif(0,1) random
    if(p_accept > np.random.uniform(	)): 
        return True
    else:
        return False


x = init("cold", N_tau, h)
    
# test run first to see how long thermalization takes.
# only 1 sweep in 100 kept as a "configuration" and used for calculating observable averages.

# store output (the x array) in a file for plotting and further analysis.
outfile = open("mcmc_ho_outfile.txt".format(N_tau,m,omega),"w")
av_x_sq_outfile = open("mcmc_defin_x2_av.txt", "w")

config_num = 0 #     counter for configuration used for observable calculations
x_squared_avg = np.zeros(math.floor(num_sweeps/N_sep))
x_avg = np.zeros(math.floor(num_sweeps/N_sep))

# acceptance rate modification 
accRate_init = 100 # number of steps to modify acceptance rate

#for sweep in range(accRate_init): 
#    h,x = MCMC_init(N_tau,x,h,idrate,m,omega)

for sweep in range(num_sweeps):
    print(sweep)
    #print(accrate)
    x = MCMC(N_tau,x,h,idrate,m,omega)
    if(sweep % N_sep == 0):
        print("accrate = {}\n".format(accrate/(N_tau*(sweep+1))))
        #<x**2> is found by averaging over all configs the average of x_i**2 over the N_tau timeslices in each config
        #x_squared_avg[config_num] = np.mean(x**2)   # store avg of each config in nparray for plotting
        #x_avg[config_num] = np.mean(x)
        av_x_sq_outfile.write("{}\n".format((1.0/N_tau)*np.sum(x**2)))                            # also output array x at each config to file for further analysis 
        config_num +=1

# plot <x^2> vs config number to determine thermalization time



#fig, (ax1,ax2) = plt.subplots(1,2) 
## ax1.xlabel("Config number")
#ax1.set_title(r'$\langle x^2 \rangle$')
#ax1.plot(x_squared_avg)
#ax2.set_title(r'$\langle x \rangle$')
##ax2.xlabel("Config number")
##ax2.ylabel(r'$\langle x \rangle')
#ax2.plot(x_avg)
#plt.show()



#<x**2> is found by averaging over all configs the average of x_i**2 over the N_tau timeslices in each config
# x_squared_ensmb_avg = x_squared_ensmb_avg/(math.floor(num_sweeps/config_freq))
        

outfile.close()
