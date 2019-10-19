import numpy as np
import random as rand
import math
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True

"""
Steps:

propose a number u which is distributed according to a uniform distro between -h and h. The value of h is adjusted to meet a predefined acceptance ratio.

"""


# initialize (thermalize) the system (cold first). The system is an array of N_tau position values corresponding to N_tau values of time.

N_tau = 120 # number of timeslices
N_sep = 100 # Every 1000th sweep is used to calculate observable quantities. This doesn't have to be so big. 100 should be fine
m = 1.0 # mass of particle
omega = 1.0
h = 1.0 # this parameter is adjusted to meet the predefined acceptance ratio idrate
num_sweeps = 1000000 # 100000 sweeps
sweepNum = 1 # for printing progress
idrate = 0.5 # this is the ideal acceptance rate. It is used to adjust h

def init(mode, N_tau, h):
    # either hot, cold or approx equilibrated
    if mode=="cold":
        return np.zeros(N_tau)
    elif mode=="hot":
        return np.random.uniform(-0.5, 0.5, N_tau) # init lattice sites to random numbers between -1/2 and 1/2

        
def MCMC(N_tau,x,h,idrate,m,omega):
    """
    Performs one sweep through the lattice.
    Returns h (the updated random step bounds, x (the updated path)
    """
    global sweepNum
    accrate = 0.0 # this is the acceptance rate.
    
    index = [0]*N_tau
    # site visiting order
    for i in range(N_tau):
        #for a random visitation of lattice sites
        index[i] = math.floor(N_tau * np.random.uniform())
        #index[i]=i # visit each lattice site in order
    #for i,x_i in enumerate(x):
    for i in range(N_tau):    
        tau = index[i]
        x_new = x[tau] + h*np.random.uniform(-1,1)     # we now propose an update to position x_i. x_i -> x_i + u
        if accept(x_new,x,tau):
            x[tau] = x_new
            accrate += 1.0/N_tau
    print("sweep: {}, h = {}".format(sweepNum,h))
    sweepNum +=1
    h = h*(accrate/idrate)
    return h,x

def accept(x_new, x, i):
    tau_plus = (i+1)%N_tau
    tau_minus = (i-1+N_tau)%N_tau
    x_plus = x[tau_plus]
    x_minus = x[tau_minus]  
    S_old = 0.5*m*(x_plus - x[i])**2 + 0.5*m*(x[i]-x_minus)**2 + 0.5*m*(omega**2)*(x[i]**2)     # we need only consider positions x_(i-1), x_i, x_i+1
    S_new = 0.5*m*(x_plus - x_new)**2 + 0.5*m*(x_new-x_minus)**2 + 0.5*m*(omega**2)*(x_new**2)   # as all others cancel in the difference delta_S
    delta_S = S_new - S_old

    p_accept =  math.exp(-delta_S)
    # accept update with probability p_accept. If the update lowers the action, then exp(-delta_S)>1, and the update is definitely made.
    # if the update doesn't lower the action, acceptance occurs with probability exp(-delta_S) as ensured by comparison with a unif(0,1) random
    if(p_accept > np.random.uniform(    )): 
        return True
    else:
        return False


x = init("cold", N_tau, h)
    
# test run first to see how long thermalization takes.
# only 1 sweep in N_sep kept as a "configuration" and used for calculating observable averages.

# store output (the x array) in a file for plotting and further analysis.
outfile = open("output_hot.txt","w")

config_num = 0 # counter for configuration used for observable calculations
x_squared_avg = np.zeros(math.floor(num_sweeps/N_sep))
x_avg = np.zeros(math.floor(num_sweeps/N_sep))
psi_0 = np.array([]) # ground state wavefunction
for sweep in range(num_sweeps):
    h,x = MCMC(N_tau,x,h,idrate,m,omega)
    if(sweep % N_sep == 0):
        #<x**2> is found by averaging over all configs the average of x_i**2 over the N_tau timeslices in each config
        x_squared_avg[config_num] = (1/N_tau)*np.mean(x**2)   # store square avg PER SITE of each config in nparray for plotting
        x_avg[config_num] = (1/N_tau)*np.mean(x)
        outfile.write("{0}\n".format(str(x)))                            # also output array x at each config to file for further analysis 
        config_num +=1
        
        # plot ground state wavefunction
        """
        At every configuration we bin x values to form a hisogram representing the ground state wavefunction.
        """
        psi_0 = np.append(psi_0,x)

# plot <x^2> and <x> vs config number to determine thermalization time
# plot also the square modulus of the groundstate wavefunction histogram |psi_0|^2
psi_0 = np.histogram(psi_0,bins='auto')


fig, axs = plt.subplots(2,2) 
# ax1.xlabel("Config number")
axs[0,0].set_title(r'$\langle x^2 \rangle$')
axs[0,0].plot(x_squared_avg)
axs[0,1].set_title(r'$\langle x \rangle$')
#ax2.xlabel("Config number")
#ax2.ylabel(r'$\langle x \rangle')
axs[0,1].plot(x_avg)
axs[1,0].plot(list(psi_0[0]))
plt.show()








# <x**2> is found by averaging over all configs the average of x_i**2 over the N_tau timeslices in each config
x_squared_ensmb_avg = np.mean(x_squared_avg)
x_ensmb_avg = np.mean(x_avg)

obs_avg_outfile = open("observable_avgs.txt", "w")
obs_avg_outfile.write("x_squared_ensmb_avg = {} \n x_ensmb_avg = {}\n".format(str(x_squared_ensmb_avg), str(x_ensmb_avg)))


outfile.close()
