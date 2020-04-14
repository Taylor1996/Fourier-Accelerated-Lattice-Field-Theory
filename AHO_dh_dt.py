import numpy as np
import matplotlib.pyplot as plt
import sys
np.set_printoptions(threshold=sys.maxsize)


"""
This program is specifically for use in obtaining the behaviour of dH vs dt
"""

# Generate an initial state. setting = 0 for now

a = 0.25
N_tau = 200 # size of the lattice (i.e., the length of the harmonic oscillator's path in time)
N_trajs = 10000 # number of trajectories (equivalent to number of measurements of observables)

t = 0
mu_sqr = -4.0
lf_steps = int(sys.argv[2])  # leapfrog steps in a trajectory
delta_t = 1.0/ float(lf_steps) # leapfrog timestep in computer time. Passed as cmd line arg for use in bash
m = 0.5 # mass of particle. Do an array divided by array element-wise if you want multiple particles with different masses.
omega = -1.1 # natural frequency of oscillator
Lambda = 1.0
f_sqr = 2.0
accrate = 0.0
idrate = 0.8 # ideal acceptance rate

def init(mode, N_tau):
    # either hot, cold or approx equilibrated
    if mode=="cold":
        return np.zeros(N_tau)
    elif mode=="hot":
        return np.random.uniform(-0.5, 0.5, N_tau) # init lattice sites to random numbers between -1/2 and 1/2



def S(q):
    # S is the action which acts like the potential in HMC: H = K + S
    # define T to be x_i+1 - x_i bit
    q_shifted  = np.array([q[(i+1)%N_tau] for i in range(len(q))])
    sqr_term = np.sum((q_shifted - q)**2)
    q_f = q**2-f_sqr
    S = (1/(2.0*a))*m*sqr_term +  Lambda*np.sum(q_f**2)  # diff form of potential for comparison with lit
    #S =  a*((1/2.0*a**2)*m*sqr_term + (1/2.0)*mu_sqr*np.sum(q**2) + Lambda*np.sum(q**4))
    return S


def grad_S(q):
    q_forward  = np.array([q[(i+1)%N_tau] for i in range(len(q))])
    q_backward  = np.array([q[(i-1)%N_tau] for i in range(len(q))])
    #grad_S = a*((m/a**2)*(-q_forward + 2*q - q_backward) + mu_sqr*q + 4*Lambda*q**3)
    grad_S = (m/a)*(-q_forward + 2*q - q_backward) + 4*Lambda*q*(q**2 - f_sqr)
    # gradient of whole action!
    return grad_S


def HMC(q_current, delta_t, lf_steps,sweep):
    global accrate 
    p = np.random.normal(0,1,N_tau)
    

    q = q_current
    p_current = p
    
    # leapfrog
    p = p - (grad_S(q) * (delta_t/2.0))
    
    for i in range(lf_steps):
        q = q + delta_t*(p) # m should not be present here
        if i != lf_steps-1:
            p = p - (grad_S(q)* delta_t)

    # final leap frog step for momentum
    p = p - (grad_S(q) * (delta_t/2.0))

    if accept(q_current,p_current,q,p):
        q_current = q
        p_current = p
        accrate +=1

    
    # should return modified delta_t so as to get accrate nearer idrate

    return q_current,p_current,accrate
    

def accept(q_current,p_current,q,p):

    # Calculation of new Hamiltonian

    U_current = S(q_current)
    K_current = (np.sum(p_current**2))/(2.0) # note that m doesn't appear here b/c the Hamiltonian dynamics is a completely artificial dynamics
    U_proposed = S(q)
    K_proposed = (np.sum(p**2))/(2.0) # actually, I'm a bit confused about the above - see Neal???

    H_current = U_current + K_current
    H_proposed = U_proposed + K_proposed
    
    delta_H = H_proposed - H_current 
    delta_H_outfile.write(str(delta_H) +" ") # write out values of delta_H to check for correct behaviour of HMC 
    # if the new state decreases the exponential of the Hamiltonian, accept the new state. Otherwise, reject with probability 1-e^(-delta_H)
    if (np.random.uniform() < np.exp(-delta_H)):
        # accept
        print("accept")
        print("accrate = {}, delta_t = {}".format(accrate/(sweep+1), str(delta_t)))
        return True 
    else:
        print("reject")
    
    print("accrate = {}, delta_t = {}".format(accrate/(sweep+1), str(delta_t)))
    return False
    

q_current = init("cold", N_tau)
#outfile = open("aho_test_m_0_N_t_{}_a_{}_N_trajs_{}".format(N_tau, a, N_trajs) + sys.argv[1],"w")
#outfile = open("aho_creutz_dh_dt_N_t_{}_a_{}_N_trajs_{}".format(N_tau, a, N_trajs) + sys.argv[1],"w")
#delta_H_outfile = open("anharmonic_creutz_dh_dt_hmc_delta_H_" + str(sys.argv[1]),"w")
delta_H_outfile = open("anharmonic_hmc_delta_H_" + str(sys.argv[1]),"w")


for sweep in range(N_trajs):
    # the current plotting code I have deals with paths, not averaged paths
    #U = S(q_current)
    #grad_S = grad_S(q_current)
    q_current,p_current,accrate = HMC(q_current, delta_t, lf_steps,sweep)
    #if ((sweep % accrate_update_freq ==0) and sweep != 0):
    #    print("accrate={}, acc_up_freq={}".format( accrate, accrate_update_freq))
    #    accrate /= accrate_update_freq
    #    delta_t = delta_t * (accrate/idrate)

    #outfile.write(str(q_current))


delta_H_outfile.close()
        
    
