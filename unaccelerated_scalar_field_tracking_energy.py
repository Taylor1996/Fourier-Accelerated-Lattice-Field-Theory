import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.stats import multivariate_normal
np.set_printoptions(threshold=sys.maxsize)


"""
run program using following syntax
python HMC.py "output_file_name" "timestep_val"


TODO: Could add equilibrate function which alters accrate until within a satisfactory neighbourhood of idrate before stopping accrate modification and 
then taking the system to equilibrium. Only after this are measurements taken. Right now, I'm not manually modifying the accrate. I just run HMC and discard the first 100 or so trajectories in post processing.
"""

if len(sys.argv) != 3:
    print("Usage: python Fourier... outfile_name num_trajectories")
    exit()

# Generate an initial state. setting = 0 for now

N_tau = 31 # size of the lattice (i.e., the length of the harmonic oscillator's path in time)
N_trajs = 1 # number of trajectories (equivalent to number of measurements of observables)
sweep = 0

t = 0

#lf_steps = int(sys.argv[2])  # leapfrog steps in a trajectory
lf_steps = 10
delta_t = 0.1
#delta_t = 1.0/lf_steps  # leapfrog timestep in computer time. Passed as cmd line arg for use in bash
m = 1.0 # mass of particle. Do an array divided by array element-wise if you want multiple particles with different masses.
omega = 1.0 # natural frequency of oscillator
accrate = 0.0
idrate = 0.8 # ideal acceptance rate

M = np.zeros((N_tau,N_tau))
for i in range(N_tau):
    for j in range(N_tau):
        if i==j:
            M[i][j] = -2
        if (j == ((i+1) % N_tau)):
            M[i][j] = 1
        if (j == ((i-1) % N_tau)):
            M[i][j] = 1
print(M) 
M*=-1.0

M += m**2*np.identity(N_tau)

#print(M)
#M_inv = np.linalg.inv(M)
#print(M_inv.dot(M))
#exit()
#print(M)
#print()
#print(M_inv)
#print( np.linalg.inv(M_inv))
#print("M_inv=")
#print(M_inv)

#A_kl = 0.0 + 0.j
FT_M = np.zeros((N_tau,N_tau), dtype=complex)

#for k in range(N_tau):
#    for l in range(N_tau):
#        A_kl = 0.0
#        for n in range(N_tau):  
#            for m in range(N_tau):
#                A_kl += M[m][n]*np.exp(-2.j*np.pi*(1/N_tau)*((m*k)+(n*l)))
#        FT_M[k][l] = A_kl
        #print("A[{}][{}] = {}".format(k,l,A_kl))       
#print(FT_M)               
#print(np.fft.fft(M_inv[1]))
#print(np.fft.fft(M_inv))


M_k_inv = np.zeros((N_tau, N_tau))
M_k= np.zeros((N_tau, N_tau))

for k in range(N_tau):
    M_k[k][k] = (m**2 + 4*(np.sin(np.pi*k/N_tau))**2)
    M_k_inv[k][k] = 1/(m**2 + 4*(np.sin(np.pi*k/N_tau))**2)
    
#print(M_k)    
#print(M_k_inv)

    
def init(mode, N_tau):
    # either hot, cold or approx equilibrated
    if mode=="cold":
        return np.zeros(N_tau)
    elif mode=="hot":
        return np.random.uniform(-0.5, 0.5, N_tau) # init lattice sites to random numbers between -1/2 and 1/2



def S(q):
    # S is the action which acts like the potential in HMC: H = K + S
    # define T to be x_i+1 - x_i bit
    global M
    q_forward  = np.array([q[(i+1)%N_tau] for i in range(len(q))])
    q_backward  = np.array([q[(i-1)%N_tau] for i in range(len(q))])
    #print(0.5*(-q_forward + 2*q - q_backward + (m**2)*q))
    S = 0.5*(-q.dot(q_forward) + 2*q.dot(q) - q_backward.dot(q) + (m**2)*q.dot(q))
    #print(S)
    #print(S)
    #S = M.dot(q)
    #print(S)
    #S = q.dot(S)
    #print(S)
    #print("S={}".format(S))
    return S

def grad_S(q):
    q_forward  = np.array([q[(i+1)%N_tau] for i in range(len(q))])
    q_backward  = np.array([q[(i-1)%N_tau] for i in range(len(q))])
    grad_S = (0.5)*(-2*q_forward + 4*q - 2*q_backward + 2*m**2*q)

    # gradient of whole action!
    return grad_S


def HMC(q_current, delta_t, lf_steps,sweep):
    global accrate 
    p = np.random.normal(0,1,N_tau)
    #p = np.random.multivariate_normal(np.zeros(N_tau), M)
    #p = multivariate_normal.rvs(mean=np.zeros(N_tau), cov=M)

    #p = np.fft.ifft(np.random.multivariate_normal(np.zeros(N_tau), N_tau*M_k) )
    # for FA, sample from multivariate normal distro with cov matrix M

    q = q_current
    p_current = p
    
    # leapfrog
    if sweep==101:
        delta_t = 0.1
        lf_steps = int(sys.argv[2])
        print("sweep 101")
        print(lf_steps)
    p = p - (grad_S(q) * (delta_t/2.0))
    p_old = p
    for i in range(lf_steps):
        #p_k = np.fft.fft(p)
        #print(i)
        #q = q + delta_t*np.fft.ifft(((N_tau)*M_k.dot(p_k)))
        q = q + delta_t*p
        
        #if sweep==0:
        #    s_traj_outfile.write("{} ".format(S(q)))
        #q = q + delta_t*(0.5*np.transpose(M_inv).dot(p) + 0.5*M_inv.dot(p) )
        
        #q = q + delta_t*p
        if i != lf_steps-1:
            p = p - (grad_S(q)* delta_t)
            p_av = (p_old + p)/2.0
            
            if sweep==101:
                H = 0.5*(p_av.dot(p_av)) + S(q)
                H_traj_outfile.write("{} ".format(H))
        p_old = p
    # final leap frog step for momentum
    p = p - (grad_S(q) * (delta_t/2.0))
    #p_av = (p_old + p)/2.0
    
    if sweep==101:
        H = 0.5*(p_av.dot(p_av)) + S(q)
        H_traj_outfile.write("{} ".format(H))
        #exit()
    if accept(q_current,p_current,q,p):
        #print("accept")
        q_current = q
        p_current = p
        accrate +=1.0
    else:
        pass
        #print("reject")

    #print("accrate = {}, delta_t = {}".format(accrate/(sweep+1), str(delta_t)))
    #outfile.write(str(q_current))

    
    # should return modified delta_t so as to get accrate nearer idrate

    return q_current,p_current,accrate
    

def HMC_single_traj(q_current, delta_t, lf_steps,sweep):
    global accrate 
    p = np.random.normal(0,1,N_tau)
    
    #p = np.random.multivariate_normal(np.zeros(N_tau), M)
    #p = multivariate_normal.rvs(mean=np.zeros(N_tau), cov=M)

    #p = np.fft.ifft(np.random.multivariate_normal(np.zeros(N_tau), N_tau*M_k) )
    # for FA, sample from multivariate normal distro with cov matrix M

    q = q_current
    p_current = p
    
    # leapfrog
    p = p - (grad_S(q) * (delta_t/2.0))
    p_old = p
    
    s_traj_outfile.write("{} ".format(S(q)))

    for i in range(lf_steps):
        #p_k = np.fft.fft(p)
        #print("sweep = {} step = {}".format(sweep, i))
        print("test")
        #q = q + delta_t*np.fft.ifft(((N_tau)*M_k.dot(p_k)))
        q = q + delta_t*p
        
        s_traj_outfile.write("{} ".format(S(q)))
        #q = q + delta_t*(0.5*np.transpose(M_inv).dot(p) + 0.5*M_inv.dot(p) )
        
        #q = q + delta_t*p
        if i != lf_steps-1:
            p = p - (grad_S(q)* delta_t)
            p_av = (p_old + p)/2.0
            H = 0.5*(p_av.dot(p_av)) + S(q)
            
            if sweep==0:
                H_traj_outfile.write("{} ".format(H))
        p_old = p
    # final leap frog step for momentum
    p = p - (grad_S(q) * (delta_t/2.0))
    p_av = (p_old + p)/2.0
    H = 0.5*(p_av.dot(p_av)) + S(q)
    if sweep==0:
        H_traj_outfile.write("{} ".format(H))

    

    #print("accrate = {}, delta_t = {}".format(accrate/(sweep+1), str(delta_t)))
    #outfile.write(str(q_current))

    
    # should return modified delta_t so as to get accrate nearer idrate

    return q_current,p_current,accrate
    

def accept(q_current,p_current,q,p):

    # Calculation of new Hamiltonian

    U_current = S(q_current)
    K_current = (1/2.0)*(np.sum(p_current**2)) # note that m doesn't appear here b/c the Hamiltonian dynamics is a completely artificial dynamics
    U_proposed = S(q)
    K_proposed = (1/2.0)*(np.sum(p**2)) # actually, I'm a bit confused about the above - see Neal???

    H_current = U_current + K_current
    H_proposed = U_proposed + K_proposed
    
    delta_H = H_proposed - H_current 
    delta_H_outfile.write(str(delta_H) +" ") # write out values of delta_H to check for correct behaviour of HMC 
    # if the new state decreases the exponential of the Hamiltonian, accept the new state. Otherwise, reject with probability 1-e^(-delta_H)
    if (np.random.uniform() < np.exp(-delta_H)):
        # accept
        s_outfile.write("{} ".format(S(q_current)))
        return True 
    
    s_outfile.write("{} ".format(S(q)))
    return False
    

q_current = init("cold", N_tau)
outfile = open("unaccelerated_1d_field_energy_tracking.dat".format(sys.argv[1]),"w")
s_outfile = open("ordered_unaccelerated_s_outfile_energy_tracking.txt", "w")
s_traj_outfile = open("ordered_var_t_s_traj_outfile_energy_tracking.txt", "w")
delta_H_outfile = open("delta_h_energy_tracking.dat","w")
H_traj_outfile = open("ordered_var_t_h_unaccel_trajenergy_tracking.dat", "w")

for sweep in range(200):
    #sweep +=1
    # the current plotting code I have deals with paths, not averaged paths
    #U = S(q_current)
    #grad_S = grad_S(q_current)
    q_current,p_current,accrate = HMC(q_current, delta_t, lf_steps,sweep)
    #if ((sweep % accrate_update_freq ==0) and sweep != 0):
    #    print("accrate={}, acc_up_freq={}".format( accrate, accrate_update_freq))
    #    accrate /= accrate_update_freq
    #    delta_t = delta_t * (accrate/idrate)

    #outfile.write(str(q_current))
 
# collect values of the action along a single trajectory for many repetitions of that same trajectory starting with different initial momenta but the same intial config.
#q_saved = q_current
q_saved = np.array([1.0]*N_tau)
delta_t = 0.001
lf_steps = 10000
for sweep in range(1):
        q_current,p_current,accrate = HMC_single_traj(q_saved, delta_t, lf_steps,sweep)



outfile.close()
delta_H_outfile.close()
s_outfile.close()
s_traj_outfile.close()
H_traj_outfile.close()
