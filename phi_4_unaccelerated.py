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

if len(sys.argv) != 2:
    print("Usage: python Fourier... num_steps")
    exit()

# Generate an initial state. setting = 0 for now

N_tau = 1001 # size of the lattice (i.e., the length of the harmonic oscillator's path in time)
N_trajs = 1000 # number of trajectories (equivalent to number of measurements of observables)
sweep = 0
Lambda = 0.5


t = 0

lf_steps = int(sys.argv[1])  # leapfrog steps in a trajectory
delta_t = 1.0/lf_steps  # leapfrog timestep in computer time. Passed as cmd line arg for use in bash
  
#delta_t = 1.0/lf_steps  # leapfrog timestep in computer time. Passed as cmd line arg for use in bash
m = 1.0 # mass of particle. Do an array divided by array element-wise if you want multiple particles with different masses.
omega = 0.1 # natural frequency of oscillator
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
#print(M) 
M*=-1.0

M += m**2*np.identity(N_tau)
#M*=0.5

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
        return np.zeros(N_tau) + np.array([1]*N_tau)
    elif mode=="hot":
        return np.random.uniform(-0.5, 0.5, N_tau) # init lattice sites to random numbers between -1/2 and 1/2



def S(q):
    # S is the action which acts like the potential in HMC: H = K + S
    # define T to be x_i+1 - x_i bit
    global M
    q_forward  = np.array([q[(i+1)%N_tau] for i in range(len(q))])
    q_backward  = np.array([q[(i-1)%N_tau] for i in range(len(q))])
    ##print(0.5*(-q_forward + 2*q - q_backward + (m**2)*q))
    S = 0.5*(-q.dot(q_forward) + 2*q.dot(q) - q_backward.dot(q) + (m**2)*q.dot(q)) + Lambda*np.sum(q**4)
    #q_shifted  = np.array([q[(i+1)%N_tau] for i in range(len(q))])
    #sqr_term = np.sum((q_shifted - q)**2)
    #S =  (1/2.0)*m*sqr_term + (1/2.0)*m*omega**2*np.sum(q**2) # total pot energy of lattice
    #print(S)
    #print(S)
    #S = M.dot(q)
    #print(S)
    #S = q.dot(S)
    #print(S)
    return S

def grad_S(q):
    q_forward  = np.array([q[(i+1)%N_tau] for i in range(len(q))])
    q_backward  = np.array([q[(i-1)%N_tau] for i in range(len(q))])
    #print(q)
    grad_S = (0.5)*(-2*q_forward + 4*q -2*q_backward + 2*m**2*q) + 4*Lambda*q**3
    # gradient of whole action!
    return grad_S


def HMC(q_current, delta_t, lf_steps,sweep):
    global accrate 
    p = np.random.normal(0,1,N_tau)
    #p = np.random.multivariate_normal(np.zeros(N_tau))
    #p = multivariate_normal.rvs(mean=np.zeros(N_tau), cov=M)

    #p = np.fft.ifft(np.random.multivariate_normal(np.zeros(N_tau), N_tau*M_k) )
    # for FA, sample from multivariate normal distro with cov matrix M
    
    q = q_current
    S(q)    
    p_current = p
    
    # leapfrog
    p = p - (grad_S(q) * (delta_t/2.0))
    for i in range(lf_steps):
        #p_k = np.fft.fft(p)
        
        #q = q + delta_t*np.fft.ifft(((N_tau)*M_k.dot(p_k)))
        q = q + delta_t*p
        
       
        #q = q + delta_t*(0.5*np.transpose(M_inv).dot(p) + 0.5*M_inv.dot(p) )
        
        #q = q + delta_t*p
        if i != lf_steps-1:
            p = p - (grad_S(q)* delta_t)
            
    # final leap frog step for momentum
    p = p - (grad_S(q) * (delta_t/2.0))
   

    if accept(q_current,p_current,q,p):
        #print("accept")
        q_current = q
        p_current = p
        accrate +=1.0
    else:
        pass
        #print("reject")
    if sweep > 0:
        
        av_sum_q_sq = (np.sum(q)**2)
        av_sum_phi_sq_outfile.write("{}\n".format(av_sum_q_sq))
        av_x_outfile.write("{}\n".format(np.mean(q_current)))
        av_x_sq_outfile.write("{}\n".format(np.mean(q_current)**2))
        x_sq_file.write("{}\n".format(np.sum(q_current**2)))
        x0_xN_file.write("{}\n".format(q[0]*q[-1]))
        print("accrate = {}, delta_t = {}".format(accrate/(sweep+1), str(delta_t)))
        outfile.write(str(q_current))

    
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
        s_outfile.write("{}\n".format(S(q_current)))
        return True 
    
    s_outfile.write("{}\n".format(S(q)))
    return False
    
#q_current = np.arange(N_tau)

q_current = init("cold", N_tau)
av_sum_phi_sq_outfile   = open("m_{}_lambda_{}_equilib_aho_av_sum_phi_sq_outfile_unaccel.txt".format(m,Lambda), "w")
outfile                 = open("m_{}_lambda_{}_equilib_aho_unaccelerated_1d_field.txt".format(m,Lambda),"w")
s_outfile               = open("m_{}_lambda_{}_equilib_aho_unaccelerated_s_outfile.txt".format(m,Lambda), "w")
delta_H_outfile         = open("m_{}_lambda_{}_equilib_aho_delta_H_unaccel.txt".format(m,Lambda),"w")
av_x_sq_outfile         = open("m_{}_lambda_{}_equilib_aho_unacc_mag_av_x_sq_outfile.txt".format(m,Lambda), "w")
av_x_outfile            = open("m_{}_lambda_{}_equilib_aho_unacc_FFT_FA_mag_av_x_outfile.txt".format(m,Lambda), "w")
x_sq_file               = open("m_{}_lambda_{}_equilib_aho_x_sq_file.txt".format(m,Lambda), "w")
x0_xN_file              = open("m_{}_lambda_{}_equilib_aho_x1_xN_sq_file.txt".format(m,Lambda), "w")

for sweep in range(N_trajs):
    sweep +=1
    # the current plotting code I have deals with paths, not averaged paths
    #U = S(q_current)
    #grad_S = grad_S(q_current)
    q_current,p_current,accrate = HMC(q_current, delta_t, lf_steps,sweep)
    #if ((sweep % accrate_update_freq ==0) and sweep != 0):
    #    print("accrate={}, acc_up_freq={}".format( accrate, accrate_update_freq))
    #    accrate /= accrate_update_freq
    #    delta_t = delta_t * (accrate/idrate)

    #outfile.write(str(q_current))

av_sum_phi_sq_outfile.close()
outfile.close()
delta_H_outfile.close()
s_outfile.close()

