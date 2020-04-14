import numpy as np
import matplotlib.pyplot as plt
import sys
np.set_printoptions(threshold=sys.maxsize)


"""
run program using following syntax
python HMC.py "output_file_name" "timestep_val"


TODO: Could add equilibrate function which alters accrate until within a satisfactory neighbourhood of idrate before stopping accrate modification and 
then taking the system to equilibrium. Only after this are measurements taken. Right now, I'm not manually modifying the accrate. I just run HMC and discard the first 100 or so trajectories in post processing.
"""

# Generate an initial state. setting = 0 for now

a=1.0
N_tau = 11   # size of the lattice in time direction 
N_x = N_tau # size of the lattice in the x position
N_trajs = 1000 # number of trajectories (equivalent to number of measurements of observables)



t = 0
lf_steps = int(sys.argv[2])  # leapfrog steps in a trajectory
delta_t = 1.0/ float(lf_steps) # leapfrog timestep in computer time. Passed as cmd line arg for use in bash
m = 1.0 # mass of particle. Do an array divided by array element-wise if you want multiple particles with different masses.
M=m
omega = 1.0 # natural frequency of oscillator
accrate = 0.0
idrate = 0.8 # ideal acceptance rate

S_out_current = 0.0  # global variables for tracking and outputting the action 
S_out = 0.0

S_outfile = open("S_avg_scalar_field.txt","w")

def init(mode, N_tau,N_x):
    # either hot, cold or approx equilibrated
    # could maybe use numpy.roll to get NNs for entire axes at a time
    if mode=="cold":
        return np.zeros((N_tau,N_x))
    elif mode=="hot":
        return np.random.uniform(-0.5, 0.5, (N_tau,N_x)) # init lattice sites to random numbers between -1/2 and 1/2

def format_phi(phi):
    phi_str = ""
    for i in range(N_tau):
        for j in range(N_x):
            phi_str += str(phi[i][j])
            phi_str += " "
    return phi_str



def S(phi):
    # S is the action which acts like the potential in HMC: H = K + S
    # define T to be x_i+1 - x_i bit
    """
    Loop over i and j. For a given lattice site, calculate (phi(rhs NN) - phi(i,j))^2 + (phi(up NN) - phi(i,j))^2 + m**2/2 * np.linalg.norm(phi(i,j))**2
    Add this to S 
    """
    S = 0.0
    for i in range(N_tau):
        for j in range(N_x):
            K = (1/2.0)*( (phi[NN_list[i][j][2][0], NN_list[i][j][2][1]] - phi[i][j])**2 + (phi[NN_list[i][j][1][0], NN_list[i][j][1][1]] - phi[i][j])**2)
            S += (K+ (a**2)*(m/2.0)*phi[i][j]**2)
    return S

def grad_S(phi):
    global NN_list
  
    # q_forward  = np.array([q[(i+1)%N_tau] for i in range(len(q))])
    #q_backward  = np.array([q[(i-1)%N_tau] for i in range(len(q))])
    #    grad_S = (m)*(-q_forward + 2*q - q_backward) + m*omega**2*q
    grad_S = np.zeros((N_tau, N_x))
    for i in range(N_tau):
        for j in range(N_x):
            grad_S[i][j] = -phi[NN_list[i][j][0][0], NN_list[i][j][0][1]] - phi[NN_list[i][j][1][0],NN_list[i][j][1][1]] - phi[NN_list[i][j][2][0],NN_list[i][j][2][1]] - phi[NN_list[i][j][3][0], NN_list[i][j][3][1]] + 4*phi[i][j] + (a**2)*(m**2)*phi[i][j]
    # gradient of whole action!
    return grad_S


def HMC(phi_current, delta_t, lf_steps,sweep):
    global accrate 
    #pi = np.random.normal(0,1/2,(N_tau,N_x))
    
    aug_cov = np.zeros((N_tau**2,N_tau**2)) # the cov matrix is (N_tau**2)x(N_tau**2) with the 1/M**2 + sin**2(...) +sin**2(...) factors along diagonal
    
    for l_t in range(N_tau): 
        for l_x in range(N_tau):
            aug_cov[l_t*N_tau + l_x][l_t*N_tau + l_x]  = (L**2)*(M**2 + 4*(np.sin(np.pi*l_t/np.pi))**2 + 4*(np.sin(np.pi*l_x/np.pi))**2)

    pi_rav = np.random.multivariate_normal(np.zeros(N_tau**2), aug_cov)
    pi = np.reshape(pi_rav, (N_tau,N_tau)) # perform mapping to 2d array

    phi = phi_current
    pi_current = pi
    
    # leapfrog
    pi = pi - (grad_S(phi) * (delta_t/2.0))
    
    for i in range(lf_steps):
        pi_k = np.fft.fft(pi)
        phi = phi + 2*delta_t*np.real(np.fft.ifft(np.reshape(np.array([pi_k[t_k][x_k]/(1) for t_k,x_k in np.ndindex(N_tau,N_tau) ]), (N_tau,N_tau))))
        prod = (N_tau**2)*np.array([pi_k[(N_tau-t_k)%N_tau][(N_tau-x_k)% N_tau]/(M**2 + 4*(np.sin(np.pi*t_k/N_tau))**2 + 4*(np.sin(np.pi*x_k/N_tau))**2) for t_k,x_k in np.ndindex(N_tau,N_tau)])
        #print(prod)
        #print("reshaped\n")
        #print(np.reshape(prod, (N_tau,N_tau)))
        #phi = phi + delta_t*np.real(np.fft.ifft(np.reshape(prod, (N_tau,N_tau)))) # prod is created as 1D array. Need to get back 2d arr indexed by t_k,x_k
        if i != lf_steps-1:
            pi = pi - (grad_S(phi)* delta_t)

    # final leap frog step for momentum
    pi = pi - (grad_S(phi) * (delta_t/2.0))

    if accept(phi_current,pi_current,phi,pi):
        print("accept")
        phi_current = phi
        pi_current = pi
        accrate +=1
        
        
    else:
        print("reject")

    print("accrate = {}, delta_t = {}".format(accrate/(sweep+1), str(delta_t)))
    outfile.write(format_phi(phi_current)) # output is now
    S_outfile.write((str(S_out_current) + " "))
    
    # should return modified delta_t so as to get accrate nearer idrate

    return phi_current,pi_current,accrate
    

def accept(phi_current,pi_current,phi,pi):
    
    global S_out
    global S_out_current 

    # Calculation of new Hamiltonian
    S_out_current = S(phi_current)
    U_current = S_out_current
    K_current = (np.sum(pi_current**2))/(2.0) # note that m doesn't appear here b/c the Hamiltonian dynamics is a completely artificial dynamics        
    S_out = S(phi)
    U_proposed = S_out
    K_proposed = (np.sum(pi**2))/(2.0) 

    H_current = U_current + K_current
    H_proposed = U_proposed + K_proposed
    
    delta_H = H_proposed - H_current 
    delta_H_outfile.write(str(delta_H) +" ") # write out values of delta_H to check for correct behaviour of HMC 
    # if the new state decreases the exponential of the Hamiltonian, accept the new state. Otherwise, reject with probability 1-e^(-delta_H)
    if (np.random.uniform() < np.exp(-delta_H)):
        # accept
        return True 
    
    return False

def getNNlist(N_tau, N_x):
    NN_list = [[[[0,0] for i in range(4)] for j in range(N_tau)] for k in range(N_x)]
    
    for i in range(N_tau):
        for j in range(N_x):
            #0,0 in bottom left corner of N_tau-N_x lattice
            # k=0 is NN to left, k=1 is NN above, k=2 is NN to right
            # Actually - I only need the RHS and above NNs
            NN_list[i][j][0] =  [(i-1) % N_tau,j]
            NN_list[i][j][1] =  [(i+1)%N_tau, j] # RHS NN
            NN_list[i][j][2] =  [i,(j+1)%N_x] # above NN
            NN_list[i][j][3] =  [i,(j-1)%N_x]
    return NN_list 


phi_current = init("cold", N_tau, N_x)
NN_list = getNNlist(N_tau, N_x) # build neighbour lists
outfile = open("2d_scalar_field_" + sys.argv[1],"w")
delta_H_outfile = open("FA_2d_scalar_field_hmc_delta_H_" + str(sys.argv[1]),"w")


for sweep in range(N_trajs):
    # the current plotting code I have deals with paths, not averaged paths
    #U = S(phi_current)
    #grad_S = grad_S(phi_current)
    phi_current,pi_current,accrate = HMC(phi_current, delta_t, lf_steps,sweep)
    #if ((sweep % accrate_update_freq ==0) and sweep != 0):
    #    print("accrate={}, acc_up_freq={}".format( accrate, accrate_update_freq))
    #    accrate /= accrate_update_freq
    #    delta_t = delta_t * (accrate/idrate)

    #outfile.write(str(phi_current))


outfile.close()
delta_H_outfile.close()
        
