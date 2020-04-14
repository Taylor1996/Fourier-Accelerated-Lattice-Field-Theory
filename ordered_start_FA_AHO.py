import numpy as np
import matplotlib.pyplot as plt
import sys
np.set_printoptions(threshold=sys.maxsize)


"""
Code for HMC simulation of anharmonic oscillator 

run program using following syntax
python HMC.py "output_file_name" "timestep_val"


TODO: Could add equilibrate function which alters accrate until within a satisfactory neighbourhood of idrate before stopping accrate modification and 
then taking the system to equilibrium. Only after this are measurements taken. Right now, I'm not manually modifying the accrate. I just run HMC and discard the first 100 or so trajectories in post processing.
"""

if(len(sys.argv) != 5):
    print("usage: python file_name lf_steps delta_t m_eff lambda")
    exit()

# Generate an initial state. setting = 0 for now

a = 1.0
N_tau = 31 # size of the lattice (i.e., the length of the harmonic oscillator's path in time)
N_trajs = 20000 # number of trajectories (equivalent to number of measurements of observables)
N_equi_trajs = 200
delta_t_equi = 0.01
lf_steps_equi = 100
t = 0
mu_sqr = -4.0
lf_steps = int(sys.argv[1])  # leapfrog steps in a trajectory
delta_t = float(sys.argv[2]) # leapfrog timestep in computer time. Passed as cmd line arg for use in bash
#m = 0.5 # mass of particle. Do an array divided by array element-wise if you want multiple particles with different masses.
m = 1.0
m_eff = float(sys.argv[3])
omega = 1.0 
#omega = -1.1 # natural frequency of oscillator
Lambda = float(sys.argv[4])
f_sqr = 2.0
accrate = 0.0
idrate = 0.8 # ideal acceptance rate

M_k_inv = np.zeros((N_tau, N_tau))
M_k= np.zeros((N_tau, N_tau))

for k in range(N_tau):
    M_k[k][k] = (m_eff**2 + 4*(np.sin(np.pi*k/N_tau))**2)
    M_k_inv[k][k] = 1/(m_eff**2 + 4*(np.sin(np.pi*k/N_tau))**2)
    

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
    #q_f = q**2-f_sqr
    #S = (1/(2.0*a))*m*sqr_term +  a*Lambda*np.sum(q_f**2)  # diff form of potential for comparison with lit
    S =  a*((1/2.0*a**2)*m*sqr_term + (1/2.0)*m**2*np.sum(q**2) + Lambda*np.sum(q**4))
    return S


def grad_S(q):
    q_forward  = np.array([q[(i+1)%N_tau] for i in range(len(q))])
    q_backward  = np.array([q[(i-1)%N_tau] for i in range(len(q))])
    grad_S = a*((m/a**2)*(-q_forward + 2*q - q_backward) + m**2*q + 4*Lambda*q**3)
    #grad_S = (m/a)*(-q_forward + 2*q - q_backward) + 4*a*Lambda*q*(q**2 - f_sqr)
    # gradient of whole action!
    return grad_S


def HMC(q_current, delta_t, lf_steps,sweep):
    global accrate 
    
     #p = np.random.normal(0,1,N_tau)
    #p = np.random.multivariate_normal(np.zeros(N_tau), M)
    #p = np.real(np.fft.ifft(np.random.multivariate_normal(np.zeros(N_tau), FT_M)))
    aug_cov = N_tau*M_k
    #aug_cov[0][0] += 0.5*m**2*N_tau
    #aug_cov = N_tau
    ##############################################################
    # aug_cov = np.zeros((N_tau, N_tau))                             #
    # for k in range(1,N_tau):                                     #
    #     aug_cov = 0.5*N_tau*m**2 + N_tau*(1/2.0)* (np.sin(k*np.pi)) #
    ##############################################################
    p_t = np.zeros(int((N_tau-1)/2), dtype='complex')
    p_t_h = np.zeros(int((N_tau-1)/2), dtype='complex')

    p_t_0 = np.random.normal(0, np.sqrt(N_tau)*m_eff)
    #print(p_t_0)
    for l in range(1, int((N_tau+1)/2)):
        #print(l)
        cov = np.sqrt(0.5*(m_eff**2 + 4*(np.sin(np.pi*l/N_tau))**2))
        p_t_r = np.random.normal(0, np.sqrt(N_tau)*cov)
        p_t_i = np.random.normal(0, np.sqrt(N_tau)*cov)
        p_t_h[l-1] = p_t_r + 1.j*p_t_i 
    p_t_end = np.conj(np.flip(p_t_h))
    #print(p_t_h)
    #print(p_t_end)
    p_t = np.concatenate((p_t_h, p_t_end))
    p_t = np.concatenate((np.array([p_t_0]), p_t))
    #print(p_t)
    p = np.real(np.fft.ifft(p_t))

   
    q = q_current
    p_current = p
    
    # leapfrog
    p = p - (grad_S(q) * (delta_t/2.0))
    
    for i in range(lf_steps):
        
        p_k = np.fft.fft(p)
    
        #q = q + delta_t*(p) # m should not be present here
        q = q + delta_t*np.real(np.fft.ifft(M_k_inv.dot(p_k)))

        if i != lf_steps-1:
            p = p - (grad_S(q)* delta_t)

    # final leap frog step for momentum
    p = p - (grad_S(q) * (delta_t/2.0))

    if accept(q_current,p_current,q,p):
        #print("accept")
        q_current = q
        p_current = p
        accrate +=1
    else:
        pass
        #print("reject")

    #print("accrate = {}, delta_t = {}".format(accrate/(sweep+1), str(delta_t)))

    
    # should return modified delta_t so as to get accrate nearer idrate

    return q_current,p_current,accrate
    

def HMC_E_track(q_current, delta_t, lf_steps,sweep):
    global accrate 
    
     #p = np.random.normal(0,1,N_tau)
    #p = np.random.multivariate_normal(np.zeros(N_tau), M)
    #p = np.real(np.fft.ifft(np.random.multivariate_normal(np.zeros(N_tau), FT_M)))
    aug_cov = N_tau*M_k
    #aug_cov[0][0] += 0.5*m**2*N_tau
    #aug_cov = N_tau
    ##############################################################
    # aug_cov = np.zeros((N_tau, N_tau))                             #
    # for k in range(1,N_tau):                                     #
    #     aug_cov = 0.5*N_tau*m**2 + N_tau*(1/2.0)* (np.sin(k*np.pi)) #
    ##############################################################
    p_t = np.zeros(int((N_tau-1)/2), dtype='complex')
    p_t_h = np.zeros(int((N_tau-1)/2), dtype='complex')

    p_t_0 = np.random.normal(0, np.sqrt(N_tau)*m_eff)
    #print(p_t_0)
    for l in range(1, int((N_tau+1)/2)):
        #print(l)
        cov = np.sqrt(0.5*(m_eff**2 + 4*(np.sin(np.pi*l/N_tau))**2))
        p_t_r = np.random.normal(0, np.sqrt(N_tau)*cov)
        p_t_i = np.random.normal(0, np.sqrt(N_tau)*cov)
        p_t_h[l-1] = p_t_r + 1.j*p_t_i 
    p_t_end = np.conj(np.flip(p_t_h))
    #print(p_t_h)
    #print(p_t_end)
    p_t = np.concatenate((p_t_h, p_t_end))
    p_t = np.concatenate((np.array([p_t_0]), p_t))
    #print(p_t)
    p = np.real(np.fft.ifft(p_t))

   
    q = q_current
    p_current = p
    
    
    U_current = S(q_current)
    ft_p_current = np.fft.fft(p_current)
    K_current = (1/2.0)*(1.0/N_tau)*np.real(np.conj(ft_p_current).dot(M_k_inv.dot(ft_p_current))) # note that m doesn't appear here b/c the Hamiltonian dynamics is a completely artificial dynamics
    H_current = U_current + K_current
    
    s_outfile.write("{}\n".format(S(q)))

    # leapfrog
    p = p - (grad_S(q) * (delta_t/2.0))
    p_old = p
    
    for i in range(lf_steps):
        
        p_k = np.fft.fft(p)
    
        #q = q + delta_t*(p) # m should not be present here
        q = q + delta_t*np.real(np.fft.ifft(M_k_inv.dot(p_k)))
        s_outfile.write("{}\n".format(S(q)))


        if i != lf_steps-1:
            p = p - (grad_S(q)* delta_t)
            p_av = (p_old + p)/2.0
            H = 0.5*(p_av.dot(p_av)) + S(q)
            H_traj_outfile.write("{}\n".format(H))
        
        p_old = p

    # final leap frog step for momentum
    p = p - (grad_S(q) * (delta_t/2.0))

    p_av = (p_old + p)/2.0
    H = 0.5*(p_av.dot(p_av)) + S(q)
    H_traj_outfile.write("{}\n".format(H))    
   

    U_proposed = S(q)
    ft_p = np.fft.fft(p)
    K_proposed = (1/2.0)*(1.0/N_tau)*np.real(np.conj(ft_p).dot(M_k_inv.dot(ft_p))) # actually, I'm a bit confused about the above - see Neal???
    #print("K_prop = {}".format(K_proposed))
    
    H_proposed = U_proposed + K_proposed
    delta_H = H_proposed - H_current 
    delta_H_outfile.write("{}\n".format(delta_H))

    
    # should return modified delta_t so as to get accrate nearer idrate

    return q_current,p_current,accrate
    

def accept(q_current,p_current,q,p):

    # Calculation of new Hamiltonian
    
    U_current = S(q_current)
    ft_p_current = np.fft.fft(p_current)
    #K_current = (1/2.0)*p_current.dot(M_inv.dot(p_current)) # note that m doesn't appear here b/c the Hamiltonian dynamics is a completely artificial dynamics
    #print("Should have K_curr = {}".format(K_current))

    K_current = (1/2.0)*(1.0/N_tau)*np.real(np.conj(ft_p_current).dot(M_k_inv.dot(ft_p_current))) # note that m doesn't appear here b/c the Hamiltonian dynamics is a completely artificial dynamics
    #print("K_curr = {}".format(K_current))
    #K_proposed = (1/2.0)*p.dot(M_inv.dot(p))
    U_proposed = S(q)
    ft_p = np.fft.fft(p)
    K_proposed = (1/2.0)*(1.0/N_tau)*np.real(np.conj(ft_p).dot(M_k_inv.dot(ft_p))) # actually, I'm a bit confused about the above - see Neal???
    #print("K_prop = {}".format(K_proposed))
    H_current = U_current + K_current
    H_proposed = U_proposed + K_proposed
    
    delta_H = H_proposed - H_current 
    #delta_H_outfile.write(str(delta_H) +" ") # write out values of delta_H to check for correct behaviour of HMC 
    # if the new state decreases the exponential of the Hamiltonian, accept the new state. Otherwise, reject with probability 1-e^(-delta_H)
    if (np.random.uniform() < np.exp(-delta_H)):
        # accept
        return True 
    
    return False
    

q_current = init("cold", N_tau)


s_outfile       = open("{}_lambda_ordered_start_m_eff_{}_dt_{}_lf_steps_{}_av_s_track_FA_AHO.txt".format(Lambda, m_eff, delta_t, lf_steps), "w")
delta_H_outfile = open("{}_lambda_ordered_start_m_eff_{}_dt_{}_lf_steps_{}_delta_H_FA_AHO.txt".format   (Lambda, m_eff, delta_t, lf_steps),"w")
H_traj_outfile  = open("{}_lambda_ordered_start_m_eff_{}_dt_{}_lf_steps_{}_h_tracking_FA_AHO.txt".format(Lambda, m_eff, delta_t, lf_steps), "w")

    


for sweep in range(100):
    # the current plotting code I have deals with paths, not averaged paths
    #U = S(q_current)
    #grad_S = grad_S(q_current)
    print(sweep)    
    HMC_E_track(q_current, delta_t, lf_steps,sweep) # track S vs t and H vs t 
    #if ((sweep % accrate_update_freq ==0) and sweep != 0):
    #    print("accrate={}, acc_up_freq={}".format( accrate, accrate_update_freq))
    #    accrate /= accrate_update_freq
    #    delta_t = delta_t * (accrate/idrate)

    #outfile.write(str(q_current))


delta_H_outfile.close()
s_outfile.close()
H_traj_outfile.close()