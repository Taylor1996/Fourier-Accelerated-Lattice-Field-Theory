import numpy as np
import matplotlib.pyplot as plt


def S(q):
    # S is the action which acts like the potential in HMC: H = K + S
    # define T to be x_i+1 - x_i bit
    q_shifted  = np.array([q[(i+1)%N_tau] for i in range(len(q))])
    sqr_term = np.sum((q_shifted - q)**2)
    U =  (1/2.0)*m*sqr_term + (1/2.0)*m*omega**2*np.sum(q**2) # total pot energy of lattice
    return U

def grad_S(q):
    q_forward  = np.array([q[(i+1)%N_tau] for i in range(len(q))])
    q_backward  = np.array([q[(i-1)%N_tau] for i in range(len(q))])
    grad_S = (m)*(-q_forward + 2*q - q_backward) + m*omega**2*q
    # gradient of whole action!
    return grad_S

def HMC(q_current, delta_t, lf_steps,sweep):
    global accrate 
#   global p_main 
#    p = p_main
    p = np.random.normal(0,1,N_tau)
    q = q_current
    p_current = p
    
    # leapfrog
    p = p - (grad_S(q) * (delta_t/2.0))
    
    for i in range(lf_steps):
        q = q + delta_t*(1/m)*(p)
        if i != lf_steps-1:
            p = p - (grad_S(q)* delta_t)

    # final leap frog step for momentum
    p = p - (grad_S(q) * (delta_t/2.0))

    
    U_current = S(q_current)
    K_current = (np.sum(p_current**2))/(2.0)
    U_proposed = S(q)
    K_proposed = (np.sum(p**2))/(2.0)

    H_current = U_current + K_current
    H_proposed = U_proposed + K_proposed
    #print("{} {}".format(H_proposed, H_current))
    delta_H = H_proposed - H_current
    # acceptance step?
    return(delta_H)
    
N_tau = 100
m = 1.0
omega = 1.0
# p_main = np.random.normal(0,1,N_tau)
timeList = list([delta_t/1000.0 for delta_t in range(10,1000,10)])
numSamples = 1000
delta_H_list = []
for delta_t in timeList:
    delta_H_avg = 0.0
    # for each timestep, average delta_H and plot 
    q_current = np.zeros(N_tau)
    # print(int(1/delta_t))
    for i in range(numSamples):
        delta_H_avg += HMC(q_current, delta_t,int(1/delta_t),i) # keep delta_t*num_steps constant
    delta_H_list.append(delta_H_avg/numSamples)


timeListSqr = [dt**2 for dt in timeList]
#timeListQuart = [dt**4 for dt in timeList]
#for i in range(len(timeListSqr)):
#    print("dt={} delta_H={}\n ".format(timeListSqr[i], delta_H_list[i]))

#plt.plot(timeListSqr, delta_H_list)
plt.plot(timeListSqr, delta_H_list)
#plt.plot(timeList, delta_H_list)
plt.xlabel("dt**2")
plt.ylabel("dH")
plt.show()
    
