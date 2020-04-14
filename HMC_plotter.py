"""
Plotter for HMC data
"""

import re
import matplotlib.pyplot as plt
import numpy as np
import itertools
import math

#delta_H_infile = open("delta_H_hmc_output_new_S.txt","r")
#for line in delta_H_infile:
#    tokensH = line.split()
#
#
#delta_H_list = [float(del_H) for del_H in tokensH]
#delta_H_list
#np.mean(np.exp(-np.array(delta_H_list)))


#infile = open("correct_mins_1000_0.25_10000.txt", "r")
infile = open("Creutz_params_303_0.25_1000000.txt", "r")
#infile = open("aho_N_t_100_a_0.25_N_trajs_100000incl_a.txt", "r")
#infile = open("aho_creutz_dh_dt_N_t_200_a_0.25_N_trajs_100000_.txt", "r")
s = infile.read()
tokens = s.split('[')[1:]

# we discard tokens[0]
for i,path in enumerate(tokens):
    new_path = re.sub('\]|\n', '', path)
    tokens[i] = new_path
#x_data = [path for path in tokens[1:]]
pathList = [[0]*len(tokens[0])]*len(tokens)
for i in range(len(tokens)):
    pathList[i] = [float(x_i) for x_i in tokens[i].split()]
    
# discard first 100 steps as thermalization
pathList = pathList[100:]

N_tau = len(pathList[0])
#print(N_tau)
#print(len(pathList))
#print(pathList)

# plotting <x> and <x^2> behaviour

N_sep = 1 # config sample frequency

##x_avg = [(1/N_tau)*np.mean(np.array(path)) for path in pathList[::N_sep]]
##x2_avg = [(1/N_tau)*np.mean(np.array(path)**2) for path in pathList[::N_sep]]
## the below are for plotting HMC, the above is for plotting ordinary standard MCMC
#x_avg = [(1/N_tau)*np.mean(np.array(path)) for path in pathList]
#x2_avg = [(1/N_tau)*np.mean(np.array(path)**2) for path in pathList]
#print(np.mean(np.array(x2_avg)))
#fig, (ax1,ax2) = plt.subplots(1,2) 
## ax1.xlabel("Config number")
#ax1.set_title(r'$\langle x^2 \rangle$')
## ax1.plot([N_sep*i for i in range(int(len(pathList)/N_sep))], x2_avg) #per site
#ax1.plot([i for i in range(len(x_avg))], x2_avg) #per site
#ax2.set_title(r'$\langle x \rangle$')
##ax2.xlabel("Config number")
##ax2.ylabel(r'$\langle x \rangle')
#ax2.plot([i for i in range(int(len(pathList)/N_sep))], x_avg) # per site 
#plt.show()

#print(len(x_avg))


#psi_0 = np.histogram([np.mean(path) for path in pathList],bins='auto', density=True)
# REMEMBER THAT WE DON'T HAVE TO SQUARE psi_0 BECAUSE IN FORMING A HISTOGRAM, WE GET A PROBABILITY DENSITY FOR x, WHICH IS 
# PRECISELY WHAT |psi_0|^2 is. 
psi_0 = np.histogram([x_i for path in pathList for x_i in path],bins="auto", density=True)
# we obtain the average of each path to ploAssertionErrort psi_0 and obtain the other observables
 


plt.title('Ground state wavefunction of the anharmonic oscillator')
plt.xlabel('x')
plt.ylabel('$|\psi_0(x)|^2$')
# I should normalize btw
# plt.scatter(list(psi_0[1])[1:], list(psi_0[0]**2))
plt.scatter(list(psi_0[1])[1:], list(psi_0[0]))
#x = np.linspace(-4,4,100)
#y = (1/(np.sqrt(np.pi)))*np.exp(-x**2)
#plt.plot(x,y)
plt.show()

# Obtain Jacknife estimates of the errors on <x>, <x^2> 



################################################################################
# N = len(O)                                                                   #
# B = int(0.01*N)                                                              #
#                                                                              #
# ############################################################################ #
# # if(N%B != 0.0):                                                          # #
# #     # B doesn't divide N. ditch first few O_i until we have divisibility # #
# #     O = O[(N/B - math.floor(N/B))* B):]                                  # #
# ############################################################################ #
#                                                                              #
# N_B = int(N/B)                                                               #
#                                                                              #
# # obtain <x>, the average over the gathered configs                          #
# tot = np.sum(O)                                                              #
# avg = tot/N                                                                  #
#                                                                              #
#                                                                              #
# o = np.zeros(N_B)                                                            #
# jack_estim = np.zeros(N_B)                                                   #
# for i in range(N_B):                                                         #
#     o[i] = np.sum(O[i*B:(i+1)*B])                                            #
#     jack_estim[i] = (1/(N-B)) * (tot - o[i])                                 #
#                                                                              #
# #np.linalg.norm(o-avg)**2                                                    #
# jack_var = ((N_B-1)/(N-B))*(np.linalg.norm(jack_estim-avg)**2)               #
#                                                                              #
# error = math.sqrt(jack_var)                                                  #
#                                                                              #
# print(avg)                                                                   #
# error                                                                        #
################################################################################
