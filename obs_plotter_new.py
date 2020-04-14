import re
import matplotlib.pyplot as plt
import numpy as np
import itertools
import math
import sys 
plt.rcParams['text.usetex'] = True


N_end = int(sys.argv[1])
#infile = open("m_0.001_equilib_phi4_lambda_0.001_FA_HO_mag_av_x_sq_outfile.txt", "r")
infile = open("m_0.001_lambda_0.001_equilib_aho_unacc_mag_av_x_sq_outfile.txt ", "r") #m_0.01_equilib_av_sum_phi_sq_outfile_unaccel.txt
#infile = open("trace_hmc_defin_x_av.txt", "r")
#infile = open("trace_s_outfile_HO_no_FFA.txt", "r")
#infile = open("trace_av_sum_phi_sq.txt", "r")
x2arr = []
for line in infile:
    x2arr.append(float(line))
mean_val = np.mean(x2arr[100:])
print(mean_val)    

import math

O = x2arr[100:]


#def Jacknife_magnetisation(M, chi):
#    
#    M = np.array(M)
#    m_arr = []
#    sum_diff_sq = 0.0
#    for i in range(len(M)):
#        resamp_i = np.delete(M,i)
#        m_i = np.mean(resamp_i)
#        m_i_sqr = np.mean(resamp_i**2)
#        chi_i = (1/N)*(m_i_sqr - m_i**2)
#        sum_diff_sq += (chi_i - chi)**2
#    
#    sigma = np.sqrt(sum_diff_sq)
#    return sigma         
    

    
#N = len(O)
#B = int(0.1*N)

###########################################################################
# if(N%B != 0.0):                                                          #
     # B doesn't divide N. ditch first few O_i until we have divisibility #
     #O = O[(N/B - math.floor(N/B))* B):]                                  #
###########################################################################


#N_B = int(N/B)
#
## obtain <x>, the average over the gathered configs
#tot = np.sum(O)
#avg = tot/N
#
#
#o = np.zeros(N_B)
#jack_estim = np.zeros(N_B)
#for i in range(N_B):
#    o[i] = (1.0/B)*np.sum(O[i*B:(i+1)*B])
#    jack_estim[i] = (1/(N-B)) * (tot - o[i])
#    
##np.linalg.norm(o-avg)**2
#jack_var = ((N_B-1)/(N_B))*(np.linalg.norm(jack_estim-avg)**2)
#
#error = (1/np.sqrt(N))*math.sqrt(jack_var)
#
#print(avg)
#print(error)


plt.xlabel("Trajectory number")
plt.ylabel(r"$\frac{1}{N_{\tau}} \left( \sum \limits_{n=0} ^{N _{\tau} -1} \phi_n^2 \right)$")
plt.title(r"$\phi^4$ theory:Plot of $\frac{1}{N_{\tau}} \left(\sum \limits_{n=0} ^{N _{\tau} -1} \phi_n^2 \right)$ vs HMC iteration (Unaccelerated)")
#plt.figure(figsize=(20,20))
plt.scatter(np.linspace(0,N_end-1,N_end), np.array(x2arr[0:N_end]), s=1.5)

plt.show()