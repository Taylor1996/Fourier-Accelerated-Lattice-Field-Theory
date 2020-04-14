import re
import matplotlib.pyplot as plt
import numpy as np
import itertools
import sys
np.set_printoptions(threshold=sys.maxsize)


"""
Make the lattice size N_tau an input
"""


#infile = open("hmc_output_new_S.txt", "r")
infile = open("aho_N_t_200_a_0.25_N_trajs_1000000_.txt", "r")
outfile = open("AHO_output_for_gamma.txt", "w")
# infile = open("output_more_runs_N_tau=1200_m=0.1_omega=0.1", "r")
s = infile.read()
tokens = s.split('[')[1:]
# we discard tokens[0]
counter = 0 
therm_len = 100


# we discard tokens[0]

new_paths = []

for i in range(1000):
    new_path = re.sub('\]|\n', '', tokens[i])
    new_paths.append(new_path)
    
#for i,path in enumerate(tokens):
#    new_path = re.sub('\]|\n', '', path)
#    tokens[i] = new_path
    
pathList = [[0]*len(new_paths[0])]*len(new_paths)
for i in range(len(new_paths)):
    pathList[i] = [float(x_i) for x_i in new_paths[i].split()]

#print(pathList)

# discard first 100 steps as thermalization
pathList = pathList[therm_len:]


N_tau = len(pathList[0])
print(N_tau)
print(len(pathList))


x_avg = [(1/N_tau)*np.mean(np.array(path)) for path in pathList]
x2_avg = [(1/N_tau)*np.mean(np.array(path)**2) for path in pathList]        


outfile.write(" ".join(map(str, x2_avg)))
