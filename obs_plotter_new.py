import re
import matplotlib.pyplot as plt
import numpy as np
import itertools
import math
import sys 


N = int(sys.argv[1])
#infile = open("mcmc_defin_x2_av.txt", "r")
infile = open("trace_hmc_defin_x_av.txt", "r")
x2arr = []
for line in infile:
    x2arr.append(float(line))

print(np.mean(x2arr[100:]))    
plt.scatter(np.linspace(0,N-1,N), x2arr[0:N])
plt.show()