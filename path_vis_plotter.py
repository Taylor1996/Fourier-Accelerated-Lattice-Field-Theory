import re
import matplotlib.pyplot as plt
import numpy as np
import itertools
import math
import sys 


N = 10
infile = open("path_vis_outfile.txt", "r")
patharr = []
for line in infile:
    patharr.append(float(line))


plt.plot( patharr[0:N])
plt.plot( patharr[N:2*N])
plt.plot( patharr[2*N:3*N])
plt.show()