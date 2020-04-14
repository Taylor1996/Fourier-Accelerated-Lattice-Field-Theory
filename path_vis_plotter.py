import re
import matplotlib.pyplot as plt
import numpy as np
import itertools
import math
import sys 


N = 11
infile = open("path_vis_outfile_m_0.1.txt", "r")
#infile = open("fa_path_vis.txt", "r")
patharr = []
for line in infile:
    patharr.append(float(line))


plt.plot(patharr[0:N]    , np.arange(N), color='blue' )
plt.plot(patharr[N:2*N]  , np.arange(N), color = 'red' )
plt.plot(patharr[2*N:3*N], np.arange(N), color = 'green')
for n in np.arange(N):
    plt.plot(np.linspace(-5,5),n*np.ones(len(np.linspace(-5,5))), color='black')
plt.plot(np.zeros(N), np.arange(N), color='black')
plt.title('Paths in harmonic oscillator Markov chain (unaccelerated)')
plt.xlabel('$x$')
plt.ylabel('t')
plt.show()