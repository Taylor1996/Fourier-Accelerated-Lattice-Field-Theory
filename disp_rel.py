import re
import matplotlib.pyplot as plt
import numpy as np
import itertools
import math
import sys 
plt.rcParams['text.usetex'] = True


t_arr = np.arange(-10001,10001)
y = np.sqrt(np.sin(0.001*t_arr)**2)

plt.plot(t_arr, y)
plt.xlabel("Hamiltonian dynamics time")
plt.ylabel(r"$ \langle S \rangle $")
plt.title(r"Intra-trajectory evolution of $\langle S \rangle $")
plt.show()