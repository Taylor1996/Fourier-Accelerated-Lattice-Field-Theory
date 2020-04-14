import matplotlib.pyplot as plt
import numpy as np 

infile = open("m_1_fourier_acc_s_1d_scalar_field_new_sampling.txt","r")
data = infile.read().split()
s_arr = np.array([float(s) for s in data])[500:]
print(np.mean(s_arr))

#plt.show()
#plt.plot(s_arr)