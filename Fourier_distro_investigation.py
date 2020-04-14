import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal
from mpl_toolkits.mplot3d import Axes3D

#Parameters to set
N_tau = 2
m=1.0

M = np.zeros((N_tau,N_tau))
for i in range(N_tau):
    for j in range(N_tau):
        if i==j:
            M[i][j] = -2
        if (j == ((i+1) % N_tau)):
            M[i][j] = 1
        if (j == ((i-1) % N_tau)):
            M[i][j] = 1
M*=-1.0
M += m**2*np.identity(N_tau)
M*=0.5
M_inv = np.linalg.inv(M)

#Create grid and multivariate normal
x = np.linspace(-10,10,500)
y = np.linspace(-10,10,500)
X, Y = np.meshgrid(x,y)
pos = np.empty(X.shape + (2,))
pos[:, :, 0] = X; pos[:, :, 1] = Y
rv = multivariate_normal(np.zeros(N_tau), np.array([[1,0],[0,1]]))
print(rv.rvs())
exit()
#Make a 3D plot
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, rv.pdf(pos),cmap='viridis',linewidth=0)
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.view_init(30, 10)
#ax.margins(1e-5,1e-5,1e-5)
plt.show()

## rotate the axes and update
#for angle in range(0, 360):
#    ax.view_init(30, angle)
#    plt.draw()
#    plt.pause(.001)