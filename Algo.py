# Heat Equation: implicit difference method

from matplotlib import pylab
import numpy as np
from scipy import linalg

#from scipy documentation for matrix computation scipy is more recommened than numpy as faster feature, unless make sure
#numpy is installed "correctly"

def initialCondition(x):
    return 4.0 * (1.0 - x) * x

N = 500  # x axis: number of steps
M = 500  # t axis: number of steps

T = 1.0
X = 1.0

nodes_u = np.zeros([N + 1, M + 1])  # initialize the U value matrix
xArray = np.linspace(0,1,N+1)
nodes_x_ini = np.linspace(0, 1, N + 1)  # initilise the x_initial values for initial condiitons
nodes_u[:, 0] = map(initialCondition, nodes_x_ini)

dx = X / N   #space step
dt = T / M   #time step
kappa = 1.0  #heat coefficient

rho = kappa * dt / (dx * dx) #parameter rho


# implicit method using tridiagonal matrix System
# Python Class Trigonal Matrix System can be utilized to sovle this problem

for k in range(0, M, 1):  # k only reachs M - 1, coz need to stop at t = T which is at index M
    # initilise the trigonal matrix
    mat_dig = np.zeros([N+1,N+1])
    for i in range(1,N,1):
        mat_dig[i][i] = 1.0 + 2. * rho
        mat_dig[i][i-1] = - rho
        mat_dig[i][i+1] = - rho
    mat_dig[0][0]=1.
    mat_dig[N][N]=1.

    #print mat_dig

    # initilise the right - hand side matrix
    rhs = nodes_u[:,k]  # nodes_u[:,k] = kth column's, index 1 to M-1 entries

    # solve the matrix equation Ac = rhs  , c = rhs/A, .solve(rhs)
    x = linalg.solve(mat_dig, rhs)
    nodes_u[:, k + 1] = x

#plot the graph
pylab.figure(figsize = (12,6))
pylab.plot(xArray, nodes_u[:,0])
pylab.plot(xArray, nodes_u[:,int(0.10/ dt)])
pylab.plot(xArray, nodes_u[:,int(0.20/ dt)])
pylab.plot(xArray, nodes_u[:,int(0.50/ dt)])
pylab.xlabel('$x$', fontsize = 15)
pylab.ylabel(r'$U(\dot, \tau)$', fontsize = 15)
pylab.title('one dimensional heat equation')
pylab.legend([r'$\tau = 0.$', r'$\tau = 0.10$', r'$\tau = 0.20$', r'$\tau = 0.50$'], fontsize = 15)
pylab.show()


#may use scipy to further improve algorithm from https://uqer.io/community/share/5534ad3ff9f06c8f33904689

#import scipy as sp
#from scipy.linalg import solve_banded

#for k in range(0, M):
#    udiag = - np.ones(N-1) * rho
#    ldiag =  - np.ones(N-1) * rho
#    cdiag =  np.ones(N-1) * (1.0 + 2. * rho)
#    mat = np.zeros((3,N-1))
#    mat[0,:] = udiag
#    mat[1,:] = cdiag
#    mat[2,:] = ldiag
#    rhs = U[1:N,k]
#    x = solve_banded ((1,1), mat,rhs)
#    U[1:N, k+1] = x
#    U[0][k+1] = 0.
#    U[N][k+1] = 0.
