#!/usr/bin/python
########################################################################################################################
# Differential Equations - Initial Value Problems: 1st order differential equation with 4th order Runge-Kutta method
# NUMERICAL METHODS: Using Matlab, Fourth Edition by John H. Mathews and Kurtis D. Fink rk4.m
# Fortran ve Python ile SayÄ±sal Fizik by Bekir KaraoÄŸlu
#
from math import *
def myRK4(fun,xb,h,N):
    """myRK4  Use 4th orderRunge-Kutta method to solve Differential Equations - Initial Value Problems
    %Input: fun    = (string) name of function
    %       xb     = x0 and y0 are the initial conditions
    %       h      = the step length
    %       N      = the number of points
    %Output: xd,yd = xd(yd) is the vector of abscissas(ordinates)
    """
    x=x0=xb[0]; y=y0=xb[1] # Assign initial values
    xd=[];yd=[] # Create arrays
    xd.append(x0) # Populate first array element for x with initial value
    yd.append(y0) # Populate first array element for y with initial value
    for i in range(N+1):
        k1 = h * f(x, y)
        k2 = h * f(x + 0.5*h, y + 0.5*k1)
        k3 = h * f(x + 0.5*h, y + 0.5*k2)
        k4 = h * f(x + h, y + k3)
        y= y + (k1 + 2*(k2 + k3) + k4)/6.0 # Calculate new y value at point i
        yd.append(y) # Populate next array element with current y value
        x=x+h        # Increment x with step length
        xd.append(x) # Populate next array element with current x value
    return(xd,yd)
# Start Modify Parameters
f = lambda x, y: x + y ; x0=0.0; y0=1.0; pfun = r'$\frac{dy}{dx}=x+y$' # ODE
n=10; h=0.1 # Vary the values for n and h to see their effects
verbose=1   # Producing detailed output for diagnostic purposes (1=yes/0=no)
# End Modify Parameters
x,yapprox=myRK4(f,[x0,y0],h,n) # Solved by 4th order Runge-Kutta method
yexact=[2*exp(x[i])-x[i]-1.0 for i in range(n+2)] # Exact Solution (Analytic)
#
from scipy.integrate import solve_ivp
sol= solve_ivp(f, [0, 1.1], [1.0], t_eval=x, rtol = 1e-8, atol = 1e-8) # Default is â€˜RK45â€™, which is the explicit Runge-Kutta method of order 5(4)
# print(sol.t)
# print(sol.y)
#
if verbose:
    print(' Step            RK4                 Exact                RK4-Exact             SciPy         ')
    print('  x               y                    y                    Error                 y           ')
    print('----- ---------------------- --------------------- --------------------- ---------------------')
    for i in range(n+1):
        print('%5.2f %21.16f %21.16f %21.16f  %21.16f' % (x[i], yapprox[i], yexact[i], abs(yapprox[i]-yexact[i]),sol.y[0,i]))
#
import numpy as np
import matplotlib.pyplot as plt
plt.title('Approximate and Exact Solution for Simple ODE: '+ pfun)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(x, yapprox, 'bo--', label='Approximate')
plt.plot(x, 2*np.exp(x)-x-1.0, 'g', label='Exact')
plt.plot(sol.t, sol.y[0], '-x', label='SciPy')
plt.grid()
plt.legend(loc='lower right')
plt.savefig('RK4.eps')
plt.show()
#
#  Step            RK4                 Exact                RK4-Exact             SciPy
#   x               y                    y                    Error                 y
# ----- ---------------------- --------------------- --------------------- ---------------------
#  0.00    1.0000000000000000    1.0000000000000000    0.0000000000000000     1.0000000000000000
#  0.10    1.1103416666666668    1.1103418361512953    0.0000001694846286     1.1103418365038888
#  0.20    1.2428051417013890    1.2428055163203395    0.0000003746189505     1.2428055171581294
#  0.30    1.3997169941250756    1.3997176151520065    0.0000006210269310     1.3997176170100418
#  0.40    1.5836484801613715    1.5836493952825408    0.0000009151211693     1.5836493990278593
#  0.50    1.7974412771936765    1.7974425414002564    0.0000012642065799     1.7974425476900568
#  0.60    2.0442359241838663    2.0442376007810177    0.0000016765971513     2.0442376098866673
#  0.70    2.3275032531935538    2.3275054149409531    0.0000021617473993     2.3275054266863835
#  0.80    2.6510791265846310    2.6510818569849350    0.0000027304003041     2.6510818708124289
#  0.90    3.0192028275601421    3.0192062223138993    0.0000033947537572     3.0192062374603896
#  1.00    3.4365594882703321    3.4365636569180902    0.0000041686477581     3.4365636726612259