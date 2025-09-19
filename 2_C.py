#!/usr/bin/python
#
from math import *
import sys
import numpy as np
import matplotlib.pyplot as plt
def mynewton(fun,dfundx,x0,xtol,ftol,verbose):
    """mynewton  Use Newton's method to find a root of the scalar equation f(x) = 0
    % Synopsis:  r = mynewton(fun,dfundx,xb)
    %            r = mynewton(fun,dfundx,xb,xtol)
    %            r = mynewton(fun,dfundx,xb,xtol,ftol)
    %            r = mynewton(fun,dfundx,xb,xtol,ftol,verbose)
    %
    % Input: fun     = (string) name of function for which roots are sought
    %        dfundx  = (string) name of derivative of function
    %        x0      = initial guess
    %        xtol    = (optional) relative x tolerance.    Default:  xtol=5*eps
    %        ftol    = (optional) relative f(x) tolerance. Default:  ftol=5*eps
    %        verbose = (optional) print switch. Default: verbose=0, no printing
    %
    % Output:  r = root (or singularity) of the function in xb(1) <= x <= xb(2)
     """
    # eps = sys.float_info.epsilon; xtol = 5 * eps; ftol = 5 * eps; verbose = 0
    xeps = max(xtol,5*eps) # Smallest tolerances are 5*eps
    feps = max(ftol,5*eps) # Smallest tolerances are 5*eps
    print("Tolerance values in x and f(x)  are %e and %e." % (xeps,feps))
    x = x0 # Assign initial  guess
    if verbose:
        print('Newton iterations for function: %s' % pfun)
        print('   k             x                    xm                    fm                dfdx')
    k = 0;  maxit = 55 # Current and max number of iterations
    converged = 0; root=0
    while k < maxit:
        k = k + 1
        xm = x-fun(x)/dfundx(x) # Computing the next point as xm=x-f(x)/f'(x)
        dx = xm-x # Difference betwwen x_{n+1}-x_n
        fm = fun(xm)      # Calculate function value at point xm
        dfm = dfundx(xm)  # Calculate function derivative value at point xm
        if verbose:
            print('%4d  %21.16f %21.16f %18.10e %18.10e' % (k,x,xm,fm,dfm))
        if (abs(fm) < feps) or (abs(dx) < xeps): # True when root is found. "and" is also acceptable
            root = xm   # Assign root variable
            converged=1 # Assign "converged" as true
            break       # Break the "while" loop
        x=xm # Set x_n=x_{n+1}
    if converged:
        print("Mynewton: Root is found as %22.16f within tolerance after %d iterations" % (root,k))
    else:
        print("Mynewton: Root is not found within tolerance after %d iterations" % k)
# Start Modify Parameters
f = lambda x: 3*x+sin(x)-exp(x); df = lambda x: 3+cos(x)-exp(x); x0=0.0; fstart=0.0; fend=2.0; pfun = '3x+sin(x)-e^x'
# f = lambda x: x-x**(1/3)-2; df= lambda x: 1-1/(3*x**(2/3)); x0=3.0; fstart=0; fend=4.0; pfun = 'x-x^(1/3)-2'
# # f = lambda x: x-np.copysign(np.abs(x)**(1./3), x)-2; df= lambda x: 1-1/(3*x**(2/3)); x0=-2.0; fstart=-4.0; fend=4.0; pfun = 'x-x^(1/3)-2'
# f = lambda x: sin(x); df= lambda x: cos(x); x0=3.0; fstart=-2*np.pi; fend=2.0*np.pi; pfun = 'sin(x)'
# End Modify Parameters
xval = []
funval = []
kfirst = fstart; klast = fend; kincrement = 0.01
for j in np.arange(kfirst, klast + kincrement, kincrement):
    xval.append(j)
    funval.append(f(j))
plt.title(pfun)
plt.xlabel('X Value')
plt.ylabel('Function Value')
plt.plot(xval, funval, '-')
plt.grid()
# Begin search for brackets
nb = 0; lb=0
for k in range(1,len(xval)):
  if np.sign(funval[lb]) != np.sign(funval[lb+1]):  # True if sign of f(x) changes in the interval
      nb = nb + 1
      circle = plt.Circle((xval[k], 0), 0.04, color='r')
      ax = plt.gca()
      ax.add_patch(circle)
  lb = lb + 1
if nb==0:
    print('No brackets found. Check [kfirst,klast] or decrease kincrement!')
# End search for brackets
plt.savefig('function_plot.eps')
plt.show()
eps = sys.float_info.epsilon
mynewton(f,df, x0,5*eps,5*eps,1)
#
from scipy import optimize
root = optimize.newton(f, x0, df)
print("SciPy Newton: Root is      %22.16f. Accept that value as exact! " % root)
#
# Tolerance values in x and f(x)  are 1.110223e-15 and 1.110223e-15.
# Newton iterations for function: 3x+sin(x)-e^x
#    k             x                    xm                    fm                dfdx
#    1     0.0000000000000000    0.3333333333333333  -6.8417728290e-02   2.5493445212e+00
#    2     0.3333333333333333    0.3601707135776337  -6.2798507057e-04   2.5022625478e+00
#    3     0.3601707135776337    0.3604216804760197  -5.6251553193e-08   2.5018142443e+00
#    4     0.3604216804760197    0.3604217029603242  -4.4408920985e-16   2.5018142041e+00
# Mynewton: Root is found as     0.3604217029603242 within tolerance after 4 iterations
# SciPy Newton: Root is          0.3604217029603244. Accept that value as exact!