import sys
from math import *


def mynewton(fun, dfundr, r0, rtol, ftol, verbose):
    """
    Use Newton's method to find a root of the scalar equation f(x) = 0
    """
    reps = max(rtol, 5*eps)  # Smallest tolerances are 5*eps
    feps = max(ftol, 5*eps)  # Smallest tolerances are 5*eps
    print("Tolerance values in x and f(x)  are %e and %e." % (reps, feps))
    r = r0  # Assign initial guess
    if verbose:
        print(f"Newton iterations for function: {fun.__name__}")
        print('   k             r                    rm                    fm                dfdr')
    k = 0  # Current number of iterations
    maxit = 10  # Maximum number of iterations
    converged = False
    root = 0
    while k < maxit:
        k += 1
        rm = r - fun(r) / dfundr(r)  # Computing the next point as rm=x-f(x)/f'(x)
        dr = rm - r  # Difference between x_{n+1}-x_n
        fm = fun(rm)  # Calculate function value at point rm
        dfm = dfundr(rm)  # Calculate function derivative value at point rm
        if verbose:
            print(f'{k:4d}  {r:21.16f} {rm:21.16f} {fm:18.10e} {dfm:18.10e}')
        if abs(fm) < feps or abs(dr) < reps:  # True when root is found. "and" is also acceptable
            root = rm  # Assign root variable
            converged = 1  # Assign "converged" as true
            break  # Break the "while" loop
        r = rm  # Set x_n=x_{n+1}
    if converged:
        print(f"Mynewton: Root is found as {root:22.16f} within tolerance after {k} iterations")
    else:
        print(f"Mynewton: Root is not found within tolerance after {k} iterations")


# Start Modify Parameters
f = lambda r: -sqrt(14.4)**2/r + 1.09*1000*(sqrt(14.4)**(-r/0.33))
df = lambda r: sqrt(14.4)**2/r**2 - (1.09*1000/0.33) * sqrt(14.4)**(-r/0.33)
r0 = 1.0
fstart = 0.0
fend = 2.0
pfun = '-sqrt(14.4)**2/r + 1.09*1000*(sqrt(14.4)**(-r/0.33))'
eps = sys.float_info.epsilon
mynewton(f, df, r0, 5*eps, 5*eps, 1)
