import numpy as np
import matplotlib.pyplot as plt

# define the function to differentiate
def f(x):
    return np.log(x)

# define the exact second derivative
def f_2(x):
    return -1 / x**2

# define the point to differentiate at
x = 1

# define the step lengths to use
hs = [0.1, 0.01, 0.001]

# loop over the step lengths
for h in hs:
    # compute the central difference approximation
    f_2_c = (f(x+h) - 2*f(x) + f(x-h)) / h**2

    # compute the forward difference approximation
    f_2_f = (f(x+2*h) - 2*f(x+h) + f(x)) / h**2

    # compute the backward difference approximation
    f_2_b = (f(x) - 2*f(x-h) + f(x-2*h)) / h**2

    # compute the exact value of the second derivative
    f_2_exact = f_2(x)

    # print the results
    print("h = {:.3f}".format(h))
    print("Central difference approximation: {:.6f}, error: {:.6f}".format(f_2_c, abs(f_2_c - f_2_exact)))
    print("Forward difference approximation: {:.6f}, error: {:.6f}".format(f_2_f, abs(f_2_f - f_2_exact)))
    print("Backward difference approximation: {:.6f}, error: {:.6f}".format(f_2_b, abs(f_2_b - f_2_exact)))
    print()
