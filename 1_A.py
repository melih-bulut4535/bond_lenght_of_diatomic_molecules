# This code defines the functions f1(x) and f2(x) to represent the two given functions.
# It then defines the interval_halving, secant, and
# newton_raphson functions to implement the corresponding root-finding algorithms.
import math

def f1(x):
    return x*math.tan(x) - 1

def f2(x):
    return math.exp(x) - x**2 + 3*x - 2

def interval_halving(f, a, b, tol=1e-6, max_iter=100):
    fa, fb = f(a), f(b)
    if fa*fb > 0:
        raise ValueError("Function has same sign at endpoints.")
    for i in range(max_iter):
        c = (a+b)/2
        fc = f(c)
        if fc == 0 or (b-a)/2 < tol:
            return c
        if fa*fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    raise ValueError("Method did not converge.")

def secant(f, a, b, tol=1e-6, max_iter=100):
    fa, fb = f(a), f(b)
    for i in range(max_iter):
        c = b - fb*(b-a)/(fb-fa)
        fc = f(c)
        if fc == 0 or abs(c-b) < tol:
            return c
        a, fa, b, fb = b, fb, c, fc
    raise ValueError("Method did not converge.")

def newton_raphson(f, dfdx, x0, tol=1e-6, max_iter=100):
    for i in range(max_iter):
        fx = f(x0)
        dfdx_val = dfdx(x0)
        if abs(fx) < tol:
            return x0
        if dfdx_val == 0:
            raise ValueError("Derivative is zero.")
        x0 = x0 - fx/dfdx_val
    raise ValueError("Method did not converge.")

# Find and write the roots of f1(x) using interval halving, secant, and Newton-Raphson methods
print("Root of f1(x) using interval halving:", interval_halving(f1, 0, 1))
print("Root of f1(x) using secant method:", secant(f1, 0, 1))
print("Root of f1(x) using Newton-Raphson method:", newton_raphson(f1, lambda x: math.tan(x) + x*(1/math.cos(x))**2, 0.5))

# Find and write the roots of f2(x) using interval halving, secant, and Newton-Raphson methods
print("Root of f2(x) using interval halving:", interval_halving(f2, 0, 1))
print("Root of f2(x) using secant method:", secant(f2, 0, 1))
print("Root of f2(x) using Newton-Raphson method:", newton_raphson(f2, lambda x: math.exp(x) - 2*x + 3, 0.5))
