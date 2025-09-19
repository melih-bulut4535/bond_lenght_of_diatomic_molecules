# This code defines two functions f1 and f2 which represent the right-hand side of the differential equations y" = yy'
# and y" - 2y + y = e^2, respectively.
# It then defines the rk4 function which solves a system of first-order differential equations
# using the fourth-order Runge-Kutta method.
import numpy as np

def f1(x, y):
    return np.array([y[1], y[0] * y[1]])

def f2(x, y):
    return np.array([y[1], 2 * y[0] - y[1] + np.exp(2*x)])

def rk4(f, x0, y0, h, n):
    x = np.linspace(x0, x0 + h * n, n + 1)
    y = np.zeros((n + 1, len(y0)))
    y[0] = y0
    for i in range(n):
        k1 = h * f(x[i], y[i])
        k2 = h * f(x[i] + h/2, y[i] + k1/2)
        k3 = h * f(x[i] + h/2, y[i] + k2/2)
        k4 = h * f(x[i] + h, y[i] + k3)
        y[i+1] = y[i] + 1/6 * (k1 + 2*k2 + 2*k3 + k4)
    return x, y

# Solving y" = yy', y(0) = 1, y'(0) = -1
x, y = rk4(f1, 0, [1, -1], 0.1, 10)
print("Solution to y'' = yy', y(0) = 1, y'(0) = -1:")
print("x        y")
for i in range(len(x)):
    print("{:.1f}   {:.5f}".format(x[i], y[i, 0]))

# Solving y" - 2y + y = e^2, y(0) = 0, y'(0) = 0
x, y = rk4(f2, 0, [0, 0], 0.1, 10)
print("\nSolution to y'' - 2y + y = e^2, y(0) = y'(0) = 0:")
print("x        y")
for i in range(len(x)):
    print("{:.1f}   {:.5f}".format(x[i], y[i, 0]))
