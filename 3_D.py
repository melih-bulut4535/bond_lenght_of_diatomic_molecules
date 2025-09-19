# This program defines a damped_harmonic_motion function that returns the right side of the system as a NumPy array.
# It also defines a runge_kutta4 function that implements the quaternary Runge-Kutta method and
# takes the damper_harmonic_motion function as an argument. Finally, it sets the b, w, h, t0, tf,
# and y0 parameters and uses the runge_kutta4 function to solve the system and plot the solution.
import numpy as np

def damped_harmonic_motion(t, y, b, w):
    return np.array([y[1], -b*y[1] - w**2*y[0]])

def runge_kutta4(f, y0, t0, tf, h, *args):
    t = np.arange(t0, tf+h, h)
    n = len(t)
    y = np.zeros((n, len(y0)))
    y[0] = y0
    for i in range(n-1):
        k1 = f(t[i], y[i], *args)
        k2 = f(t[i] + h/2, y[i] + h/2*k1, *args)
        k3 = f(t[i] + h/2, y[i] + h/2*k2, *args)
        k4 = f(t[i] + h, y[i] + h*k3, *args)
        y[i+1] = y[i] + h/6*(k1 + 2*k2 + 2*k3 + k4)
    return t, y

# Set parameters
b = 0.5
w = 1
h = 0.01
t0 = 0
tf = 10
y0 = np.array([1, 1])

# Solve the system using Runge-Kutta method
t, y = runge_kutta4(damped_harmonic_motion, y0, t0, tf, h, b, w)

# Print the first few rows of the data table for part 1
print('t       y(t)    v(t)')
for i in range(10):
    print(f'{t[i]:.2f}  {y[i,0]:.4f}  {y[i,1]:.4f}')

# Plot the solution for part 1
import matplotlib.pyplot as plt
plt.title("PART 1")
plt.plot(t, y[:,0], label='y(t)')
plt.plot(t, y[:,1], label='v(t)')
plt.xlabel('Time')
plt.legend()
plt.show()

# Solve the system for part 2
b = 0
t, y = runge_kutta4(damped_harmonic_motion, y0, t0, tf, h, b, w)

# Plot the solution for part 2
plt.title("PART 2")
plt.plot(t, y[:,0], label='y(t)')
plt.plot(t, y[:,1], label='v(t)')
plt.xlabel('Time')
plt.legend()
plt.show()
