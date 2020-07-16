
import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import matplotlib.pyplot as plt
from math import e

# Constants
p = 2.5
A = 1
C = 1
D = 1
INFINITY = 1000
gamma_min = 1

# Synchrotron function.
def F(x):
    # special.kv() is the modified Bessel function of the second kind.
    # integrate.quad performs general integration.
    return x * integrate.quad(lambda y: special.kv(5/3, y), x, INFINITY)[0]

# Total emitted synchrotron power per frequency for a single electron
def P(w, gamma):
    return A * F(w/D*pow(gamma, 2))

def N(gamma):
    if (gamma < gamma_min):
        return 0
    else:
        return pow(gamma, -p)

gamma_1 = 1
gamma_2 = 2
def Ptot(w):
    return integrate.quad(lambda gamma: P(w, gamma)*N(gamma), gamma_1, gamma_2)[0]

w = np.arange(0, 5, 0.05)
'''
Want to plot in this way but gives error
#y = Ptot(w)
#y = F(w)
'''

"""
Limits chosen so that the behaviour of the function is clear while detail
can still be resolved.
"""
lower_limit = 0
upper_limit = 3
# Gives sufficient plotting precision.
num_of_points = 500

# Define an empty arrays to store F(x) values.
y0 = [0]*(num_of_points)
y1 = [0]*(num_of_points)

max_P = 0
max_w = 0
# Calculate and store y for each point.
for i in range(0, num_of_points):
    # Calculates x between upper and lower limits.
    w = i * (upper_limit - lower_limit)/num_of_points
    y0[i] = N(w)
    y1[i] = Ptot(w)
    if y1[i] > max_P:
        max_P = y1[i]
        max_w = w

print(f"Max at {max_w, max_P}")

w = np.linspace(lower_limit, upper_limit, num=num_of_points)

max_P = max(y1)

for i in range(0, num_of_points):
    y1[i] = y1[i]/max_P


fig, axs = plt.subplots(3)

# Electron power law distribution
axs[0].plot(w, y0)
axs[0].set_xbound(0, 8)
axs[0].set_ybound(0,1.1)
axs[0].set_ylabel(r'N($\gamma$)')

# Ptot(w) compared with w^1/3 for small w
axs[1].plot(w, y1)
axs[1].plot(w, pow(w-0.01, 1/3)+0.5)
axs[1].set_xbound(0, 0.3)
axs[1].set_ybound(0,1.1)
axs[1].set_ylabel(r'Ptot($\omega$)')

# Ptot(w) compared to e^x for larger w
axs[2].plot(w, y1)
axs[2].plot(w, pow(e, -(w - 0.15)))
axs[2].set_xbound(0, 1.2)
axs[2].set_ybound(0,1.1)
axs[2].set_ylabel(r'Ptot($\omega$)')


plt.show()










