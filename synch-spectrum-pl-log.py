
import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import matplotlib.pyplot as plt
from math import e

# Constants
p = 2.5
# Capital constants can be extrapolated from Rybicki & Lightman (1979).
# In this code all constants are set to 1 and the spectrum is normalized to
# peak at 1.
A = 1
C = 1
D = 1
INFINITY = 100

# Synchrotron function.
def F(x):
    # special.kv() is the modified Bessel function of the second kind.
    # integrate.quad performs general integration.
    return x * integrate.quad(lambda y: special.kv(5/3, y), x, INFINITY)[0]

# Total emitted synchrotron power per frequency for a single electron
def P(w, gamma):
    return A * F(w/(D*pow(gamma, 2)))

# Number of electrons per unit energy.
# The electron distribution used here is a simple power law with a sharp cutoff
# at gamma_min.
def N(gamma):
    gamma_min = 1
    if (gamma < gamma_min):
        return 0
    else:
        return pow(gamma, -p)

# Total power emitted per unit frequency for a distribution of electrons.
def Ptot(w):
    gamma_1 = 1
    gamma_2 = 2
    return integrate.quad(lambda gamma: P(w, gamma)*N(gamma), gamma_1, gamma_2)[0]

"""
Limits chosen so that the behaviour of the function is clear while detail
can still be resolved.
"""
lower_limit = 0
upper_limit = 100
num_of_points = 500

# Define an empty arrays to store F(x) values.
# Numbers correspond to their figure indexes.
# 0 -> Electron power dist.
# 1 -> Ptot(w) for small w
# 2 -> Ptot(w) for larger w
y0 = [0]*(num_of_points)
y1 = [0]*(num_of_points)

max_P = 0
max_w = 0
# Calculate and store y for each point.
for i in range(0, num_of_points):
    # Calculates w between upper and lower limits.
    w = i * (upper_limit - lower_limit)/num_of_points
    y0[i] = N(w)
    y1[i] = Ptot(w)
    # Calculate spectral peak for normalization
    if y1[i] > max_P:
        max_P = y1[i]
        max_w = w

w = np.linspace(lower_limit, upper_limit, num=num_of_points)

# Normalize the spectrum.
max_P = max(y1)
for i in range(0, num_of_points):
    y1[i] = y1[i]/max_P

fig, axs = plt.subplots(2, constrained_layout=True)

# Electron power law distribution
axs[0].set_title('Electron power distribution')
axs[0].plot(w, y0)
axs[0].set_xlabel(r'$\gamma$')
axs[0].set_ylabel(r'N($\gamma$)')
axs[0].set_xscale('log')
axs[0].set_yscale('log')

# Synchrotron spectrum compared with fits.
# The w^1/3 lower tail and exponential upper tail are expected from summing 
# F(x) over a power law distribution, as F(x) is approximately proportional to 
# these for limitting values of x as described in R&L (1979). 
axs[1].set_title('Synchrotron spectrum')
axs[1].plot(w, y1, label=r'$P_{tot}(\omega)$')
# TODO: Fits need shifting.
axs[1].plot(w, pow(w, 1/3), '--', label=r'$\omega^\frac{1}{3}$ fit')
axs[1].plot(w, pow(w, -(p-1)/2), '--', label=r'$\omega^{\frac{-(p-1)}{2}}$ fit')
axs[1].plot(w, pow(e, -w, '--', label=r'e$^{-\omega}$ fit')
axs[1].set_ylabel(r'N($\gamma$)')
axs[1].set_ylabel(r'$P_{tot}(\omega)$')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].set_ylim(0.1, 1.1)
axs[1].legend()

plt.show()










