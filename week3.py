# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math
from functions1 import *

# Here we call the metropolis algorithm one time
c = 0.35
iterations = 30000
dim = 6 # total dimensions = (spatial dimensions) x (number of particles)
sp_dim = dim//2
gamma = 1
walkers = 20

positions_tracker, energy_tracker, variation_expectation_1, variation_expectation_2, accept_ratio = metropolis_helium(c, iterations, walkers, dim)

# We print the observed values for this specific c
print(np.mean(energy_tracker))
print(accept_ratio)

# +
# We plot the positions of the walkers to see if they make sense

xdata = np.linspace(1, iterations, iterations)

plt.plot(xdata, positions_tracker[:,0,0], 'b-', label='x1 for walker 1')
plt.title('Displacement of x1 of walker1 for d=0.2')
plt.xlabel('Iterations')
plt.ylabel('Position')
plt.legend()
plt.show()

plt.plot(xdata, positions_tracker[:,1,3], 'r-', label='x2 for walker 2')
plt.title('Displacement of x2 of walker2 for d=0.2')
plt.xlabel('Iterations')
plt.ylabel('Position')
plt.legend()
plt.show()

plt.plot(xdata, positions_tracker[:,2,1], 'b-', label='y1 for walker 3')
plt.title('Displacement of y1 of walker3 for d=0.2')
plt.xlabel('Iterations')
plt.ylabel('Position')
plt.legend()
plt.show()

plt.plot(xdata, positions_tracker[:,3,5], 'r-', label='z2 for walker 4')
plt.title('Displacement of z2 of walker4 for d=0.2')
plt.xlabel('Iterations')
plt.ylabel('Position')
plt.legend()
plt.show()
# -

# Now we compute the values via the fucntion (it accepts any value of c and it returns the metropolis integration for it)
c1 = 0.5 # now we choose c = 0.5
E, variance1, variance2, accept_ratio, var1, var2 = mean_energy_helium1(c1, iterations, walkers, dim)

# We print the important values
print(E)
print(variance1, variance2)

# This time the function, given a c initial parameter calls the metropolis functiona nd finds another value of c closer to the optimal.
# This is repeated a certain amount of times
c = 1 # we choose now this c
c_tracker, c, mean_energies, energy_tracker, variance1_tracker, variance2_tracker, accept_ratio, n = helium_opt1(c, walkers, dim)

print(mean_energies, variance1_tracker[100])
print(n) # n is the number of times we call metropolis for a different value of c

# +
# Now we can plot the evolution of the parameter c during every loop of the function "helium_opt1":
# The minimum of the helium energy can be found searching for the minimum over the mean values

loops = np.arange(0, n+1)

plt.plot(loops, c_tracker[:n+1], 'b-', label='value of the parameter')
plt.title('Value of the parameter c obtained in every loop of the optimization')
plt.xlabel('loop nr')
plt.ylabel('parameter c')
plt.legend()
plt.show()

plt.plot(loops, energy_tracker[:n+1], 'r-', label='value of the mean energy')
plt.title('Value of the mean energy E obtained in every loop of the optimization')
plt.xlabel('loop nr')
plt.ylabel('Mean energy')
plt.legend()
plt.show()

plt.plot(c_tracker[:n+1], variance1_tracker[:n+1], 'g-', label='Variance1')
plt.title('Variance over all configurations')
plt.xlabel('parameter c')
plt.ylabel('Variance')
plt.legend()
plt.show()

plt.plot(c_tracker[1:n+1], variance2_tracker[1:n+1], 'y-', label='Variance2')
plt.title('Variance over the mean energies of the walkers')
plt.xlabel('parameter c')
plt.ylabel('Variance')
plt.legend()
plt.show()

plt.plot(c_tracker[:n+1], variance1_tracker[:n+1], 'g-', label='Variance1')
plt.plot(c_tracker[1:n+1], variance2_tracker[1:n+1], 'y-', label='Variance2')
plt.title('Variance over all configurations and variance over the mean energies of each walker')
plt.xlabel('parameter c')
plt.ylabel('Variance')
plt.legend()
plt.show()

plt.plot(c_tracker[:n+1], energy_tracker[:n+1], 'b-', label='Mean energy')
plt.title('Mean energy over the parameter c')
plt.xlabel('parameter c')
plt.ylabel('Mean energy')
plt.legend()
plt.show()

# +
"""Now we repeat the same steps but by using the autocorrelation2 function to compute the means of all observables and their variances. 
There are two major points that differentiate to our previous approaches both for Helium and Hydrogen: 
A) Now we do not treat the different configurations as statistically independent but rather as sets of correlated data (autocorrelation)
B) We use autocorrelation2, in which, unlike in autocorrelation1, we do not compute the total variance as the mean of each walker's variance 
but as the square root of the mean of their squares --> More correct statistically"""

# We begin by testing how the autocorrelation1 function works on Helium before we proceed to the autocorrelation2 (mean_energy_helium2)
c = 0.05
E, variance1, variance2, accept_ratio, var1, var2 = mean_energy_helium3(c, iterations, walkers, dim)
print(E)
print(variance2)

# +
# Works well so we move to the autocorrelation2 (mean_energy_helium3)

c = 0.05
E, variance1, variance2, accept_ratio, var1, var2 = mean_energy_helium3(c, iterations, walkers, dim)
print(E)
print(variance2)
# -

# no matter how we define the variance it is always too small when using the autocorrelation function. We yet to see if it follows the patterns found in literature:
print(variance1, variance2)
print(accept_ratio)

# We finnaly do the minimization of c taking into account correlated sets of configurations (helium_opt2)
c = 1
c_tracker, c, mean_energies, energy_tracker, variance1_tracker, variance2_tracker, accept_ratio, n = helium_opt2(c, walkers, dim)

# +
# We plot the values obtainedfor every c going towards the minimum

loops = np.arange(0, n+1)

plt.plot(loops, c_tracker[:n+1], 'b-', label='value of the parameter')
plt.title('Value of the parameter c obtained in every loop of the optimization')
plt.xlabel('loop nr')
plt.ylabel('parameter c')
plt.legend()
plt.show()

plt.plot(loops, energy_tracker[:n+1], 'r-', label='value of the mean energy')
plt.title('Value of the mean energy E obtained in every loop of the optimization')
plt.xlabel('loop nr')
plt.ylabel('Mean energy')
plt.legend()
plt.show()

plt.plot(loops, variance1_tracker[:n+1], 'g-', label='Variance1')
plt.title('Variance over all configurations')
plt.xlabel('loop number')
plt.ylabel('Variance')
plt.legend()
plt.show()

plt.plot(loops[1:], variance2_tracker[1:n+1], 'y-', label='Variance2')
plt.title('Variance over the mean energies of the walkers')
plt.xlabel('loop number')
plt.ylabel('Variance')
plt.legend()
plt.show()

plt.plot(c_tracker[:n+1], variance1_tracker[:n+1], 'g-', label='Variance1')
plt.plot(c_tracker[1:n+1], variance2_tracker[1:n+1], 'y-', label='Variance2')
plt.title('Variance over all configurations and variance over the mean energies of each walker')
plt.xlabel('parameter c')
plt.ylabel('Variance')
plt.legend()
plt.show()

plt.plot(c_tracker[:n+1], energy_tracker[:n+1], 'b-', label='Mean energy')
plt.title('Mean energy over the parameter c')
plt.xlabel('parameter c')
plt.ylabel('Mean energy')
plt.legend()
plt.show()

# +
"""Once we have obtained the optimal value of the parameter c, we can plug it into the metropolis algorithm
and observe the evolution of the energy around the expectation value of a walkers-
This way we can go through one final metropolis with the minimized value of c and obtain the minimum value of the energy"""

walkers = 5
positions_tracker, energy_tracker, variation_expectation_1, variation_expectation_2, accept_ratio = metropolis_helium(c, iterations, walkers, dim)
E, variance1, variance2, accept_ratio, var1, var2 = mean_energy_helium3(c, iterations, walkers, dim)
# -

print(E)
print(variance1, variance2)
print(accept_ratio)


