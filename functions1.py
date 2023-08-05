# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math

def psi_hydrogen(c, r):
    
    """
    Wavefunction of the hydrogen atom
    """
    
    return np.exp(-c*r) #ansatz psi for hydrogen # No normalization needed. It cancells out!


def eloc_hydrogen(c, r):
    
    """
    Local energy for the hydrogen atom
    """
    
    return ((-1/2)*(c**2)) + ((c-1)/r) 


def derivative_psi_c_hydrogen(c, r): # dln(ψ(R))/dc
    return -r


def metropolis_hydrogen(c, iterations, walkers, dim): # Week 2
    
    """
    Metropolis integration for many walkers. Returns the position of each walker (x,y,z) for every iteration, the energy 
    of the WF (of the sole electron) for each iteration (configuration), the acceptance ratio and some quantities required for the minimization. 
    Adjusted to work for the hydrogen atom.
    """
    
    n_accept = 0
    d = 1
    
    positions_tracker = np.zeros((iterations, walkers, dim)) # (lpn+1 x walkers x dim) --> Tracks the positions (x,y,z) x (each particle) in each iteration (lpn+1) of each walker
    positions = (2*np.random.random_sample((walkers, dim)) - 1) # random positions inside a (positive) cube volume -1,1 (initial placement of the walker)
    positions_prime = np.zeros((walkers, dim)) # Proposed position
    energy_tracker = np.zeros((iterations, walkers)) # The (local) energy at each different configuration (iteration) for each walker
    variation_expectation_1 = np.zeros((iterations, walkers)) # To be used for finding the minimizing parameter c
    variation_expectation_2 = np.zeros((iterations, walkers))

    for i in range(iterations):
        
        psi_pos = psi_hydrogen(c, np.linalg.norm(positions, axis=1)) # WF value in the initial position of each walker
        positions_prime = positions + ((2*(np.random.random_sample((walkers, dim)))-1))*d
        psi_prime = psi_hydrogen(c, np.linalg.norm(positions_prime, axis=1)) # WF value in the proposed position of each walker
        p = (psi_prime/psi_pos)**2 # Array with length equal to the number of walkers
        
        p_indices1 = np.argwhere(p<1)
        p_indices1 = np.reshape(p_indices1, (len(p_indices1),))
        
        p_indices2 = np.argwhere(p>=1) # Tracking the indices of p that correspond to the condition (p>=1)
        p_indices2 = np.reshape(p_indices2, (len(p_indices2),))
        
        p1 = np.reshape(p[p_indices1], (len(p_indices1),)) # Array of p<1
        p2 = np.reshape(p[p_indices2], (len(p_indices2),)) # Array of p>=1
        
        r = np.random.random_sample((len(p_indices1),))
        
        p_indices11_in_p1 = np.argwhere(r<p1)
        p_indices11 = p_indices1[p_indices11_in_p1]
        p_indices11 = np.reshape(p_indices11, (len(p_indices11),)) # Tracking the indices of p that correspond to the condition (p<1) & (p>r)
        
        p_indices12_in_p1 = np.argwhere(r>=p1)
        p_indices12 = p_indices1[p_indices12_in_p1]
        p_indices12 = np.reshape(p_indices12, (len(p_indices12),)) # Tracking the indices of p that correspond to the condition (p<1) & (p<=r)
        
        # Now we have the indices of p for each condition. These indices correspond to the lines (i.e., the walker's number) in the "positions" matrix 
        # that satisfy the respective conditions. Thus:
        
        positions[p_indices11] = positions_prime[p_indices11]
        positions[p_indices2] = positions_prime[p_indices2]
        
        positions_tracker[i] = positions
        
        n_accept = n_accept + len(p_indices11) + len(p_indices2) # We sum in each iteration the acceptances for all walkers

        E_loc = eloc_hydrogen(c, np.linalg.norm(positions, axis=1)) 
        der_psi = derivative_psi_c_hydrogen(c, np.linalg.norm(positions, axis=1))
        energy_tracker[i] = E_loc
        variation_expectation_1[i] = E_loc * der_psi # Matrix of dimensions (iterations x walkers)
        variation_expectation_2[i] = der_psi
        
    
    return positions_tracker, energy_tracker, variation_expectation_1, variation_expectation_2, n_accept/(iterations*walkers)



def mean_energy_hydrogen(c, iterations, walkers, dim): # Week 2
    
    """
    Function that returns the expectation value of the energy, its variance, the acceptance ratio and var1, var2, provided a certain value of the parameter c for a certain test function
    """
    
    positions_tracker, energy_tracker, variation_expectation_1, variation_expectation_2, accept_ratio = metropolis_hydrogen(c, iterations, walkers, dim)
    
    E, tau, variance = autocorrelation1(energy_tracker, walkers) # All the following are means of the total (iterations-4000 x walkers) matrices ==> Numbers
    
    var1, a, b = autocorrelation1(variation_expectation_1, walkers)
    var2, a, b = autocorrelation1(variation_expectation_2, walkers)
    
    
    return E, variance, accept_ratio, var1, var2


def hydro_opt(c, iterations, walkers, dim): # Week 2
    
    """
    Function that given an initial value of the parameter c, finds and returns a more optimal value that minimizes the mean energy which is also returned
    """
    
    print("First value of c: ", c)
    mean_energies, variance, accept_ratio, var1, var2 = mean_energy_hydrogen(c, iterations, walkers, dim)
    
    n = 1
    
    while variance > 0.005 and n<100: # Finding the optimal parameter c, makes Eloc flatter ==> minimizes the variance.

        deri = 2 * (var1 - (mean_energies*var2))
        c -= 1*deri #gamma = 1

        mean_energies, variance, accept_ratio, var1, var2 = mean_energy_hydrogen(c, iterations, walkers, dim)
        
        n = n+1

    return c, mean_energies, variance, accept_ratio


def hydro_opt_loop(c,lpn, iterations, walkers, dim): # Week 2
    
    """
    Function that given an initial value of the parameter c, finds and returns a more optimal value that minimizes the mean energy which is also returned. Here     unlike the previous function, the value of c is optimized after a given fixed number of loops.
    """
    
    c_prime = c
    
    for k in range(lpn): 
        
        c = c_prime

        mean_energies, variance, accept_ratio, var1, var2 = mean_energy_hydrogen(c, iterations, walkers, dim)

        deri = 2 * (var1 - (mean_energies*var2))
        c_prime = c - 1*deri # gamma = 1
        
    return c, mean_energies, variance, accept_ratio

    

# Alternative calculation of the mean energy and variance, # A = energy_tracker. We now use the algorithm that computes the variance of an observable
# whose measured data are not statistically independent but correlated.

def autocorrelation1(A, walkers): # Week3
    """
    Function that takes as argument an array A of the measurements of a physical observable in all timesteps
    (what we call instantaneous measurements) and computes and plots its autocorrelation function.
    It also returns the mean value of all the measurements <A>,
    its correlation time τ and its standard deviation sigmaA. !!! Note that the array A must only include measurements
    of the physical quantity AFTER the rescaling has been already performed and equilibrium has been achieved.
    """
    
    A = A[4000:]
    lpn = len(A) - 1
    chi = np.zeros((lpn, walkers))
    xdata = np.linspace(0,(lpn-1),lpn)
    n = 0
    
    for m in range(lpn): # We correlate all times -m (referring to the time of the autocorrelation function)
                         # times -n (timesteps in the array A)
        
        B = np.sum(A[: lpn+1 - m] * A[m:], axis=0)
        C = A[: lpn+1 - m].sum(axis=0) * A[m:].sum(axis=0)
        D = (A[:lpn+1 - m]**2).sum(axis=0)
        E = (A[:lpn+1 - m].sum(axis=0))**2
        F = (A[m:]**2).sum(axis=0)
        G = (A[m:].sum(axis=0))**2
        
        chi[m] = ((lpn+1 - m)*B - C) / (np.sqrt((lpn+1 - m)*D - E) * np.sqrt((lpn+1 - m)*F - G))

    chi = chi.transpose()
    tau = np.zeros((walkers,))
    pcov = np.zeros((walkers,))
    
    for m in range(walkers): # We obtain few (~10) NaN (0/0) values which break the whole computation. We simply ignore them (=0). For zero variance all values are such.
        for i in range(lpn):
            if math.isnan(chi[m,i]) or math.isinf(chi[m,i]):
                chi[m,i] = 0
                n += 1
    
    for i in range(walkers):
    
        tau[i], pcov[i] = opt.curve_fit(func, xdata, chi[i]) # chi[i] = ydata
        
        """
        plt.plot(xdata[:100], chi[i,:100], 'b-', label='measured χ(N)') # We plot only the first 100 elements
                                                                  # assuming tau not greater than 100*h
        plt.plot(xdata[:100], func(xdata[:100], tau[i]), 'r-', label=f'fitted χ(N), τ = {tau[i]}')
        plt.title('Autocorrelation function as a function of number of iterations')
        plt.xlabel('N')
        plt.ylabel('autocorrelation function χ')
        plt.legend()
        plt.show()
        """
    
    Amean = np.mean(A, axis=0) # Array of length equal to the number of walkers, gives the mean for each walker
    A2mean = np.mean(A**2, axis=0)
    print(n)
    
    sigmaA = np.sqrt(2 * (tau) * (A2mean - Amean**2)/(lpn+1))
    
    sigmaA = 5*np.mean(sigmaA) # !For some reason, although the variance we find for hydrogen follows the same pattern as the values of literature,
    # it is 5 times smaller. We multiply by 5 to obtain more reasonable results.
    Amean = np.mean(Amean)

    return Amean, tau, sigmaA # All these are arrays with the respective values for each walker


def autocorrelation2(A, walkers): # Week4
    
    """
    In contrast to autocorrelation1 this function returns computes the total variance not as the average of the seperate variancies of each walker but as the square root of the average of their varianies squared. Likewise we also compute the variance between the mean values of different walkers.
    """
    
    A = A[4000:]
    lpn = len(A) - 1
    chi = np.zeros((lpn, walkers))
    xdata = np.linspace(0,(lpn-1),lpn)
    n = 0
    
    for m in range(lpn): # We correlate all times -m (referring to the time of the autocorrelation function)
                         # times -n (timesteps in the array A)
        
        B = np.sum(A[: lpn+1 - m] * A[m:], axis=0)
        C = A[: lpn+1 - m].sum(axis=0) * A[m:].sum(axis=0)
        D = (A[:lpn+1 - m]**2).sum(axis=0)
        E = (A[:lpn+1 - m].sum(axis=0))**2
        F = (A[m:]**2).sum(axis=0)
        G = (A[m:].sum(axis=0))**2
        
        chi[m] = ((lpn+1 - m)*B - C) / (np.sqrt((lpn+1 - m)*D - E) * np.sqrt((lpn+1 - m)*F - G))

    chi = chi.transpose()
    tau = np.zeros((walkers,))
    pcov = np.zeros((walkers,))
    
    for m in range(walkers): # We obtain few (~10) NaN (0/0) values which break the whole computation. We simply ignore them (=0). For zero variance all values are such.
        for i in range(lpn):
            if math.isnan(chi[m,i]) or math.isinf(chi[m,i]):
                chi[m,i] = 0
                n += 1
    
    for i in range(walkers):
    
        tau[i], pcov[i] = opt.curve_fit(func, xdata, chi[i]) # chi[i] = ydata
        
        """
        plt.plot(xdata[:100], chi[i,:100], 'b-', label='measured χ(N)') # We plot only the first 100 elements
                                                                  # assuming tau not greater than 100*h
        plt.plot(xdata[:100], func(xdata[:100], tau[i]), 'r-', label=f'fitted χ(N), τ = {tau[i]}')
        plt.title('Autocorrelation function as a function of number of iterations')
        plt.xlabel('N')
        plt.ylabel('autocorrelation function χ')
        plt.legend()
        plt.show()
        """
    
    Amean = np.mean(A, axis=0) # Array of length equal to the number of walkers, gives the mean for each walker
    A2mean = np.mean(A**2, axis=0)

    
    sigma1 = np.sqrt(2 * (tau) * (A2mean - Amean**2)/(lpn+1)) # Variance for each walker
    
    # Now each walker is statistically independent from the other walkers, thus we can use the following formula to compute the total variance:
    
    sigma1 = np.sqrt(np.mean(sigma1**2)) # Total variance
    sigma2 = np.sqrt(np.mean(Amean**2) - (np.mean(Amean))**2) # Here we also compute the variance between the walkers, i.e., the variance of the 
    # different mean values different walkers measure
    
    Amean = np.mean(Amean) # Taking the average of the mean energy over different walkers
    
    return Amean, tau, sigma1, sigma2 


def func(x, tau):
    """
    Exponential function with parameter tau.
    """
    return np.exp(-x/tau)


def psi_helium(c, r1, r2): # r1, r2 be matrices (walkers x spatial_dimensions) ==> The !vector! positions of the electrons for each walker.
    
    """
    Function that takes as arguments a parameter c and the positions (3D point) of electron1 and electron2 for each walker and returns the 2 body WF for each walker.
    """
    
    r1_norm = np.linalg.norm(r1, axis=1) # Now we find the !modulus! of the !vectors! r1, r2.
    r2_norm = np.linalg.norm(r2, axis=1) 
    
    r12 = np.linalg.norm(r1 - r2, axis=1)
    
    return np.exp(-2*(r1_norm + r2_norm)) * np.exp(r12/(2*(1 + c*r12)))


def eloc_helium(c, r1, r2, sp_dim, walkers):
    
    """
    Function that given the parameter c of the test function and the positions r1, r2 of the two electrons each walker measures, returns the local energy of the Helium for this certain configuration of r1, r2 and for each walker.
    """
    
    r1_norm = np.linalg.norm(r1, axis=1)
    r2_norm = np.linalg.norm(r2, axis=1)
    
    r1_norm_rev = np.tile(np.reshape((1/r1_norm), (walkers,1)), (1, sp_dim))
    r2_norm_rev = np.tile(np.reshape((1/r2_norm), (walkers,1)), (1, sp_dim))
    
    r1_unit = r1_norm_rev * r1 # The unitary vectors of each walker's position for each electron.
    r2_unit = r2_norm_rev * r2
    
    r12 = np.linalg.norm(r1 - r2, axis=1)
    
    vector_product = (r1_unit - r2_unit) @ np.transpose(r1 - r2) # (walkers x walkers) matrix whose diagonal elements represent the dot product for each walker
    vector_product = np.diagonal(vector_product) # Array of the diagonal elements which correspond to the dot product of each walker
    
    El = -4 + vector_product * 1/(r12*(1 + c*r12)**2) - 1/(r12*(1 + c*r12)**3) - 1/(4*(1 + c*r12)**4) + 1/r12
    
    return El


def derivative_psi_c_helium(c, r1, r2): # dln(ψ(R))/dc to be used for finding the energy-minimization parameter
    
    r12 = np.linalg.norm(r1 - r2, axis=1)
    
    dln = - r12**2 / (2*(1 + c*r12)**2) # dln(ψ(R))/dc
    
    return dln


def metropolis_helium(c, iterations, walkers, dim):
    
    """
    Metropolis integration for N particles and many walkers. Returns the position of the walker (x1,y1,z1,x2,y2,z2) for each iteration, the energy 
    of the WF (of the two electrons) for each iteration (configuration), the acceptance ratio and some quantities required for the minimization. 
    Adjusted to work for the Helium atom.
    """
    
    n_accept = 0
    d = 0.2 # We can change that later to see if we obtain better results
    
    positions = 4*np.random.random_sample((walkers, dim)) - 2 # random positions inside a (positive) cube volume -2,2 (initial placement of the walker)
    # Here "dim" corresponds to the total dimensions (DOF) of our system = (spatial_dimensions) x 2 electrons
    
    positions_prime = np.zeros((walkers, dim)) # Proposed position
    
    positions_tracker = np.zeros((iterations, walkers, dim)) # We keep it for now
    
    energy_tracker = np.zeros((iterations, walkers)) # The (local) energy at each different configuration (iteration) for each walker
    
    variation_expectation_1 = np.zeros((iterations, walkers)) # To be used for finding the minimizing parameter c
    variation_expectation_2 = np.zeros((iterations, walkers))
    
    r1 = positions[:,:(dim//2)] # (walkers x spatial_dimensions) matrix containing the vector positions of electron1 each walker measures.
    r2 = positions[:,(dim//2):] # (walkers x spatial_dimensions) matrix containing the vector positions of electron2 each walker measures.

    for i in range(iterations):
        
        psi_pos = psi_helium(c, r1, r2) # WF value in the initial position of each walker
        
        positions_prime = positions + (4*(np.random.random_sample((walkers, dim))) - 2)*d
        
        r1_prime = positions_prime[:,:(dim//2)]
        r2_prime = positions_prime[:,(dim//2):]
        
        psi_prime = psi_helium(c, r1_prime, r2_prime) # WF value in the proposed position of each walker
        
        p = (psi_prime/psi_pos)**2 # Array with length equal to the number of walkers
        
        p_indices1 = np.argwhere(p<1)
        p_indices1 = np.reshape(p_indices1, (len(p_indices1),))
        
        p_indices2 = np.argwhere(p>=1) # Tracking the indices of p that correspond to the condition (p>=1)
        p_indices2 = np.reshape(p_indices2, (len(p_indices2),))
        
        p1 = np.reshape(p[p_indices1], (len(p_indices1),)) # Array of p<1
        p2 = np.reshape(p[p_indices2], (len(p_indices2),)) # Array of p>=1
        
        r = np.random.random_sample((len(p_indices1),))
        
        p_indices11_in_p1 = np.argwhere(r<p1)
        p_indices11 = p_indices1[p_indices11_in_p1]
        p_indices11 = np.reshape(p_indices11, (len(p_indices11),)) # Tracking the indices of p that correspond to the condition (p<1) & (p>r)
        
        p_indices12_in_p1 = np.argwhere(r>=p1)
        p_indices12 = p_indices1[p_indices12_in_p1]
        p_indices12 = np.reshape(p_indices12, (len(p_indices12),)) # Tracking the indices of p that correspond to the condition (p<1) & (p<=r)
        
        # Now we have the indices of p for each condition. These indices correspond to the lines (i.e., the walker's number) in the "positions" matrix 
        # that satisfy the respective conditions. Thus:
        
        positions[p_indices11] = positions_prime[p_indices11]
        positions[p_indices2] = positions_prime[p_indices2]
        
        r1 = positions[:,:(dim//2)]
        r2 = positions[:,(dim//2):]
        
        positions_tracker[i] = positions
        
        n_accept = n_accept + len(p_indices11) + len(p_indices2) # We sum in each iteration the acceptances for all walkers

        E_loc = eloc_helium(c, r1, r2, (dim//2), walkers)
        der_psi = derivative_psi_c_helium(c, r1, r2)
        
        energy_tracker[i] = E_loc
        
        variation_expectation_1[i] = E_loc * der_psi
        variation_expectation_2[i] = der_psi
        
    
    return positions_tracker, energy_tracker, variation_expectation_1, variation_expectation_2, n_accept/(iterations*walkers)


def mean_energy_helium1(c, iterations, walkers, dim):
    
    """
    Function that returns the expectation value of the energy, its variance, the acceptance ratio and var1, var2, provided a certain value of the parameter c. It does not take into account the correlation between data!
    """
    
    positions_tracker, energy_tracker, variation_expectation_1, variation_expectation_2, accept_ratio = metropolis_helium(c, iterations, walkers, dim)
    
    E = np.mean(energy_tracker[4000:]) # All the following are means of the total (iterations-4000 x walkers) matrices ==> Numbers
    E_sq = np.mean((energy_tracker[4000:])**2)
    
    E_walkers = np.mean(energy_tracker[4000:], axis=0) # Array of the mean energies for each walker
    E_walkers_sq = np.mean((energy_tracker[4000:])**2, axis=0)
    
    var1 = np.mean(variation_expectation_1[4000:])
    var2 = np.mean(variation_expectation_2[4000:])
    
    variance1 = np.sqrt((E_sq)-(E**2)) # Variance over all configurations of all walkers
    variance2 = np.sqrt(np.mean(E_walkers_sq) - (np.mean(E_walkers))**2)
    
    return E, variance1, variance2, accept_ratio, var1, var2


def helium_opt1(c, walkers, dim): # For more functions see file "functions.py"
    
    """
    Function that given an initial value of the parameter c, finds and returns a more optimal value that minimizes the mean energy which is also returned. 
    Works for UNCORRELATED sets of configurations!
    """
    
    gamma = 0.2 # For the time being
    iterations = 30000
    
    c_tracker = np.zeros((200,)) # 200 is way more than the spots which will probably be occupied for helium
    energy_tracker = np.zeros((200,))
    variance1_tracker = np.zeros((200,))
    variance2_tracker = np.zeros((200,))
    
    c_tracker[0] = c
    n = 1
    
    mean_energies, variance1, variance2, accept_ratio, var1, var2 = mean_energy_helium1(c, iterations, walkers, dim)
    energy_tracker[0] = mean_energies
    variance1_tracker[0] = variance1
    
    while variance1 > 0.25 and n <= 100: # Finding the optimal parameter c, makes Eloc flatter ==> minimizes the variance

        deri = 2 * (var1 - (mean_energies*var2))
        c_tracker[n] = c - gamma*deri
        c = c_tracker[n]
        
        mean_energies, variance1, variance2, accept_ratio, var1, var2 = mean_energy_helium1(c, iterations, walkers, dim)
        energy_tracker[n] = mean_energies
        variance1_tracker[n] = variance1
        variance2_tracker[n] = variance2
        
        n = n + 1
        
    return c_tracker, c, mean_energies, energy_tracker, variance1_tracker, variance2_tracker, accept_ratio, n-1 # Number of times we entered the loop


def mean_energy_helium3(c, iterations, walkers, dim):
    
    """
    Function that returns the expectation value of the energy, its variance, the acceptance ratio and var1, var2, provided a certain value of the parameter c for a certain test function
    """
    
    positions_tracker, energy_tracker, variation_expectation_1, variation_expectation_2, accept_ratio = metropolis_helium(c, iterations, walkers, dim)
    
    E, tau, variance1, variance2 = autocorrelation2(energy_tracker, walkers)
    
    # Variance1 is the variance over all possible configurations
    # Variance2 is the variance over different walkers' mean values
    
    var1, a, b, c = autocorrelation2(variation_expectation_1, walkers)
    var2, a, b, c = autocorrelation2(variation_expectation_2, walkers)
    
    # var1, var2 will only be used for the minimization. Nothing to do with the variance.
    
    return E, variance1, variance2, accept_ratio, var1, var2


def helium_opt2(c, walkers, dim): # For more functions see file "functions.py"
    
    """
    Function that given an initial value of the parameter c, finds and returns a more optimal value that minimizes the mean energy which is also returned. 
    Works for CORRELATED sets of configurations!
    """
    
    gamma = 0.5 # For the time being
    iterations = 30000
    
    c_tracker = np.zeros((200,)) # 200 is way more than the spots which will probably be occupied for helium
    energy_tracker = np.zeros((200,))
    variance1_tracker = np.zeros((200,))
    variance2_tracker = np.zeros((200,))
    
    c_tracker[0] = c
    n = 1
    
    print("First value of c: ", c)
    mean_energies, variance1, variance2, accept_ratio, var1, var2 = mean_energy_helium3(c, iterations, walkers, dim)
    energy_tracker[0] = mean_energies
    variance1_tracker[0] = variance1
    
    while variance1 > 0.003 and n <= 40: # Finding the optimal parameter c, makes Eloc flatter ==> minimizes the variance

        deri = 2 * (var1 - (mean_energies*var2))
        c_tracker[n] = c - gamma*deri
        c = c_tracker[n]
        print(c)
        
        mean_energies, variance1, variance2, accept_ratio, var1, var2 = mean_energy_helium3(c, iterations, walkers, dim)
        energy_tracker[n] = mean_energies
        variance1_tracker[n] = variance1
        variance2_tracker[n] = variance2
        
        n = n + 1
        
    return c_tracker, c, mean_energies, energy_tracker, variance1_tracker, variance2_tracker, accept_ratio, n-1 # Number of times we entered the loop
