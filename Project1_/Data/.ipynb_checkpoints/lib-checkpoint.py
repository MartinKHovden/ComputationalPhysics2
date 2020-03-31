import numpy as np
import matplotlib.pyplot as plt

def extractInfo(file, alpha, num_dims, num_particles):
    """Function for extracting the information of interest. Assumes that the 
    file contains one last line with computational time of the code. 
    """
    
    file = open(file, "r")
    
    temp = np.array([i for i in file])
    local_energies = np.array([float(j) for j in temp[:-1]])
    
    time = temp[-1].split(" ")[-1]
    number_of_metropolis_steps = len(local_energies)
    
    total_energy = np.mean(local_energies)
    
    total_enegy_bootstrap = tsboot(local_energies, 2**12, 2**10)
    energy_per_particle = total_energy/float(num_particles)
    
    blocking_total_energy, blocking_variance = block(local_energies)
    
    total_energy_development = np.cumsum(local_energies)
    total_energy_development = np.array([energy/float(num+1) for num, energy in enumerate(total_energy_development)])
    
    accepted_steps = 0
    for i in range(1, len(local_energies)):
        if(local_energies[i] != local_energies[i-1]):
            accepted_steps += 1
    
    acceptance_ratio = accepted_steps/number_of_metropolis_steps
    
    
    return total_energy, energy_per_particle, blocking_variance, total_energy_development, time, acceptance_ratio, total_energy_bootstrap


# from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, sqrt
from numpy.linalg import inv

def block(x):
    # preliminaries
    n = len(x)
    d = int(np.log2(n))
    s, gamma = np.zeros(d), np.zeros(d)
    mu = np.mean(x)

    # estimate the auto-covariance and variances 
    # for each blocking transformation
    for i in np.arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*np.sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = np.var(x)
        # perform blocking transformation
        x = 0.5*(x[0:-1:2] + x[1::2])
   
    # generate the test observator M_k from the theorem
    M = (np.cumsum( ((gamma/s)**2*2**np.arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =np.array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in np.arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    return mu, s[k]/2**(d-k)

def exactValue(alpha, num_dims, num_particles):
    """Returns the exact energy for the non-interacting case in a spherical trap. 
    """
    return (alpha/2. + 1./(8*alpha))*num_particles*num_dims

def tsboot(data,statistic,R,l):
    t = zeros(R); n = len(data); k = int(ceil(float(n)/l));
    inds = arange(n); t0 = time()
    
    # time series bootstrap
    for i in range(R):
        # construct bootstrap sample from
        # k chunks of data. The chunksize is l
        _data = concatenate([data[j:j+l] for j in randint(0,n-l,k)])[0:n];
        t[i] = statistic(_data)

    # analysis
    print ("Runtime: %g sec" % (time()-t0)); print ("Bootstrap Statistics :")
    print ("original           bias      std. error")
    print ("%8g %14g %15g" % (statistic(data), \
                             mean(t) - statistic(data), \
                             std(t) ))
    return t

def stat(data):
    return mean(data)
