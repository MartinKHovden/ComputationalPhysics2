import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
import os 

def plot_energy(method, filedir, n_particles, n_dims, iterations):
    """Function for plotting the energy as a function of optimization iteration
    
    Parameters:
    -----------
    method: string 
        String for which method. "bf" = brute-force, "is" = importance sampling and "gs" = Gibbs sampling
        
    filedir: string
        The directory of where the files are located. 
    
    """

    filedir += method + "/"
    
    interacting = filedir.split("/")[1]
    print(interacting)

    plt.rcParams["figure.figsize"] = 20,10
    plt.rcParams["lines.linewidth"] = 2
    plt.rcParams["xtick.labelsize"] =15
    plt.rcParams["ytick.labelsize"] =15
    plt.rcParams["axes.labelsize"] = 20
    plt.rcParams["legend.fontsize"] = 20
    plt.rcParams["legend.title_fontsize"] = 15
    
    for filename in os.listdir(filedir):
        
        temp1 = filename.split("_")
        
        if int(temp1[2]) == n_particles and int(temp1[5]) == n_dims and int(temp1[-1].split(".")[0]) == iterations:
            temp=temp1
            
            
            local_energies = np.loadtxt(filedir + filename, skiprows=1)
            label = "lr = " + str(temp[7]) + " , hidden nodes = " + str(temp[9])
            plt.plot(local_energies, label= label,linestyle="-", marker = ".", markersize = 10)

    method_full_name = ""

    if method == "bf":
        method_full_name += "Metropolis Brute-force"
        
    if method == "is":
        method_full_name += "Metropolis Importance Sampling"
        
    if method == "gs":
        method_full_name += "Gibbs Sampling"

    title = "Energy as a function of gradient descent iteration using " + method_full_name +  " \n with " + str(iterations) + " mc-iterations" + " for " + str(temp[2]) + " particles in " + str(temp[5]) +  " dimensions. " + interacting  


    plt.xlabel("Iteration")
    plt.ylabel("Energy")
    plt.legend(title="Parameters")
    plt.title(title, fontsize=20)
    
    filename = "plots/" + interacting + "_grid_search_" + method +"_energy_n_particles_" + str(temp[2]) + "_n_dims_" + str(temp[5]) + ".png"
    
    plt.savefig(filename)
        
    plt.show()
    

def create_table_non_interacting(method, filedir, n_particles, n_dims, iterations):
    
    optimal_energies = []
    time_values = []
    network_architectures = []
    
    filedir += (method + "/")
    
    interacting = filedir.split("/")[1]

    for filename in os.listdir(filedir):
        
        temp1 = filename.split("_")
        
#         print(temp1)
                        
        if int(temp1[2]) == n_particles and int(temp1[5]) == n_dims and int(temp1[-1].split(".")[0]) == iterations:
            
            temp = temp1
            info = open(filedir + filename).readlines()[0].split(" ")
            
            local_energies = np.loadtxt(filedir + filename, skiprows=1)
            optimal_energies.append(float(f"{np.mean(local_energies[-10:]):.4f}"))
            
            time = float(info[1].split("=")[1])
            time_values.append(float(f"{time:.4f}"))
            
            network_architectures.append("lr = " + str(temp[7]) + " , hidden nodes = " + str(temp[9]))
            
    exact_value = 0.5*n_particles*n_dims

    errors = np.array(optimal_energies) - exact_value
    relative_errors = abs(errors)/exact_value
    
    method_full_name = ""

    if method == "bf":
        method_full_name += "Metropolis Brute-force"
        cap = "Optimal energy, error and time for different network architectures using " + method_full_name + " on a non-interacting system with " + str(n_particles) + " particles in " + str(n_dims) + " dimensions. Using $\sigma^2 = $" + temp[11] + ", MC-iterations = " + temp[-1].split(".")[0] + " and MC step-length = " + temp[15] + " and " + str(len(local_energies)) + " optimization steps."
        
    if method == "is":
        method_full_name += "Metropolis Importance Sampling"
        cap = "Optimal energy, error and time for different network architectures using " + method_full_name + " on a non-interacting system with " + str(n_particles) + " particles in " + str(n_dims) + " dimensions. Using $\sigma^2 = $" + temp[11] + ", MC-iterations = " + temp[-1].split(".")[0] + " and MC step-length = " + temp[15] + " and " + str(len(local_energies)) + " optimization steps."
        
    if method == "gs":
        method_full_name += "Gibbs Sampling"
        cap = "Optimal energy, error and time for different network architectures using " + method_full_name + " on a non-interacting system with " + str(n_particles) + " particles in " + str(n_dims) + " dimensions. Using $\sigma^2 = $" + temp[11] + " and MC-iterations = " + temp[-1].split(".")[0]  + " and " + str(len(local_energies)) + " optimization steps."


    df = pd.DataFrame({"Network Architecture": network_architectures, "Optimal Energy": optimal_energies, "Abs error": abs(errors), "Rel error": relative_errors, "Time [s]": time_values}) 
    df.style.set_properties(**{"text_align": "right"})
    print(df.to_latex(index=False, caption = cap))

    
def create_table_interacting(method, filedir, n_particles, n_dims, iterations):
    
    optimal_energies = []
    time_values = []
    network_architectures = []
    
    filedir += (method + "/")
    
    interacting = filedir.split("/")[1]

    for filename in os.listdir(filedir):
        
        temp1 = filename.split("_")
        
#         print(temp1)
                        
        if int(temp1[2]) == n_particles and int(temp1[5]) == n_dims and int(temp1[-1].split(".")[0]) == iterations:
            
            temp = temp1
            info = open(filedir + filename).readlines()[0].split(" ")
            
            local_energies = np.loadtxt(filedir + filename, skiprows=1)
            optimal_energies.append(float(f"{np.mean(local_energies[-10:]):.4f}"))
            
            time = float(info[1].split("=")[1])
            time_values.append(float(f"{time:.4f}"))
            
            network_architectures.append("lr = " + str(temp[7]) + " , hidden nodes = " + str(temp[9]))
            
    exact_value = 3.0

    errors = np.array(optimal_energies) - exact_value
    relative_errors = abs(errors)/exact_value
    
    method_full_name = ""

    if method == "bf":
        method_full_name += "Metropolis Brute-force"
        cap = "Optimal energy, error and time for different network architectures using " + method_full_name + " on a interacting system with " + str(n_particles) + " particles in " + str(n_dims) + " dimensions. Using $\sigma^2 = $" + temp[11] + ", MC-iterations = " + temp[-1].split(".")[0] + " and MC step-length = " + temp[15] + " and " + str(len(local_energies)) + " optimization steps."
        
    if method == "is":
        method_full_name += "Metropolis Importance Sampling"
        cap = "Optimal energy, error and time for different network architectures using " + method_full_name + " on a interacting system with " + str(n_particles) + " particles in " + str(n_dims) + " dimensions. Using $\sigma^2 = $" + temp[11] + ", MC-iterations = " + temp[-1].split(".")[0] + " and MC step-length = " + temp[15]  + " and " + str(len(local_energies)) + " optimization steps."
        
    if method == "gs":
        method_full_name += "Gibbs Sampling"
        cap = "Optimal energy, error and time for different network architectures using " + method_full_name + " on a interacting system with " + str(n_particles) + " particles in " + str(n_dims) + " dimensions. Using $\sigma^2 = $" + temp[11] + " and MC-iterations = " + temp[-1].split(".")[0]   + " and " + str(len(local_energies)) + " optimization steps."


    df = pd.DataFrame({"Network Architecture": network_architectures, "Optimal Energy": optimal_energies, "Abs error": abs(errors), "Rel error": relative_errors, "Time [s]": time_values}) 
    df.style.set_properties(**{"text_align": "right"})
    print(df.to_latex(index=False, caption = cap))
    