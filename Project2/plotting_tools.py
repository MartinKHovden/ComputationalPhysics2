import numpy as np 
import matplotlib.pyplot as plt 
import os 

def plot_energy(method, filedir):
    """Function for plotting the energy as a function of optimization iteration
    
    Parameters:
    -----------
    method: string 
        String for which method. "bf" = brute-force, "is" = importance sampling and "gs" = Gibbs sampling
        
    filedir: string
        The directory of where the files are located. 
    
    """

    filedir += method + "/"

    plt.rcParams["figure.figsize"] = 20,10
    plt.rcParams["lines.linewidth"] = 2
    plt.rcParams["xtick.labelsize"] =15
    plt.rcParams["ytick.labelsize"] =15
    plt.rcParams["axes.labelsize"] = 20
    plt.rcParams["legend.fontsize"] = 20
    plt.rcParams["legend.title_fontsize"] = 15


    for filename in os.listdir(filedir):
        local_energies = np.loadtxt(filedir + filename, skiprows=1)
        temp = filename.split("_")
        label = "lr = " + str(temp[1]) + " , hidden nodes = " + str(temp[3])
        plt.plot(local_energies, label= label)

    method_full_name = ""

    if method == "bf":
        method_full_name += "Metropolis Brute-force"
    if method == "is":
        method_full_name += "Metropolis Importance Sampling"
    if method == "gs":
        method_full_name += "Gibbs Sampling"

    title = "Energy as a function of gradient descent iteration using " + method_full_name +  " with " + str(temp[-1].split(".")[0]) + " mc-iterations" 


    plt.xlabel("Iteration")
    plt.ylabel("Energy")
    plt.legend(title="Parameters")
    plt.title(title, fontsize=15)
    plt.show()