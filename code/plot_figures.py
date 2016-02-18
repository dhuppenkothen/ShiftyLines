import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np


def plot_test_noshift1(datadir="../data/"):
    """
    Plot results for the data set test_noshift1.dat. 
    The initial run here had one fixed redshift (which could have been zero), 
    which means that the redshift parameter is in column 9 in the posterior_samples.txt 
    file. Note: This is hard-coded below! If the model changes, especially the 
    number of hyperparameters, then those columns should change, too!

    """ 

    ndshift_column = 9 # column with the correct Doppler shift distribution

    data = np.loadtxt(datadir+"test_noshift1.dat")
    sample = np.atleast_2d(np.loadtxt("test_noshift1_posterior_sample.txt"))

    plt.figure(figsize=(16,6))
    plt.plot(data[:,0], data[:,2], linestyle="steps-mid")
    for s in sample[:20]:
        plt.plot(data[:,0], s[-len(data[:,0]):])

    plt.xlabel("Wavelength [Angstrom]")
    plt.ylabel("Flux")
    plt.title("First Test: No Doppler shift, all lines same width and amplitude")
    plt.savefig(datadir+"test_noshift1_samples.png", format="png")
    plt.close()

    plt.figure(figsize=(12,9))
    plt.hist(sample[:,ndshift_column], bins=30, histtype="stepfilled")
    plt.xlabel("Doppler shift $v/c$")
    plt.ylabel("p(Doppler shift)")
    plt.savefig(datadir+"test_noshift1_dshift.png", format="png")
    plt.close()

    return




