import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np


def plot_results(data_file, posterior_file, output_root):
    """
    Plot results for a data set. 
    The initial run here had one fixed redshift (which could have been zero), 
    which means that the redshift parameter is in column 9 in the posterior_samples.txt 
    file. Note: This is hard-coded below! If the model changes, especially the 
    number of hyperparameters, then those columns should change, too!

    """ 

    ndshift_column = 9 # column with the correct Doppler shift distribution

    data = np.loadtxt(data_file)
    sample = np.atleast_2d(np.loadtxt(posterior_file))

    plt.figure(figsize=(16,6))
    plt.plot(data[:,0], data[:,2], linestyle="steps-mid")
    for s in sample[:20]:
        plt.plot(data[:,0], s[-len(data[:,0]):])

    plt.xlabel("Wavelength [Angstrom]")
    plt.ylabel("Flux")
    plt.title("First Test: No Doppler shift, all lines same width and amplitude")
    plt.savefig(output_root+"_samples.png", format="png")
    plt.close()

    plt.figure(figsize=(12,9))
    plt.hist(sample[:,ndshift_column], bins=30, histtype="stepfilled")
    plt.xlabel("Doppler shift $v/c$")
    plt.ylabel("p(Doppler shift)")
    plt.savefig(output_root+"_dshift.png", format="png")
    plt.close()

    return


def plot_test_noshift1(datadir="../data/"):
    data_file = datadir+"test_noshift1.dat"
    posterior_file = datadir+"test_noshift1_posterior_sample.txt"
    output_root = datadir+"test_noshift1"
    plot_results(data_file, posterior_file, output_root)
    return


def plot_test_noshift2(datadir="../data/"):
    data_file = datadir+"test_noshift2.dat"
    posterior_file = datadir+"test_noshift2_posterior_sample.txt"
    output_root = datadir+"test_noshift2"
    plot_results(data_file, posterior_file, output_root)
    return

