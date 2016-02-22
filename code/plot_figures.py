import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import cPickle as pickle

def plot_samples(data, sample, fout, close=True):
    """
    FIRST PLOT: SPECTRUM + SAMPLES FROM THE POSTERIOR
    """
    
    # number of posterior samples
    nsamples = sample.shape[0]

    # randomly pick some samples from the posterior to plot
    s_ind = np.random.choice(np.arange(nsamples, dtype=int), size=20)

    # the middle of the wavelength bins for plotting
    wmid = data["wavelength_left"] + (data["wavelength_right"] - data["wavelength_left"])/2.

    # the error on the data
    yerr = np.zeros_like(wmid) + data['err']

    plt.figure(figsize=(14,6))
    plt.errorbar(wmid, data["fake_flux"], yerr=yerr, fmt="o")
    for i in s_ind:
        plt.plot(wmid, sample[i,-wmid.shape[0]:], lw=2, alpha=0.7)

    plt.xlim(wmid[0], wmid[-1])
    plt.gca().invert_xaxis()
    plt.xlabel("Wavelength [Angstrom]")
    plt.ylabel("Normalized Flux")
    plt.tight_layout()
    plt.savefig(fout+"_samples.png", format="png")
    if close:
        plt.close()
    
    return

def plot_bkg(data, sample, fout, close=True):
    """
    PLOT THE BACKGROUND POSTERIOR
    """
    fig = plt.figure(figsize=(12,9))
    # Plot a histogram and kernel density estimate
    ax = fig.add_subplot(111)
    sns.distplot(sample[:,0], hist_kws={"histtype":"stepfilled"}, ax=ax)
    plt.xticks(rotation=45)

    _, ymax = ax.get_ylim()
    ax.vlines(np.exp(data["logbkg"]), 0, ymax, lw=3, color="black")

    ax.set_xlabel("Normalized Background Flux")
    ax.set_ylabel("$N_{\mathrm{samples}}$")
    plt.tight_layout()
    plt.savefig(fout+"_bkg.png", format="png")
    if close:
        plt.close()
    return
    
def plot_ou_bkg(sample, fout, close=True):
    """
    PLOT THE POSTERIOR FOR THE OU PROCESS
    """
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(14,6))
    sns.distplot(sample[:,1], hist_kws={"histtype":"stepfilled"}, ax=ax1)
    #_, ymax = ax.get_ylim()
    #ax.vlines(np.exp(data["logbkg"]), 0, ymax, lw=3, color="black")

    ax1.set_xlabel(r"OU time scale $\tau$")
    ax1.set_ylabel(r"$N_{\mathrm{samples}}$")

    sns.distplot(sample[:,2], hist_kws={"histtype":"stepfilled"}, ax=ax2)
    plt.xticks(rotation=45)

    #_, ymax = ax.get_ylim()
    #ax.vlines(np.exp(data["logbkg"]), 0, ymax, lw=3, color="black")

    ax2.set_xlabel(r"OU amplitude")
    ax2.set_ylabel(r"$N_{\mathrm{samples}}$")

    fig.tight_layout()

    plt.savefig(fout+"_ou.png", format="png")
    if close:
        plt.close()
    return

def plot_logamp_hyper(data, sample, fout, close=True):
    """
    PLOT THE POSTERIOR FOR THE LOG-AMP HYPERPARAMETERS
    """    
    
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(14,6))
    sns.distplot(np.log(sample[:,5]), hist_kws={"histtype":"stepfilled"}, ax=ax1)
    _, ymax = ax1.get_ylim()
    ax1.vlines(data["logamp"], 0, ymax, lw=3, color="black")

    ax1.set_xlabel(r"$\mu_{\mathrm{\log{A}}}$")
    ax1.set_ylabel(r"$N_{\mathrm{samples}}$")
    ax1.set_title("Location parameter of the amplitude prior")


    sns.distplot(sample[:,6], hist_kws={"histtype":"stepfilled"}, ax=ax2)
    plt.xticks(rotation=45)

    #_, ymax = ax.get_ylim()
    #ax.vlines(np.exp(data["logbkg"]), 0, ymax, lw=3, color="black")

    ax2.set_xlabel(r"$\sigma_{\mathrm{\log{A}}}$")
    ax2.set_ylabel(r"$N_{\mathrm{samples}}$")
    ax2.set_title("Scale parameter of the amplitude prior")

    fig.tight_layout()

    plt.savefig(fout+"_logamp_prior.png", format="png")
    if close:
        plt.close()
    return

def plot_logq_hyper(data, sample, fout, close=True):
    """
    PLOT THE POSTERIOR FOR THE LOG-Q HYPERPARAMETERS
    """
    
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(14,6))
    sns.distplot(sample[:,7], hist_kws={"histtype":"stepfilled"}, ax=ax1)
    _, ymax = ax1.get_ylim()
    ax1.vlines(data["logq"], 0, ymax, lw=3, color="black")

    ax1.set_xlabel(r"$\mu_{\mathrm{\log{q}}}$")
    ax1.set_ylabel(r"$N_{\mathrm{samples}}$")
    ax1.set_title(r"Location parameter of the $q$ prior")


    sns.distplot(sample[:,8], hist_kws={"histtype":"stepfilled"}, ax=ax2)
    plt.xticks(rotation=45)

    #_, ymax = ax.get_ylim()
    #ax.vlines(np.exp(data["logbkg"]), 0, ymax, lw=3, color="black")

    ax2.set_xlabel(r"$\sigma_{\mathrm{\log{q}}}$")
    ax2.set_ylabel(r"$N_{\mathrm{samples}}$")
    ax2.set_title(r"Scale parameter of the $q$ prior")

    fig.tight_layout()

    plt.savefig(fout+"_logq_prior.png", format="png")
    if close:
        plt.close()
    return

def plot_threshold(sample, fout, close=True):
    """
    PLOT THE POSTERIOR FOR THE THRESHOLD PARAMETER
    """
    
    fig = plt.figure(figsize=(12,9))
    # Plot a historgram and kernel density estimate
    ax = fig.add_subplot(111)
    sns.distplot(sample[:,9], hist_kws={"histtype":"stepfilled"}, ax=ax)
    plt.xticks(rotation=45)

    ax.set_xlabel("Threshold parameter")
    ax.set_ylabel("$N_{\mathrm{samples}}$")
    plt.tight_layout()

    plt.savefig(fout+"_pp.png", format="png")
    if close:
        plt.close()
    return

def plot_dshift(data, sample, fout, dshift_ind=0, close=True):
    """ 
    PLOT THE POSTERIOR FOR THE DOPPLER SHIFT
    """
    
    fig = plt.figure(figsize=(12,9))
    plt.locator_params(axis = 'x', nbins = 6)

    # Plot a historgram and kernel density estimate
    ax = fig.add_subplot(111)
    sns.distplot(sample[:,11+dshift_ind], hist_kws={"histtype":"stepfilled"}, ax=ax)
    plt.xticks(rotation=45)


    _, ymax = ax.get_ylim()
    ax.vlines(data["dshift"], 0, ymax, lw=3, color="black")

    ax.set_xlabel(r"Doppler shift $d=v/c$")
    ax.set_ylabel("$N_{\mathrm{samples}}$")
    plt.tight_layout()

    plt.savefig(fout+"_dshift%i.png"%dshift_ind, format="png")
    if close:
        plt.close()
    return

def plot_logamp(data, sample, fout, ndshift, nlines, ncolumns=3, 
                dshift_ind=0, close=True):
    """
    PLOT THE POSTERIOR FOR THE LINE LOG-AMPLITUDES
    """
    
    nrows = int(nlines/ncolumns)+1

    fig = plt.figure(figsize=(nrows*4,ncolumns*4))
    plt.locator_params(axis = 'x', nbins = 6)

    # index of column where the log-amplitudes start:
    start_ind = 11 + ndshift + dshift_ind*nlines
    
    # log-amplitudes
    for i in range(nlines):
        ax = plt.subplot(nrows, ncolumns, i+1)    
#        ax.hist(sample[:,start_ind+i], histtype="stepfilled", alpha=0.7)
        sns.distplot(sample[:,start_ind+i], hist_kws={"histtype":"stepfilled"}, ax=ax)

        plt.locator_params(axis = 'x', nbins = 6)
        xlabels = ax.get_xticklabels() 
        for l in xlabels: 
            l.set_rotation(45) 

        _, ymax = ax.get_ylim()
        ax.vlines(data["logamp"][i], 0, ymax, lw=3, color="black")

        ax.set_xlabel(r"$\log{A}$")
        ax.set_ylabel("$N_{\mathrm{samples}}$")

    plt.tight_layout()
    if dshift_ind == 0:
        plt.savefig(fout+"_logamp.png", format="png")
    else:
        plt.savefig(fout+"_logamp%i.png"%dshift_ind, format="png")

    if close:
        plt.close()
    
    return

def plot_logq(data, sample, fout, ndshift, nlines, ncolumns=3, 
              dshift_ind=0, close=True):
    """
    PLOT THE POSTERIOR FOR THE LINE LOG-Q
    """

    nrows = int(nlines/ncolumns)+1

    fig = plt.figure(figsize=(nrows*4,ncolumns*4))
    plt.locator_params(axis = 'x', nbins = 6)

    # set starting index for the logq values:
    start_ind = 11 + ndshift + nlines*(dshift_ind + 1)
    
    # log-amplitudes
    for i in range(nlines):
        ax = plt.subplot(nrows, ncolumns, i+1)    
        #ax.hist(sample[:,start_ind+i], histtype="stepfilled", alpha=0.7)
        sns.distplot(sample[:,start_ind+i], hist_kws={"histtype":"stepfilled"}, ax=ax)
        
        plt.locator_params(axis = 'x', nbins = 6)
        xlabels = ax.get_xticklabels() 
        for l in xlabels: 
            l.set_rotation(45) 

        _, ymax = ax.get_ylim()
        ax.vlines(data["logq"][i], 0, ymax, lw=3, color="black")

        ax.set_xlabel(r"$\log{q}$")
        ax.set_ylabel("$N_{\mathrm{samples}}$")

    plt.tight_layout()
    if dshift_ind == 0:
        plt.savefig(fout+"_logq.png", format="png")
    else:
        plt.savefig(fout+"_logq%i.png"%dshift_ind, format="png")
    
    if close:
        plt.close()
    return


def plot_signs(data, sample, fout, ndshift, nlines, ncolumns=3, 
               dshift_ind=0, close=True):
    """
    PLOT THE POSTERIOR FOR THE LINE AMPLITUDE SIGNS
    """

    nrows = int(nlines/ncolumns)+1
    
    fig = plt.figure(figsize=(nrows*4,ncolumns*4))
    plt.locator_params(axis = 'x', nbins = 6)

    # set starting index for the logq values:
    start_ind = 11 + ndshift + nlines*(dshift_ind + 2)
    
    # log-amplitudes
    for i in range(nlines):
        ax = plt.subplot(nrows, ncolumns, i+1)    
#        ax.hist(sample[:,start_ind+i], histtype="stepfilled", alpha=0.7)
        sns.distplot(sample[:,start_ind+i], hist_kws={"histtype":"stepfilled"}, ax=ax)

        plt.locator_params(axis = 'x', nbins = 6)
        xlabels = ax.get_xticklabels() 
        for l in xlabels: 
            l.set_rotation(45) 

        ax.set_xlabel(r"Emission/absorption line sign")
        ax.set_ylabel("$N_{\mathrm{samples}}$")

    plt.tight_layout()
    if dshift_ind == 0:
        plt.savefig(fout+"_signs.png", format="png")
    else:
        plt.savefig(fout+"_signs%i.png"%dshift_ind, format="png")
    
    if close:
        plt.close()
    return


def plot_posterior_summary(froot, datadir="../data/", nlines=8, 
                           ndshift=1, ncolumns=3, close=True):
    """
    Plot summeries of the posterior distribution. Mostly histograms.
    Plots a bunch of Figures to png files.
    
    Parameters
    ----------
    froot: str
        The root string of the data file and file with posterior samples to be loaded.
        
    datadir: str
        The directory with the data.
        Default: "../data/"
        
    nlines: int
        The number of lines in the model.
        Default: 8
        
    ndshift: int
        The number of (possible) Doppler shifts in the model.
        Default: 1
        
    ncolumns: int
        The number of columns in multi-panel plots. Default: 3
    
    close: bool
        Close plots at the end of the plotting? Default: True
    """
    
    
    # the pickle file with the data + parameters:
    f = open(datadir+froot+"_data.pkl")
    data = pickle.load(f)
    f.close()

    print("Keys in data dictionary: " + str(data.keys()))

    # the posterior samples
    sample = np.atleast_2d(np.loadtxt(datadir+froot+"_posterior_sample.txt"))
    nsamples = sample.shape[0]
    print("We have %i samples from the posterior."%nsamples)
    
    # set the directory path and file name for output files:
    fout = datadir+froot
    
    # plot the spectrum with some draws from the posterior
    plot_samples(data, sample, fout, close=close)

    # plot a histogram of the background parameter
    plot_bkg(data, sample, fout, close=close)
    
    # plot histograms of the OU parameters
    plot_ou_bkg(sample, fout, close=close)
    
    # plot the hyper parameters of the log-amp prior
    plot_logamp_hyper(data, sample, fout, close=close)
    
    # plot the hyper parameters of the log-q prior
    plot_logq_hyper(data, sample, fout, close=close)
    
    # plot the threshold for the amplitude sign
    plot_threshold(sample, fout, close=close)
    
    # for the varying number of Doppler shifts, plot their posterior
    if ndshift == 1:
        plot_dshift(data, sample, fout, dshift_ind=0, close=close)
        plot_logamp(data, sample, fout, ndshift, nlines, 
                    ncolumns=ncolumns, dshift_ind=0, close=close)
        plot_logq(data, sample, fout, ndshift, nlines, 
                  ncolumns=ncolumns, dshift_ind=0, close=close)
        plot_signs(data, sample, fout, ndshift, nlines, 
                   ncolumns=ncolumns, dshift_ind=0, close=close)        
    else:
        for i in range(ndshift):
            plot_dshift(data, sample, fout, dshift_ind=i, close=close)
            plot_logamp(data, sample, fout, ndshift, nlines, 
                        ncolumns=ncolumns, dshift_ind=i, close=close)
            plot_logq(data, sample, fout, ndshift, nlines, 
                      ncolumns=ncolumns, dshift_ind=i, close=close)
            plot_signs(data, sample, fout, ndshift, nlines, 
                       ncolumns=ncolumns, dshift_ind=i, close=close)        

    return


    
