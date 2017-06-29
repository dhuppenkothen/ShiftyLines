
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd

from clarsach import XSpectrum, ARF, RMF

import astropy.constants
import astropy.units as u

def load_data(filenames, names="None"):
    """
    Load one or more data files from a FITS file.

    **Warning: Currently only works for Chandra/HETG!**

    Parameters
    ----------
    filenames : str or iterable of str
        A single file name or a list of file names

    names : str or list of str, optional
        Some designation of the data files in filenames;
        Useful for distinguishing multiple files.

    Returns
    -------
    data : dict
        A dictionary of {"name": XSpectrum} pairs
    """

    # If names are not given, then just label them
    # 1, 2, ... etc.
    if names is None:
        if isinstance(filenames, str):
            names = ["1"]
            filenames = [filenames]
        else:
            names = ["%i"%i for i in range(len(filenames))]

    # make an empty dictionary for the data
    data = {}

    # loop over all file names, extract a spectrum from
    # each and store it in the dictionary
    for n,f in zip(names, filenames):
        spec = XSpectrum(filename=f)
        data[n] = spec

    return data


def plot_ndoppler(samples, idx=14):
    """
    Plot the number of Doppler shifts that the model infers.

    Parameters
    ----------
    samples : numpy.ndarray of shape (nsamples, nparameters)
        An array with samples from a DNest4 run (e.g. by importing
        posterior_sample.txt)

    idx : int
        The index of the column where the Doppler shift lives in the
        posterior sample file.

    Returns:
    fig, ax : matplotlib objects
        The plotting objects for further use.

    """


    # count the number of Doppler shifts for all
    # samples
    s = pd.Series(samples[:, idx]).value_counts()

    fig, ax = plt.subplots(1, 1, figsize=(8,5))

    sns.barplot(s.index, s.values/s.values.sum(), ax=ax)
    ax.set_title("Posterior distribution of the number of Doppler shifts")
    ax.set_xlabel("Number of Doppler shifts")
    ax.set_ylabel("Probability")

    return fig, ax


def convert_doppler_to_velocity(doppler):
    """
    Convert a Doppler shift to velocity.

    The Doppler shift d is defined such that:

        $\lambda_d = \lambda (1 + d)$

    Or for energies rather than wavelengths:

        $E_d = E / (1 + d)$

    We convert to velocity in km/s by multiplying this Doppler shift
    with the speed of light c.

        $v = c * d$

    Parameters
    ----------
    doppler : float
        Doppler shift

    Returns
    -------
    vel : float
        Velocity in km/s

    """
    vel = astropy.constants.c*doppler
    return vel.to(u.km/u.s).value


def plot_doppler_shifts(samples, ndoppler, idx=15, velocity=False):
    """
    Plot the Doppler shifts found in a spectrum.

    Parameters
    ----------
    samples : numpy.ndarray
        The array of samples from a DNest4 run. The array has the shape
        (nsamples, npars), where nsamples is the number of samples in the
        run, and npars is the number of parameters (which also includes the
        model values if those had been printed out).

    ndoppler : int
        The maximum number of Doppler shifts allowed by the model.

    idx : int
        The index of the column containing the first Doppler shift

    velocity : bool, default False
        Should the code plot the distribution of Doppler shifts in units of
        dimensionless Doppler shift, or in units of velocity (i.e. km/s)?


    Returns
    -------
    fig, axes : matplotlib Figure and Axes objects

    """

    if ndoppler <= 2:
        fig, axes = plt.subplots(1, ndoppler, figsize=(3 * ndoppler, 3))

    else:
        nrows = int(np.ceil(ndoppler / 4.))
        fig, axes = plt.subplots(nrows, 4, figsize=(12, nrows * 3))
        axes = np.hstack(axes)

    for n in range(ndoppler):
        s = samples[:, idx + n]

        if velocity:
            vel = list(map(convert_doppler_to_velocity, s))
            xlabel = "Doppler shift in km/s"
        else:
            vel = s
            xlabel = "Doppler shift"

        if ndoppler == 1:
            sns.distplot(np.array(vel), ax=axes)
            axes.set_xlabel(xlabel)
        else:
            sns.distplot(np.array(vel), ax=axes[n])
            axes[n].set_xlabel(xlabel)

    plt.tight_layout()
    return fig, axes


def make_flux_spectra(spec, counts=None):
    """
    Making flux spectra out of a a single data set.

    Parameters
    ----------
    spec : clarsach.XSpectrum object
        An XSpectrum object containing the X-ray spectrum

    counts : iterable
        An array with counts to be converted to the flux
        If this is None, use the `counts` attribute of the `spec`
        object.

    Returns
    -------
    flux_spectrum:
        The flux spectrum

    """

    # make a flat spectrum so that I can integrate
    # ARF and RMF only
    flat_model = np.ones_like(spec.counts) * (
    spec.rmf.energ_hi - spec.rmf.energ_lo)

    # now apply ARF to the flat model
    m_arf = spec.arf.apply_arf(flat_model)

    # apply RMF to the flat model
    m_rmf = spec.rmf.apply_rmf(m_arf)

    # divide the observed counts by the flat model with the ARF/RMF applied
    if counts is None:
        flux_spec = spec.counts / m_rmf / spec.exposure
    else:
        flux_spec = counts / m_rmf / spec.exposure

    return flux_spec


def plot_data_model(data, samples, flux=False,
                    min_kev=1.95, max_kev=2.05, ncols=2, lines=None, dshift=0.0):

    """
    Plot the data (i.e. the spectrum or spectra) together with
    posterior realizations of the model.

    Parameters
    ----------
    data : dictionary
        A dictionary with the spectrum or spectra

    samples : numpy.ndarray
        The array of samples from a DNest4 run. The array has the shape
        (nsamples, npars), where nsamples is the number of samples in the
        run, and npars is the number of parameters (which also includes the
        model values if those had been printed out).


    flux : bool, default False
        Plot the spectra in units of counts or units of flux?

    min_kev : float
        Lower edge of the part of the spectrum that is of interest. Used
        if the model actually only computes a physical model for part of the
        spectrum, and also to set x-axis limits for plotting

    max_kev : float
        Upper edge of the part of the spectrum that is of interest. Used
        if the model actually only computes a physical model for part of the
        spectrum, and also to set x-axis limits for plotting

    ncols : int
        The number of columns used to make the grid of subplots. The number
        of rows will be automatically computed.

    lines : iterable, default None
        A list with line rest energies in keV for plotting the
        location of the lines on top of the spectrum. If None, no lines will
        be plotted.

    dshift : float, default 0
        A Doppler shift to be applied to the the lines before plotting them.
        Can be used to show where the centroid of the line is given a Doppler
        shift either assumed or obsered.

    Returns
    -------
    fig, axes : matplotlib Figure and Axes objects

    """
    nrows = int(np.ceil(len(data.keys()) / 2.))

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * nrows, 3 * ncols))
    axes = np.hstack(axes)

    len_all = 0

    keys = ["MEG M1", "MEG P1", "HEG M1", "HEG P1"]
    for i, k in enumerate(keys):

        d = data[k]

        print("key: " + keys[i])
        idx = np.argwhere(
            ((d.bin_lo >= min_kev) & (d.bin_hi <= max_kev))).flatten()

        # note: np.max(idx) is the *last* element of the relevant array, which will
        # be excluded if it is used as the last index in an array slice. So we'll use
        # max(idx)+1 as the correct value
        min_idx = np.min(idx)
        max_idx = np.max(idx) + 1

        len_data = len(idx)
        print("len(data): " + str(len_data))

        if flux:
            d.flux = make_flux_spectra(d)

            axes[-1 - i].plot(d.bin_lo[min_idx:max_idx],
                              d.flux[min_idx:max_idx], color="black", lw=2,
                              linestyle="steps-mid")
        else:
            axes[-1 - i].plot(d.bin_lo[min_idx:max_idx],
                              d.counts[min_idx:max_idx], color="black", lw=2,
                              linestyle="steps-mid")

        for s in samples[-50:]:

            if len_all == 0:
                model_counts = s[-len_all - len_data:]
            else:
                model_counts = s[-len_all - len_data:-len_all]

            if flux:
                model_pad = np.zeros_like(d.bin_lo)
                model_pad[min_idx:max_idx] = model_counts

                model_flux = make_flux_spectra(d, model_pad)
                axes[-1 - i].plot(d.bin_lo[min_idx:max_idx],
                                  model_flux[min_idx:max_idx], c="red",
                                  alpha=0.3)

            else:
                axes[-1 - i].plot(d.bin_lo[min_idx:max_idx], model_counts,
                                  c="red", alpha=0.3)

            axes[-1 - i].set_xlim(min_kev, max_kev)

        if not lines is None:
            ax_ylim = axes[-1 - i].get_ylim()

            for l in lines:
                ax_ylim = axes[-1 - i].get_ylim()

                axes[-1 - i].vlines(l / (1. + dshift), ax_ylim[0],
                                    0.8 * ax_ylim[1], lw=3,
                                    color=sns.color_palette()[1])

        len_all += len_data

        axes[i].set_xlabel("Wavelength in keV")

        if flux:
            axes[i].set_ylabel("Photon flux")
        else:
            axes[i].set_ylabel("Counts")

        axes[i].set_title(k)

    plt.tight_layout()

    return fig, axes


def extract_parameters(samples, ndoppler, nlines, istart=17, pp_ind=13):

    """
    Extract the line parameters from a DNest4 run. Current parameters are:
        * line amplitude
        * line width
        * sign parameter determining whether we have an emission or
          absorption line. This parameter in the code is between 0 and 1
        * numerical translation of the sign parameter, either -1 (absorption)
          or 1 (emission)
        * presence parameter determining whether the line exists in the data.
          This parameter in the code is between 0 and 1
        * numerical translation of the presence parameter, either 0 (line absent)
          or 1 (line present)

    Parameters
    ----------
    samples : numpy.ndarray
        The array of samples from a DNest4 run. The array has the shape
        (nsamples, npars), where nsamples is the number of samples in the
        run, and npars is the number of parameters (which also includes the
        model values if those had been printed out).

    ndoppler : int
        The maximum number of Doppler shifts allowed by the model.

    nlines : int
        The number of possible lines in the data

    istart : int, default 17
        The index where the columns start that contain the parameters of the
        lines themselves. Depends on the number of interim parameters sampled
        for the priors as well as the number of allowed Doppler shifts

    pp_ind : int, default 13
        The index of the interim parameter for the sign parameter, which sets
        whether the sign parameter translates to either -1 for an absorption
        line or +1 for an emission line. Depends on the number of other interim
        parameters in the model.

    Returns
    -------
    par_dict : dictionary
        A dictionary with the parameters. Each element in the dictionary has a
        numerical key ranging from 0 to `ndoppler`. Each element is itself
        a dictionary with the following keys:
            * `amp`: line amplitude
            * `width`: line width
            * `sign`: absorption/emission line model parameter
            * `sign_num`: binary (-1/1) absorption/emission line parameter
            * `presence`: line presence model parameter
            * `presence_num`: binary (0,1) line presence model parameter
        Each of these contains an array of shape (nsamples, nlines), so each
        column corresponds to the posterior sample of a single line.

    """

    par_dict = {}

    for n in range(ndoppler):
        amp = np.zeros((samples.shape[0], nlines))
        width = np.zeros((samples.shape[0], nlines))
        sign = np.zeros((samples.shape[0], nlines))
        presence = np.zeros((samples.shape[0], nlines))

        for i, s in enumerate(samples):
            amp[i, :] = np.exp(
                [s[istart + n + (ndoppler * i)] for i in range(nlines)])
            width[i, :] = np.exp(
                [s[istart + (nlines * ndoppler) + n + (ndoppler * i)] for i in
                 range(nlines)])
            sign[i, :] = [
                s[istart + 2 * (nlines * ndoppler) + n + (ndoppler * i)] for i
                in range(nlines)]
            presence[i, :] = [
                s[istart + 3 * (nlines * ndoppler) + n + (ndoppler * i)] for i
                in range(nlines)]

        pp = samples[:, pp_ind]
        pp_presence = samples[:, pp_ind + 1]

        sign_num = np.zeros_like(sign.T)
        # make a mask: -1 for sign < pp, +1 for sign > pp
        sign_num[sign.T < pp] = -1.0
        sign_num[sign.T >= pp] = 1.0

        presence_num = np.zeros_like(presence.T)
        # make a mask: -1 for sign < pp, +1 for sign > pp
        presence_num[presence.T < pp_presence] = 0.0
        presence_num[presence.T >= pp_presence] = 1.0

        doppler_dict = {"amp": amp,
                        "width": width,
                        "sign": sign,
                        "sign_num": sign_num.T,
                        "presence": presence,
                        "presence_num": presence_num.T}

        par_dict["%i" % n] = doppler_dict

    return par_dict


def plot_line_parameter(par_dict, lines, parname, ncols=4, bins=20, kde=True):

    """
    Plot the marginal posterior distributions of a single line parameter
    for all lines in the model.

    Parameters
    ----------
    par_dict : dictionary
        A dictionary of parameters as returned by the `extract_parameters()`
        function.

    lines : iterable
        The list of possible lines in the model.

    parname : str
        The name of the parameter to plot. One of
            * `amp`
            * `width`
            * `sign`
            * `sign_num`
            * `presence`
            * `presence_num`

    ncols : int, default 4
        The number of columns in the subplot matrix.

    bins : {int|iterable}, default 20
        The number of bins or an array of bins for the histogram. Takes any
        valid input to `matplotlib.pyplot.hist`.

    kde : bool, default True
        Should the plot include a KDE estimate?

    Returns
    -------
    fig, axes : matplotlib Figure and Axes objects

    """
    nlines = len(lines)

    par = par_dict[parname]

    if nlines <= ncols:
        nrows = 1
        ncols = nlines

    else:
        nrows = int(np.ceil(nlines / 4.))

    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * nrows, 3 * ncols))

    if nlines != 1:
        axes = np.hstack(axes)
    print("bins: " + str(bins))
    for i, l in enumerate(lines):
        if nlines == 1:
            sns.distplot(np.array(par[:, i]), ax=axes, bins=bins, kde=kde)
            axes.set_xlabel("%s" % parname)
            axes.set_title("%.3f keV" % l)
        else:
            sns.distplot(np.array(par[:, i]), ax=axes[i], bins=bins, kde=kde)
            axes[i].set_xlabel("%s" % parname)
            axes[i].set_title("%.3f keV" % l)

    return fig, axes


def plot_amplitude(par_dict, lines, doppler_ind=0, ncols=4, bins=20):
    """
    Plot the marginal posterior distributions of the line amplitude
    for all lines in the model.

    Parameters
    ----------
    par_dict : dictionary
        A dictionary of parameters as returned by the `extract_parameters()`
        function.

    lines : iterable
        The list of possible lines in the model.

    doppler_ind : int, default 0
        If the model allows multiple Doppler shifts, use this to pick the
        Doppler shift of interest

    ncols : int, default 4
        The number of columns in the subplot matrix.

    bins : {int|iterable}, default 20
        The number of bins or an array of bins for the histogram. Takes any
        valid input to `matplotlib.pyplot.hist`.

    Returns
    -------
    fig, axes : matplotlib Figure and Axes objects

    """
    fig, axes = plot_line_parameter(par_dict["%i" % doppler_ind], lines, "amp",
                                ncols, bins=bins)
    plt.tight_layout()
    return fig, axes


def plot_width(par_dict, lines, doppler_ind=0, ncols=4, bins=20):
    """
    Plot the marginal posterior distributions of the line width
    for all lines in the model.

    Parameters
    ----------
    par_dict : dictionary
        A dictionary of parameters as returned by the `extract_parameters()`
        function.

    lines : iterable
        The list of possible lines in the model.

    doppler_ind : int, default 0
        If the model allows multiple Doppler shifts, use this to pick the
        Doppler shift of interest

    ncols : int, default 4
        The number of columns in the subplot matrix.

    bins : {int|iterable}, default 20
        The number of bins or an array of bins for the histogram. Takes any
        valid input to `matplotlib.pyplot.hist`.

    Returns
    -------
    fig, axes : matplotlib Figure and Axes objects

    """
    fig, axes = plot_line_parameter(par_dict["%i" % doppler_ind], lines, "width",
                                    ncols, bins=bins)
    plt.tight_layout()
    return fig, axes


def plot_sign(par_dict, lines, doppler_ind=0, ncols=4):
    """
    Plot the marginal posterior distributions of the line sign (-1, 1)
    for all lines in the model.

    Parameters
    ----------
    par_dict : dictionary
        A dictionary of parameters as returned by the `extract_parameters()`
        function.

    lines : iterable
        The list of possible lines in the model.

    doppler_ind : int, default 0
        If the model allows multiple Doppler shifts, use this to pick the
        Doppler shift of interest

    ncols : int, default 4
        The number of columns in the subplot matrix.

    Returns
    -------
    fig, axes : matplotlib Figure and Axes objects

    """
    bins = [-1.5, 0.0, 1.5]
    fig, axes = plot_line_parameter(par_dict["%i" % doppler_ind], lines, "sign_num",
                                    ncols, bins, kde=False)
    plt.tight_layout()
    return fig, axes


def plot_presence(par_dict, lines, doppler_ind=0, ncols=4):
    """
    Plot the marginal posterior distributions of the line presence parameter
    for all lines in the model.

    Parameters
    ----------
    par_dict : dictionary
        A dictionary of parameters as returned by the `extract_parameters()`
        function.

    lines : iterable
        The list of possible lines in the model.

    doppler_ind : int, default 0
        If the model allows multiple Doppler shifts, use this to pick the
        Doppler shift of interest

    ncols : int, default 4
        The number of columns in the subplot matrix.

    Returns
    -------
    fig, axes : matplotlib Figure and Axes objects

    """
    bins = [-0.5, 0.5, 1.5]
    fig, axes = plot_line_parameter(par_dict["%i" % doppler_ind], lines,
                               "presence_num", ncols, bins, kde=False)
    plt.tight_layout()
    return fig, axes


def plot_unfolded_amplitude(par_dict, lines, doppler_ind=0, ncols=4, bins=20):
    """
    Plot the marginal posterior distributions of the unfolded amplitude
    (amplitude * sign) for all lines in the model.

    Parameters
    ----------
    par_dict : dictionary
        A dictionary of parameters as returned by the `extract_parameters()`
        function.

    lines : iterable
        The list of possible lines in the model.

    doppler_ind : int, default 0
        If the model allows multiple Doppler shifts, use this to pick the
        Doppler shift of interest

    ncols : int, default 4
        The number of columns in the subplot matrix.

    bins : {int|iterable}, default 20
        The number of bins or an array of bins for the histogram. Takes any
        valid input to `matplotlib.pyplot.hist`.

    Returns
    -------
    fig, axes : matplotlib Figure and Axes objects

    """
    nlines = len(lines)

    amp = par_dict["%i" % doppler_ind]["amp"]
    sign = par_dict["%i" % doppler_ind]["sign_num"]

    if nlines <= ncols:
        nrows = 1
        ncols = nlines

    else:
        nrows = int(np.ceil(nlines / 4.))

    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * nrows, 3 * ncols))

    if nlines != 1:
        axes = np.hstack(axes)

    for i, l in enumerate(lines):
        uf_amp = amp[:, i] * sign[:, i]

        if nlines == 1:
            sns.distplot(np.array(uf_amp), ax=axes, bins=bins)
            axes.set_xlabel("Unfolded amplitude")
            axes.set_title("%.3f keV" % l)
        else:
            sns.distplot(np.array(uf_amp), ax=axes[i], bins=bins)
            axes[i].set_xlabel("Unfolded amplitude")
            axes[i].set_title("%.3f keV" % l)

    plt.tight_layout()
    return fig, axes


def plot_marginal_posteriors(data, samples, lines, ndoppler=1, idx=17, ncols=4,
                             dshift=0.0, doppler_ind=0, bins=20):
    """
    Plot marginal posteriors for the Doppler shifts, line parameters and
    realizations of the model.

    Parameters
    ----------
    data : dictionary
        A dictionary with the observed spectra, as returned by `load_data`

    samples : numpy.ndarray
        The array of samples from a DNest4 run. The array has the shape
        (nsamples, npars), where nsamples is the number of samples in the
        run, and npars is the number of parameters (which also includes the
        model values if those had been printed out).

    lines : iterable
        The list of possible lines in the model.

    ndoppler : int
        The maximum number of Doppler shifts allowed by the model.

    idx : int
        The index of the column containing the first Doppler shift

    ncols : int, default 4
        The number of columns in the subplot matrix.

    dshift : float, default 0
        A Doppler shift to apply to the rest energies of the lines used when
        overplotting these positions onto the posterior realizations of the model

    doppler_ind : int, default 0
        If the model allows multiple Doppler shifts, use this to pick the
        Doppler shift of interest

    bins : {int|iterable}, default 20
        The number of bins or an array of bins for the histogram. Takes any
        valid input to `matplotlib.pyplot.hist`.

    Returns
    -------
    par_dict : dictionary
        A dictionary with the parameters. Each element in the dictionary has a
        numerical key ranging from 0 to `ndoppler`. Each element is itself
        a dictionary with the following keys:
            * `amp`: line amplitude
            * `width`: line width
            * `sign`: absorption/emission line model parameter
            * `sign_num`: binary (-1/1) absorption/emission line parameter
            * `presence`: line presence model parameter
            * `presence_num`: binary (0,1) line presence model parameter
        Each of these contains an array of shape (nsamples, nlines), so each
        column corresponds to the posterior sample of a single line.
    """

    nlines = len(lines)

    plot_doppler_shifts(samples, 1, idx=16, velocity=False)
    plot_data_model(data, samples, lines=lines, dshift=dshift)
    plot_data_model(data, samples, lines=lines, dshift=dshift, flux=True)

    par_dict = extract_parameters(samples, ndoppler=ndoppler, nlines=nlines)

    _, _ = plot_amplitude(par_dict, lines, doppler_ind, ncols, bins)
    _, _ = plot_width(par_dict, lines, doppler_ind, ncols, bins)
    _, _ = plot_sign(par_dict, lines, doppler_ind, ncols)
    _, _ = plot_presence(par_dict, lines, doppler_ind, ncols)
    _, _ = plot_unfolded_amplitude(par_dict, lines, doppler_ind, ncols, bins)

    return par_dict


def one_line_viz(datadir, model_ind, lines):
    """
    Visualize the results of the runs


    """
    filenames = [datadir + "oneline_%i_heg_p1_src.pha" % model_ind,
                 datadir + "oneline_%i_heg_m1_src.pha" % model_ind,
                 datadir + "oneline_%i_meg_p1_src.pha" % model_ind,
                 datadir + "oneline_%i_meg_m1_src.pha" % model_ind]

    names = ["HEG P1", "HEG M1", "MEG P1", "MEG M1"]

    data = load_data(filenames, names=names)

    samples = np.loadtxt(
        datadir + "oneline_%i_posterior_sample.txt" % model_ind)
    print("There are %i samples in the posterior." % len(samples))

    par_dict = plot_all(data, samples, lines)

    return data, samples, par_dict
