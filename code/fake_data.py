import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import cPickle as pickle
import scipy.special


def gaussian_cdf(x, w0, width):
    return 0.5*(1. + scipy.special.erf((x-w0)/(width*np.sqrt(2.))))

def spectral_line(wleft, wright, w0, amplitude, width):
    """
    Use the CDF of a Gaussian distribution to define spectral 
    lines. We use the CDF to integrate over the energy bins, 
    rather than taking the mid-bin energy.
    
    Parameters
    ----------
    wleft: array
        Left edges of the energy bins
        
    wright: array
        Right edges of the energy bins
        
    w0: float
        The centroid of the line
        
    amplitude: float
        The amplitude of the line
        
    width: float
        The width of the line
        
    Returns
    -------
    line_flux: array
        The array of line fluxes integrated over each bin
    """
    line_flux = amplitude*(gaussian_cdf(wright, w0, width)-
                           gaussian_cdf(wleft, w0, width))
    return line_flux

def fake_spectrum(wleft, wright, line_pos, logbkg=np.log(0.09), err=0.007, dshift=0.0,
                  sample_logamp=False, sample_logq=False, sample_signs=False,
                  logamp_hypermean=None, logamp_hypersigma=np.log(0.08), nzero=0, 
                  logq_hypermean=np.log(500), logq_hypersigma=np.log(50)):
    """
    Make a fake spectrum with emission/absorption lines.
    The background is constant, though that should later become an OU process or 
    something similar.
    
    NOTE: The amplitude *must not* fall below zero! I'm not entirely sure how to deal 
    with that yet!
    
    Parameters
    ----------
    wleft: np.ndarray
        Left edges of the energy bins
        
    wright: np.ndarray
        Right edges of the energy bins
    
    line_pos: np.ndarray
        The positions of the line centroids
    
    bkg: float
        The value of the constant background
        
    err: float
        The width of the Gaussian error distribution
        
    dshift: float, default 0.0
        The Doppler shift of the spectral lines.
        
    sample_amp: bool, default False
        Sample all amplitudes? If not, whatever value is set in 
        `amp_hypermean` will be set as collective amplitude for all 
        lines
        
    sample_width: bool, default False
        Sample all line widths?  If not, whatever value is set in 
        `width_hypersigma` will be set as collective amplitude for all 
        lines
    
    sample_signs: bool, default False
        Sample the sign of the line amplitude (i.e. whether the line is an 
        absorption or emission line)? If False, all lines will be absorption 
        lines
        
    logamp_hypermean: {float | None}, default None
        The mean of the Gaussian prior distribution on the line amplitude. If None, 
        it is set to the same value as `bkg`.
        
    logamp_hypersigma: float, default 0.08
        The width of the Gaussian prior distribution on the line amplitude. 
        
    nzero: int, default 0
        The number of lines to set to zero amplitude
        
    logq_hypermean: float, default 0.01
        The mean of the Gaussian prior distribution on the 
        q-factor, q=(line centroid wavelength)/(line width)
        
    logq_hypersigma: float, default 0.01
        The width of the Gaussian prior distribution on the
        q-factor, q=(line centroid wavelength)/(line width)
    
    Returns
    -------
    model_flux: np.ndarray
        The array of model line fluxes for each bin
        
    fake_flux: np.ndarray
        The array of fake fluxes (with errors) for each bin
        
    """
    
    # number of lines
    nlines = line_pos.shape[0]
    
    # shift spectral lines
    line_pos_shifted = line_pos*(1. + dshift)
   
    # if sampling the amplitudes
    if sample_logamp:  
        # sample line amplitudes
        logamps = np.random.normal(logamp_hypermean, logamp_hypersigma, size=nlines)
    else:
        logamps = np.zeros(nlines)+logamp_hypermean
        
    
    amps = np.exp(logamps)
    
    if nzero > 0:
        zero_ind = np.random.choice(np.arange(nlines), size=nzero)
        for z in zero_ind:
            amps[int(z)] = 0.0
    
    if sample_signs:
        # sample sign of the amplitudes
        signs = np.random.choice([-1., 1.], p=[0.5, 0.5], size=nlines)
    else:
        # all lines are absorption lines
        signs = -1.*np.ones(nlines)
        
    # include signs in the amplitudes
    amps *= signs
    
    if sample_logq:
        # widths of the lines
        logq = np.random.normal(logq_hypermean, logq_hypersigma, size=nlines)
    else:
        logq = np.ones(nlines)*logq_hypermean
        
        
    widths = line_pos_shifted/np.exp(logq)
    
    model_flux = np.zeros_like(wleft) + np.exp(logbkg)
    
    for si, a, w in zip(line_pos_shifted, amps, widths):
        model_flux += spectral_line(wleft, wright, si, a, w)
      
    fake_flux = model_flux + np.random.normal(0.0, 0.007, size=model_flux.shape[0])
    
    pars = {"wavelength_left": wleft, "wavelength_right": wright, "err":err,
            "model_flux": model_flux, "fake_flux": fake_flux, "logbkg":logbkg,
            "dshift": dshift, "line_pos": line_pos_shifted, "logamp": logamps,
           "signs": signs, "logq": logq }
    
    
    return pars
    
    
