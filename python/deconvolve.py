
import sherpa.astro.ui as ui
import numpy as np


def get_spectrum(filename, data_id="1"):
    ui.load_data(id=data_id, filename=filename)
    d = ui.get_data(data_id)
    arf = d.get_arf()
    rmf = d.get_rmf()

    return d, arf, rmf

def deconvolve(spec, arf, rmf):
    """
    Attempt a rough deconvolution of the specrum by applying the 
    ARF and RMF to a simple array filled with 1s, and then dividing 
    out the result.

    Parameters
    ----------
    spec : numpy.ndarray
        The array with the spectrum to be deconvolved. Can be 
        a data or model array.

    arf : a sherpa.astro.data.DataARF object
        The object with the ARF

    rmf : a sherpa.astro.data.DataRMF object
        The object with the RMF


    Returns
    -------
    ratio : numpy.ndarray
         The deconvolved spectrum

    """

    # make an array of ones that we'll convolve
    # with the responses
    m = np.ones_like(spec)

    # multiply array with ARF
    m_arf = arf.apply_arf(m)
    # convolve output of previous operaton with RMF
    m_rmf = rmf.apply_rmf(m_arf)
    
    # divide out the effect of the responses
    ratio = spec/m_rmf

    return ratio
    
