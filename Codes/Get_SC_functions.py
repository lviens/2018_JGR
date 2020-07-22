"""
L. Viens 02/29/2018
This function computes Single-station Cross-correlation (SC) functions from the data recorded by the vertical and horizontal components at one seismic station. See Viens et al. (2018, JGR) for more details.
The SC functions are computed in the frequency domain.
"""

import numpy as np
import obspy.signal as obssig

def Comupte_SC_functions(dat_Z, dat_HZ, delta, len_corr):
    """
    Function to compute Single-station Cross-correlation functions (SCFs).
    INPUT: 
        - dat_Z  = data recoded by the vertical component in time domain
        - dat_HZ = data recoded by one horizontal component (either NS or EW) in time domain
        - delta  = Sampling rate (in Hz)
        - len_corr = duration of the SC function in second (as the causal and acausal parts of the SC function are computed, the total duration in 2 times len_corr)
    OUTPUT:
        - corr = SC function in time domain.
    """
    n = len(dat_Z) 
    dt = 1/delta
    
    # 1-bit normalization
    dat_Z = np.sign(dat_Z)
    dat_HZ = np.sign(dat_HZ)

    # Compute FFT
    fft_Z = np.fft.fft(dat_Z, n*4)
    fft_HZ = np.fft.fft(dat_HZ, n*4)
    sj2 = obssig.util.smooth(np.absolute(fft_Z), 20)

    # Compute SC functions in Fourier domain and IFFT to return in the time domain
    cc_t = np.real(np.fft.ifft( (fft_HZ * np.conj(fft_Z))/ (sj2**2) ))  
    # Rearrange causal and acausal parts
    corr_tot = np.concatenate((cc_t[len(cc_t)/2:], cc_t[:len(cc_t)/2]))
    # Select the data between -len_corr and len_corr second around the zero time-lag
    corr  = corr_tot[len(corr_tot)/2-len_corr*dt:len(corr_tot)/2+len_corr*dt]

    return corr
