#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This file defines a python implementation of X-Pipeline's xoptimalsnr
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys
import numpy as np
from scipy.interpolate import interp1d
import pylab

def optimalSNR(h,Fs,t0,S,f,Fmin,Fmax):

    # OPTIMALSNR is a python implementation of xoptimalsnr.m
    #
    # https://trac.ligo.caltech.edu/xpipeline/browser/trunk/share/xoptimalsnr.m
    #
    #
    #   h     is an array of timeseries data.  For now assume single column.
    #
    #   Fs    is the sample rate, in Hz, of h(t)
    #
    #   t0    is the start time of the timeseries h
    #
    #   S     is the **one-sided power** spectrum of the noise
    #
    #   f     is the frequencies corresponding to the power spectral densities in S
    #
    #   Fmin  is the lower limit for the frequency-domain calculations.
    #
    #   Fmax  is the upper limit for the frequency-domain calculations.
    #
    # S will be interpolated to match the frequencies of FFT[h(t)].
    #
    # Note: for white Gaussian noise with variance s^2, sampled at Fs, the corresponding 
    # one-sided PSD is a constant value (2 s^2 / Fs).
    # 
    # E.g., a noise time-series n = s * numpy.random.randn(Fs*T) will have one-sided PSD of
    # S(f) = (2 * s**2 / Fs) for f = numpy.arange(0,Fs/2,1/T).
    #
    #
    # returns [SNR, h_rss, h_peak, Fchar, bw, Tchar, dur]
    #
    # Note that h_rss is the square root of the total signal power.

    N = len(h)
    T = N/Fs
    dF = 1/T

    t = np.arange(0,T,1/Fs)

    h_rss = np.sqrt(np.sum(h**2)/Fs)

    # Calculate FFT and normalize by sampling rate to get power density per frequency bin

    hf = np.fft.fft(h)/Fs
    
    F = np.fft.fftfreq(len(h),1/Fs)

    if F[N/2-1] < Fmax:
        print 'Warning: Nyquist frequency of DFT is below requested Fmax!'
        print 'Calculating SNR in frequency band', Fmin, 'to', Fs/2
        Fmax = F[N/2-1]

    # hf will be in usual screwy FFT order: [0 1 ... N/2 -N/2 ... -2 -1]
    # find the indices of the samples closest to Fmin, Fmax (positive frequencies only!)

    idx_min = abs(F[0:N/2-1]-Fmin).argmin()
    idx_max = abs(F[0:N/2-1]-Fmax).argmin()

    freq = F[idx_min:idx_max]

    # Generate the one-sided power spectrum in the interval [Fmin,Fmax]
    # Only keep positive frequency information, include a factor of two to account for negative frequencies
    # Strip off imaginary part (which is zero after multiplying by the conjugate)
    #
    # Note we have to include a factor of two here (even though there is another factor of two
    # in the integral over frequency below!) because we have assumed that the noise is provided in 
    # a one-sided PSD.  So we need an extra factor of two here to balance the factor of two that is
    # assumed to be included in the noise PSD.

    Pf = 2*np.real(hf[idx_min:idx_max]*np.conj(hf[idx_min:idx_max]))

    # If we are including 0Hz and the Nyquist, count those bins half
    if freq[0]==0.0:
        Pf[0] = 0.5*Pf[0]
    if freq[-1]==Fs/2:
        Pf[-1] = 0.5*Pf[-1]

    # We cannot interpolate the noise beyond where it is defined
    
    if Fmin < f[0] or Fmax > f[-1]:
        print f[0], Fmin
        print f[-1], Fmax
        print 'Noise spectrum does not cover desired frequency range.'
        sys.exit()

    noise_interp = interp1d(f,S,kind='linear')
    Sf = noise_interp(freq)

    # SNR^2 vs frequency on interval [Fmin:Fmax]
    SNR2f = Pf/Sf

    # Integrate over frequency to get SNR, factor of two for negative frequencies
    SNR = np.sqrt(2*dF*sum(SNR2f))

    Fchar = np.sum(freq*SNR2f)/np.sum(SNR2f)

    bw = np.sqrt( np.sum( (SNR2f*(freq-Fchar)**2) / np.sum(SNR2f) ))

    # Power in the time domain

    Pt = h**2

    h_peak = np.sqrt(max(Pt))

    Tchar = t0 + np.sum(t*Pt) / np.sum(Pt)

    dur = np.sqrt( np.sum( Pt*(t-Tchar)**2) / np.sum(Pt) )

    return SNR, h_rss, h_peak, Fchar, bw, Tchar, dur
