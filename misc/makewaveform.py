#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A Python implementation of xmakewaveform.py
#
# WARNING: waveforms have not been checked with any rigor -- due diligence required before 
# using them for quantitative results!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import sys

def makewaveform(type,params,T,T0,fs):

    #  type: string indentifying waveform type
    #
    #  params: tilde-delimited string containing parameters for the specified waveform type 
    #
    #  T: Scalar. Duration of the waveform.
    #  T0: Scalar.  Desired peak time of the waveform.  Must be less than T.
    #  fs: Scalar.  Sampling rate (in Hz). T*fs must be an integer
    #
    #  
    #  t: Times at which the waveform is sampled, starting from zero.
    #  hp: Plus polarization waveform [strain].
    #  hc: Cross polarization waveform [strain].
    #
    #  Example:
    #
    #  t, hp, hc = makewaveform('chirplet','2.5e-21~0.1~910~200~0~0',2.0,1.0,4096.0)
    #
    #  t, hp, hc = makewaveform('inspiral','1.4~1.4~0~20',2.0,1.0,4096.0)

    G = 6.673*1e-11             #-- the gravitational constant in m^3/(kg*s^2)
    c = 299792458               #-- the speed of light in m/s
    pc = 3.0856775807*1e16      #-- a parsec in m
    Mpc = pc*1e6                #-- a mega-parsec in m
    ms = 1.98892*1e30*G/c**3    #-- in geometrized unit 1 solar mass = 4.92e-6 seconds

    if T0 >= T:
        print 'Peak of waveform must be at a time smaller than the waveform duration!'
        sys.exit()

    fs = float(fs)

    t = arange(0,T-1/fs,1/fs)

    parameters = params.split('~')

    if type=='inspiral':    

        ##############################################################################
        #
        # Inspiral function taken from X-pipeline implementation of 2PN waveform, 
        # https://trac.ligo.caltech.edu/xpipeline/browser/trunk/share/inspiral2pn.m
        # See references / caveats therein.
        # 
        # Also see Chapter 2 of D Brown's thesis, arXiv:0705.1514.
        #
        # 
        # Format: makewaveform('insp2pn','mass1~mass2~iota~d',T,T0,Fs)
        #
        # mass1, mass2 are masses in the system [in units Msun]
        # 
        # iota is inclination angle [rad]
        #
        # d is distance from source to detector [Mpc]
        #
        ##############################################################################


        # ---- Required parameters.
        mass1 = float(parameters[0])
        mass2 = float(parameters[1])
        iota = float(parameters[2])
        distance = float(parameters[3])

        m1 = mass1*ms
        m2 = mass2*ms

        # total mass
        M = m1+m2

        # reduced mass
        mu = m1*m2/M

        # ---- mass difference of the system [s].
        dm = m1 - m2

        # ---- Mass ratio.
        eta = mu/M

        # ---- Chirp mass.
        Mchirp = mu**(3.0/5.0) * M**(2.0/5.0)

        # ---- Distance from the source to the detector [s].
        r = distance*Mpc/c

        # ---- Time before collapse when separation reaches that of a test particle
        #      in the innermost stable circular orbit of Schwarzschild mass M.
        tplunge = (405*M/(16*eta))
        # ---- ISCO frequency?
        fisco = 1/(6**(3/2)*pi*M)


        # We only want to calculate the waveform BEFORE the coalescence;
        # after the coalescence all these factors will diverge.
        idx = where(t<T0)
        t_insp = t[idx]

        # Rescale the time parameter into dimensionless time variable
        # See D Brown's thesis, equation 2.99
        tau = (eta/(5*M) * (T0-t_insp))**(-1.0/8.0);

        # ---- GW phase.
        Phi0 = 0  #-- phase at coalescence time

        phic0 = -1/eta * tau**(-5.0)
        phic2 = -1/eta * (3715.0/8064.0 + 55.0*eta/96.0) * tau**(-3.0)
        phic3 = -1/eta * (-3.0*pi/4.0) * tau**(-2)
        phic4 = -1/eta * (9275495.0/14450688.0 + 284875.0*eta/258048.0 + 1855.0*eta**2/2048.0) * tau**(-1)

        phi = Phi0 + phic0 + phic2 + phic3 + phic4

        omega0 = tau**3
        omega2 = (743.0/2688.0 + 11.0/32.0 * eta) * tau**5
        omega3 = -3.0*pi/10.0 * tau**6
        omega4 = (1855099.0/14450688.0 + 56975.0/258048.0 * eta + 371.0/2048.0 * eta**2) * tau**7

        omega = (1.0/8.0) * (omega0 + omega2 + omega3 + omega4)

        freq = ediff1d(phi)/(2.0*pi/fs)
        [badEvoln] = where(ediff1d(freq)<0)

        if badEvoln.size:
            stop_idx = badEvoln[0]
            t_insp = t_insp[0:stop_idx+1]
            tau = tau[0:stop_idx+1]
            omega = omega[0:stop_idx+1]
            phi = phi[0:stop_idx+1]
            freq = freq[0:stop_idx+1]

        # ---- Inspiral GW amplitude.
        x = omega**(2.0/3.0)
        amp = 2.0*mu*x/r

        # ---- Construct plus and cross waveforms.
        cc=cos(iota)
        ss=sin(iota)

        hplus = amp * (1+cc**2) * cos(2*phi)
        hcross = 2 * amp * cc * sin(2*phi)

        hp = zeros(shape(t))
        hc = zeros(shape(t))

        hp[:len(hplus)] += hplus
        hc[:len(hcross)] += hcross

        return t,hp,hc

    if type=='chirplet':

        # ---- Chirplet - Gaussian-modulated sinusoid with frequency
        #      changing linearly with time.  Put chirping cosine-Gaussian
        #      in plus polarization, chirping sine-Gaussian in cross.

        # parameters are hrss~tau~freq~alpha~delta~iota
        # tau is duration (seconds)
        # alpha is the chirp parameter (>0 chirps upward in freq)
        # delta is the phase at the peak of the envelope
        # iota is the inclination angle (of the source?)
        # Q = f*tau*pi*sqrt(8)

        # ---- Required parameters.
        h_rss = float(parameters[0])
        tau = float(parameters[1])
        f0 = float(parameters[2])

        # ---- Optional parameters.
        alpha = 0
        delta = 0
        iota = 0
        if len(parameters)>=4:
            alpha = float(parameters[3])
        elif len(parameters)>=5:
            delta = float(parameters[4])
        elif len(parameters)>=6:
            iota = float(parameters[5])
        
        # ---- Waveform.
        
        sinusoid = exp(2j*pi*f0*(t-T0))
        envelope = exp(-1*(t-T0)**2/(4*tau**2))
        chirp = exp((1j*alpha)*(t-T0)**2 / (4*tau**2))
        phase = exp(1j*delta)
        normalization = 1/(2*pi*tau**2)**(1/4.0)

        h = h_rss * sinusoid * envelope * chirp * phase * normalization

        hp = 1/2.0*(1+(cos(iota))**2) * real(h)
        hc = cos(iota) * imag(h)

        return t,hp,hc
