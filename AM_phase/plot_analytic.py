#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A script to fit the data collected for the cavity width measurement
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
import matplotlib.gridspec as gridspec
#from scipy.optimize import leastsq
from lmfit import minimize, Parameters, fit_report

## Parametric function: 'v' is the vector of parameters, 'x' the independent variable
def Irf(params,f_RF):

    f_FSR = params['FSR'].value
    T = params['T'].value
    A0 = params['A0'].value
    f0 = params['f0'].value
    Q = params['Q'].value

    AF = params['audio_freq'].value

    if f0==0 and Q==0:
        p0 = 2*pi * A0  # this expression includes the constant terms df, Gamma, P_in
    else:
        p0 = 2*pi * A0 * sqrt( (f0/(2*Q))**2 / ( (f_RF-f0)**2 + (f0/(2*Q))**2 ) )

    phi_a = 2*pi*AF / f_FSR
    phi_f = 2*pi*f_RF / f_FSR

    r = sqrt(1-T)
    R = r**2
    g_tr = 1   # assuming a critically coupled cavity

    # The RF parts - original
    S_Irf = sin(phi_f/2)*(1 + 2*R - (1+2*cos(phi_f))*R**2)
    S_Qrf = cos(phi_f/2)*(1 - 2*R - (1-2*cos(phi_f))*R**2)

    S_tr = p0 * g_tr**2 * (-2*R / (1 - 2*R*cos(phi_f) + R**2)**2) * sin(phi_f/2) * ( S_Irf + 1j*S_Qrf )

    audioI = ( (exp(2j*phi_a)-1)    / (2j*phi_a) ) * ( (1-R)*exp(1j*phi_f) - R*(1-2*R*cos(phi_a) + R) )
    audioQ = ( (exp(1j*phi_a)-1)**2 / (2*phi_a)  ) * ( (1+R)*exp(1j*phi_f) + R*(1-2*R*cos(phi_a) - R) )

    S_tr_Iaf = S_tr * ( (1-R)*(exp(1j*phi_f) - R) / ( (exp(1j*phi_f + 1j*phi_a) - R)*(exp(1j*phi_f) - R*exp(1j*phi_a)) * (1-2*R*cos(phi_a)+R**2) ) ) * audioI
    S_tr_Qaf = S_tr * ( (1-R)*(exp(1j*phi_f) - R) / ( (exp(1j*phi_f + 1j*phi_a) - R)*(exp(1j*phi_f) - R*exp(1j*phi_a)) * (1-2*R*cos(phi_a)+R**2) ) ) * audioQ

    # the real part is the RF in-phase signal
    return real(S_tr_Iaf), real(S_tr_Qaf), imag(S_tr_Iaf), imag(S_tr_Qaf)


## Error function we seek to minimize
def e2(params,f_RF,rII,rQI,rQQ,rIQ):

    VII, VIQ, VQI, VQQ = Irf(params,f_RF)

    FSR = params['FSR'].value
    f0 = params['f0'].value
    Q = params['Q'].value

    if f0==0 and Q==0:    
        theta = params['theta0'].value
    else:
        theta = params['theta0'].value + arctan((f_RF-f0) / (f0/(2*Q)) )

    phiAF = params['phi0'].value

    # rotate audio and RF I & Q

    yII = cos(phiAF)*(cos(theta)*rII - sin(theta)*rQI) - sin(phiAF)*(cos(theta)*rIQ - sin(theta)*rQQ)
    yQI = cos(phiAF)*(sin(theta)*rII + cos(theta)*rQI) - sin(phiAF)*(sin(theta)*rIQ + cos(theta)*rQQ)
    yIQ = sin(phiAF)*(cos(theta)*rII - sin(theta)*rQI) + cos(phiAF)*(cos(theta)*rIQ - sin(theta)*rQQ)
    yQQ = sin(phiAF)*(sin(theta)*rII + cos(theta)*rQI) + cos(phiAF)*(sin(theta)*rIQ + cos(theta)*rQQ)

    return concatenate([(VQI-yQI),(VII-yII),(VQQ-yQQ),(VIQ-yIQ)])




fontP = FontProperties()
fontP.set_size('medium')

fig_width_pt = 600  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
#golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio, = 0.618
golden_mean = 0.8         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size = [fig_width,fig_height]
matplotlib.rcParams.update(
        {'axes.labelsize': 16, 
            'font.size':   16, 
            'legend.fontsize': 16, 
            'xtick.labelsize': 16, 
            'ytick.labelsize': 16, 
            'text.usetex': False,
            'figure.figsize': fig_size,
            #'font.family': "serif",
            #'font.serif': ["Times New Roman"],
            'font.family': "sans-serif",
            #'font.serif': ["Times New Roman"],
            'savefig.dpi': 250,
            'xtick.major.size':8,
            'xtick.minor.size':4,
            'ytick.major.size':8,
            'ytick.minor.size':4
            })




c = 299792458.0
lamd = 1064.0e-9

audio_freq = 110.0

freq = 6*linspace(1.00e+6, 1.1e+6, 300000)
L = 143.4259
T = 2500e-6
n = 6

#freq = linspace(8.95e+6, 9.25e+6, 10240)
#L = 16.4717
#T = 6000e-6
#n = 1

# starting parameters for fit
params_FSR = Parameters()

params_FSR.add('audio_freq', value=audio_freq, vary=False)
params_FSR.add('FSR', value=c/(2*L), vary=False)
params_FSR.add('T', value=T, vary=False)
params_FSR.add('A0',value=1, vary=False)
params_FSR.add('f0',value=0.0, vary=False)
params_FSR.add('theta0',value=0.0, vary=False)
params_FSR.add('Q',value=0.0, vary=False)
params_FSR.add('phi0',value=0.0, vary=False)


VII, VIQ, VQI, VQQ = Irf(params_FSR,freq)

FSR = params_FSR['FSR']
R = 1 - params_FSR['T']
finesse = pi*sqrt(R) / (1- R)
f_pole = FSR / (2*finesse)
f_p = f_pole / sqrt(R) + n*FSR
f_m = -1*f_pole / sqrt(R) + n*FSR

fignum=0

from matplotlib.font_manager import FontProperties

fontP = FontProperties()
fontP.set_size('medium')

fignum=fignum+1
pylab.figure(fignum)

pylab.plot(freq/1e6,VII,'b-',linewidth=1.2)
pylab.plot(freq/1e6,VIQ,'y-',linewidth=0.8)
pylab.plot(freq/1e6,VQI,'r-',linewidth=1.2)
pylab.plot(freq/1e6,VQQ,'g-',linewidth=0.8)

#pylab.plot(f_p/1e6,0,'k*',markersize=10)
#pylab.plot(f_m/1e6,0,'k*',markersize=10)

print n*FSR, f_m, n*FSR-f_m

#pylab.xlim(9.08, 9.118)
pylab.xlim(6.2695, 6.2718)

pylab.ylim(-1200,800)

pylab.grid(True)
pylab.ylabel(r'$S_{tr}$ (arb)', fontsize=14)
pylab.xlabel('Frequency (MHz)', fontsize=14)
pylab.xticks(fontsize=12)
pylab.yticks(fontsize=12)

#pylab.savefig('plots/analytic_test_ligo.png',bbox_inches='tight')
pylab.savefig('plots/analytic_test_virgo_lof.png',bbox_inches='tight')

pylab.close()
