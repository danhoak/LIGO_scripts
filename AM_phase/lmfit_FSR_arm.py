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
    T_i = params['T_i'].value
    T_e = params['T_e'].value
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

    r_i = sqrt(1-T_i)
    r_e = sqrt(1-T_e)
    R = r_i * r_e
    g_refl = -1*(1-r_i**2) / (1 - R)**2

    # The RF parts
    S_Irf = sin(phi_f/2)*(1 + 2*R - (1+2*cos(phi_f))*(R**2 + r_e**2) + R*r_e**2 * (2+R))
    S_Qrf = cos(phi_f/2)*(1 - 2*R - (1-2*cos(phi_f))*(R**2 - r_e**2) + R*r_e**2 * (2-R))

    S_refl = p0 * g_refl**2 * (-2*R / (1 - 2*R*cos(phi_f) + R**2)**2) * sin(phi_f/2) * ( S_Irf + 1j*S_Qrf )

    #S_refl = p0 * g_refl**2 * (-2*R / (1 - 2*R*cos(phi_f) + R**2)**2) * sin(phi_f/2) * ( S_Irf + 1j*S_Qrf )

    audioI = sin(phi_a)/phi_a * ( (exp(2j*phi_f) + R*r_e**2)*(1-R) - exp(1j*phi_f)*((1+R)*(R-r_e**2)-2*cos(phi_a)*(R**2-r_e**2)) )
    audioQ = (cos(phi_a)-1)/phi_a * ( (exp(2j*phi_f) + R*r_e**2)*(1+R) + exp(1j*phi_f)*((1-R)*(R-r_e**2)-2*cos(phi_a)*(R**2+r_e**2)) )

    S_pre = exp(1j*phi_a)*(1-R)*(exp(1j*phi_f) - R) / ((exp(1j*phi_f + 1j*phi_a) - R)*(exp(1j*phi_f) - R*exp(1j*phi_a))*(1-2*R*cos(phi_a)+R**2)*(exp(1j*phi_f) - r_e**2))

    S_refl_Iaf = S_refl * S_pre * audioI
    S_refl_Qaf = S_refl * S_pre * audioQ

    # the real part is the RF in-phase signal
    return real(S_refl_Iaf), real(S_refl_Qaf), imag(S_refl_Iaf), imag(S_refl_Qaf)


## Error function we seek to minimize
def e2(params,f_RF,rII,rQI,rQQ,rIQ):

    VII, VIQ, VQI, VQQ = Irf(params,f_RF)

    FSR = params['FSR'].value
    f0 = params['f0'].value
    Q = params['Q'].value

    if f0==0 and Q==0:    
        theta = params['theta0'].value + params['theta1'].value*(f_RF - FSR)
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
            'savefig.dpi': 200,
            'xtick.major.size':8,
            'xtick.minor.size':4,
            'ytick.major.size':8,
            'ytick.minor.size':4
            })




c = 299792458.0
lamd = 532.0e-9

audio_freq = 303.0
data = genfromtxt('data/alsY_303_17Jan.txt')
print shape(data)

#### Begin fit of cavity length

# First pick out narrow region around central zero crossing to fit FSR
#idx1 = 940
#idx2 = 990

FSR_freq = data[:,1]/666.0

#FSR_freq = arange(-600,600,0.1)

# sign flip on RF-Q phase signals -- due to cosine/sine convention in demod board?
rII = data[:,2]
rIQ = data[:,3]
rQI = data[:,4]
rQQ = data[:,5]

# starting parameters for fit
params_FSR = Parameters()

# green trans params
T_etm = 0.31
T_itm = 0.01

params_FSR.add('audio_freq', value=audio_freq, vary=False)
params_FSR.add('FSR', value=37.5259e3, vary=True)
params_FSR.add('T_e', value=0.01, vary=False)
params_FSR.add('T_i', value=0.31, vary=False)
params_FSR.add('A0',value=0.2, vary=True)
params_FSR.add('f0',value=0.0, vary=False)
params_FSR.add('theta0',value=0.0, vary=True)
params_FSR.add('theta1',value=0.0, vary=True)
params_FSR.add('Q',value=0.0, vary=False, min=0, max=200)
params_FSR.add('phi0',value=0.0, vary=True)

# First use the Nelder-Mead method to get close to final solution (it is more robust to local minima)...
out = minimize(e2, params_FSR, args=(FSR_freq,rII,rQI,rQQ,rIQ), method='nelder')

# ...then finish with Levenberg-Marquardt
out = minimize(e2, params_FSR, args=(FSR_freq,rII,rQI,rQQ,rIQ), method='leastsq')

FII, FIQ, FQI, FQQ = Irf(params_FSR,FSR_freq)


FSR = params_FSR['FSR'].value
theta = params_FSR['theta0'].value + params_FSR['theta1'].value*(FSR_freq-FSR)
phiAF = params_FSR['phi0'].value

cavity_length = c / (2*FSR)
cavity_length_error = cavity_length * (params_FSR['FSR'].stderr / params_FSR['FSR'].value)

#r = sqrt(1-params_FSR['T'].value)
#finesse = pi*sqrt(r**2)/(1-r**2)
#cavity_pole = FSR/(2*finesse)

fII = cos(phiAF)*(cos(theta)*rII - sin(theta)*rQI) - sin(phiAF)*(cos(theta)*rIQ - sin(theta)*rQQ)
fQI = cos(phiAF)*(sin(theta)*rII + cos(theta)*rQI) - sin(phiAF)*(sin(theta)*rIQ + cos(theta)*rQQ)
fIQ = sin(phiAF)*(cos(theta)*rII - sin(theta)*rQI) + cos(phiAF)*(cos(theta)*rIQ - sin(theta)*rQQ)
fQQ = sin(phiAF)*(sin(theta)*rII + cos(theta)*rQI) + cos(phiAF)*(sin(theta)*rIQ + cos(theta)*rQQ)

print 
print 'Results of Cavity Length Fit'
print
print 'FSR [Hz] =', FSR
print FSR*666, params_FSR['FSR'].stderr*666
print
print 'Cavity length [m] =', cavity_length, '+\-', cavity_length_error
#print
#print 'Finesse =', finesse
#print 'Cavity pole [Hz] =', cavity_pole
print
print fit_report(params_FSR)
print

chisq = sum(((FQI-fQI)**2 + (FII-fII)**2 + (FIQ-fIQ)**2 + (FQQ-fQQ)**2))
#print 'chisq =', chisq
reduced_chisq = chisq/(4*len(FSR_freq)-7-1)
print 'Reduced chi sq is', reduced_chisq




fignum=0

from matplotlib.font_manager import FontProperties

fontP = FontProperties()
fontP.set_size('medium')

gs = gridspec.GridSpec(2,1,height_ratios=[3,1])

fignum=fignum+1
pylab.figure(fignum)

pylab.subplot(gs[0])

phiF = 2*pi*(FSR_freq/FSR - 1)

"""
#pylab.plot(FSR_freq-FSR,fII,'c+',markersize=5)
pylab.plot(FSR_freq-FSR,fIQ,'y+',markersize=5)
#pylab.plot(FSR_freq-FSR,fQI,'r+',markersize=5)
pylab.plot(FSR_freq-FSR,fQQ,'g+',markersize=5)

#pylab.plot(FSR_freq-FSR,FII,'b--',linewidth=1)
pylab.plot(FSR_freq-FSR,FIQ,'b--',linewidth=1)
#pylab.plot(FSR_freq-FSR,FQI,'r--',linewidth=1)
pylab.plot(FSR_freq-FSR,FQQ,'r--',linewidth=1)
"""

pylab.plot(phiF,fII,'c+',markersize=5)
pylab.plot(phiF,fIQ,'y+',markersize=5)
pylab.plot(phiF,fQI,'r+',markersize=5)
pylab.plot(phiF,fQQ,'g+',markersize=5)

pylab.plot(phiF,FII,'k--',linewidth=1)
pylab.plot(phiF,FIQ,'k--',linewidth=1)
pylab.plot(phiF,FQI,'k--',linewidth=1)
pylab.plot(phiF,FQQ,'k--',linewidth=1)

#pylab.xlim(9.0988, 9.1012)
#pylab.ylim(-0.04, 0.07)

pylab.grid(True)
pylab.ylabel('Signal [arb]')
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)

pylab.title('Cavity Length Fit - Data: solid; Fit: dashed, FSR = ' + str(FSR) + ' Hz', fontsize=10)

#pylab.legend(('RF I AF I','RF I AF Q','RF Q AF I','RF Q AF Q'),fancybox=True,fontsize=11,bbox_to_anchor=[1.07, 1.07])

pylab.subplot(gs[1])


pylab.plot(phiF,fII-FII,'c+',markersize=5)
pylab.plot(phiF,fIQ-FIQ,'y+',markersize=5)
pylab.plot(phiF,fQI-FQI,'r+',markersize=5)
pylab.plot(phiF,fQQ-FQQ,'g+',markersize=5)

#pylab.plot(FSR_freq-FSR,fII-FII,'c+',markersize=5)
#pylab.plot(FSR_freq-FSR,fIQ-FIQ,'y+',markersize=5)
#pylab.plot(FSR_freq-FSR,fQI-FQI,'r+',markersize=5)
#pylab.plot(FSR_freq-FSR,fQQ-FQQ,'g+',markersize=5)

#pylab.xlim(9.0988, 9.1012)

pylab.grid(True)
#pylab.xlabel('Frequency [MHz]')
#pylab.xlabel('Frequency [Hz]')
pylab.xlabel(r'$\phi_f$ (rad)')
pylab.ylabel('Residuals')

pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)

pylab.savefig('plots/global_analytic_arm.pdf')
