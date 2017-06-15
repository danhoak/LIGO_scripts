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
            'savefig.dpi': 200,
            'xtick.major.size':8,
            'xtick.minor.size':4,
            'ytick.major.size':8,
            'ytick.minor.size':4
            })




c = 299792458.0
lamd = 1064.0e-9

audio_freq = 1000.0
data = genfromtxt('data/REFL9_1000_18Dec.txt')
print shape(data)


#### Begin fit of cavity length

# First pick out narrow region around central zero crossing to fit FSR
#idx1 = 940
#idx2 = 990

idx1 = 470
idx2 = 530

FSR_freq = data[idx1:idx2,1]

# sign flip on RF-Q phase signals -- due to cosine/sine convention in demod board?
rII = data[idx1:idx2,2]
rIQ = data[idx1:idx2,3]
rQI = -1*data[idx1:idx2,4]
rQQ = -1*data[idx1:idx2,5]

# starting parameters for fit
params_FSR = Parameters()

params_FSR.add('audio_freq', value=audio_freq, vary=False)
params_FSR.add('FSR', value=9.1e6, vary=True)
params_FSR.add('T', value=6000e-6, vary=False)
params_FSR.add('A0',value=0.0002, vary=True)
params_FSR.add('f0',value=0.0, vary=False)
params_FSR.add('theta0',value=0.0, vary=False)
params_FSR.add('Q',value=0.0, vary=False)
params_FSR.add('phi0',value=0.0, vary=True)

# First use the Nelder-Mead method to get close to final solution (it is more robust to local minima)...
out = minimize(e2, params_FSR, args=(FSR_freq,rII,rQI,rQQ,rIQ), method='nelder')

# ...then finish with Levenberg-Marquardt
out = minimize(e2, params_FSR, args=(FSR_freq,rII,rQI,rQQ,rIQ), method='leastsq')

FII, FIQ, FQI, FQQ = Irf(params_FSR,FSR_freq)

FSR = params_FSR['FSR'].value
theta = 0.0
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




#### Begin fit for cavity width


# Pick a subset of the data to fit - use this to ignore edge effects and noise
#idx1 = 0
#idx2 = 1960

idx1 = 0
idx2 = 935


freq = data[idx1:idx2,1]

# sign flip on RF-Q phase signals -- due to cosine/sine convention in demod board?
rII = data[idx1:idx2,2]
rIQ = data[idx1:idx2,3]
rQI = -1*data[idx1:idx2,4]
rQQ = -1*data[idx1:idx2,5]

"""
idxa = 330
idxb = 440
freq = hstack((data[idx1:idxa,1],data[idxb:idx2,1]))

# sign flip on RF-Q phase signals -- due to cosine/sine convention in demod board?
rII = hstack((data[idx1:idxa,2],data[idxb:idx2,2]))
rIQ = hstack((data[idx1:idxa,3],data[idxb:idx2,3]))
rQI = -1*hstack((data[idx1:idxa,4],data[idxb:idx2,4]))
rQQ = -1*hstack((data[idx1:idxa,5],data[idxb:idx2,5]))
"""

# starting parameters for fit
# Use FSR from cavity length fit
params = Parameters()

params.add('audio_freq', value=audio_freq, vary=False)
params.add('FSR', value=params_FSR['FSR'].value, vary=False)
params.add('T', value=params_FSR['T'].value, vary=True)
params.add('A0',value=params_FSR['A0'].value, vary=True)
params.add('f0',value=9.11e6,vary=True, min=9e6,max=12.0e6)
params.add('theta0',value=0.0, vary=True, min=0.0, max=2*pi)
params.add('Q',value=61.0,vary=True,min=10, max=100)
params.add('phi0',value=params_FSR['phi0'].value, vary=True, min=0.0, max=2*pi)
#params.add('phi0',value=0.78, vary=False, min=0.0, max=2*pi)

# First use the Nelder-Mead method to get close to final solution (it is more robust to local minima)...
out = minimize(e2, params, args=(freq,rII,rQI,rQQ,rIQ), method='nelder')

# ...then finish with Levenberg-Marquardt
out = minimize(e2, params, args=(freq,rII,rQI,rQQ,rIQ), method='leastsq')

VII, VIQ, VQI, VQQ = Irf(params,freq)

FSR = params['FSR'].value
f0 = params['f0'].value
Q = params['Q'].value
Q_err = params['Q'].stderr

if f0==0 and Q==0:    
    theta = params['theta0'].value
else:
    theta = params['theta0'].value + arctan((freq-f0) / (f0/(2*Q)) )

if f0==0 and Q==0:    
    theta_up = params['theta0'].value
else:
    theta_up = params['theta0'].value + arctan((freq-f0) / (f0/(2*(Q+Q_err))) )

if f0==0 and Q==0:    
    theta_down = params['theta0'].value
else:
    theta_down = params['theta0'].value + arctan((freq-f0) / (f0/(2*(Q-Q_err))) )



phiAF = params['phi0'].value

r = sqrt(1-params['T'].value)
finesse = pi*sqrt(r**2)/(1-r**2)
cavity_pole = FSR/(2*finesse)

cavity_pole_error = cavity_pole * sqrt( (params_FSR['FSR'].stderr/FSR)**2 + (2*pi*params['T'].stderr/params['T'].value)**2 )

#print min(theta)*180/pi
#print max(theta)*180/pi
theta_shift = (max(theta) - min(theta))*180/pi
theta_shift_up = (max(theta_up) - min(theta_up))*180/pi
theta_shift_down = (max(theta_down) - min(theta_down))*180/pi

yII = cos(phiAF)*(cos(theta)*rII - sin(theta)*rQI) - sin(phiAF)*(cos(theta)*rIQ - sin(theta)*rQQ)
yQI = cos(phiAF)*(sin(theta)*rII + cos(theta)*rQI) - sin(phiAF)*(sin(theta)*rIQ + cos(theta)*rQQ)
yIQ = sin(phiAF)*(cos(theta)*rII - sin(theta)*rQI) + cos(phiAF)*(cos(theta)*rIQ - sin(theta)*rQQ)
yQQ = sin(phiAF)*(sin(theta)*rII + cos(theta)*rQI) + cos(phiAF)*(sin(theta)*rIQ + cos(theta)*rQQ)

print 
print 'Results from Cavity Width Fit'
print
print 'FSR [Hz] =', FSR
print '--> Change in RF phase across frequency range [deg] =', theta_shift, theta_shift_up, theta_shift_down
print
print 'Audio phase shift [deg] =', phiAF*180/pi
print
print 'Finesse =', finesse
print 'Cavity pole [Hz] =', cavity_pole, '+\-', cavity_pole_error
print
print 'EOM central frequency', f0
print
print fit_report(params)
print


chisq = sum(((VQI-yQI)**2 + (VII-yII)**2 + (VIQ-yIQ)**2 + (VQQ-yQQ)**2))
#print 'chisq =', chisq
reduced_chisq = chisq/(4*len(freq)-7-1)
print 'Reduced chi sq is', reduced_chisq




# the reduced chisq is our estimate of sigma
#print '1-sigma uncertainty is', sqrt(reduced_chisq)
#print


fignum=0

from matplotlib.font_manager import FontProperties

fontP = FontProperties()
fontP.set_size('medium')

gs = gridspec.GridSpec(2,1,height_ratios=[3,1])

fignum=fignum+1
pylab.figure(fignum)

ax1 = pylab.subplot(gs[0])

mm = 5

phiF = 2*pi*(freq/FSR - 1)

"""
pylab.plot(freq/1e6,yII,'c+-',markersize=mm)
pylab.plot(freq/1e6,yIQ,'y+-',markersize=mm)
pylab.plot(freq/1e6,yQI,'r+-',markersize=mm)
pylab.plot(freq/1e6,yQQ,'g+-',markersize=mm)

pylab.plot(freq/1e6,VII,'k--',linewidth=1.8)
pylab.plot(freq/1e6,VIQ,'k--',linewidth=1.8)
pylab.plot(freq/1e6,VQI,'k--',linewidth=1.8)
pylab.plot(freq/1e6,VQQ,'k--',linewidth=1.8)
"""

pylab.plot(phiF,yII,'c+-',markersize=mm)
pylab.plot(phiF,yIQ,'y+-',markersize=mm)
pylab.plot(phiF,yQI,'r+-',markersize=mm)
pylab.plot(phiF,yQQ,'g+-',markersize=mm)

pylab.plot(phiF,VII,'k--',linewidth=1.8)
pylab.plot(phiF,VIQ,'k--',linewidth=1.8)
pylab.plot(phiF,VQI,'k--',linewidth=1.8)
pylab.plot(phiF,VQQ,'k--',linewidth=1.8)

#a=0
#pylab.plot([(FSR-cavity_pole-a)*1e-6, (FSR-cavity_pole-a)*1e-6], [-0.3,0.2],'m--')
#pylab.plot([(FSR+cavity_pole+a)*1e-6, (FSR+cavity_pole+a)*1e-6], [-0.3,0.2],'m--')

pylab.xlim(-0.01, 0.01)
#pylab.xlim(9.08, 9.118)
#pylab.xlim(9.0905, 9.0925)
#pylab.xlim(9.1075, 9.1105)
#pylab.xlim(45.48, 45.53)
#pylab.ylim(-0.35, 0.2)
#pylab.ylim(-0.02, 0.02)

pylab.grid(True)
#pylab.ylabel('Signal (arb)', fontsize=20)
pylab.ylabel(r'$S_{tr}$ (arb)', fontsize=20)
pylab.xticks(fontsize=18)
pylab.yticks(fontsize=18)
ax1.get_yaxis().set_label_coords(-0.13,0.5)

#pylab.xticks(visible=False)
#pylab.title('Cavity Width Fit - Data: solid; Fit: dashed', fontsize=16)

#pylab.legend(('RF I AF I','RF I AF Q','RF Q AF I','RF Q AF Q'),fancybox=True,fontsize=16,bbox_to_anchor=[1.07, 0.52])

pylab.legend(('RF in-phase, audio in-phase','RF in-phase, audio quad-phase','RF quad-phase, audio in-phase','RF quad-phase, audio quad-phase'),
             fancybox=True,loc='upper center', bbox_to_anchor=(0.47, 1.3),ncol=2, fontsize=14)

#pylab.title('REFL9 1000Hz')

ax2 = pylab.subplot(gs[1])

#pylab.plot(freq/1e6,yII-VII,'c+',markersize=mm)
#pylab.plot(freq/1e6,(yIQ-VIQ),'y+',markersize=mm)
#pylab.plot(freq/1e6,yQI-VQI,'r+',markersize=mm)
#pylab.plot(freq/1e6,(yQQ-VQQ),'g+',markersize=mm)

pylab.plot(phiF,yII-VII,'c+',markersize=mm)
pylab.plot(phiF,(yIQ-VIQ),'y+',markersize=mm)
pylab.plot(phiF,yQI-VQI,'r+',markersize=mm)
pylab.plot(phiF,(yQQ-VQQ),'g+',markersize=mm)

pylab.xlim(-0.01, 0.01)
#pylab.xlim(9.0905, 9.0925)
#pylab.xlim(9.1075, 9.1105)
#pylab.xlim(9.08, 9.118)

pylab.grid(True)
#pylab.xlabel('Frequency (MHz)', fontsize=20,labelpad=10)
pylab.xlabel(r'$\phi_f$ (rad)', fontsize=20,labelpad=10)
pylab.ylabel('Residuals', fontsize=20)
ax2.get_yaxis().set_label_coords(-0.13,0.5)

pylab.yticks(array([-0.010, 0.0, 0.010]))

pylab.xticks(fontsize=18)
pylab.yticks(fontsize=18)

#pylab.savefig('plots/REFL9_1000_full.pdf',bbox_inches='tight')
pylab.savefig('plots/REFL9_1000_full_phi.pdf',bbox_inches='tight')

#pylab.savefig('/Users/dhoak/Desktop/LIGO/aLIGO/IMC_demod/AMLengthv4/global_analytic_lmfit.pdf',bbox_inches='tight')







fignum=fignum+1
pylab.figure(fignum)

ax1 = pylab.subplot(gs[0])

#pylab.plot(FSR_freq/1e6,fII,'c+',markersize=5)
#pylab.plot(FSR_freq/1e6,fIQ,'y+',markersize=5)
#pylab.plot(FSR_freq/1e6,fQI,'r+',markersize=5)
#pylab.plot(FSR_freq/1e6,fQQ,'g+',markersize=5)

#pylab.plot(FSR_freq/1e6,FII,'k--',linewidth=1)
#pylab.plot(FSR_freq/1e6,FIQ,'k--',linewidth=1)
#pylab.plot(FSR_freq/1e6,FQI,'k--',linewidth=1)
#pylab.plot(FSR_freq/1e6,FQQ,'k--',linewidth=1)

pylab.plot(FSR_freq-FSR,fII,'c+',markersize=8)
pylab.plot(FSR_freq-FSR,fIQ,'y+',markersize=8)
pylab.plot(FSR_freq-FSR,fQI,'r+',markersize=8)
pylab.plot(FSR_freq-FSR,fQQ,'g+',markersize=8)

pylab.plot(FSR_freq-FSR,FII,'k--',linewidth=1)
pylab.plot(FSR_freq-FSR,FIQ,'k--',linewidth=1)
pylab.plot(FSR_freq-FSR,FQI,'k--',linewidth=1)
pylab.plot(FSR_freq-FSR,FQQ,'k--',linewidth=1)

pylab.xlim(-1100,1100)
#pylab.xlim(9.0988, 9.1012)
#pylab.ylim(-0.04, 0.07)

pylab.grid(True)
pylab.ylabel('Signal [arb]')
pylab.xticks(fontsize=12)
pylab.yticks(fontsize=12)
ax1.get_yaxis().set_label_coords(-0.10,0.5)

pylab.title('Cavity Length Fit - Data: solid; Fit: dashed, FSR = ' + str(FSR) + ' Hz', fontsize=12)

#pylab.legend(('RF I AF I','RF I AF Q','RF Q AF I','RF Q AF Q'),fancybox=True,fontsize=11,bbox_to_anchor=[1.07, 1.07])

ax2 = pylab.subplot(gs[1])

#pylab.plot(FSR_freq/1e6,fII-FII,'c+',markersize=5)
#pylab.plot(FSR_freq/1e6,fIQ-FIQ,'y+',markersize=5)
#pylab.plot(FSR_freq/1e6,fQI-FQI,'r+',markersize=5)
#pylab.plot(FSR_freq/1e6,fQQ-FQQ,'g+',markersize=5)

pylab.plot(FSR_freq-FSR,fII-FII,'c+',markersize=8)
pylab.plot(FSR_freq-FSR,fIQ-FIQ,'y+',markersize=8)
pylab.plot(FSR_freq-FSR,fQI-FQI,'r+',markersize=8)
pylab.plot(FSR_freq-FSR,fQQ-FQQ,'g+',markersize=8)

pylab.xlim(-1100,1100)
#pylab.xlim(9.0988, 9.1012)

pylab.grid(True)
#pylab.xlabel('Frequency [MHz]')
pylab.xlabel('RF Frequency - FSR [Hz]')
pylab.ylabel('residuals')

pylab.xticks(fontsize=12)
pylab.yticks(fontsize=12)
ax2.get_yaxis().set_label_coords(-0.10,0.5)

pylab.savefig('plots/global_analytic_zoom.pdf',bbox_inches='tight')





OSA_data = genfromtxt('data/kiwamu_OSA_sideband.txt')
f_MHz = OSA_data[:,0]
OSA = OSA_data[:,1]

fignum=fignum+1
pylab.figure(fignum)

f = arange(8.4e6, 10e6, 1e3)

A = 0.4
f0 = params['f0'].value
Q = params['Q'].value

#A = 0.5
#f0 = 9.137707e6
#Q = 100

pylab.subplot(2,1,1)

pylab.plot(f_MHz,OSA,'ro')

lorentz = A * sqrt( (f0/(2*Q))**2 / ( (f-f0)**2 + (f0/(2*Q))**2 ) )
pylab.plot(f*1e-6, lorentz, 'k')

lorentz = A * sqrt( (f0/(2*(Q+Q_err)))**2 / ( (f-f0)**2 + (f0/(2*(Q+Q_err)))**2 ) )
pylab.plot(f*1e-6, lorentz, 'b')

lorentz = A * sqrt( (f0/(2*(Q-Q_err)))**2 / ( (f-f0)**2 + (f0/(2*(Q-Q_err)))**2 ) )
pylab.plot(f*1e-6, lorentz, 'r')

pylab.plot([FSR*1e-6, FSR*1e-6],[0.0, 1.0],'g-')

pylab.plot([min(freq)*1e-6, min(freq)*1e-6],[0.0, 1.0],'k--',linewidth=0.8)
pylab.plot([max(freq)*1e-6, max(freq)*1e-6],[0.0, 1.0],'k--',linewidth=0.8)

pylab.grid(True)
pylab.ylim(0.0,0.6)
#pylab.ylim(0.3,0.4)
#pylab.xlim(min(freq)/1e6,max(freq)/1e6)
pylab.ylabel('Amplitude [arb]',fontsize=10)
pylab.xticks(fontsize=9)
pylab.yticks(fontsize=9)

#pylab.legend(('Kiwamu data','Resonance fit from IMC width data'),fancybox=True,loc=1,fontsize=9)

pylab.title('9MHz EOM Resonance - OSA measurement and fit parameters', fontsize=10)

pylab.subplot(2,1,2)

phi = arctan((f-f0) / (f0/(2*Q)) )
pylab.plot(f*1e-6, phi*180/pi, 'k')

phi = arctan((f-f0) / (f0/(2*(Q+Q_err))) )
pylab.plot(f*1e-6, phi*180/pi, 'b')

phi = arctan((f-f0) / (f0/(2*(Q-Q_err))) )
pylab.plot(f*1e-6, phi*180/pi, 'r')

pylab.plot([min(freq)*1e-6, min(freq)*1e-6],[-100.0, 100],'k--',linewidth=0.8)
pylab.plot([max(freq)*1e-6, max(freq)*1e-6],[-100.0, 100],'k--',linewidth=0.8)

pylab.grid(True)
#pylab.xlim(8.8,9.4)
#pylab.ylim(-30,10)
#pylab.xlim(min(freq)/1e6,max(freq)/1e6)
pylab.xlabel('Frequency [MHz]',fontsize=10)
pylab.ylabel('Phase shift [deg]',fontsize=10)
pylab.xticks(fontsize=9)
pylab.yticks(fontsize=9)

pylab.savefig('plots/EOM_resonance.pdf')


