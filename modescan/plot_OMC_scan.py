#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script to plot frames for a movie of OMC sweeps over time
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
from scipy import signal
import matplotlib.gridspec as gridspec
from scipy.special import jn
#from scipy.signal import find_peaks_cwt, argrelmin
#from extrema import minima
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters, fit_report

def lorentzian(f,params):

    A = params['A'].value
    f0 = params['f0'].value
    fc = params['fc'].value
    return A/(1.0+((f-f0)/fc)**2)

def error_function(params,volts,data):

    fit = lorentzian(volts, params)
    return data - fit

def fit_peak(volts,data,A0,f0):

    params = Parameters()
    params.add('A', A0, vary=True)
    params.add('f0', f0, vary=True)
    params.add('fc', 0.05, vary=True)

    out = minimize(error_function, params, args=(volts,data), method='nelder')
    out = minimize(error_function, params, args=(volts,data), method='leastsq')

    fit_vals = [params['A'].value, params['f0'].value, params['fc'].value]
    fit_errors = [params['A'].stderr, params['f0'].stderr, params['fc'].stderr]

    return fit_vals, fit_errors


def freqfit(v,params):

    a0 = params['a0'].value
    a1 = params['a1'].value
    a2 = params['a2'].value
    a3 = params['a3'].value
    a4 = params['a4'].value
    a5 = params['a5'].value

    return a0 + a1*v + a2*v**2 + a3*v**3 + a4*v**4 + a5*v**5

def error_function2(params,volts,freqs):
    fit = freqfit(volts,params)
    return freqs - fit

def invfreqfit(v,params):

    a0 = params['a0'].value
    a1 = params['a1'].value
    a2 = params['a2'].value
    a3 = params['a3'].value

    return a0 + a1*v + a2*v**2 + a3*v**3

def error_function3(params,volts,freqs):
    fit = invfreqfit(freqs,params)
    return volts - fit

def fit_pzt_volts(volts,freqs,FSR_slope):

    params = Parameters()
    params.add('a0', 350.0e6, vary=True)
    params.add('a1', -1/FSR_slope, vary=True)
    params.add('a2', 0.0, vary=True)
    params.add('a3', 0.0, vary=True)
    params.add('a4', 0.0, vary=False)
    params.add('a5', 0.0, vary=False)

    out = minimize(error_function2, params, args=(volts,freqs), method='nelder')
    out = minimize(error_function2, params, args=(volts,freqs), method='leastsq')

    fit_vals = [params['a0'].value, params['a1'].value, params['a2'].value, params['a3'].value, params['a4'].value, params['a5'].value]
    fit_errors = [params['a0'].stderr, params['a1'].stderr, params['a2'].stderr, params['a3'].stderr, params['a4'].stderr, params['a5'].stderr]

    return fit_vals, fit_errors


def fit_pzt_freq(freqs,volts,FSR_slope):

    params = Parameters()
    params.add('a0', 1.0, vary=True)
    params.add('a1', FSR_slope, vary=True)
    params.add('a2', 0.0, vary=True)
    params.add('a3', 0.0, vary=True)
    params.add('a4', 0.0, vary=False)
    params.add('a5', 0.0, vary=False)

    out = minimize(error_function3, params, args=(volts,freqs), method='nelder')
    out = minimize(error_function3, params, args=(volts,freqs), method='leastsq')

    fit_vals = [params['a0'].value, params['a1'].value, params['a2'].value, params['a3'].value]
    fit_errors = [params['a0'].stderr, params['a1'].stderr, params['a2'].stderr, params['a3'].stderr]

    return fit_vals, fit_errors


fontP = FontProperties()
fontP.set_size('medium')

fig_width_pt = 600  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
#golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio, = 0.618
golden_mean = 0.7         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size = [fig_width,fig_height]
matplotlib.rcParams.update(
        {'axes.labelsize': 14, 
            'font.size':   14, 
            'legend.fontsize': 14, 
            'xtick.labelsize': 12, 
            'ytick.labelsize': 12, 
            'text.usetex': False,
            'figure.figsize': fig_size,
            'font.family': "serif",
            'font.serif': ["Times New Roman"],
            'savefig.dpi': 200,
            'xtick.major.size':8,
            'xtick.minor.size':4,
            'ytick.major.size':8,
            'ytick.minor.size':4
            })


def decimate(x, q, n=None, ftype='iir', axis=-1):
    if not isinstance(q, int):
        raise TypeError("q must be an integer")
    if n is None:
        if ftype == 'fir':
            n = 30
        else:
            n = 3
    if ftype == 'fir':
        b = signal.firwin(n + 1, 1. / q, window='hamming')
        a = 1.
    else:
        b, a = signal.cheby1(n, 0.05, 0.8 / q)

    y = signal.filtfilt(b, a, x)
    sl = [slice(None)] * y.ndim
    sl[axis] = slice(None, None, q)
    return y[sl]

data_string = 'on'

data = genfromtxt('data/OMC_modescan_full_lock_offset_' + data_string + '_4Apr2015.txt')

xmin = 30
xmax = 90
ymin = 1e-3
ymax = 50

FSR = 261.72
HOMs = 57.33

f1 = 9.1002e6
f2 = 5*f1

times1 = decimate(data[:,0],4)
pzt1 = decimate(data[:,1],4)
dcpd1 = decimate(data[:,2],4)

times = data[:,0]
pzt = data[:,1]
dcpd = data[:,2]

tstart = min(times)
tstop = max(times)

dt = (tstop-tstart)/float(len(times))

peak_threshold = max(dcpd)*0.8

print dt, 1/dt

# Now find the sections of data where the voltage was increasing
positive_segs = []
start_flag = False
for i in range(1,len(pzt1)):

    if not(start_flag) and (pzt1[i]-pzt1[i-1])>0:
        start_flag = True
        start_idx = i

    elif start_flag and (pzt1[i]-pzt1[i-1])<0:
        start_flag = False
        positive_segs.append([start_idx, i])

# If the loop ends inside of an upgoing segment, set the last sample to the end
if start_flag:
    positive_segs.append([start_idx, i])

seg_peaks = []
i=0
# Now find the peaks in the segments
for seg in positive_segs:

    seg_length = seg[1]-seg[0]
    if len(seg)==2 and seg_length>100:

        # pad the segment by a few samples on either end to avoid transients
        dcpd_sum = dcpd1[seg[0]+10:seg[1]-10]
        peak_idx = signal.find_peaks_cwt(dcpd_sum,array([int(floor(seg_length*0.02))]))

        carrier_peaks = []
        for jj in range(len(peak_idx)):

            # Sometimes it doesn't find the peak???
            peakmax = argmax(dcpd_sum[peak_idx[jj]-20:peak_idx[jj]+20])

            if dcpd_sum[peak_idx[jj]+peakmax-20] > peak_threshold:
                carrier_peaks.append(peak_idx[jj]+peakmax-20+seg[0]+10)
                
        seg_peaks.append(carrier_peaks)

# Now we should have an array of peak indices
# These should include at least two sideband peaks
# Grab the times
peak_times = []
for seg in seg_peaks:
    if len(seg)>=4:

        # We are looking for the pairs of the 45MHz sidebands.  These should be tall, similar heights, and separated by ~12V
        group_count = 0
        peak_group = []
        for j in range(len(seg)):

            idx1 = seg[j]
            for k in range(j+1,len(seg)):

                idx2 = seg[k]
                if (pzt1[idx2]-pzt1[idx1])>10 and (pzt1[idx2]-pzt1[idx1])<15 and (dcpd1[idx1]/dcpd1[idx2])<1.1 and (dcpd1[idx1]/dcpd1[idx2])>0.9:

                    # We have found a pair of sidebands
                    peak_group.append([idx1,idx2])
                    group_count+=1

        if group_count==2:
            peak_times.append([times1[peak_group[0][0]],times1[peak_group[0][1]],times1[peak_group[1][0]],times1[peak_group[1][1]]])


# Now connect the times to the indices of the full data array
# Get the frequency as a linear fit to the voltage, using the FSR

a, b = shape(peak_times)

print
print 'We found ' + str(a) + ' sweeps.'
print
print 'Starting fits.'

print
print '45MHz Sidebands'

sb45_idx = zeros((a,b))
FSR_calib1 = zeros(a)
sb45_fits = zeros((a,b,3))
sb45_errors = zeros((a,b,3))
for i in range(a):

    sb45_idx[i,0] = argmin(abs(times - peak_times[i][0]))
    sb45_idx[i,1] = argmin(abs(times - peak_times[i][1]))
    sb45_idx[i,2] = argmin(abs(times - peak_times[i][2]))
    sb45_idx[i,3] = argmin(abs(times - peak_times[i][3]))

    s = 200
    # Fit the 45MHz peaks
    for j in range(4):

        volts = pzt[sb45_idx[i,j]-s:sb45_idx[i,j]+s]
        data = dcpd[sb45_idx[i,j]-s:sb45_idx[i,j]+s]
        A0 = dcpd[sb45_idx[i,j]]
        f0 = pzt[sb45_idx[i,j]]
        sb45_fits[i,j,:], sb45_errors[i,j,:] = fit_peak(volts,data,A0,f0)

    # Get the frequency calibration in terms of the PZT voltage
    FSR_calib1[i] = (sb45_fits[i,3,1] - sb45_fits[i,1,1])/(FSR*1e6)
    

print
print 'Carrier'

fit_volts = zeros((a,6))
freq_fitting = zeros((a,4))
invfreq_fitting = zeros((a,4))

tem00_idx = zeros((a,2))
tem00_fits = zeros((a,2,3))
tem00_errors = zeros((a,2,3))
tem00_guess = zeros((a,2))
for i in range(a):
    LSB2 = pzt[sb45_idx[i,3]]
    full_volts = pzt[sb45_idx[i,0]:sb45_idx[i,3]+500]

    sb1 = sb45_fits[i,0,0]/(1.0+((full_volts-sb45_fits[i,0,1])/sb45_fits[i,0,2])**2)
    sb2 = sb45_fits[i,1,0]/(1.0+((full_volts-sb45_fits[i,1,1])/sb45_fits[i,1,2])**2)
    sb3 = sb45_fits[i,2,0]/(1.0+((full_volts-sb45_fits[i,2,1])/sb45_fits[i,2,2])**2)
    sb4 = sb45_fits[i,3,0]/(1.0+((full_volts-sb45_fits[i,3,1])/sb45_fits[i,3,2])**2)

    tem_model = sb1 + sb2 + sb3 + sb4

    fit_vals = invfreq_fitting[i,:]

    r = 100
    for j in range(2):

        # Based on the frequency, guess the PZT voltage
        #freq = (j*FSR*1e6)
        #tem_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2

        freq = (j*FSR*1e6 + f2)
        tem_guess = LSB2 - freq*FSR_calib1[i]
        #print tem_guess, tem_guess1

        guess_idx = argmin(abs(pzt[sb45_idx[i,0]:sb45_idx[i,3]] - tem_guess)) + sb45_idx[i,0]
        tem00_guess[i,j] = guess_idx

        peak_idx = signal.find_peaks_cwt(dcpd[guess_idx-r:guess_idx+r],array([r/10])) + guess_idx-r
        peak_heights = []
        for peak in peak_idx:
            peak_heights.append(dcpd[peak] - tem_model[peak-sb45_idx[i,0]])

        tem00_idx[i,j] = peak_idx[argmax(peak_heights)]

        s = 35
        volts = pzt[tem00_idx[i,j]-s:tem00_idx[i,j]+s]
        data = dcpd[tem00_idx[i,j]-s:tem00_idx[i,j]+s] - tem_model[tem00_idx[i,j]-sb45_idx[i,0]-s:tem00_idx[i,j]-sb45_idx[i,0]+s]
        A0 = dcpd[tem00_idx[i,j]] - tem_model[tem00_idx[i,j]-sb45_idx[i,0]]
        f0 = pzt[tem00_idx[i,j]]
        tem00_fit, tem00_errors = fit_peak(volts,data,A0,f0)
        tem00_fits[i,j,:] = tem00_fit


    fit_volts[i,:] = array([sb45_fits[i,0,1], tem00_fits[i,1,1], sb45_fits[i,1,1], sb45_fits[i,2,1], tem00_fits[i,0,1], sb45_fits[i,3,1]])
    fit_freqs = array([FSR*1e6+f2, FSR*1e6, FSR*1e6-f2, f2, 0.0, -1*f2])

    #print FSR_calib1[i], 39/(FSR*1e6)
    fit_vals, fit_errors = fit_pzt_volts(fit_volts[i,:],fit_freqs,FSR_calib1[i])
    freq_fitting[i,:] = fit_vals[0:4]

    fit_vals, fit_errors = fit_pzt_freq(fit_freqs,fit_volts[i,:],FSR_calib1[i])
    invfreq_fitting[i,:] = fit_vals[0:4]


print
print 'Carrier HOMs'

hom_text = ['CR1','CR2','CR3','CR4','CR5','CR6','CR7','CR8','CR9','CR10','CR11','CR12','CR13','CR14','CR15']

Hnum = 9
hom_idx = zeros((a,Hnum))
hom_fits = zeros((a,Hnum,3))

for i in range(a):
    USB1 = pzt[sb45_idx[i,0]]
    LSB2 = pzt[sb45_idx[i,3]]
    full_volts = pzt[sb45_idx[i,0]:sb45_idx[i,3]+500]

    sb1 = sb45_fits[i,0,0]/(1.0+((full_volts-sb45_fits[i,0,1])/sb45_fits[i,0,2])**2)
    sb2 = sb45_fits[i,1,0]/(1.0+((full_volts-sb45_fits[i,1,1])/sb45_fits[i,1,2])**2)
    sb3 = sb45_fits[i,2,0]/(1.0+((full_volts-sb45_fits[i,2,1])/sb45_fits[i,2,2])**2)
    sb4 = sb45_fits[i,3,0]/(1.0+((full_volts-sb45_fits[i,3,1])/sb45_fits[i,3,2])**2)

    car1 = tem00_fits[i,0,0]/(1.0+((full_volts-tem00_fits[i,0,1])/tem00_fits[i,0,2])**2)
    car2 = tem00_fits[i,1,0]/(1.0+((full_volts-tem00_fits[i,1,1])/tem00_fits[i,1,2])**2)

    tem_model = sb1 + sb2 + sb3 + sb4 + car1 + car2

    fit_vals = invfreq_fitting[i,:]

    for j in range(Hnum):
        freq = remainder((j+1)*HOMs*1e6,FSR*1e6)
        hom_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

        #hom_guess1 = LSB2 - remainder(f2+(j+1)*HOMs*1e6,FSR*1e6)*FSR_calib1[i]

        #print j+1, freq*1e-6, hom_guess
        guess_idx = argmin(abs(pzt[sb45_idx[i,0]:sb45_idx[i,3]] - hom_guess)) + sb45_idx[i,0]

        r = 200
        peak_idx = signal.find_peaks_cwt(dcpd[guess_idx-r:guess_idx+r],array([r/8])) + guess_idx-r
        peak_heights = []
        for peak in peak_idx:
            peak_heights.append(dcpd[peak] - tem_model[peak-sb45_idx[i,0]])

        hom_idx[i,j] = peak_idx[argmax(peak_heights)]

        s = 35
        volts = pzt[hom_idx[i,j]-s:hom_idx[i,j]+s]
        data = dcpd[hom_idx[i,j]-s:hom_idx[i,j]+s] - tem_model[hom_idx[i,j]-sb45_idx[i,0]-s:hom_idx[i,j]-sb45_idx[i,0]+s]
        A0 = dcpd[hom_idx[i,j]] - tem_model[hom_idx[i,j]-sb45_idx[i,0]]
        f0 = pzt[hom_idx[i,j]]
        hom_fit, hom_errors = fit_peak(volts,data,A0,f0)
        hom_fits[i,j,:] = hom_fit



ehom_text = ['CR4','CR5','CR9']
print
print 'Extra Carrier HOMs'

eHnum_locs = array([4*HOMs*1e6-FSR*1e6, 5*HOMs*1e6, 9*HOMs*1e6-2*FSR*1e6])
eHnum = [3, 4, 9]  # index needs to be one smaller, keeps code consistent
ehom_idx = zeros((a,len(eHnum)))
ehom_fits = zeros((a,len(eHnum),3))

for i in range(a):
    USB1 = pzt[sb45_idx[i,0]]
    LSB2 = pzt[sb45_idx[i,3]]
    full_volts = pzt[sb45_idx[i,0]:sb45_idx[i,3]+500]

    sb1 = sb45_fits[i,0,0]/(1.0+((full_volts-sb45_fits[i,0,1])/sb45_fits[i,0,2])**2)
    sb2 = sb45_fits[i,1,0]/(1.0+((full_volts-sb45_fits[i,1,1])/sb45_fits[i,1,2])**2)
    sb3 = sb45_fits[i,2,0]/(1.0+((full_volts-sb45_fits[i,2,1])/sb45_fits[i,2,2])**2)
    sb4 = sb45_fits[i,3,0]/(1.0+((full_volts-sb45_fits[i,3,1])/sb45_fits[i,3,2])**2)

    car1 = tem00_fits[i,0,0]/(1.0+((full_volts-tem00_fits[i,0,1])/tem00_fits[i,0,2])**2)
    car2 = tem00_fits[i,1,0]/(1.0+((full_volts-tem00_fits[i,1,1])/tem00_fits[i,1,2])**2)

    tem_model = sb1 + sb2 + sb3 + sb4 + car1 + car2

    for j in range(Hnum):
        tem_model += hom_fits[i,j,0]/(1.0+((full_volts-hom_fits[i,j,1])/hom_fits[i,j,2])**2)        

    fit_vals = invfreq_fitting[i,:]

    for j in range(len(eHnum)):
        freq = eHnum_locs[j]
        ehom_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

        ehom_guess1 = LSB2 - (FSR*1e6+remainder(f2+(eHnum[j]+1)*HOMs*1e6,FSR*1e6))*FSR_calib1[i]

        #print j+1, freq*1e-6, ehom_guess, ehom_guess1
        guess_idx = argmin(abs(pzt[sb45_idx[i,0]:sb45_idx[i,3]] - ehom_guess)) + sb45_idx[i,0]

        r = 200
        peak_idx = signal.find_peaks_cwt(dcpd[guess_idx-r:guess_idx+r],array([r/6])) + guess_idx-r
        peak_heights = []
        for peak in peak_idx:
            peak_heights.append(dcpd[peak] - tem_model[peak-sb45_idx[i,0]])

        ehom_idx[i,j] = peak_idx[argmax(peak_heights)]

        s = 35
        volts = pzt[ehom_idx[i,j]-s:ehom_idx[i,j]+s]
        data = dcpd[ehom_idx[i,j]-s:ehom_idx[i,j]+s] - tem_model[ehom_idx[i,j]-sb45_idx[i,0]-s:ehom_idx[i,j]-sb45_idx[i,0]+s]
        A0 = dcpd[ehom_idx[i,j]] - tem_model[ehom_idx[i,j]-sb45_idx[i,0]]
        f0 = pzt[ehom_idx[i,j]]
        ehom_fit, ehom_errors = fit_peak(volts,data,A0,f0)
        ehom_fits[i,j,:] = ehom_fit



#print invfreq_fitting

print
print 'Fitting the PZT value to optical frequency using Carrier HOMs'


# Update the fitting with the carrier HOMs
fit_volts = zeros((a,6+Hnum+len(eHnum_locs)))
freq_fitting = zeros((a,4))
invfreq_fitting = zeros((a,4))

for i in range(a):

    fit_volts1 = array([sb45_fits[i,0,1], tem00_fits[i,1,1], sb45_fits[i,1,1], sb45_fits[i,2,1], tem00_fits[i,0,1], sb45_fits[i,3,1]])

    fit_volts[i,:] = hstack((fit_volts1, hom_fits[i,:,1], ehom_fits[i,:,1]))

    fit_freqs1 = array([FSR*1e6+f2, FSR*1e6, FSR*1e6-f2, f2, 0.0, -1*f2])
    HOM_freqs = remainder(arange(1,Hnum+1)*HOMs*1e6,FSR*1e6)

    fit_freqs = hstack((fit_freqs1, HOM_freqs, eHnum_locs))

    #print FSR_calib1[i], 39/(FSR*1e6)
    fit_vals, fit_errors = fit_pzt_volts(fit_volts[i,:],fit_freqs,FSR_calib1[i])
    freq_fitting[i,:] = fit_vals[0:4]

    fit_vals, fit_errors = fit_pzt_freq(fit_freqs,fit_volts[i,:],FSR_calib1[i])
    invfreq_fitting[i,:] = fit_vals[0:4]

print
print 'Inverse frequency fitting parameters (MHz --> Volts):'
print invfreq_fitting

print
print 'Frequency fitting parameters (Volts --> MHz):'
print freq_fitting




print
print 'LSB 45MHz HOMs'

lsb_text = ['LSB1', 'LSB2', 'LSB1', 'LSB3', 'LSB4']
lsb_locs = array([HOMs*1e6-f2, 2*HOMs*1e6-f2, HOMs*1e6+FSR*1e6-f2, 3*HOMs*1e6-f2, 4*HOMs*1e6-f2])

lsb_idx = zeros((a,len(lsb_locs)))
lsb_fits = zeros((a,len(lsb_locs),3))
for i in range(a):
    USB1 = pzt[sb45_idx[i,0]]
    LSB2 = pzt[sb45_idx[i,3]]
    full_volts = pzt[sb45_idx[i,0]:sb45_idx[i,3]+500]

    sb1 = sb45_fits[i,0,0]/(1.0+((full_volts-sb45_fits[i,0,1])/sb45_fits[i,0,2])**2)
    sb2 = sb45_fits[i,1,0]/(1.0+((full_volts-sb45_fits[i,1,1])/sb45_fits[i,1,2])**2)
    sb3 = sb45_fits[i,2,0]/(1.0+((full_volts-sb45_fits[i,2,1])/sb45_fits[i,2,2])**2)
    sb4 = sb45_fits[i,3,0]/(1.0+((full_volts-sb45_fits[i,3,1])/sb45_fits[i,3,2])**2)

    car1 = tem00_fits[i,0,0]/(1.0+((full_volts-tem00_fits[i,0,1])/tem00_fits[i,0,2])**2)
    car2 = tem00_fits[i,1,0]/(1.0+((full_volts-tem00_fits[i,1,1])/tem00_fits[i,1,2])**2)

    tem_model = sb1 + sb2 + sb3 + sb4 + car1 + car2

    fit_vals = invfreq_fitting[i,:]

    for j in range(len(lsb_locs)):
        freq = lsb_locs[j]
        lsb_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

        lsb_guess1 = LSB2 - lsb_locs[j]*FSR_calib1[i]
        guess_idx = argmin(abs(pzt[sb45_idx[i,0]:sb45_idx[i,3]] - lsb_guess)) + sb45_idx[i,0]

        r = 80
        peak_idx = signal.find_peaks_cwt(dcpd[guess_idx-r:guess_idx+r],array([r/10])) + guess_idx-r
        peak_heights = []
        for peak in peak_idx:
            peak_heights.append(dcpd[peak] - tem_model[peak-sb45_idx[i,0]])

        lsb_idx[i,j] = peak_idx[argmax(peak_heights)]

        s = 30
        volts = pzt[lsb_idx[i,j]-s:lsb_idx[i,j]+s]
        data = dcpd[lsb_idx[i,j]-s:lsb_idx[i,j]+s] - tem_model[lsb_idx[i,j]-sb45_idx[i,0]-s:lsb_idx[i,j]-sb45_idx[i,0]+s]
        A0 = dcpd[lsb_idx[i,j]] - tem_model[lsb_idx[i,j]-sb45_idx[i,0]]
        f0 = pzt[lsb_idx[i,j]]
        lsb_fit, lsb_errors = fit_peak(volts,data,A0,f0)
        lsb_fits[i,j,:] = lsb_fit


print
print 'USB 45MHz HOMs'

usb_text = ['USB1', 'USB2']
usb_locs = array([HOMs*1e6+f2, 2*HOMs*1e6+f2])

usb_idx = zeros((a,len(usb_locs)))
usb_fits = zeros((a,len(usb_locs),3))
for i in range(a):
    USB1 = pzt[sb45_idx[i,0]]
    LSB2 = pzt[sb45_idx[i,3]]
    full_volts = pzt[sb45_idx[i,0]:sb45_idx[i,3]+500]

    sb1 = sb45_fits[i,0,0]/(1.0+((full_volts-sb45_fits[i,0,1])/sb45_fits[i,0,2])**2)
    sb2 = sb45_fits[i,1,0]/(1.0+((full_volts-sb45_fits[i,1,1])/sb45_fits[i,1,2])**2)
    sb3 = sb45_fits[i,2,0]/(1.0+((full_volts-sb45_fits[i,2,1])/sb45_fits[i,2,2])**2)
    sb4 = sb45_fits[i,3,0]/(1.0+((full_volts-sb45_fits[i,3,1])/sb45_fits[i,3,2])**2)

    car1 = tem00_fits[i,0,0]/(1.0+((full_volts-tem00_fits[i,0,1])/tem00_fits[i,0,2])**2)
    car2 = tem00_fits[i,1,0]/(1.0+((full_volts-tem00_fits[i,1,1])/tem00_fits[i,1,2])**2)

    tem_model = sb1 + sb2 + sb3 + sb4 + car1 + car2

    fit_vals = invfreq_fitting[i,:]

    for j in range(len(usb_locs)):
        freq = usb_locs[j]
        usb_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

        usb_guess1 = LSB2 - usb_locs[j]*FSR_calib1[i]
        guess_idx = argmin(abs(pzt[sb45_idx[i,0]:sb45_idx[i,3]] - usb_guess)) + sb45_idx[i,0] 

        #print freq*1e-6, usb_guess, usb_guess1

        r = 80
        peak_idx = signal.find_peaks_cwt(dcpd[guess_idx-r:guess_idx+r],array([r/8])) + guess_idx-r
        peak_heights = []
        for peak in peak_idx:
            peak_heights.append(dcpd[peak] - tem_model[peak-sb45_idx[i,0]])

        usb_idx[i,j] = peak_idx[argmax(peak_heights)]

        s = 30
        volts = pzt[usb_idx[i,j]-s:usb_idx[i,j]+s]
        data = dcpd[usb_idx[i,j]-s:usb_idx[i,j]+s] - tem_model[usb_idx[i,j]-sb45_idx[i,0]-s:usb_idx[i,j]-sb45_idx[i,0]+s]
        A0 = dcpd[usb_idx[i,j]] - tem_model[usb_idx[i,j]-sb45_idx[i,0]]
        f0 = pzt[usb_idx[i,j]]
        usb_fit, usb_errors = fit_peak(volts,data,A0,f0)
        usb_fits[i,j,:] = usb_fit



print
print 'LSB 9MHz HOMs'

lsb9_text = ['lsb2', 'lsb3', 'lsb4', 'lsb4', 'lsb5', 'lsb5', 'lsb6', 'lsb7', 'lsb8']
lsb9_locs = array([2*HOMs*1e6-f1, 3*HOMs*1e6-f1, 4*HOMs*1e6-f1, 4*HOMs*1e6-FSR*1e6-f1, 5*HOMs*1e6-f1, 5*HOMs*1e6-FSR*1e6-f1, 6*HOMs*1e6-FSR*1e6-f1, 7*HOMs*1e6-FSR*1e6-f1, 8*HOMs*1e6-FSR*1e6-f1])

lsb9_idx = zeros((a,len(lsb9_locs)))
lsb9_fits = zeros((a,len(lsb9_locs),3))
for i in range(a):
    USB1 = pzt[sb45_idx[i,0]]
    LSB2 = pzt[sb45_idx[i,3]]
    full_volts = pzt[sb45_idx[i,0]:sb45_idx[i,3]+500]

    sb1 = sb45_fits[i,0,0]/(1.0+((full_volts-sb45_fits[i,0,1])/sb45_fits[i,0,2])**2)
    sb2 = sb45_fits[i,1,0]/(1.0+((full_volts-sb45_fits[i,1,1])/sb45_fits[i,1,2])**2)
    sb3 = sb45_fits[i,2,0]/(1.0+((full_volts-sb45_fits[i,2,1])/sb45_fits[i,2,2])**2)
    sb4 = sb45_fits[i,3,0]/(1.0+((full_volts-sb45_fits[i,3,1])/sb45_fits[i,3,2])**2)

    car1 = tem00_fits[i,0,0]/(1.0+((full_volts-tem00_fits[i,0,1])/tem00_fits[i,0,2])**2)
    car2 = tem00_fits[i,1,0]/(1.0+((full_volts-tem00_fits[i,1,1])/tem00_fits[i,1,2])**2)

    tem_model = sb1 + sb2 + sb3 + sb4 + car1 + car2

    tem_model += lsb_fits[i,0,0]/(1.0+((full_volts-lsb_fits[i,0,1])/lsb_fits[i,0,2])**2)
    tem_model += lsb_fits[i,1,0]/(1.0+((full_volts-lsb_fits[i,1,1])/lsb_fits[i,1,2])**2)
    tem_model += lsb_fits[i,2,0]/(1.0+((full_volts-lsb_fits[i,2,1])/lsb_fits[i,2,2])**2)
    tem_model += lsb_fits[i,3,0]/(1.0+((full_volts-lsb_fits[i,3,1])/lsb_fits[i,3,2])**2)
    tem_model += lsb_fits[i,4,0]/(1.0+((full_volts-lsb_fits[i,4,1])/lsb_fits[i,4,2])**2)

    tem_model += usb_fits[i,0,0]/(1.0+((full_volts-usb_fits[i,0,1])/usb_fits[i,0,2])**2)
    tem_model += usb_fits[i,1,0]/(1.0+((full_volts-usb_fits[i,1,1])/usb_fits[i,1,2])**2)

    fit_vals = invfreq_fitting[i,:]

    for j in range(len(lsb9_locs)):
        freq = lsb9_locs[j]
        lsb9_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

        guess_idx = argmin(abs(pzt[sb45_idx[i,0]:sb45_idx[i,3]] - lsb9_guess)) + sb45_idx[i,0]

        r = 80
        peak_idx = signal.find_peaks_cwt(dcpd[guess_idx-r:guess_idx+r],array([r/10])) + guess_idx-r
        peak_heights = []
        for peak in peak_idx:
            peak_heights.append(dcpd[peak] - tem_model[peak-sb45_idx[i,0]])

        lsb9_idx[i,j] = peak_idx[argmax(peak_heights)]

        s = 30
        volts = pzt[lsb9_idx[i,j]-s:lsb9_idx[i,j]+s]
        data = dcpd[lsb9_idx[i,j]-s:lsb9_idx[i,j]+s] - tem_model[lsb9_idx[i,j]-sb45_idx[i,0]-s:lsb9_idx[i,j]-sb45_idx[i,0]+s]
        A0 = dcpd[lsb9_idx[i,j]] - tem_model[lsb9_idx[i,j]-sb45_idx[i,0]]
        f0 = pzt[lsb9_idx[i,j]]
        lsb9_fit, lsb9_errors = fit_peak(volts,data,A0,f0)
        lsb9_fits[i,j,:] = lsb9_fit


print
print 'USB 9MHz HOMs'

usb9_text = ['usb0','usb0','usb1', 'usb2','usb3','usb4','usb4','usb5','usb5','usb6','usb7','usb8']
usb9_locs = array([f1, FSR*1e6+f1,HOMs*1e6+f1, 2*HOMs*1e6+f1, 3*HOMs*1e6+f1, 4*HOMs*1e6+f1, 4*HOMs*1e6-FSR*1e6+f1, 5*HOMs*1e6+f1, 5*HOMs*1e6-FSR*1e6+f1, 6*HOMs*1e6-FSR*1e6+f1, 7*HOMs*1e6-FSR*1e6+f1, 8*HOMs*1e6-FSR*1e6+f1])

usb9_idx = zeros((a,len(usb9_locs)))
usb9_fits = zeros((a,len(usb9_locs),3))
for i in range(a):
    USB1 = pzt[sb45_idx[i,0]]
    LSB2 = pzt[sb45_idx[i,3]]
    full_volts = pzt[sb45_idx[i,0]:sb45_idx[i,3]+500]

    sb1 = sb45_fits[i,0,0]/(1.0+((full_volts-sb45_fits[i,0,1])/sb45_fits[i,0,2])**2)
    sb2 = sb45_fits[i,1,0]/(1.0+((full_volts-sb45_fits[i,1,1])/sb45_fits[i,1,2])**2)
    sb3 = sb45_fits[i,2,0]/(1.0+((full_volts-sb45_fits[i,2,1])/sb45_fits[i,2,2])**2)
    sb4 = sb45_fits[i,3,0]/(1.0+((full_volts-sb45_fits[i,3,1])/sb45_fits[i,3,2])**2)

    car1 = tem00_fits[i,0,0]/(1.0+((full_volts-tem00_fits[i,0,1])/tem00_fits[i,0,2])**2)
    car2 = tem00_fits[i,1,0]/(1.0+((full_volts-tem00_fits[i,1,1])/tem00_fits[i,1,2])**2)

    tem_model = sb1 + sb2 + sb3 + sb4 + car1 + car2

    tem_model += lsb_fits[i,0,0]/(1.0+((full_volts-lsb_fits[i,0,1])/lsb_fits[i,0,2])**2)
    tem_model += lsb_fits[i,1,0]/(1.0+((full_volts-lsb_fits[i,1,1])/lsb_fits[i,1,2])**2)
    tem_model += lsb_fits[i,2,0]/(1.0+((full_volts-lsb_fits[i,2,1])/lsb_fits[i,2,2])**2)
    tem_model += lsb_fits[i,3,0]/(1.0+((full_volts-lsb_fits[i,3,1])/lsb_fits[i,3,2])**2)
    tem_model += lsb_fits[i,4,0]/(1.0+((full_volts-lsb_fits[i,4,1])/lsb_fits[i,4,2])**2)

    tem_model += usb_fits[i,0,0]/(1.0+((full_volts-usb_fits[i,0,1])/usb_fits[i,0,2])**2)
    tem_model += usb_fits[i,1,0]/(1.0+((full_volts-usb_fits[i,1,1])/usb_fits[i,1,2])**2)

    fit_vals = invfreq_fitting[i,:]

    for j in range(len(usb9_locs)):
        freq = usb9_locs[j]
        usb9_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

        guess_idx = argmin(abs(pzt[sb45_idx[i,0]:sb45_idx[i,3]] - usb9_guess)) + sb45_idx[i,0] 

        r = 80
        peak_idx = signal.find_peaks_cwt(dcpd[guess_idx-r:guess_idx+r],array([r/8])) + guess_idx-r
        peak_heights = []
        for peak in peak_idx:
            peak_heights.append(dcpd[peak] - tem_model[peak-sb45_idx[i,0]])

        usb9_idx[i,j] = peak_idx[argmax(peak_heights)]

        s = 30
        volts = pzt[usb9_idx[i,j]-s:usb9_idx[i,j]+s]
        data = dcpd[usb9_idx[i,j]-s:usb9_idx[i,j]+s] - tem_model[usb9_idx[i,j]-sb45_idx[i,0]-s:usb9_idx[i,j]-sb45_idx[i,0]+s]
        A0 = dcpd[usb9_idx[i,j]] - tem_model[usb9_idx[i,j]-sb45_idx[i,0]]
        f0 = pzt[usb9_idx[i,j]]
        usb9_fit, usb9_errors = fit_peak(volts,data,A0,f0)
        usb9_fits[i,j,:] = usb9_fit






# Run through the sweeps and find the best fit

fit_chisq = zeros(a)

for f in range(a):

    # Grab data between the 45MHz sideband peaks, extending slightly more than 1 FSR
    idx0 = sb45_idx[f,0]-1600
    idx1 = sb45_idx[f,3]+1600

    # grab the fit parameters for this sweep
    fit_vals = freq_fitting[f,:]
    
    # grab the pzt values at the location of the big peaks (45MHz sidebands and carrier) - have to do this by hand because the format is different
    fit_volts1 = array([sb45_fits[f,0,1], tem00_fits[f,1,1], sb45_fits[f,1,1], sb45_fits[f,2,1], tem00_fits[f,0,1],sb45_fits[f,3,1]])

    # gather all the voltages together
    fit_volts = hstack((fit_volts1, lsb_fits[f,:,1], usb_fits[f,:,1], hom_fits[f,:,1], ehom_fits[f,:,1], lsb9_fits[f,:,1], usb9_fits[f,:,1]))

    # set the frequency for the big peaks
    fit_freqs1 = array([FSR*1e6+f2, FSR*1e6, FSR*1e6-f2, f2, 0.0, -1*f2])

    # set the frequency for the HOMs
    HOM_freqs = remainder(arange(1,Hnum+1)*HOMs*1e6,FSR*1e6)

    # gather all the frequencies together
    fit_freqs = hstack((fit_freqs1, lsb_locs, usb_locs, HOM_freqs, eHnum_locs, lsb9_locs, usb9_locs))

    # convert the peak voltages into frequency
    v = fit_volts
    freq_result = fit_vals[0] + fit_vals[1]*v + fit_vals[2]*v**2 + fit_vals[3]*v**3

    # compare fit to prediction
    fit_chisq[f] = sum( (fit_freqs*1e-6 - freq_result*1e-6)**2 )

    print f, fit_chisq[f]


f = argmin(fit_chisq)

print
print f

fignum=0

fignum=fignum+1
pylab.figure(fignum)

gs = gridspec.GridSpec(2,1,height_ratios=[3,1])

pylab.subplot(gs[0])

idx0 = sb45_idx[f,0]-1600
idx1 = sb45_idx[f,3]+1600

v = pzt[idx0:idx1]

fit_vals = freq_fitting[f,:]
freq_result = fit_vals[0] + fit_vals[1]*v + fit_vals[2]*v**2 + fit_vals[3]*v**3

fit_volts1 = array([sb45_fits[f,0,1], tem00_fits[f,1,1], sb45_fits[f,1,1], sb45_fits[f,2,1], tem00_fits[f,0,1],sb45_fits[f,3,1]])
fit_volts = hstack((fit_volts1, lsb_fits[f,:,1], usb_fits[f,:,1], hom_fits[f,:,1], ehom_fits[f,:,1], lsb9_fits[f,:,1], usb9_fits[f,:,1]))
#fit_volts = hstack((fit_volts1, lsb_fits[f,:,1], usb_fits[f,:,1], hom_fits[f,:,1], ehom_fits[f,:,1], usb9_fits[f,:,1]))

fit_freqs1 = array([FSR*1e6+f2, FSR*1e6, FSR*1e6-f2, f2, 0.0, -1*f2])
HOM_freqs = remainder(arange(1,Hnum+1)*HOMs*1e6,FSR*1e6)

fit_freqs = hstack((fit_freqs1, lsb_locs, usb_locs, HOM_freqs, eHnum_locs, lsb9_locs, usb9_locs))
#fit_freqs = hstack((fit_freqs1, lsb_locs, usb_locs, HOM_freqs, eHnum_locs, usb9_locs))

#print fit_freqs

#print lsb9_fits[f,:,1]

pylab.plot(v,freq_result*1e-6,'r-')
pylab.plot(fit_volts,fit_freqs*1e-6,'bo')

pylab.grid(True, which='both', linestyle=':', alpha=0.4)
pylab.ylabel('Optical Frequency [MHz]')
pylab.xlim(xmin,xmax)
pylab.title('OMC Mode Scan - Frequency Calibration - DARM Offset ' + data_string,fontsize=14)

pylab.subplot(gs[1])

v = fit_volts
freq_result = fit_vals[0] + fit_vals[1]*v + fit_vals[2]*v**2 + fit_vals[3]*v**3

pylab.plot(v, fit_freqs*1e-6 - freq_result*1e-6,'bo')

pylab.xlim(xmin,xmax)
#pylab.ylim(-2,2)
pylab.grid(True, linestyle=':', alpha=0.4)
pylab.xlabel('PZT2 Output [V]')
pylab.ylabel('residuals [MHz]')

pylab.savefig('full_freq_fit.pdf')



"""
lsb_text = ['LSB1', 'LSB2', 'LSB1', 'LSB3', 'LSB4','LSB5','LSB5','LSB6','LSB6','LSB7','LSB8']
lsb_locs = array([HOMs*1e6-f2, 2*HOMs*1e6-f2, HOMs*1e6+FSR*1e6-f2, 3*HOMs*1e6-f2, 4*HOMs*1e6-f2, 5*HOMs*1e6-f2, 5*HOMs*1e6-FSR*1e6-f2, 6*HOMs*1e6-f2, 6*HOMs*1e6-FSR*1e6-f2, 7*HOMs*1e6-FSR*1e6-f2,8*HOMs*1e6-FSR*1e6-f2])

lsb9_text = ['lsb0','lsb1', 'lsb2', 'lsb1', 'lsb3', 'lsb4','lsb5','lsb5','lsb6','lsb6','lsb7','lsb8']
lsb9_locs = array([-1*f1,HOMs*1e6-f1, 2*HOMs*1e6-f1, HOMs*1e6+FSR*1e6-f1, 3*HOMs*1e6-f1, 4*HOMs*1e6-f1, 5*HOMs*1e6-f1, 5*HOMs*1e6-FSR*1e6-f1, 6*HOMs*1e6-f1, 6*HOMs*1e6-FSR*1e6-f1, 7*HOMs*1e6-FSR*1e6-f1,8*HOMs*1e6-FSR*1e6-f1])

usb9_text = ['usb0','usb1', 'usb2', 'usb3','usb3','usb4','usb4','usb5','usb6','usb7','usb8','usb8']
usb9_locs = array([f1, HOMs*1e6+f1, 2*HOMs*1e6+f1, 3*HOMs*1e6+f1, 3*HOMs*1e6+f1-FSR*1e6, 4*HOMs*1e6+f1, 4*HOMs*1e6+f1-FSR*1e6, 5*HOMs*1e6+f1-FSR*1e6, 6*HOMs*1e6+f1-FSR*1e6, 7*HOMs*1e6+f1-FSR*1e6,8*HOMs*1e6+f1-FSR*1e6,8*HOMs*1e6+f1-2*FSR*1e6])

usb_text = ['USB1', 'USB2', 'USB3','USB3','USB4','USB4','USB5','USB6','USB7','USB8','USB8']
usb_locs = array([HOMs*1e6+f2, 2*HOMs*1e6+f2, 3*HOMs*1e6+f2, 3*HOMs*1e6+f2-FSR*1e6, 4*HOMs*1e6+f2, 4*HOMs*1e6+f2-FSR*1e6, 5*HOMs*1e6+f2-FSR*1e6, 6*HOMs*1e6+f2-FSR*1e6, 7*HOMs*1e6+f2-FSR*1e6,8*HOMs*1e6+f2-FSR*1e6,8*HOMs*1e6+f2-2*FSR*1e6])

mhom_text = ['CR10','CR11','CR12','CR14','CR15','CR14','CR13','CR16','CR17']
mhom_locs = array([remainder(10*HOMs*1e6,FSR*1e6), remainder(11*HOMs*1e6,FSR*1e6), remainder(12*HOMs*1e6,FSR*1e6), remainder(14*HOMs*1e6,FSR*1e6), remainder(15*HOMs*1e6,FSR*1e6), remainder(14*HOMs*1e6,FSR*1e6)+FSR*1e6, remainder(13*HOMs*1e6,FSR*1e6), remainder(16*HOMs*1e6,FSR*1e6), remainder(17*HOMs*1e6,FSR*1e6)])
"""


fignum=fignum+1
pylab.figure(fignum)

idx0 = sb45_idx[f,0]-1600
idx1 = sb45_idx[f,3]+1600

volts = pzt[idx0:idx1]
power = dcpd[idx0:idx1]

model = zeros(shape(volts))

fit_vals = invfreq_fitting[f,:]

#print 'Carrier peaks'
for i in range(2):

    pfit = tem00_fits[f,i,0]/(1.0+((volts-tem00_fits[f,i,1])/tem00_fits[f,i,2])**2)
    model += pfit
    pylab.plot([tem00_fits[f,i,1], tem00_fits[f,i,1]],[ymin,ymax],'k-',linewidth=0.4)
    pylab.annotate('CR0', xy = (tem00_fits[f,i,1], ymax), xytext = (1, -4), textcoords = 'offset points', ha = 'left', va = 'top',fontsize=10)

    #print 'CR0', tem00_fits[f,i,2]

#print
sb_text = ['USB0','LSB0','USB0','LSB0']
freq = array([FSR*1e6+f2, FSR*1e6-f2, f2, -1*f2])
for i in range(4):
    sb_guess = fit_vals[0] + fit_vals[1]*freq[i] + fit_vals[2]*freq[i]**2 + fit_vals[3]*freq[i]**3
    pfit = sb45_fits[f,i,0]/(1.0+((volts-sb45_fits[f,i,1])/sb45_fits[f,i,2])**2)
    model += pfit

    #print sb_text[i], sb45_fits[f,i,2]
    
    pylab.plot([sb_guess, sb_guess],[ymin,ymax],'r-',linewidth=0.4)
    pylab.annotate(sb_text[i], xy = (sb_guess, ymax), xytext = (1, -10), textcoords = 'offset points', ha = 'left', va = 'top',fontsize=10)

#print
for i in range(len(lsb_locs)):
    freq = lsb_locs[i]
    lsb_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

    pfit = lsb_fits[f,i,0]/(1.0+((volts-lsb_fits[f,i,1])/lsb_fits[f,i,2])**2)
    model += pfit

    #print lsb_text[i], lsb_fits[f,i,2]

    pylab.plot([lsb_guess, lsb_guess],[ymin,5],'g-',linewidth=0.4)
    pylab.annotate(lsb_text[i], xy = (lsb_guess, 5), xytext = (1, 0), textcoords = 'offset points', ha = 'left', va = 'bottom',fontsize=7)

#print
for i in range(len(lsb9_locs)):
    freq = lsb9_locs[i]
    lsb9_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

    pfit = lsb9_fits[f,i,0]/(1.0+((volts-lsb9_fits[f,i,1])/lsb9_fits[f,i,2])**2)
    model += pfit

    #print lsb9_text[i], lsb9_fits[f,i,0]

    pylab.plot([lsb9_guess, lsb9_guess],[ymin,10],'g-',linewidth=0.4)
    pylab.annotate(lsb9_text[i], xy = (lsb9_guess, 10), xytext = (1, 0), textcoords = 'offset points', ha = 'left', va = 'bottom',fontsize=7)

#print
for i in range(len(usb9_locs)):
    freq = usb9_locs[i]
    usb9_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

    pfit = usb9_fits[f,i,0]/(1.0+((volts-usb9_fits[f,i,1])/usb9_fits[f,i,2])**2)
    model += pfit

    #print usb9_text[i], usb9_fits[f,i,2]

    pylab.plot([usb9_guess, usb9_guess],[ymin,15],'g-',linewidth=0.4)
    pylab.annotate(usb9_text[i], xy = (usb9_guess, 15), xytext = (1, 0), textcoords = 'offset points', ha = 'left', va = 'bottom',fontsize=7)

#print
for i in range(len(usb_locs)):
    freq = usb_locs[i]
    usb_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

    pfit = usb_fits[f,i,0]/(1.0+((volts-usb_fits[f,i,1])/usb_fits[f,i,2])**2)
    model += pfit

    #print usb_text[i], usb_fits[f,i,2]

    pylab.plot([usb_guess, usb_guess],[ymin,3],'m-',linewidth=0.4)
    pylab.annotate(usb_text[i], xy = (usb_guess, 5), xytext = (1, -10), textcoords = 'offset points', ha = 'left', va = 'top',fontsize=7)

print
for i in range(Hnum):
    freq = remainder((i+1)*HOMs*1e6,FSR*1e6)
    hom_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

    pfit = hom_fits[f,i,0]/(1.0+((volts-hom_fits[f,i,1])/hom_fits[f,i,2])**2)
    model += pfit

    #print hom_text[i], hom_fits[f,i,2]

    pylab.plot([hom_guess, hom_guess],[ymin,2],'b-',linewidth=0.4)
    pylab.annotate(hom_text[i], xy = (hom_guess, 5), xytext = (1, -22), textcoords = 'offset points', ha = 'left', va = 'top',fontsize=7)

print
for i in range(len(eHnum_locs)):
    freq = eHnum_locs[i]
    ehom_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

    pfit = ehom_fits[f,i,0]/(1.0+((volts-ehom_fits[f,i,1])/ehom_fits[f,i,2])**2)
    model += pfit

    #print ehom_text[i], ehom_fits[f,i,0]

    pylab.plot([ehom_guess, ehom_guess],[ymin,2],'b-',linewidth=0.4)
    pylab.annotate(ehom_text[i], xy = (ehom_guess, 5), xytext = (1, -22), textcoords = 'offset points', ha = 'left', va = 'top',fontsize=7)

"""
for i in range(len(mhom_locs)):
    freq = mhom_locs[i]
    mhom_guess = fit_vals[0] + fit_vals[1]*freq + fit_vals[2]*freq**2 + fit_vals[3]*freq**3

    pylab.plot([mhom_guess, mhom_guess],[ymin,1.4],'b-',linewidth=0.4)
    pylab.annotate(mhom_text[i], xy = (mhom_guess, 5), xytext = (1, -34), textcoords = 'offset points', ha = 'left', va = 'top',fontsize=7)
"""

pylab.semilogy(pzt[idx0:idx1],dcpd[idx0:idx1],'k',linestyle='-',linewidth = 0.8)
pylab.semilogy(volts,model,'orange',linestyle='-',linewidth = 1.8,alpha=0.8)

pylab.grid(True, which='both', linestyle=':', alpha=0.4)

pylab.ylabel('OMC DCPD Sum [mA]')
pylab.xlabel('PZT2 Output [V]')

pylab.xlim(xmin,xmax)
pylab.ylim(ymin,ymax)

#pylab.xlim(35,38)

pylab.ylim(1e-2,3e1)

pylab.title('OMC Mode Scan - Full Lock - DARM Offset ' + data_string,fontsize=14)

pylab.savefig('full_lock_mode_fit_offset_' + data_string + '.pdf')


