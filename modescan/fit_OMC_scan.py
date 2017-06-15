#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Code to fit peaks in OMC full-lock modescan
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from scipy import signal
from lmfit import minimize, Parameters, fit_report

def lorentzian_error_function(params,x,data):

    A = params['A'].value
    f0 = params['f0'].value
    fc = params['fc'].value
    fit = A/(1.0+((x-f0)/fc)**2)
    return data - fit

def fit_peak(volts,data,A0,f0):

    params = Parameters()
    params.add('A', A0, vary=True)
    params.add('f0', f0, vary=True)
    params.add('fc', 0.05, vary=True)

    out = minimize(lorentzian_error_function, params, args=(volts,data), method='nelder')
    out = minimize(lorentzian_error_function, params, args=(volts,data), method='leastsq')

    fit_vals = [params['A'].value, params['f0'].value, params['fc'].value]
    fit_errors = [params['A'].stderr, params['f0'].stderr, params['fc'].stderr]

    return fit_vals, fit_errors


def poly_error_function(params,x,y):

    a0 = params['a0'].value
    a1 = params['a1'].value
    a2 = params['a2'].value
    a3 = params['a3'].value
    a4 = params['a4'].value
    a5 = params['a5'].value
    fit = a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 + a5*x**5
    return y - fit


def fit_pzt_volts(volts,freqs,FSR_slope):

    params = Parameters()
    params.add('a0', 350.0e6, vary=True)
    params.add('a1', -1/FSR_slope, vary=True)
    params.add('a2', 0.0, vary=True)
    params.add('a3', 0.0, vary=True)
    params.add('a4', 0.0, vary=False)
    params.add('a5', 0.0, vary=False)

    # Note order of arguments: we're fitting frequencies to volts.  This is used for predicting the location of peaks.
    out = minimize(poly_error_function, params, args=(volts,freqs), method='nelder')
    out = minimize(poly_error_function, params, args=(volts,freqs), method='leastsq')

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

    # Note order of arguments: we're fitting volts to frequencies.  This is used for converting voltage to optical frequency.
    out = minimize(poly_error_function, params, args=(freqs,volts), method='nelder')
    out = minimize(poly_error_function, params, args=(freqs,volts), method='leastsq')

    fit_vals = [params['a0'].value, params['a1'].value, params['a2'].value, params['a3'].value]
    fit_errors = [params['a0'].stderr, params['a1'].stderr, params['a2'].stderr, params['a3'].stderr]

    return fit_vals, fit_errors


# Function that does the peak identification and fitting
# Inputs are dcpd data, pzt data, an array of mode optical frequencies to fit, fitting parameters to 
# convert mode freq to pzt voltage, two scalars that define the width of data to fit, and the sum of previous fits
def mode_fitting(dcpd, pzt, mode_locs, freq_fit, r, s, sum_model):

    mode_fits = zeros((len(mode_locs),3))
    mode_errors = zeros((len(mode_locs),3))

    for j in range(len(mode_locs)):
        freq = mode_locs[j]
        mode_guess = freq_fit[0] + freq_fit[1]*freq + freq_fit[2]*freq**2 + freq_fit[3]*freq**3
        gidx = argmin(abs(pzt - mode_guess))

        dcpd_new = dcpd-sum_model

        peak_idx = signal.find_peaks_cwt(dcpd_new[gidx-r:gidx+r],array([r/5])) + gidx-r

        # Associate the highest peak in the [-r,+r] span as the mode we're looking for
        peak_heights = []
        for peak in peak_idx:
            peak_heights.append(dcpd_new[peak])
        pidx = peak_idx[argmax(peak_heights)]

        mode_fits[j,:], mode_errors[j,:] = fit_peak(pzt[pidx-s:pidx+s],dcpd_new[pidx-s:pidx+s],dcpd_new[pidx],pzt[pidx])

    return mode_fits, mode_errors




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


# inputs are pzt data, dcpd data, and sample rates for both
# outputs are peak parameters and pzt to optical frequency calibration parameters
def fit_OMC_scan(pzt_raw, dcpd, pzt_fs, dcpd_fs, FSR, HOMs, f1, f2):

    # if the sample rates are different, interpolate the pzt data up to the dcpd sample rate
    # deal with this later
    pzt = pzt_raw

    seg_length = len(dcpd)

    # peaks are about 0.1V wide - use this to set a fundamental length for peak fitting
    scan_velocity = len(pzt)/(max(pzt)-min(pzt))  # samples per volt
    peak_width = int(0.1*scan_velocity)

    r = 2*peak_width  # width for finding peaks
    s = peak_width    # width for fitting peaks

    # first job is to identify the 45MHz sideband peaks
    # set a threshold on how high we expect them to be
    # we expect they are the tallest thing in the scan - if Gamma45 changes this may no longer be true
    peak_threshold = max(dcpd)*0.8

    seg_peaks = []

    # decimate dcpd before searching for sidebands, otherwise it takes too long
    dcpd_short = decimate(dcpd,4)
    pzt_short = decimate(pzt,4)

    short_length = len(pzt_short)

    # find peaks with characteristic width around 3% of the scan
    peak_idx = signal.find_peaks_cwt(dcpd_short,array([int(floor(short_length*0.02))]))

    # Now go back to the full-sample-rate data and find the max value around each peak
    for jj in range(len(peak_idx)):

        pidx = argmin(abs(pzt-pzt_short[peak_idx[jj]]))

        peakmax = argmax(dcpd[pidx-r:pidx+r])

        #print pzt[pidx-r+peakmax], dcpd[pidx-r+peakmax]

        if dcpd[pidx-r+peakmax] > peak_threshold:
            seg_peaks.append(pidx-r+peakmax)

    # Now we should have an array of peak indices - these should include two pairs of sideband peaks
    sb45_idx = zeros(4)
    if len(seg_peaks)>=4:

        # We are looking for the pairs of the 45MHz sidebands.  These should be tall, similar heights, and separated by ~12V
        group_count = 0
        peak_group = []
        for j in range(len(seg_peaks)):

            idx1 = seg_peaks[j]
            for k in range(j+1,len(seg_peaks)):

                idx2 = seg_peaks[k]

                if (pzt[idx2]-pzt[idx1])>10 and (pzt[idx2]-pzt[idx1])<15 and (dcpd[idx1]/dcpd[idx2])<1.1 and (dcpd[idx1]/dcpd[idx2])>0.9:

                    # We have found a pair of sidebands
                    peak_group.append([idx1,idx2])
                    group_count+=1
                    
        # If the sweep includes two FSRs, use the lower one, to keep things consistent.
        if group_count==2 or group_count==3:
            sb45_idx = array([peak_group[0][0], peak_group[0][1], peak_group[1][0], peak_group[1][1]])
            #sb45_idx = array([peak_group[1][0], peak_group[1][1], peak_group[2][0], peak_group[2][1]])
        else:
            print 'Did not find a set of sideband pairs!'
            print 'Number of sideband pairs:', group_count

    else:
        print 'Not enough peaks! Cannot find sidebands.'



    print
    print 'Starting fits.'

    print '45MHz Sidebands'

    sb45_text = ['USB0','LSB0','USB0','LSB0']
    sb45_fits = zeros((4,3))
    sb45_errors = zeros((4,3))

    for i in range(4):

        # The DCPD sum saturates at about 40mA.  The 45MHz sidebands exceed this when the IFO is locked at high power.
        # Need to work around the side of the peak to avoid the flat top.

        # First, grab a chunk of data around the peak
        volts_full = pzt[sb45_idx[i]-1.3*r:sb45_idx[i]+1.3*r]
        data_full = dcpd[sb45_idx[i]-1.3*r:sb45_idx[i]+1.3*r]

        # If the peak is saturated, only use data less than 30mA photocurrent
        if dcpd[sb45_idx[i]] > 30:
            sidx = where(data_full < 30)[0]
            volts = volts_full[sidx]
            data = data_full[sidx]
        else:
            volts = volts_full
            data = data_full            
        print data

        sb45_fits[i,:], sb45_errors[i,:] = fit_peak(volts,data,dcpd[sb45_idx[i]],pzt[sb45_idx[i]])

    # Get a preliminary frequency calibration in terms of the PZT voltage
    FSR_calib1 = (sb45_fits[3,1] - sb45_fits[1,1])/FSR
    #FSR_calib1 = (sb45_fits[2,1] - sb45_fits[0,1])/FSR
    
    # add up 45MHz sideband peaks to subtract out later
    sum_model = zeros(shape(pzt))
    for i in range(4):
        sum_model += sb45_fits[i,0]/(1.0+((pzt-sb45_fits[i,1])/sb45_fits[i,2])**2)


    print 'Carrier TEM00 Modes'

    tem00_text = ['CR0','CR0']
    tem00_fits = zeros((2,3))
    tem00_errors = zeros((2,3))
    for j in range(2):

        # Use the location of the 45MHz sidebands and the preliminary frequency
        # calibration to make a guess for the location of the carrier
        CR_freq = (j*FSR + f2) # this is the difference between the CR peak and the second LSB45
        tem_guess = pzt[sb45_idx[3]] - CR_freq*FSR_calib1
        gidx = argmin(abs(pzt - tem_guess))

        # since the calibration is pretty rough at this point, broaden the search for peaks
        # is this going to screw things up when the DARM offet is small?? --> yes, need to swap commented code below
        w = 3*r

        # find the peaks near the guess
        # sometimes find_peaks_cwt doesn't return the maximum, make sure we find it
        peak_idx = signal.find_peaks_cwt(dcpd[gidx-w:gidx+w],array([w/5])) + gidx-w
        peak_heights = []
        peak_maxes = []
        for peak in peak_idx:
            peak_max = argmax(dcpd[peak-r:peak+r])+(peak-r)
            peak_heights.append(dcpd[peak_max] - sum_model[peak_max])
            peak_maxes.append(peak_max)

            #peak_maxes.append(peak)
            #peak_heights.append(dcpd[peak] - sum_model[peak])

        pidx = peak_maxes[argmax(peak_heights)]
        volts = pzt[pidx-s:pidx+s]
        data = (dcpd - sum_model)[pidx-s:pidx+s]
        tem00_fits[j,:], tem00_errors[j,:] = fit_peak(volts,data,data[s],pzt[pidx])

    # add up CR0 peaks to subtract out later
    for i in range(2):
        sum_model += tem00_fits[i,0]/(1.0+((pzt-tem00_fits[i,1])/tem00_fits[i,2])**2)


    # fit pzt voltage to optical frequency with six data points
    fit_volts = hstack((tem00_fits[:,1], sb45_fits[:,1]))
    fit_freqs = array([0.0, FSR, FSR+f2, FSR-f2, f2, -1*f2])

    fit_vals, fit_errors = fit_pzt_volts(fit_volts,fit_freqs,FSR_calib1)
    freq_fitting = fit_vals[0:4]

    # fit_vals is used after this for MHz --> volts conversion
    fit_vals, fit_errors = fit_pzt_freq(fit_freqs,fit_volts,FSR_calib1)
    invfreq_fitting = fit_vals[0:4]


    print 'LSB 45MHz HOMs'

    lsb_text = ['LSB1', 'LSB2', 'LSB1', 'LSB3', 'LSB5', 'LSB5']
    lsb_locs = array([HOMs-f2, 2*HOMs-f2, HOMs+FSR-f2, 3*HOMs-f2, 5*HOMs-f2, 5*HOMs-FSR-f2])

    #lsb_text = ['LSB1', 'LSB2', 'LSB1', 'LSB3', 'LSB4']
    #lsb_locs = array([HOMs-f2, 2*HOMs-f2, HOMs+FSR-f2, 3*HOMs-f2, 4*HOMs-f2])

    lsb_fits, lsb_errors = mode_fitting(dcpd, pzt, lsb_locs, fit_vals, r, s, sum_model)

    for i in range(len(lsb_locs)):
        sum_model += lsb_fits[i,0]/(1.0+((pzt-lsb_fits[i,1])/lsb_fits[i,2])**2)


    print 'USB 45MHz HOMs'

    usb_text = ['USB1', 'USB2']
    usb_locs = array([HOMs+f2, 2*HOMs+f2])

    usb_fits, usb_errors = mode_fitting(dcpd, pzt, usb_locs, fit_vals, r, s, sum_model)

    for i in range(len(usb_locs)):
        sum_model += usb_fits[i,0]/(1.0+((pzt-usb_fits[i,1])/usb_fits[i,2])**2)



    print 'Carrier HOMs'

    hom_text = ['CR1','CR2','CR3','CR4','CR5','CR6','CR7','CR8','CR9','CR10','CR11','CR12','CR13','CR14','CR15']
    Hnum = 9
    hom_locs = []
    for j in range(Hnum):
        hom_locs.append(remainder((j+1)*HOMs,FSR))

    # Have to narrow the range to keep the CR1 fit from getting spoiled by junk nearby
    hom_fits, hom_errors = mode_fitting(dcpd, pzt, hom_locs, fit_vals, int(1*r), int(0.9*s), sum_model)

    #hom_fits, hom_errors = mode_fitting(dcpd, pzt, hom_locs, fit_vals, r, s, sum_model)

    for i in range(len(hom_locs)):
        sum_model += hom_fits[i,0]/(1.0+((pzt-hom_fits[i,1])/hom_fits[i,2])**2)



    # Update the fitting with the carrier HOMs
    # previous fitting arrays were [cr0, sb45]

    start_fit_freqs = fit_freqs
    fit_volts = hstack((fit_volts, hom_fits[:,1]))
    fit_freqs = hstack((fit_freqs, hom_locs))

    fit_vals, fit_errors = fit_pzt_volts(fit_volts,fit_freqs,FSR_calib1)
    freq_fitting = fit_vals[0:4]

    fit_vals, fit_errors = fit_pzt_freq(fit_freqs,fit_volts,FSR_calib1)
    invfreq_fitting = fit_vals[0:4]


    print 'Extra Carrier HOMs'

    #ehom_text = ['CR4','CR5','CR9','CR11','CR12']
    #ehom_locs = array([4*HOMs-FSR, 5*HOMs, 9*HOMs-2*FSR, 11*HOMs-2*FSR, 12*HOMs-2*FSR])

    #ehom_text = ['CR4','CR5','CR9','CR12']
    #ehom_locs = array([4*HOMs-FSR, 5*HOMs, 9*HOMs-2*FSR, 12*HOMs-2*FSR])

    #ehom_text = ['CR4','CR5','CR9']
    #ehom_locs = array([4*HOMs-FSR, 5*HOMs, 9*HOMs-2*FSR])

    ehom_text = ['CR4','CR5','CR9']
    ehom_locs = array([4*HOMs-FSR, 5*HOMs, 9*HOMs-2*FSR])

    #ehom_text = ['CR4','CR5']
    #ehom_locs = array([4*HOMs-FSR, 5*HOMs])

    ehom_fits, ehom_errors = mode_fitting(dcpd, pzt, ehom_locs, fit_vals, r, s, sum_model)

    for i in range(len(ehom_locs)):
        sum_model += ehom_fits[i,0]/(1.0+((pzt-ehom_fits[i,1])/ehom_fits[i,2])**2)



    print 'LSB 9MHz HOMs'

    lsb9_text = ['lsb0', 'lsb2', 'lsb3', 'lsb5', 'lsb5', 'lsb6', 'lsb7', 'lsb4', 'lsb4', 'lsb9', 'lsb9', 'lsb8']
    lsb9_locs = array([FSR-f1, 2*HOMs-f1, 3*HOMs-f1, 5*HOMs-f1, 5*HOMs-FSR-f1, 6*HOMs-FSR-f1, 7*HOMs-FSR-f1, 4*HOMs-FSR-f1, 4*HOMs-f1, 9*HOMs-FSR-f1, 
                       9*HOMs-2*FSR-f1, 8*HOMs-FSR-f1])

    #lsb9_text = ['lsb2', 'lsb3', 'lsb5', 'lsb5', 'lsb6', 'lsb7', 'lsb0', 'lsb4', 'lsb4']
    #lsb9_locs = array([2*HOMs-f1, 3*HOMs-f1, 5*HOMs-f1, 5*HOMs-FSR-f1, 6*HOMs-FSR-f1, 7*HOMs-FSR-f1, FSR-f1, 4*HOMs-FSR-f1, 4*HOMs-f1])

    lsb9_fits, lsb9_errors = mode_fitting(dcpd, pzt, lsb9_locs, fit_vals, r, s, sum_model)

    for i in range(len(lsb9_locs)):
        sum_model += lsb9_fits[i,0]/(1.0+((pzt-lsb9_fits[i,1])/lsb9_fits[i,2])**2)


    print 'USB 9MHz HOMs'

    usb9_text = ['usb0','usb0','usb1', 'usb2','usb3','usb4','usb4','usb5','usb5','usb6', 'usb7', 'usb8']
    usb9_locs = array([f1, FSR+f1,HOMs+f1, 2*HOMs+f1, 3*HOMs+f1, 4*HOMs+f1, 4*HOMs-FSR+f1, 5*HOMs+f1, 5*HOMs-FSR+f1, 6*HOMs-FSR+f1, 7*HOMs-FSR+f1, 
                       8*HOMs-FSR+f1])

    #usb9_text = ['usb0','usb0','usb1', 'usb2','usb3','usb4','usb4','usb5','usb5','usb6','usb7','usb8']
    #usb9_locs = array([f1, FSR+f1,HOMs+f1, 2*HOMs+f1, 3*HOMs+f1, 4*HOMs+f1, 4*HOMs-FSR+f1, 5*HOMs+f1, 5*HOMs-FSR+f1, 6*HOMs-FSR+f1, 7*HOMs-FSR+f1, 8*HOMs-FSR+f1])

    #usb9_text = ['usb0','usb0','usb1', 'usb2','usb3','usb4','usb4','usb5','usb5','usb6', 'usb7']
    #usb9_locs = array([f1, FSR+f1,HOMs+f1, 2*HOMs+f1, 3*HOMs+f1, 4*HOMs+f1, 4*HOMs-FSR+f1, 5*HOMs+f1, 5*HOMs-FSR+f1, 6*HOMs-FSR+f1, 7*HOMs-FSR+f1])

    usb9_fits, usb9_errors = mode_fitting(dcpd, pzt, usb9_locs, fit_vals, r, s, sum_model)

    for i in range(len(usb9_locs)):
        sum_model += usb9_fits[i,0]/(1.0+((pzt-usb9_fits[i,1])/usb9_fits[i,2])**2)






    print 'Fitting the PZT value to optical frequency'

    # Update the fitting with the latest modes that were identified
    # previous fitting arrays were [cr0, sb45, hom]

    fit_volts = hstack((fit_volts, ehom_fits[:,1], lsb9_fits[:,1], usb9_fits[:,1]))

    HOM_freqs = remainder(arange(1,Hnum+1)*HOMs,FSR)
    fit_freqs = hstack((fit_freqs, ehom_locs, lsb9_locs, usb9_locs))

    fit_vals, fit_errors = fit_pzt_volts(fit_volts,fit_freqs,FSR_calib1)
    freq_fitting = fit_vals[0:4]

    fit_vals, fit_errors = fit_pzt_freq(fit_freqs,fit_volts,FSR_calib1)
    invfreq_fitting = fit_vals[0:4]

    #print 'Inverse frequency fitting parameters (MHz --> Volts):'
    #print invfreq_fitting
    
    #print 'Frequency fitting parameters (Volts --> MHz):'
    #print freq_fitting

    # all output arrays are given in the following order:
    # TEM00, SB45, LSB, USB, HOM, eHOM, lsb9, usb9

    # results of the fits
    mode_fits = vstack((tem00_fits, sb45_fits, lsb_fits, usb_fits, hom_fits, ehom_fits, lsb9_fits, usb9_fits))
    mode_errors = vstack((tem00_errors, sb45_errors, lsb_errors, usb_errors, hom_errors, ehom_errors, lsb9_errors, usb9_errors))

    # predicted optical frequencies for the mode fits
    mode_freqs = hstack((start_fit_freqs, lsb_locs, usb_locs, hom_locs, ehom_locs, lsb9_locs, usb9_locs))

    mode_text = tem00_text + sb45_text + lsb_text + usb_text + hom_text[0:Hnum] + ehom_text + lsb9_text + usb9_text


    return mode_fits, mode_errors, freq_fitting, invfreq_fitting, mode_freqs, mode_text
