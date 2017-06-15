#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to perform some analysis of mode fitting results 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from scipy import signal
from scipy.special import jn
from uncertainties import ufloat
from uncertainties import unumpy


# Calculates contrast defect
# input: photocurrent in carrier HOMs (in mA), with errors
# output: contrast defect
def contrast_defect(pcar):

    # These numbers and their uncertainty may change
    pin = ufloat(2.8,0.1)
    #pin = ufloat(23.1,0.1)  # high power

    J9 = ufloat(0.9909,0.0007)  # for Gamma9 = 0.191 +/- 0.005
    J45 = ufloat(0.9799,0.0007) # for Gamma45 = 0.284 +/- 0.005

    eIO = ufloat(0.88,0.02)
    #gcr = ufloat(33,2)
    gcr = ufloat(38,2)  # better alignment
    tSRM = ufloat(0.37,0.001)

    pAS = pin * J9**2 * J45**2 * eIO * gcr * tSRM

    tOFI = ufloat(0.95,0.02)
    #tOM1 = ufloat(0.95,0.001)
    tOM1 = ufloat(0.9992,20e-6)   # new OM1
    tOM3 = ufloat(0.99,0.002)
    PDq = ufloat(0.75,0.02)

    ploss = tOFI * tOM1 * tOM3 * PDq

    return (pcar / (pAS * ploss)) * 1e6


def calc_modescan(mode_fits, mode_errors, volts_to_MHz, mode_freqs, mode_text, FSR):

    fit_frequencies = volts_to_MHz[0] + volts_to_MHz[1]*mode_fits[:,1] + volts_to_MHz[2]*mode_fits[:,1]**2 + volts_to_MHz[3]*mode_fits[:,1]**3
    frequency_errors = mode_freqs - fit_frequencies

    # We calculate the fit errors on mode frequency using the linear portion of the PZT calibration
    # It's pretty linear to begin with, so this should be ok
    # This is clearly a wild underestimate, the peak positions can vary from sweep to seep, but the fits are very precise
    mode_error_frequency = abs(mode_errors[:,1]*volts_to_MHz[1])
    freq_error_chisq = sum((frequency_errors)**2 / mode_error_frequency**2)
    print 'Frequency error chisq/NDOF:', freq_error_chisq/len(mode_text)

    pcar = 0
    USB = 0
    LSB = 0
    usb = 0
    lsb = 0
    USBn = 0
    LSBn = 0
    for i in range(len(mode_text)):

        # use the first CR0 mode (at f=0), not the second (at f=FSR)
        if mode_text[i]=='CR0' and mode_freqs[i]/1e6<1:
            CR0 = ufloat(mode_fits[i,0],mode_errors[i,0])

        if mode_text[i]=='CR1':
            CR1 = ufloat(mode_fits[i,0],mode_errors[i,0])

        if mode_text[i]=='CR2':
            CR2 = ufloat(mode_fits[i,0],mode_errors[i,0])

        # use the 45MHz peaks internal to the span of the FSR
        if mode_text[i]=='USB0' and mode_freqs[i]<FSR and mode_freqs[i]/1e6>1:
            #print mode_text[i], mode_fits[i,0], mode_freqs[i]/1e6
            USB += ufloat(mode_fits[i,0],mode_errors[i,0])

        if mode_text[i]=='LSB0' and mode_freqs[i]<FSR and mode_freqs[i]/1e6>1:
            #print mode_text[i], mode_fits[i,0], mode_freqs[i]/1e6
            LSB += ufloat(mode_fits[i,0],mode_errors[i,0])            

        # use the 9MHz peaks internal to the span of the FSR
        if mode_text[i]=='usb0' and mode_freqs[i]<FSR and mode_freqs[i]/1e6>1:
            #print mode_text[i], mode_fits[i,0], mode_freqs[i]/1e6
            usb += ufloat(mode_fits[i,0],mode_errors[i,0])

        if mode_text[i]=='lsb0' and mode_freqs[i]<FSR and mode_freqs[i]/1e6>1:
            #print mode_text[i], mode_fits[i,0], mode_freqs[i]/1e6
            lsb += ufloat(mode_fits[i,0],mode_errors[i,0])      
            
        # Sum of carrier higher order modes within the FSR for contrast defect calculation
        if mode_text[i][0:2]=='CR' and mode_freqs[i]<(FSR-1) and mode_freqs[i]/1e6>1:
            #print mode_text[i], mode_fits[i,0]
            pcar += ufloat(mode_fits[i,0],mode_errors[i,0])

        # Sum of 45MHz upper sideband higher order modes
        if mode_text[i][0:3]=='USB' and int(mode_text[i][3])>0:# and mode_freqs[i]/1e6<260 and mode_freqs[i]/1e6>50:
            USBn += ufloat(mode_fits[i,0],mode_errors[i,0])

        # Sum of 45MHz upper sideband higher order modes
        if mode_text[i][0:3]=='LSB' and int(mode_text[i][3])>0:#mode_freqs[i]/1e6<(260-46) and mode_freqs[i]/1e6>0:
            LSBn += ufloat(mode_fits[i,0],mode_errors[i,0])


    # calculate contrast defect using carrier higher order modes
    # this should be an overestimate since it includes some misalignment & imperfect mode-matching
    CD = contrast_defect(pcar*1e-3)

    # carrier mode-matching and alignment are only interesting when DARM offset is on
    CR_MM = 1 - CR2/CR0
    CR_align = 1 - CR1/CR0  # is this correct??

    SB45_bal = USB/LSB
    SB9_bal = usb/lsb

    #pin = ufloat(2.8,0.1)
    #Gamma45 = 0.284
    #SB45_input = pin*jn(0,Gamma45)*jn(1,Gamma45)
    #SB45_trans = (USB+LSB)/SB45_input


    print 'Contrast defect is:', CD
    print 'Carrier mode matching is:', CR_MM
    print 'Carrier alignment is:', CR_align
    print '45MHz power is:', USB+LSB+USBn+LSBn
    print '45MHz sideband balance is:', SB45_bal, USB, LSB
    print '45MHz HOM content is:', (LSBn+USBn)/(USB+LSB)
    print '9MHz sideband balance is:', SB9_bal, usb, lsb
    #print 'USBn:', USBn, 'LSBn:', LSBn
    #print '45MHz sideband transmission is:', SB45_trans


    return fit_frequencies, unumpy.nominal_values(CD)
