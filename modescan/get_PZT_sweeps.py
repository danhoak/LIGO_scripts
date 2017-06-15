#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Code that scans a data file and finds indices where PZT voltage was ramping
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from scipy import signal


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


# we will decimate the data before finding the positive-voltage sweep segments
# assume that the ramp rate is less than 1V/sec (?)
# in this case probably 16Hz sample rate is fine

# expect ramps that go from less than 5V to more than 95V

# input is a time vector (tscan) and a pzt vector, with the same sample rate (assume this is the full sample rate)
# output is a 2D array of indices that correspond to the start and stop of a ramp, with some padding
# hard ceiling at 95V on the PZT
def get_PZT_sweeps(tscan,pzt):

    dt = tscan[1] - tscan[0]
    fs = 1/dt

    ds_factor = int(fs/16)

    times1 = decimate(tscan,ds_factor)
    pzt1 = decimate(pzt,ds_factor)

    # Now find the sections of data where the voltage was increasing
    positive_segs = []
    start_flag = False
    num_scans = 0
    for i in range(1,len(pzt1)):

        # start a segment if we aren't in a segment, and the slope of the pzt is positive, and the voltage IS greater than...5V?
        if not(start_flag) and (pzt1[i]-pzt1[i-1])>0 and pzt1[i] > 5:
            start_flag = True
            start_idx = i

        # stop a segment if we aren't in a segment, and the slope of the pzt is positive, and the voltage WAS greater than...85V?
        # bail out if we reach 95V
        #elif start_flag and ( ( (pzt1[i]-pzt1[i-1])<0 and pzt1[i-1] > 85 ) or pzt1[i] > 95 ):
        elif start_flag and ( ( (pzt1[i]-pzt1[i-1])<0 and pzt1[i-1] > 65 ) or pzt1[i] > 70 ):
            start_flag = False
            positive_segs.append([start_idx, i])
            num_scans += 1


    ramp_indices = []
    for seg in positive_segs:

        seg_length = seg[1]-seg[0]

        # If the ramp length is less than 16Hz and 2V per second for 60 volts, return an error
        if seg_length < 16*60*0.5:
            #print 'Scan segment is too short!'
            continue
        else:
            # get the start times of the ramp
            seg_start = argmin(abs(tscan - times1[seg[0]]))
            seg_stop = argmin(abs(tscan - times1[seg[1]]))

            ramp_indices.append([seg_start, seg_stop])


    return ramp_indices
