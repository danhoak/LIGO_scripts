#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This routine plots a spectrogram of the input data using numpy rfft (instead of the crappy psd or specgram)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np
import matplotlib
from matplotlib import mlab
matplotlib.use("Agg")
import pylab
from makewaveform import *
#from optimalSNR import *

Fs = 4096.0
duration = 3.0

t = np.arange(0,duration-1.0/Fs,1/Fs)

noise_amp = 7.0
noise = noise_amp*np.random.randn(len(t))

S_freqs = np.arange(0.0,Fs/2+1,1.0)
S = (2*noise_amp**2/Fs)*np.ones(np.shape(S_freqs))

# make an inspiral
t1, hp1, hc1 = makewaveform('inspiral','1.4~1.4~0~20',duration,duration/2.0,Fs)

A = 3.0
tau = 0.25
f = 500.0
params = str(A) + '~' + str(tau) + '~' + str(f) + '~0~0~0'
t2, hp2, hc2 = makewaveform('chirplet',params,duration,duration/2.0,Fs)

#[SNR, h_rss, h_peak, Fchar, bw, Tchar, dur] = optimalSNR(hp2,Fs,0,S,S_freqs,Fmin=0.0,Fmax=Fs/2-1)
#print SNR

A = 2.0
tau = 0.10
f = 1000.0
params = str(A) + '~' + str(tau) + '~' + str(f) + '~0~0~0'
t3, hp3, hc3 = makewaveform('chirplet',params,duration,duration/2.0,Fs)

#[SNR, h_rss, h_peak, Fchar, bw, Tchar, dur] = optimalSNR(hp3,Fs,0,S,S_freqs,Fmin=0.0,Fmax=Fs/2-1)
#print SNR

A = 1.0
tau = 0.01
f = 1660.0
params = str(A) + '~' + str(tau) + '~' + str(f) + '~0~0~0'
t4, hp4, hc4 = makewaveform('chirplet',params,duration,duration/2.0,Fs)

#[SNR, h_rss, h_peak, Fchar, bw, Tchar, dur] = optimalSNR(hp4,Fs,0,S,S_freqs,Fmin=0.0,Fmax=Fs/2-1)
#print SNR

y = noise + hp2 + hp3 + hp4

# 16 Hz resolution, 1/32 sec time steps

step_size = Fs/16.0
overlap = 0.5*step_size

nyquist = Fs/2
dt = overlap/Fs
df = Fs/step_size

Pxx, freq, ts = mlab.specgram(y, NFFT=int(step_size), Fs=int(Fs), noverlap=int(overlap))

Axx = sqrt(Pxx)

SNRxx = sqrt(Pxx*Fs/(noise_amp**2))

print amax(SNRxx), mean(SNRxx)

SNR_total = 0
pixels = 0
for row in SNRxx:
    for i in row:
        if i > 5.5:
            pixels += 1
            SNR_total += i
            #print i

print pixels, SNR_total

fixed_powspec = ma.fix_invalid(Pxx,fill_value = 0)
cluster_mask = ones(shape(Pxx))
cluster_mask_threshold = ma.masked_less(fixed_powspec.data,0.1)

# Find the indices where we want to set the mask to zero
mdx1 = where(cluster_mask_threshold.mask)

# Set the mask to zero at those places
cluster_mask[mdx1] = 0


fignum=0

fignum=fignum+1
pylab.figure(fignum)

#pylab.imshow(Pxx,aspect='auto',origin='lower',cmap=pylab.cm.jet, interpolation='nearest')
pylab.imshow(cluster_mask,aspect='auto',origin='lower',cmap=pylab.cm.gray_r, interpolation='nearest')

#pylab.colorbar()

ytick_locs = np.array(range(0,int(nyquist/df)+1,int(256/df)))
ytick_labels = freq[ytick_locs].astype(int)
pylab.yticks(ytick_locs, ytick_labels)

xtick_locs = np.array(range(-1,int(duration/dt)-1,int(0.5/dt)))
xtick_labels = ts[xtick_locs]
xtick_labels[0] = 0
pylab.xticks(xtick_locs, xtick_labels)

#pylab.ylim(0,40)
#pylab.xlim(120,200)

pylab.xlabel('Time (seconds)')
pylab.ylabel('Frequency (Hz)')
#pylab.title('Spectrogram')

pylab.savefig("SG_spectrogram.png")
