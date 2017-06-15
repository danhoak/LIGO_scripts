#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This program compares the chisq and Rayleigh distributions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
from scipy import stats

### Generate two zero mean unit vriance normal distributions
### root sum of squares will be Rayleigh distributed

mu = 0.0
sigma = 1.0

x = random.normal(mu,sigma,1000000)
y = random.normal(mu,sigma,1000000)
r = sqrt(x**2 + y**2)
print 'Mean, median of r = sqrt(x**2 + y**2):', mean(r), median(r)


### Generate two zero mean unit variance normal distributions
### root sum of squares should be rician distributed

mu = 2.0
sigma = 1.0

u = random.normal(mu,sigma,1000000)
v = random.normal(mu,sigma,1000000)
#rice = sqrt(u**2 + v**2)
rice = sqrt((u-mean(u))**2 + (v-mean(v))**2)
rice_corr = rice/median(rice)
print 'Mean, median of mean-subtracted r = sqrt(u**2 + v**2):', mean(rice), median(rice)

# Generate a normal distribution for comparison
norm = random.normal(mean(rice),std(rice),1000000)


### numerical sanity checks - how to normalize such that median == 1?

r_med = (r/median(r))
r_med_corr = (r/median(r))*sqrt(log(4))

p50 = percentile(r,50)
p_med_50 = percentile(r_med,50)
p_med_corr_50 = percentile(r_med_corr,50)

print
print 'Median of R = sqrt(X**2 + Y**2):', median(r)
#print '50th percentile of R:', p50
print 'Median of R/median(R):', p_med_50
#print 'Median of R / sqrt(ln(4)):', median(r) / sqrt(log(4))
#print 'Median of R/median(R) * sqrt(ln(4)):', p_med_corr_50

# Note that median of corrected data is not 1!  But this data matches the rayleigh distribution.



"""
p975 = percentile(r,97.5)
p_med_975 = percentile(r_med,97.5)
p_med_corr_975 = percentile(r_med_corr,97.5)

print p975, p_med_975, p_med_corr_975
print p975/p50, p_med_975/p_med_50, p_med_corr_975/p_med_corr_50
"""

### Generate histograms, get the bin centers so the plotting isn't screwed up

bins = arange(0,8,0.07)

# Histogram of r = sqrt(X**2 + Y**2)
r_hist, r_bins = histogram(r,bins,normed=True)
r_bin_center = range(len(r_bins)-1)
for i in range(len(r_bins)-1):
    r_bin_center[i] = r_bins[i]+(r_bins[i+1]-r_bins[i])/2.0

# Histogram of nonzero mean data
rice_hist, rice_bins = histogram(rice,bins,normed=True)
rice_bin_center = range(len(rice_bins)-1)
for i in range(len(rice_bins)-1):
    rice_bin_center[i] = rice_bins[i]+(rice_bins[i+1]-rice_bins[i])/2.0

# Histogram of nonzero mean data, normalized to median
rice_corr_hist, rice_corr_bins = histogram(rice_corr,bins,normed=True)
rice_corr_bin_center = range(len(rice_corr_bins)-1)
for i in range(len(rice_corr_bins)-1):
    rice_corr_bin_center[i] = rice_corr_bins[i]+(rice_corr_bins[i+1]-rice_corr_bins[i])/2.0

# Normally-distributed data with 
norm_hist, norm_bins = histogram(norm,bins,normed=True)
norm_bin_center = range(len(norm_bins)-1)
for i in range(len(norm_bins)-1):
    norm_bin_center[i] = norm_bins[i]+(norm_bins[i+1]-norm_bins[i])/2.0

# Histogram of r = sqrt(X**2 + Y**2), normalized by median
r_med_hist, r_med_bins = histogram(r_med,bins,normed=True)
r_med_bin_center = range(len(r_med_bins)-1)
for i in range(len(r_med_bins)-1):
    r_med_bin_center[i] = r_med_bins[i]+(r_med_bins[i+1]-r_med_bins[i])/2.0

# Histogram of r = sqrt(X**2 + Y**2), normalized by median, with correction factor for median bias
r_med_corr_hist, r_med_corr_bins = histogram(r_med_corr,bins,normed=True)
r_med_corr_bin_center = range(len(r_med_corr_bins)-1)
for i in range(len(r_med_corr_bins)-1):
    r_med_corr_bin_center[i] = r_med_corr_bins[i]+(r_med_corr_bins[i+1]-r_med_corr_bins[i])/2.0


### Now generate a rayleigh distribution for comparison

r_pdf = random.rayleigh(1.0,1000000)

rice_pdf = stats.rice.rvs(1.0,size=1000000)
rice_pdf_corr = rice_pdf / median(rice_pdf) * sqrt(log(4))
#chisq = random.chisquare(2,1000000)

print
print 'Median of Rice distribution with b=1.0:', median(rice_pdf)
print 'Median of normalized Rice distribution:', median(rice_pdf/median(rice_pdf))
print 'Median of normalized, corrected Rice distribution:', median(rice_pdf_corr)


r_pdf_hist, r_pdf_bins = histogram(r_pdf,bins,normed=True)
r_pdf_bin_center = range(len(r_pdf_bins)-1)
for i in range(len(r_pdf_bins)-1):
    r_pdf_bin_center[i] = r_pdf_bins[i]+(r_pdf_bins[i+1]-r_pdf_bins[i])/2.0

rice_pdf_hist, rice_pdf_bins = histogram(rice_pdf,bins,normed=True)
rice_pdf_bin_center = range(len(rice_pdf_bins)-1)
for i in range(len(rice_pdf_bins)-1):
    rice_pdf_bin_center[i] = rice_pdf_bins[i]+(rice_pdf_bins[i+1]-rice_pdf_bins[i])/2.0

rice_pdf_corr_hist, rice_pdf_corr_bins = histogram(rice_pdf_corr,bins,normed=True)
rice_pdf_corr_bin_center = range(len(rice_pdf_corr_bins)-1)
for i in range(len(rice_pdf_corr_bins)-1):
    rice_pdf_corr_bin_center[i] = rice_pdf_corr_bins[i]+(rice_pdf_corr_bins[i+1]-rice_pdf_corr_bins[i])/2.0

"""
chisq_hist, chisq_bins = histogram(chisq,bins,normed=True)
chisq_bin_center = range(len(chisq_bins)-1)
for i in range(len(chisq_bins)-1):
    chisq_bin_center[i] = chisq_bins[i]+(chisq_bins[i+1]-chisq_bins[i])/2.0
"""


fontP = FontProperties()
fontP.set_size('small')

fignum=0
fignum=fignum+1
pylab.figure(fignum)

pylab.plot(r_bin_center[0:len(r_hist)],r_hist,'r--',linewidth=3.0,label='R = sqrt[x^2 + y^2]')

pylab.plot(r_med_bin_center[0:len(r_med_hist)],r_med_hist,'b--',linewidth=1.0,label='R/median(R)')
pylab.plot(r_med_corr_bin_center[0:len(r_med_corr_hist)],r_med_corr_hist,'ko',markersize=4.0,label='R/median(R) * sqrt(ln(4))')

pylab.plot(r_pdf_bin_center[0:len(r_pdf_hist)],r_pdf_hist,'g-',linewidth=1.0,label='Rayleigh PDF')

pylab.plot(rice_pdf_bin_center[0:len(rice_pdf_hist)],rice_pdf_hist,'ro',markersize=4.0,label=r'Rice PDF, $\nu$=1')
#pylab.plot(rice_pdf_corr_bin_center[0:len(rice_pdf_corr_hist)],rice_pdf_corr_hist,'ro',markersize=4.0,label='Rice PDF / median * sqrt(ln(4))')

#pylab.plot(rice_bin_center[0:len(rice_hist)],rice_hist,'ro',markersize=4.0,label='Mean-subtracted')
#pylab.plot(rice_corr_bin_center[0:len(rice_corr_hist)],rice_corr_hist,'m--',linewidth=1.4,label='Nonzero mean, median normalized')
#pylab.plot(norm_bin_center[0:len(norm_hist)],norm_hist,'m--',linewidth=1.4,label='Normal')

#pylab.axis([0.0, 8.0, 0.0, 0.8])

pylab.xlim(0,8)
#pylab.ylim(0,1)

pylab.grid(True)
pylab.xlabel('R',fontsize=12)
pylab.ylabel('PDF(R)',fontsize=12)
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)

#pylab.title(ifo + ' ' + epoch + ' -- Distribution of Detector Noise by Frequency Band',fontsize=12)
#pylab.legend(('R = sqrt[x^2 + y^2]','R/median(R)','R/median(R) * sqrt(log(4))','Rayleigh PDF','Nonzero mean','Nonzero mean, median normalized'),loc=1,prop=fontP)

pylab.legend(loc=1,prop=fontP)

pylab.savefig('rayleigh_median_compare.pdf')
