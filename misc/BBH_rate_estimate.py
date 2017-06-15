#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version of Steve Fairhurst's rate estimates for GW150914
#
# https://geo2.arcca.cf.ac.uk/~spxsf2/LVC/notebooks/Rates.html
#
# Takeaway is this: ER8 is equivalent to S5/S6
# We didn't observe anything in S5/S6, thus an observation in ER8 could not have been probable 
# (i.e., we got lucky, maybe only a little lucky, but still)
# This reduces the expectation for the rest of O1.
#
#
# Furthermore, not seeing an event in O1 never reduces the probability of seeing an event in ER8
# to below 2-sigma
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
from scipy.stats import poisson

# S5, S6 from IMR papers

time = {}
time["s5"] = 0.7  # after cat3 vetoes - seems pretty low...
time["s6"] = 0.35
time["er8"] = 0.044  # 16 days
time["o1"] = 0.123  # 3 months at 50% duty cycle

bns_range = {}
bns_range["s5"] = 13.3
bns_range["s6"] = 17.7
bns_range["er8"] = 70.
bns_range["o1"] = 70.

# Correct for ratio of 30-30 BBH to BNS range, divide by 10x to keep things small
# S5/S6 comes from IMR paper
#
# S5: BNS is 13.3, BBH with Mchirp = 28 is 148
# S6: BNS is 17.6, BBH with Mchirp = 28 is 194
#
# ER8 sensitivity comes from pyCBC software injections, 700Mpc for 30-30 BBH
# Surprisingly difficult to get a straight answer abour range / horizon distance
# O1 assumed the same as ER8
range_fac = {}
range_fac["s5"] = 1.11  # did S6 improve sensitivity to BNS but sensitivity to BBH did not change as much?
range_fac["s6"] = 1.1
range_fac["er8"] = 1.
range_fac["o1"] = 1.


vt = {}
for run in time.keys():
    vt[run] = time[run] * bns_range[run]**3 * range_fac[run]**3
print vt

print time.keys()

fignum=0

matplotlib.rcParams.update({'savefig.dpi':250})

fignum=fignum+1
pylab.figure(fignum)

i=1


colors = {'er8': 'k', 'o1': 'r', 's5': 'g', 's6': 'b'}

# Log prior on rate
def prob_log(N, vt, vt0):
    return vt**N * vt0 / (vt + vt0)**(N+1)

fudge_factor = 0.5
print 'Probability to observe exactly one event in ER8 is: %.2f' % poisson.pmf(1, fudge_factor)

for run in vt.keys():

    pylab.subplot(2,2,i)

    # set the probability to observe 1 event in ER8
    # for poisson distribution mean of 0.7 returns a 50% chance of >0 events
    """
    mu = 0.7*vt[run]/vt["er8"]
    print mu
    x = range(13)
    pylab.bar(x, poisson.pmf(x, mu), color= colors[run], alpha = 0.7, label=run)
    """

    x = arange(50)
    p = prob_log(x, vt[run], vt["er8"]/fudge_factor)
    pylab.bar(x, p, color= colors[run], alpha = 0.7, label=run)
    print("In run %s, expected number of events %.2f" % (run, sum(x*p)))

    pylab.legend()
    pylab.grid()
    pylab.xlim([0,20])

    i+=1
#pylab.xlabel("Number of events")
#pylab.ylabel("Probability")

#pylab.savefig('rates_no_uncertainty.png',bbox_inches='tight')
#pylab.savefig('rates_no_uncertainty.png')
pylab.savefig('rates_log_prior.png')
