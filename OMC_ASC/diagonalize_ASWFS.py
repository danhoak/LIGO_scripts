#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *

#########

# data from Sep 28 2014
# 3.9Hz drive to tip tilts
# phase of 25deg taken as positive, -155deg is negative

# pitch (in dB):

#            OM1      OM2
# AS_A     -63.8     -52.1
# AS_B   -1*-52.3   -1*-75.3


# yaw

#           OM1      OM2
# AS_A     -66.3   -1*-54.7
# AS_B   -1*-52.3   -71.8

pitch_sense_dB = array([[-63.8, -52.1],[-52.3, -75.3]])
pitch_mask = array([[1., 1.], [-1., -1.]])

yaw_sense_dB = array([[-66.3, -54.7],[-52.3, -71.8]])
yaw_mask = array([[1., -1.], [-1., 1.]])

pitch_sense = 10**(pitch_sense_dB/20) * pitch_mask

yaw_sense = 10**(yaw_sense_dB/20) * yaw_mask

print
print 'Pitch Sensing Matrix (dB):'
print pitch_sense_dB

print
print 'Yaw Sensing Matrix (dB):'
print yaw_sense_dB


print
print 'Pitch Sensing Matrix:'
print pitch_sense

print
print 'Yaw Sensing Matrix:'
print yaw_sense

print
print linalg.inv(pitch_sense)

print
print linalg.inv(yaw_sense)


"""
# Use the simple one...
dither2DOF = linalg.inv(sense2/300)

print dither2DOF

posx_inputs = array([4.89-p1, -4.75-p2, -0.31-y1, 0.11-y2])

print
print posx_inputs
print dot(posx_inputs,dither2DOF.T)

test = array([1, 0, 0, 0])
print
print dot(test, dither2DOF)
"""
