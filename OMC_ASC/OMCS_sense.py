#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *

#########

# data from Jan 31 2015
# 3.9Hz drive to OM3 and OMCS
# choose +/- 180 deg to be negative

# pitch (in dB):

#           OM3      OMCS
# AS_A   -1*-7.9    -48.5
# AS_B   -1*-1.16   -31.2


# yaw

#           OM3      OMCS
# AS_A   -1*-11.3   -1*-49.7
# AS_B   -1*-2.78   -1*-30.6

pitch_sense_dB = array([[-7.9, -48.5],[-1.16, -31.2]])
pitch_mask = array([[-1., 1.], [-1., 1.]])

 yaw_sense_dB = array([[-11.3, -49.7],[-2.78, -30.6]])
yaw_mask = array([[-1., -1.], [-1., -1.]])

pitch_sense = 10**(pitch_sense_dB/20) * pitch_mask

yaw_sense = 10**(yaw_sense_dB/20) * yaw_mask

print pitch_sense_dB
print yaw_sense_dB

print
print 'Pitch Sensing Matrix:'
print pitch_sense

print
print 'Yaw Sensing Matrix:'
print yaw_sense

print
print 'Pitch Control Matrix'
print linalg.inv(pitch_sense)

print
print 'Yaw Control Matrix'
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
