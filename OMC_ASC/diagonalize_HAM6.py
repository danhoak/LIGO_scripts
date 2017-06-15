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

#            OM1      OM2       OM3
# AS_A     -63.8     -52.1      --
# AS_B   -1*-52.3   -1*-75.3    --
# OMCA   -1*-57.6   -1*-54.2  -1*-62.5
# OMCB   -1*-62.5   -1*-52.8  -1*-58.5

# yaw

#           OM1      OM2       OM3
# AS_A     -66.3   -1*-54.7    --
# AS_B   -1*-52.3   -71.8      --
# OMCA   -1*-58.5   -55.5    -1*-62.4
# OMCB     -63.8   -1*-54.8   -59.2

pitch_sense_dB = array([[-63.8, -52.1, 1.0],[-52.3, -75.3, 1.0],[-57.6, -54.2, -62.5],[-62.5, -52.8, -58.5]])
pitch_mask = array([[1., 1., 0.0], [-1., -1., 0.0],[-1.0, -1.0, -1.0],[-1.0, -1.0, -1.0]])

yaw_sense_dB = array([[-66.3, -54.7, 1.0],[-52.3, -71.8, 1.0],[-58.5, -55.5, -62.4],[-63.8, -54.8, -59.2]])
yaw_mask = array([[1., -1., 0.0], [-1., 1., 0.0],[-1., 1., -1.],[1., -1., 1.]])

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


yaw_sense_DOF = yaw_sense.copy()
pitch_sense_DOF = pitch_sense.copy()

# POS X
yaw_sense_DOF[2,0] = 381*yaw_sense[2,0] + 46.3*yaw_sense[3,0]
yaw_sense_DOF[2,1] = 381*yaw_sense[2,1] + 46.3*yaw_sense[3,1]
yaw_sense_DOF[2,2] = 381*yaw_sense[2,2] + 46.3*yaw_sense[3,2]

# ANG X
yaw_sense_DOF[3,0] = -787.4*yaw_sense[2,0] + -1067.5*yaw_sense[3,0]
yaw_sense_DOF[3,1] = -787.4*yaw_sense[2,1] + -1067.5*yaw_sense[3,1]
yaw_sense_DOF[3,2] = -787.4*yaw_sense[2,2] + -1067.5*yaw_sense[3,2]

# POS Y
pitch_sense_DOF[2,0] = 518*pitch_sense[2,0] + -56.2*pitch_sense[3,0]
pitch_sense_DOF[2,1] = 518*pitch_sense[2,1] + -56.2*pitch_sense[3,1]
pitch_sense_DOF[2,2] = 518*pitch_sense[2,2] + -56.2*pitch_sense[3,2]

# ANG Y
pitch_sense_DOF[3,0] = -1070*pitch_sense[2,0] + 1294*pitch_sense[3,0]
pitch_sense_DOF[3,1] = -1070*pitch_sense[2,1] + 1294*pitch_sense[3,1]
pitch_sense_DOF[3,2] = -1070*pitch_sense[2,2] + 1294*pitch_sense[3,2]



# We only want AS_B, POS, ANG
# Rows 1,2,3

#yaw_sense_3by3 = yaw_sense_DOF[1:,:]
#pitch_sense_3by3 = pitch_sense_DOF[1:,:]

yaw_sense_3by3 = yaw_sense_DOF[[0,2,3],:]
pitch_sense_3by3 = pitch_sense_DOF[[0,2,3],:]

print
print 'Yaw 3-by-3:'
print yaw_sense_3by3
print
print 'Pitch 3-by-3:'
print pitch_sense_3by3

print
print 'YAW control:'
print linalg.inv(yaw_sense_3by3)

print
print 'PITCH control:'
print linalg.inv(pitch_sense_3by3)



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
