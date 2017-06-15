#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# This script replicates the calculations to find the correct OMC QDP -> POS,ANG basis,
# and the output matrix, for the OMC centering loops.
#
# The calculations are described by Koji in T1400585
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *

# First we need to use the pitch & yaw transfer functions of the OMs to
# get the calibration of the drives into radians

# DC response:

R_OM1_DC_P = 0.890
R_OM2_DC_P = 1.137
R_OM3_DC_P = 1.087

R_OM1_DC_Y = 0.936
R_OM2_DC_Y = 1.331
R_OM3_DC_Y = 0.962

# AC response (pit at 3.9Hz, yaw at 2.9Hz)

R_OM1_AC_P = -0.249
R_OM2_AC_P = -0.261
R_OM3_AC_P = -0.283

R_OM1_AC_Y = -0.446
R_OM2_AC_Y = -0.481
R_OM3_AC_Y = -0.450

"""
# some guesses for yaw at 3.9Hz
R_OM1_AC_Y = -0.225
R_OM2_AC_Y = -0.225
R_OM3_AC_Y = -0.225
"""

# Calibration matrix:

C = array([[R_OM1_DC_P/R_OM1_AC_P, 0, 0, 0, 0, 0],
           [0, R_OM2_DC_P/R_OM2_AC_P, 0, 0, 0, 0],
           [0, 0, R_OM3_DC_P/R_OM3_AC_P, 0, 0, 0],
           [0, 0, 0, R_OM1_DC_Y/R_OM1_AC_Y, 0, 0],
           [0, 0, 0, 0, R_OM2_DC_Y/R_OM2_AC_Y, 0],
           [0, 0, 0, 0, 0, R_OM3_DC_Y/R_OM3_AC_Y]])

print
print 'Calibration matrix, AC drives to DC response.  From SUS transfer functions.'
print C


# Now we input the transfer function of the OMs to the QPDs

# First in dB units, and using -150deg as negative phase:

T_dB = array([[-59.6, -55.7, -64.1, 0, 0, 0],
              [-63.9, -54.6, -60.4, 0, 0, 0],
              [0, 0, 0, -53.1, -50.3, -57.4],
              [0, 0, 0, -59.1, -50.1, -54.7]])

T_sign = array([[-1, -1, -1, 0, 0, 0],
                [-1, -1, -1, 0, 0, 0],
                [0, 0, 0, -1, 1, -1],
                [0, 0, 0, 1, -1, 1]])

T = 10**(T_dB/20.0) * T_sign

print
print 'Transfer coefficients, OM --> QPD'
print T


# Define an actuation matrix; this picks out which OMs we use for the actuation
# Our choice of OM1 and OM3 was driven by the nominal Guoy phase between the mirrors (70deg)

"""
# For OM1,3
A = array([[1, 0, 0, 0],
           [0, 0, 0, 0],
           [0, 1, 0, 0],
           [0, 0, 1, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 1]])
"""

# For OM1,2
A = array([[1, 0, 0, 0],
           [0, 1, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 1, 0],
           [0, 0, 0, 1],
           [0, 0, 0, 0]])


"""
# For OM2,3
A = array([[0, 0, 0, 0],
           [1, 0, 0, 0],
           [0, 1, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 1, 0],
           [0, 0, 0, 1]])
"""

# Calculate the conversion from the QPDA/B basis to the cavity waist POS, ANG basis
# This is done using the geometry of the OMC breadboard and our calculation of the waist size on the QPDs,
# and the location of OM3 relative to the QPDs.  We use the OM3 -> QPD transfer function numbers to calculate the waist size

# Note that the cavity waist is fixed, but the waist of the input beam to the OMC may vary based on the PRC/SRC, etc.
# But, we only care about the cavity waist; we can't change the input beam waist size and location with the OMs

# Distance from OM3 to QPDA, B

D_OM3_QPDA = 0.520
D_OM3_QPDB = 0.962

T_OM3_QPDA_P = T[0,2]
T_OM3_QPDA_Y = T[2,5]

T_OM3_QPDB_P = T[1,2]
T_OM3_QPDB_Y = T[3,5]

#print T_OM3_QPDA_P
#print T_OM3_QPDA_Y

#print T_OM3_QPDB_P
#print T_OM3_QPDB_Y

#print R_OM3_AC_P
#print R_OM3_AC_Y

w_AP =  2*sqrt(8/pi) * D_OM3_QPDA * R_OM3_AC_P*1e-6 / T_OM3_QPDA_P
w_AY =  2*sqrt(8/pi) * D_OM3_QPDA * R_OM3_AC_Y*1e-6 / T_OM3_QPDA_Y

w_BP =  2*sqrt(8/pi) * D_OM3_QPDB * R_OM3_AC_P*1e-6 / T_OM3_QPDB_P
w_BY = -2*sqrt(8/pi) * D_OM3_QPDB * R_OM3_AC_Y*1e-6 / T_OM3_QPDB_Y

# The minus sign in the last term comes from the odd number of mirrors between OM3 and QPDB

#print
#print w_AP
#print w_AY
#print w_BP
#print w_BY

#print w_BP/w_AP
#print w_BY/w_AY

pit_waists = array([[w_AP*sqrt(pi/8)],[w_BP*sqrt(pi/8)]])
yaw_waists = array([[w_AY*sqrt(pi/8)],[-1*w_BY*sqrt(pi/8)]])

#pit_waists = array([[0.66e-3*sqrt(pi/8)],[0.80e-3*sqrt(pi/8)]])
#yaw_waists = array([[0.55e-3*sqrt(pi/8)],[-1*0.75e-3*sqrt(pi/8)]])

# Use the spot sizes to calculate the matrix that relates QPDA/B signals to the POS, ANG basis

# Distances from QPDA/B to OMC waist position

L_QPDA_waist = 0.0434
L_QPDB_waist = 0.4840

waist_prop = array([[1.0, L_QPDA_waist],
                    [1.0, L_QPDB_waist]])

pitch_input = pit_waists.T*linalg.inv(waist_prop)*1e6

yaw_input = yaw_waists.T*linalg.inv(waist_prop)*1e6

#print pitch_input
#print yaw_input

input_matrix = array([[0, 0, yaw_input[0,0], yaw_input[0,1]],
                      [pitch_input[0,0], pitch_input[0,1], 0, 0],
                      [0, 0, yaw_input[1,0], yaw_input[1,1]],
                      [pitch_input[1,0], pitch_input[1,1], 0, 0]])

print
print 'Input matrix, QDP --> DOF'
print input_matrix

print
print 'Output matrix, DOF --> OMs'
print -1*dot(A,linalg.inv(dot(input_matrix,dot(T,dot(C,A)))))*1e3
