#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *

#########

# matrices are in the following format: 
#
#  POSX->P1   POSY->P1   ANGX->P1   ANGY->P1
#  POSX->P2   POSY->P2   ANGX->P2   ANGY->P2
#  POSX->Y1   POSY->Y1   ANGX->Y1   ANGY->Y1
#  POSX->Y2   POSY->Y2   ANGX->Y2   ANGY->Y2
#
#  ...and include the subtraction of the pre-existing offset
#

#########

"""
# Measurements using 200-ct offsets

# pre-existing offsets
p1 = -0.49
p2 = -0.72
y1 = -0.68
y2 = 0.25

sense1 = array([[-0.21-p1, -3.16-p1, -1.08-p1,  3.13-p1],
                [-0.54-p2, -2.38-p2, -0.42-p2, -3.33-p2],
                [ 3.12-y1, -0.60-y1, -4.21-y1, -0.64-y1],
                [ 0.97-y2,  0.21-y2,  2.11-y2,  0.19-y2]])
"""

#########

# Measurements using 300-ct offsets

# pre-existing offsets
p1 = -0.50
p2 = -0.78
y1 = -0.60
y2 = 0.26

sense2 = array([[-0.03-p1, -4.68-p1, -1.37-p1,  4.89-p1],
                [-0.34-p2, -3.21-p2, -0.05-p2, -4.75-p2],
                [ 5.95-y1, -0.56-y1, -6.19-y1, -0.31-y1],
                [ 1.38-y2,  0.18-y2,  3.08-y2,  0.11-y2]])



#### Construct the pitch and yaw arrays by hand from the 300-ct offset data

"""
# pitch = POSY->P1  ANGY->P1
#         POSY->P2  ANGY->P2
#
# yaw   = POSX->Y1  ANGX->Y1
#         POSX->Y2  ANGX->Y2
#

pitch = array([[-4.68-p1,  4.89-p1],
               [-3.21-p2, -4.75-p2]])

yaw   = array([[ 5.95-y1, -6.19-y1],
               [ 1.38-y2,  3.08-y2]])
"""



#### More complicated version - includes A2A coupling (?)

# pitch = POSY->P1  ANGY->P1
#         POSY->P2  ANGY->P2
#         POSY->Y1  ANGY->Y1
#         POSY->Y2  ANGY->Y2
#
# yaw   = POSX->P1  ANGX->P1
#         POSX->P2  ANGX->P2
#         POSX->Y1  ANGX->Y1
#         POSX->Y2  ANGX->Y2

pitch = array([[-4.68-p1,  4.89-p1],
               [-3.21-p2, -4.75-p2],
               [-0.56-y1, -0.31-y1],
               [ 0.18-y2,  0.11-y2]])

yaw   = array([[-0.03-p1, -1.37-p1], 
               [-0.34-p2, -0.05-p2],
               [ 5.95-y1, -6.19-y1],
               [ 1.38-y2,  3.08-y2]])

#print sense1/200
#print sense2/300
#print (sense1/200) / (sense2/300)


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
