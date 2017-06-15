#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *

#########

# matrices are in the following format: 
#
#  POSY->P1   ANGY->P1
#  POSY->P2   ANGY->P2


#  POSX->Y1   ANGX->Y1
#  POSX->Y2   ANGX->Y2

# ANGX EXC:
# Y1 = 0.000152 -123deg
# Y2 = 0.000335 40deg

# POSX EXC:
# Y1 = 0.000357 -9deg
# Y2 = 0.001091 134deg



# ANGY EXC:
# P1 = 0.000411 -22deg
# P2 = 0.000650 160deg

# POSY EXC:
# P1 = 0.000839 170deg
# P2 = 0.000776 -52deg

# phase convention: pitch, -40deg is positive
#                   yaw, 0deg is positive

"""
pit_sense = array([[-0.000839,0.000411],
                   [0.000776, -0.000650]])


yaw_sense = array([[0.000357, -0.000152],
                   [-0.001091, 0.000335]])


print linalg.inv(pit_sense)
print
print linalg.inv(yaw_sense)
"""


pit_sense = array([[0.045,-0.03],
                   [0.041-0.013, 0.05-0.013]])/30


yaw_sense = array([[0.086, -0.006],
                   [-0.03, 0.055]])/30


print linalg.inv(pit_sense)
print
print linalg.inv(yaw_sense)
