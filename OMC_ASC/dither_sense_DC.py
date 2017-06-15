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

pit_sense = array([[0.045,-0.03],
                   [0.041-0.013, 0.05-0.013]])/30


yaw_sense = array([[0.086, -0.006],
                   [-0.03, 0.055]])/30


print linalg.inv(pit_sense)
print
print linalg.inv(yaw_sense)
