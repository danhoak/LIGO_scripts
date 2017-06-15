#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script to generate movie of OMC mode scans
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys, os, commands
from numpy import *
#from scipy.interpolate import interp2d
from PIL import Image
#import pylab


# code to number the image files sensibly

files = os.listdir('sweep_images/')

cmd = 'rm sw_images/*.gif'
commands.getoutput(cmd)
print len(files)
for i in range(len(files)):

    file_in = files[i]
    print file_in
    if i < 10:
        file_out = 'sw_images/OMC_0' + str(i) + '0.gif'
    else:
        file_out = 'sw_images/OMC_' + str(i) + '0.gif'

    im = Image.open('sweep_images/' + file_in)
    im.save(file_out,"GIF")

cmd = 'convert -delay 100 -loop 0 sw_images/OMC*.gif sweep_movie.gif'
stin, out, err = os.popen3(cmd)
pid, status = os.wait()

