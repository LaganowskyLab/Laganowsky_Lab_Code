

fetch 1RC2
bg_color white
as cartoon

python

from pymol import cmd
import numpy
import sys
sys.path.insert(0, "/Users/alaganowsky/opt/anaconda3/lib/python3.7/site-packages")

import matplotlib
from matplotlib import cm

# get data
d = numpy.loadtxt("./AqpZUVPD6mjTMCDL.txt")

# normalize
raw = d[::,1]
if numpy.min(raw) == 0:
    m = numpy.max(d[::,1])
else:
    print("min not equal to zero")
    sys.exit()

# get colormap
cm = cm.get_cmap("gnuplot")

for resi, intensity in d:

    # normalize and get cmap rgba value
    color = cm(intensity/m)

    # set color
    cmd.set_color( "color_%i" % resi, color[:3])   
    cmd.color( "color_%i" % resi, "polymer and resi %i" % resi)


python end
