#!/usr/bin/env python
#
# Generates a number of different projections of periodic boxes onto the
# sky, each covering a different redshift, using box2sky.
# Stitches these together to get a single "light cone" output.
# Uses interpolation of the positions and velocities to the
# appropriate redshift.
#
from __future__ import print_function,division


import numpy  as np
import ndfilehandler as FH
import subprocess
import sys


__author__ = "Martin White"
__version__ = "1.0"
__email__  = "mjwhite@lbl.gov"



def stitch():
    """    
    Does the work.
    """
    print("This is currently a stub.")
    print("-- this functionality is being moved to the")
    print("-- individual mock sections.")
    #



if __name__=="__main__":
    if len(sys.argv)!=1:
        raise RuntimeError,\
          "Usage: %s"%sys.argv[0]
    else:
        stitch()
    #
