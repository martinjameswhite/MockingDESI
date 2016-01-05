#!/usr/bin/env python
#
# Read a FITS file containing the 3D positions and velocities of our QSOs,
# and their absolute magnitudes, and convert it to a FileHandler "file"
# that can be processed easily.
#
# Uses Erin Sheldon's FITSIO package:
#	https://github.com/esheldon/fitsio
# to read the FITS file.
#
from __future__ import print_function,division


import numpy  as np
import fitsio as F
import ndfilehandler as FH
import sys


__author__ = "Martin White"
__version__ = "1.0"
__email__  = "mjwhite@lbl.gov"



def convert(finp,fout):
    """    
    Does the work of converting the FITS file to the FileHandler format.
    finp: is a FITS filename, with the data assumed to be in ext 1.
    fout: is a FileHandler output filename.
    """
    # Read the relevant columns from the FITS file.
    fits = F.read(finp,columns=['POS','VEL','MIZ2'],ext=1)
    # Put them in a dictionary and write them to a FileHandler file.
    data = {}
    data['pos']  = fits['POS'][:,:]
    data['vel']  = fits['VEL'][:,:]
    data['flw']  = fits['VEL'][:,:]	# Assume all are centrals.
    data['Miz2'] = fits['MIZ2'][:]
    FH.write_file(fout,data)
    #



if __name__=="__main__":
    if len(sys.argv)!=3:
        raise RuntimeError,\
          "Usage: %s <FITS-file-name> <FH-file-name>"%sys.argv[0]
    else:
      finp = sys.argv[1]
      fout = sys.argv[2]
      convert(finp,fout)
    #
