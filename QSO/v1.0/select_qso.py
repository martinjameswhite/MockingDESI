#!/usr/bin/env python
#
# Given an octant of QSOs apply an approximation to the
# DESI selection.
#
#
from __future__ import print_function,division


import numpy  as np
import ndfilehandler as FH
import sys
import os


__author__ = "Martin White"
__version__ = "1.0"
__email__  = "mjwhite@lbl.gov"



def calc_selection(octant):
    """    
    Decides whether a given QSO is "in" the DESI survey based on an
    approximation to the selection function.
    Writes a field, "in_survey" containing 1 or 0 (currently an integer).
    """
    # Read the QSO data.
    data = FH.read_file(octant)
    magg = data['GMAG']
    insurv=np.zeros(magg.size,dtype='i4')
    # An approximation to the survey selection function.
    rr     = np.random.uniform(size=magg.size)
    bright = np.nonzero( (magg< 22.0)             )[0]
    medium = np.nonzero( (magg>=22.0)&(magg<22.5) )[0]
    faint  = np.nonzero( (magg>=22.5)&(magg<23.0) )[0]
    insurv[bright[np.nonzero( rr[bright]<0.83 )[0]]]=1
    insurv[medium[np.nonzero( rr[medium]<0.72 )[0]]]=1
    insurv[faint [np.nonzero( rr[faint ]<0.37 )[0]]]=1
    #
    data = {}
    data['IN_SURVEY'] = insurv
    FH.write_file(octant,data)
    #



if __name__=="__main__":
    if len(sys.argv)!=2:
        raise RuntimeError,\
          "Usage: %s <oct-file-name>"%sys.argv[0]
    else:
      foct = sys.argv[1]
      calc_selection(foct)
    #
