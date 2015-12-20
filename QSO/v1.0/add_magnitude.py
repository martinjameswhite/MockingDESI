#!/usr/bin/env python
#
# Given an octant of QSOs, generated by the box2sky code,
# add a field giving their g-band apparent magnitude.
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



def calc_magnitude(box,octant):
    """    
    Generates the apparent magnitude given the absolute magnitude,
    obtained from the "box", and the distance modulus and K-correction.
    """
    # Read the Mi(z=2) magnitudes for the box.
    miz2 = FH.read_file(box)['Miz2'][:]
    # Read the index for each QSO in the octant, and get the Mi(z=2).
    data = FH.read_file(octant)
    zz   = data['Z']
    dmod = data['DMOD']
    miz2 = miz2[data['INDX']]
    # Now convert to apparent i-band magnitude using the k-correction.
    # If a tabulated k-correction is available, use that, otherwise
    # default to a power-law continuum approximation.
    # See discussion in Ross++13, Appendix B and Section 4.
    kfile=os.getenv('MOCKINGDESI_BASE')+"/data/qso-iband-k-correction.txt"
    if os.path.exists(kfile):
        print("Using K-correction from "+kfile)
        kcorr = np.loadtxt(kfile)
        kcorr = np.interp(zz,kcorr[:,1],kcorr[:,2])
    else:
        print("Using power-law K-correction")
        alpha = -0.5
        kcorr = -2.5*(1+alpha)*np.log10( (1+zz)/(1+2.0) )
    magi = miz2 + dmod + kcorr	# e.g. Ross++13, Eq. 5
    magg = magi + 0.255		# e.g. Ross++13, Eq. B7.
    # and write the results
    data = {}
    data['GMAG'] = magg.astype('f4')
    FH.write_file(octant,data)
    #



if __name__=="__main__":
    if len(sys.argv)!=3:
        raise RuntimeError,\
          "Usage: %s <box-file-name> <oct-file-name>"%sys.argv[0]
    else:
      fbox = sys.argv[1]
      foct = sys.argv[2]
      calc_magnitude(fbox,foct)
    #