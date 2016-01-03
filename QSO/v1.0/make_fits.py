#!/usr/bin/env python
#
# Takes the results produced by the upstream codes and packages
# them into "final" FITS files which can also be moved to an
# appropriate location.
#
# Uses Erin Sheldon's FITSIO package:
#	https://github.com/esheldon/fitsio
# to read the FITS file.
#
from __future__ import print_function,division


__author__ = "Martin White"
__version__ = "1.0"
__email__  = "mjwhite@lbl.gov"


import numpy  as np
import fitsio as F
import ndfilehandler as FH
import subprocess
import sys



def make_fits(ioct=0):
    """   
    Does the work of making the FITS file by catenating the
    redshift slices for each octant.
    In future this could have more complex functionality.
    """
    zlist = np.arange(0.50,3.76,0.25)
    dat   = FH.read_file("v1.0_qso_z%4.2f_oct%d"%(zlist[0],ioct))
    for zcen in zlist[1:]:
        dat1 = FH.read_file("v1.0_qso_z%4.2f_oct%d"%(zcen,ioct))
        for f in dat.keys():
            dat[f] = np.append(dat[f],dat1[f])
    #
    fout = "v1.0_qso_lc_oct%d.fits"%ioct
    subprocess.check_call(["/bin/rm","-f",fout])
    header = {'Version':1.0,'Class':'QSO','Model':'CW13'}
    fits   = F.FITS(fout,'rw')
    fits.write(dat,header=header)
    fits.close()
    #



if __name__=="__main__":
    for ioct in range(1):
        make_fits(ioct)
    #
