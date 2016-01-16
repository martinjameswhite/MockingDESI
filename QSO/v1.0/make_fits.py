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
    redshift slices for each octant -- this needs to be coordinated
    with "driver.py", which creates those slices.

    Since the QLF is changing fairly rapidly with z and we have
    relatively coarse binning in z at present, the catenated file
    has jumps across bin boundaries that are not ideal.  To deal
    with this we make the shells wider than they need to be and
    then randomly select objects from the files based on their distance
    (in redshift) from the shell center.  This is simple, but not optimal.

    In future this code could have more complex functionality which
    helped resolve this issue.
    """
    dz    = 0.25	# The redshift spacing of the slices.
    zlist = np.arange(0.50,3.76,dz)
    # Use a "triangle" function, peaked at the central redshift of the
    # slice, to select the QSOs.  For constantly spaced slices of equal
    # width the total selection probability is 1, for more complex cases
    # better logic needs to be used.
    dat = FH.read_file("v1.0_qso_z%4.2f_oct%d"%(zlist[0],ioct))
    wt  = 1 - np.abs(dat['Z']-zlist[0])/dz
    ww  = np.nonzero( np.random.uniform(size=dat['Z'].size)<wt )[0]
    print("Keeping %8d of %10d for zcen=%.2f"%(len(ww),dat['Z'].size,zlist[0]))
    for f in dat.keys():
        dat[f] = dat[f][ww]
    for zcen in zlist[1:]:
        dat1 = FH.read_file("v1.0_qso_z%4.2f_oct%d"%(zcen,ioct))
        wt   = 1 - np.abs(dat1['Z']-zcen)/dz
        ww   = np.nonzero( np.random.uniform(size=dat1['Z'].size)<wt )[0]
        print("Keeping %8d of %10d for zcen=%.2f"%(len(ww),dat1['Z'].size,zcen))
        for f in dat.keys():
            dat[f] = np.append(dat[f],dat1[f][ww])
    #
    fout = "v1.0_qso_lc_oct%d.fits"%ioct
    subprocess.check_call(["/bin/rm","-f",fout])
    header = {'Version':1.0,'Class':'QSO','Model':'CW13'}
    fits   = F.FITS(fout,'rw')
    fits.write(dat,header=header)
    fits.close()
    #



if __name__=="__main__":
    for ioct in range(8):
        make_fits(ioct)
    #
