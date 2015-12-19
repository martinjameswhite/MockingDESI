#!/usr/bin/env python
# 
# A python script to create a Mangle mask for the DESI survey
# outline, consisting of one polygon per pointing.
# To get the DESI survey area, or the DESI footprint or angular mask,
# use by "desi-tiles.fits" on the DESI Wiki in 
#	code/desimodel/trunk/data/footprint
# Speed improvements can be obtained by using the Mangle tools to
# pixelize the mask, but this code does not do that.
#
# Uses Erin Sheldon's FITSIO package:
#	https://github.com/esheldon/fitsio
#
from __future__ import print_function,division


import numpy  as np
import fitsio as F
import sys


__author__ = "Martin White"
__version__ = "1.0"
__email__  = "mjwhite@lbl.gov"


def read_centers(fn):
    """   
    Reads the RA, DEC of the plate centers from a FITS file.
    """
    # Read the tile centers from the supplied file.  For now
    # keep everything.  We could apply cuts later, or store
    # a mask value in each polygon.
    dat = F.read(fn,1,columns=['RA','DEC','IN_DESI'])
    ww  = np.where(dat['IN_DESI']>0)
    ra  = dat['RA'][ww]
    dec = dat['DEC'][ww]
    return( (ra,dec) )
    #



def make_mask(plyfn,ra,dec,diam=3.2):
    """   
    Writes a Mangle .ply file containing polygons, each
    of one cap, centered on the supplied RA, DEC.
    Each polygon has a cap of diameter "diam".
    """
    # Do a sanity check on the input.
    if diam>180.:
        raise RuntimeError,"Weird value of diam."
    # Work out the unit vectors pointing to each RA, DEC
    # and the opening angle for each cap.
    th        = np.pi/180*(90.-dec)
    phi       = np.pi/180*ra
    nhat      = np.zeros( (ra.size,4) )
    nhat[:,0] = np.sin(th)*np.cos(phi)
    nhat[:,1] = np.sin(th)*np.sin(phi)
    nhat[:,2] = np.cos(th)
    nhat[:,3] = 1-np.cos(np.pi/360.*diam)
    area      = 2.*np.pi*nhat[0,3]
    # Now make a Mangle .ply file, with one polygon of one
    # cap for each tile position.
    ff   = open(plyfn,"w")
    ff.write("%d polygons\n"%len(ra))
    for i in range(ra.size):
        ff.write("polygon %10d ( 1 caps, 1.0 weight, %12.8f str)\n"%\
          (i,area))
        ff.write("%23.20f %23.20f %23.20f %23.20f\n"%\
          (nhat[i,0],nhat[i,1],nhat[i,2],nhat[i,3]))
    ff.close()
    #




if __name__=="__main__":
    if len(sys.argv)!=2:
        raise RuntimeError,"Usage: %s <tile.fits>"%sys.argv[0]
    else:
        diam   = 3.21	# DESI pointing diameter, in degrees.
        ra,dec = read_centers(sys.argv[1])
        make_mask("tiles.ply",ra,dec,diam)
    #
