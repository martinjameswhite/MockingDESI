#!/usr/bin/env python
#
# Given an octant of QSOs and dNdzdg 
# Select QSO to mathc the dNdgdz
#
#
from __future__ import print_function,division


import numpy  as np
import ndfilehandler as FH
import sys
import os


__author__ = "Shadab Alam"
__version__ = "1.0"
__email__  = "shadaba@andrew.cmu.edu"




def Apply_dNdgdz(dNdzdgfile,magg,redshift,Obs_Area):
   #Obsd_Area should be in deg^-2 after mask

   data=np.loadtxt(dNdzdgfile)
   nmag=np.sum(data[:,0]==data[0,0],dtype='int')
   nz=np.int(data.shape[0]/nmag)
   zbins=data[::nmag,0]
   magbins=data[:nmag,1]
   dNdzdg=data[:,2].reshape(nz,nmag)

   dmag=0.5*np.mean(magbins[1:]-magbins[:-1])
   dz=0.5*np.mean(zbins[1:]-zbins[:-1])

   zmin=np.min(redshift)
   zmax=np.max(redshift)

   magg_min=np.min(magg)
   magg_max=np.max(magg)

   print("nz ,nmag: %i %i %f %f" %(redshift.size,nmag,magg_min,magg_max))
   insurv=np.zeros(magg.size,dtype='i4')
   for zz in range(0,nz):
      zbin_min=zbins[zz]-dz
      zbin_max=zbins[zz]+dz
      if(zbin_max<zmin or zbin_min>zmax):
	 continue
      for gg in range(0,nmag):
         subsamp=np.nonzero((magg>=magbins[gg]-dmag)&(magg<magbins[gg]+dmag)
	             &(redshift>=zbin_min)&(redshift<zbin_max))[0]

	 dz_bin=min([dz,np.abs(zbins[zz]-zmin),np.abs(zbins[zz]-zmax)])
         Nsel_qso=np.int(dNdzdg[zz][gg]*Obs_Area*(dz+dz_bin)*2*dmag)
	 if(subsamp.size<Nsel_qso):
	    print("Not Enough QSOs in simulation (z=%f,magg=%f, Nsim=%i , Nsel=%i ):"%(
	                                zbins[zz] ,magbins[gg],subsamp.size,Nsel_qso))
	    continue
            raise RuntimeError,\
                   "Not Enough QSOs in simulation (z=%f,magg=%f, Nsim=%i , Nsel=%i ):"%(
			 zbins[zz] ,magbins[gg],subsamp.size,Nsel_qso)

         sel=np.random.choice(subsamp, size=Nsel_qso, replace=False)
	 print('%f %f Nsim , Nsel Ngot: %i %i %i' %( zbins[zz] ,magbins[gg],subsamp.size,Nsel_qso,sel.size))
         insurv[sel]=1

   
   print('NQSO in sim %i , NQSO selceted %i' %(redshift.size,np.sum(insurv)))
   return insurv

def calc_selection(octant,fselect,Obs_Area):
    """    
    Decides whether a given QSO is "in" the DESI survey based on
    the selection function.
    Writes a field, "in_surveyD" containing 1 or 0 (currently an integer).
    """
    # Read the QSO data.
    data = FH.read_file(octant)
    magg = data['GMAG']
    redshift=data['Z']

    # Survey selection function.
    insurv=Apply_dNdgdz(fselect,magg,redshift,Obs_Area)
    
    #
    data = {}
    data['IN_SURVEYD'] = insurv
    FH.write_file(octant,data)
    #



if __name__=="__main__":
    if len(sys.argv)!=4:
        raise RuntimeError,\
          "Usage: %s <oct-file-name, dNdgdZfile, Obs_Area>"%sys.argv[0]
    else:
      foct = sys.argv[1]
      fselect=sys.argv[2]
      Obs_Area=np.float(sys.argv[3])
      calc_selection(foct,fselect,Obs_Area)

    #
