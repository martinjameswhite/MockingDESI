#!/usr/bin/env python
#
# Wrapper code to generate octants of QSOs for each redshift
# slice given precomputed "boxes of QSOs".
#
# Writes and submit scripts to run "box to sky", then generates
# observed magnitudes and applied a mock DESI selection to the
# resulting objects.
#
# This code should be run interactively from the mock directory
# as "./driver.py".
#
#
from __future__ import print_function,division


import numpy   as np
import fits2fh as F2FH
import subprocess


__author__ = "Martin White"
__version__ = "1.0"
__email__  = "mjwhite@lbl.gov"




# Some paths at NERSC for where the mocks live.
mockbase = "/project/projectdirs/desi/mocks/qso_v1.0/"
# The SLURM batch submission command on Cori
sbatch   = "/opt/slurm/default/bin/sbatch"


# A fiducial cosmology.
OmM = 0.30
hub = 0.68
Lbox= 2873.	# Converted to Mpc/h holding Mpc fixed.



def submit_jobs():
    """    
    Does the work of writing the scripts and submitting the jobs.
    """
    # Submit a different job for each redshift slice.
    for zcen in np.arange(0.50,3.76,0.25):
        fb = "v1.0_qso_z%4.2f"%zcen
        F2FH.convert(mockbase+fb+".fits",fb)
        #
        ff = open("driver_script.sh","w")
        ff.write("""#!/bin/bash -l
#SBATCH -J Box2Sky
#SBATCH -t 0:10:00
#SBATCH -n 4
#SBATCH -o Box2Sky.out
#SBATCH -e Box2Sky.err
#SBATCH -p shared
#SBATCH -A desi
#
export OMP_NUM_THREADS=4
#\n""")
        ff.write("# The DESI footprint\n")
        ff.write("mask=../../data/tiles.pix.ply\n")
        ff.write("# The file to process.\n")
        ff.write("fb=%s\n"%fb)
        ff.write("# Cori\n")
        ff.write("srun -n 1 -c 4 ../../generic/box2sky/box2sky \\\n")
        ff.write("  %f %f %f "%(OmM,hub,Lbox))
        ff.write(" %f %f %f ${mask} ${fb}\n"%(zcen,zcen-0.25/2,zcen+0.25/2))
        ff.write("#\n")
        for ioct in range(1):	# For each octant, do ...
            ff.write("./add_magnitude.py %s %s_oct%d\n#\n"%(fb,fb,ioct))
            ff.write("./select_qso.py %s_oct%d\n#\n"%(fb,ioct))
        ff.close()
        subprocess.check_call([sbatch,"driver_script.sh"])
    subprocess.check_call(["/bin/rm","-f","driver_script.sh"])
    #




if __name__=="__main__":
    submit_jobs()
    #
