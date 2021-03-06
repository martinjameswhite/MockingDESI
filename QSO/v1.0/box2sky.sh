#!/bin/bash -l
#SBATCH -J Box2Sky
#SBATCH -t 0:10:00
#SBATCH -n 4
#SBATCH -o Box2Sky.out
#SBATCH -e Box2Sky.err
#SBATCH -p shared
#SBATCH -A desi
#
export OMP_NUM_THREADS=4
#
# The DESI footprint
mask=../../data/tiles.pix.ply
# For the OuterRim simulation
OmM=0.265
hub=0.710
Box=3000.	# Mpc/h
# A fiducial cosmology
OmM=0.30
hub=0.68
Box=2873.	# Mpc/h converted to have same Mpc value as above
# The file to process.
zcen=2.000
zmin=1.875
zmax=2.125
fb=qso_z2.00
# Cori
srun -n 1 -c 4 ../../generic/box2sky/box2sky \
  ${OmM} ${hub} ${Box} ${zcen} ${zmin} ${zmax} ${mask} ${fb}
#
