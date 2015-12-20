#!/bin/bash -l
#SBATCH -J Box2Sky
#SBATCH -t 0:20:00
#SBATCH -n 4
#SBATCH -o Box2Sky.out
#SBATCH -e Box2Sky.err
#SBATCH -p shared
#SBATCH -A desi
#
export OMP_NUM_THREADS=4
#
# The DESI footprint
mask=../data/tiles.pix.ply
# For the OuterRim simulation
OmM=0.265
hub=0.710
Box=3000.	# Mpc/h
# The file to process.
zmin=1.5
zmax=2.5
fb=qso_z2.00
# Cori
srun -n 1 -c 4 ../generic/box2sky/box2sky \
  ${OmM} ${hub} ${Box} ${zmin} ${zmax} ${mask} ${fb}
#
