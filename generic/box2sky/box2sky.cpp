#include	<cstdlib>
#include	<cmath>
#include	<iostream>
#include	<iomanip>
#include	<fstream>
#include	<sstream>
#include	<string>
#include	<list>
#include	<vector>
#include	<exception>
#include	"omp.h"

#include	"spline.h"
#include	"mangle.h"
#include	"filehandler.h"



/*

Reads a series of fields which comprise a box of mock "objects".
Replicates this box as necessary.
Applies a Mangle mask to the output and produces a file containing
ra, dec and z of the mock objects, plus the distance modulus and
the 'true' z (i.e. without peculiar velocities).

Currently a serial code, but it could be easily made OpenMP parallel
and with a little work an MPI code.

Author:		Martin White	(UCB,LBNL)
Written:	20-Dec-2015
Modified:

*/


#define		Nres	1000000	// Reserve space for this many objects.
#define		Nrepl	3	// Number of periodic replications.





struct	ObsObj {
  float	ra,dec,zz,zr,dmod,wt;
  int	indx;	// Object index in the original array.
};






void	myexit(const int flag) {
  std::cerr<<"myexit called with flag "<<flag<<std::endl;
  std::cout.flush();
  std::cerr.flush();
  exit(1);
}


void	myexception(std::exception& e) {
  std::cout<<"myexception called with "<<e.what()<<std::endl;
  std::cout.flush();
  std::cerr.flush();
  exit(1);
}


double	periodic(const double x) {
// Wraps x periodically in the range [0,1).
  double tmp = x;
  if (tmp >= 1.0) tmp = x-floor(x);
  if (tmp <  0.0) tmp = 1.0-(-x-floor(-x));
  return(tmp);
}







class Cosmology {
private:
  double Ez(double zz) const {
    // The dimensionless Hubble parameter at redshift zz.
    // This version ignores radiation, so shouldn't be used above z~100.
    double zp1  = 1.0 + zz;
    double rhom = OmM*zp1*zp1*zp1;
    double rhoK = OmK*zp1*zp1;
    double rhox = OmX*pow(zp1,3.0*(1+w0+wa))*exp(-3*wa*zz/zp1);
    double Ez2  = rhom+rhoK+rhox;
    return(sqrt(Ez2));
  }
public:
double	OmM,OmX,OmK,hh,w0,wa;
  double chi(const double zz) const {
    // The comoving distance to redshift zz, in Mpc/h.
    // Specialized for lower redshift.
    const double Lhub=2997.925;	// In Mpc/h.
    const int NN=2500;		// Must be even.
    const double dz = zz/NN;
    double sum= 1.0 + 1.0/Ez(zz);
    int i,wt;
    for (wt=4, i=1; i<NN; ++i) {
      sum += 1.0/Ez(i*dz) * wt;
      wt   = 8/wt;
    }
    sum *= dz/3.0;
    return(sum * Lhub);
  }
  Spline<double> z_of_chi() const {
    // Fit the redshift-distance relation as a spline.
    const int    Nspl=1024;
    const double zmax=10.0;
    std::vector<double> zval,cval;
    try {
      zval.resize(Nspl);
      cval.resize(Nspl);
    } catch(std::exception& e) {myexception(e);}
    zval[0]=cval[0]=0;
    // Do it this way so we can make this loop OpenMP parallel.
#pragma omp parallel for shared(zval,cval)
    for (int i=1; i<Nspl; ++i) {
      zval[i] = 0.0 + i*(zmax-0.0)/Nspl;
      cval[i] = chi(zval[i]);
    }
    Spline<double> zchi(cval,zval);
    return(zchi);
  }
  void set(const double OmM0, const double hh0) {
    OmM     = OmM0;
    OmX     = 1-OmM0;
    OmK     = 0;
    w0      = -1;
    wa      = 0;
  }
  Cosmology() {}  // Don't do anything.
  Cosmology(const double OmM0, const double hh0) {
    set(OmM0,hh0);
  }
  ~Cosmology() {} // Don't need to do anything
};














std::vector<struct ObsObj>	convert2obs(const std::vector<float>& pos,
                                            const std::vector<float>& vel,
                                            const Spline<double>& zchi,
                                            const double boxside,
                                            const double hubble,
                                            const double offset[],
                                            const int ic,
                                            const double zcen,
                                            const double zmin,
                                            const double zmax,
                                            const Mangle::MaskClass& M) {
  // Compute the RA, DEC, z (obs), z (real-space) and distance modulus fields
  // for an array of 3D positions, assuming an offset from the origin given
  // by offset[] and an observation "corner" of the box given by 0<=ic<8,
  // Returns an array of ObsObj structs.
  // Positions are interpolated, linearly, assuming constant velocity.
  // The returned values are cut on zmin<zobs<zmax.
  // This could be made a method of a template class if helpful.
  if (pos.size()%3!=0 || pos.size()!=vel.size()) {
    std::cerr<<"Do not understand length of pos or vel."<<std::endl;
    std::cerr.flush();
    myexit(1);
  }
  const int nobj = pos.size()/3; // Assume only a few objects, fits in an int.
  std::vector<struct ObsObj> Objs;
  try {Objs.resize(nobj);} catch(std::exception& e) {myexception(e);}
#pragma omp parallel for shared(Objs)
  for (int nn=0; nn<nobj; ++nn) {
    double rr=1e-10,nhat[3],cpos[3],cvel[3];
    struct ObsObj O;
    O.indx = nn;
    // Assume loop-unrolling is enabled.
    for (int idim=0; idim<3; ++idim)	// Do the coordinate transformation
      if ( (ic&(1<<idim))>0 ) {		// for viewing from a corner.
        cpos[idim] = 1.0 - pos[3*nn+idim] + offset[idim];
        cvel[idim] = 0.0 - vel[3*nn+idim];
      }
      else {
        cpos[idim] = pos[3*nn+idim] + offset[idim];
        cvel[idim] = vel[3*nn+idim];
      }
    for (int idim=0; idim<3; ++idim)
      rr += cpos[idim]*cpos[idim];
    rr   = sqrt(rr);
    O.zr = zchi(rr*boxside,-1.0);		// Real-space "redshift".
    // Now change cpos, based on cvel, to account for the
    // difference in redshift between zcen and now.
    // Assume this is small enough one iteration is enough.
    double dlna=log( (1+zcen)/(1+O.zr) );
    rr = 1e-10;
    for (int idim=0; idim<3; ++idim) {
      cpos[idim] = periodic(cpos[idim]+cvel[idim]*dlna);
      rr += cpos[idim]*cpos[idim];
    }
    rr   = sqrt(rr);
    O.zr = zchi(rr*boxside,-1.0);		// Real-space "redshift".
    for (int idim=0; idim<3; ++idim)
      nhat[idim] = cpos[idim]/rr;
    double vlos=0;
    for (int idim=0; idim<3; ++idim)
      vlos += cvel[idim]*nhat[idim];
    for (int idim=0; idim<3; ++idim)
      cpos[idim] += vlos*nhat[idim];
    rr=1e-10;
    for (int idim=0; idim<3; ++idim)
      rr += cpos[idim]*cpos[idim];
    rr   = sqrt(rr);
    O.zz = zchi(rr*boxside,-1.0);		// Obs redshift.
    if (O.zz>zmin && O.zz<zmax) {
      O.ra = 180/M_PI*atan2(cpos[1],cpos[0]);	// RA  in degrees.
      if (O.ra<0) O.ra+=360.0;			// Periodically wrap.
      O.dec = 90-180/M_PI*acos(cpos[2]/rr);	// DEC in degrees.
      O.dmod= 25+5*log10(rr*boxside/hubble);	// Distance modulus.
      // This next step can be slow if the mask isn't pixelized.
      O.wt  = M.getweight((90.-O.dec)*M_PI/180.,O.ra*M_PI/180.);
    }
    else {
      O.ra = O.dec = O.dmod = 0;
      O.wt =-1;
    }
    // Put it into the array -- do it this way rather than using
    // push_back so we can do OpenMP parallel.
    Objs[nn] = O;
  }
  return(Objs);
}



void	write_objs(const std::vector<struct ObsObj>& Objs,
                   const char fname[]) {
  // Generate a FileHandler file containing the relevant fields in Objs.
  const int nobj = Objs.size();
  std::vector<float> fvec;
  std::vector<int>   ivec;
  try {
    fvec.resize(nobj);
    ivec.resize(nobj);
  } catch(std::exception& e) {myexception(e);}
  std::vector<long> ndims(1);  ndims[0]=nobj;
  // Do RA, DEC, Z, ZTRUE, DMOD, INDX.
  for (int i=0; i<nobj; ++i) fvec[i]=Objs[i].ra;
  FileHandler::write_float(fname,  "RA",fvec,ndims);
  for (int i=0; i<nobj; ++i) fvec[i]=Objs[i].dec;
  FileHandler::write_float(fname, "DEC",fvec,ndims);
  for (int i=0; i<nobj; ++i) fvec[i]=Objs[i].zz;
  FileHandler::write_float(fname,   "Z",fvec,ndims);
  for (int i=0; i<nobj; ++i) fvec[i]=Objs[i].zr;
  FileHandler::write_float(fname,  "ZR",fvec,ndims);
  for (int i=0; i<nobj; ++i) fvec[i]=Objs[i].dmod;
  FileHandler::write_float(fname,"DMOD",fvec,ndims);
  fvec.resize(0);
  for (int i=0; i<nobj; ++i) ivec[i]=Objs[i].indx;
  FileHandler::write_int(fname,"INDX",ivec,ndims);
}






int	main(int argc, char **argv)
{
  if (argc!=9) {
    std::cout<<"Usage: box2sky "
             <<"<OmM> <hub> <boxside> <zcen> <zmin> <zmax> <mask> <fbase>"
             <<std::endl;
    std::cout<<"--Box side should be in Mpc/h."<<std::endl;
    std::cout<<"--Mask should be an ascii, Mangle ply file."<<std::endl;
    myexit(1);
  }
  double OmM = atof(argv[1]);
  double hub = atof(argv[2]);
  double box = atof(argv[3]);
  double zcen= atof(argv[4]);
  double zmin= atof(argv[5]);
  double zmax= atof(argv[6]);
  std::cout<<"# Mock making code running on "<<omp_get_max_threads()
           <<" threads."<<std::endl;
  std::cout<<"# Using OmM="<<OmM<<", hub="<<hub<<" and box="<<box<<"Mpc/h."
           <<std::endl;
  std::cout<<"# Output is taken to be at z="<<zcen<<std::endl;
  std::cout<<"# Will cut output to "<<zmin<<"<z<"<<zmax<<std::endl;
  std::cout.flush();

  // Load the cosmology and the Mangle mask.
  Cosmology         C(OmM,hub);
  Spline<double> zchi=C.z_of_chi();
  Mangle::MaskClass M(argv[7]);
  std::cout<<"# Read mask containing "<<M.npolygons()
           <<" polygons from "<<argv[7]<<std::endl;
  std::cout.flush();

  // Read the phase space data for the objects.
  std::vector<long> ndims;
  std::vector<float> pos=FileHandler::read_float(argv[8],"pos",ndims);
  if (ndims.size()!=2 || ndims[1]!=3) {
    std::cout<<"Expecting to read (Nobj,3) array from "<<argv[8]<<"[pos]"
             <<std::endl;
    std::cout<<"Instead got dimensions: "<<std::endl;
    for (int i=0; i<ndims.size(); ++i)
      std::cout<<"  ndims["<<i<<"]="<<ndims[i]<<std::endl;
  }
  std::vector<float> vel=FileHandler::read_float(argv[8],"vel",ndims);
  if (ndims.size()!=2 || ndims[1]!=3) {
    std::cout<<"Expecting to read (Nobj,3) array from "<<argv[8]<<"[vel]"
             <<std::endl;
    std::cout<<"Instead got dimensions: "<<std::endl;
    for (int i=0; i<ndims.size(); ++i)
      std::cout<<"  ndims["<<i<<"]="<<ndims[i]<<std::endl;
  }
  std::cout<<"# Read pos/vel information for "<<ndims[0]
           <<" objects from "<<argv[8]<<std::endl;
  std::cout.flush();

  // Now we can subsample these galaxies by the above ratio and apply the mask.
  // The mask checking step can be slow, especially if the mask is not
  // pixelized, but by doing it this way we keep data size manageable.
  for (int ic=0; ic<8; ++ic) {	// For each octant...
    std::vector<struct ObsObj> Objs;
    try {Objs.reserve(Nres);} catch(std::exception& e) {myexception(e);}
    // Accumulate objects from all of the replications.
    double offset[3];
    for (int ix=-Nrepl; ix<=Nrepl; ++ix) {
      offset[0] = (double)ix;
      for (int iy=-Nrepl; iy<=Nrepl; ++iy) {
        offset[1] = (double)iy;
        for (int iz=-Nrepl; iz<=Nrepl; ++iz) {
          offset[2] = (double)iz;
          try {
            std::vector<struct ObsObj> O;
            O = convert2obs(pos,vel,zchi,box,hub,offset,ic,zcen,zmin,zmax,M);
            for (int i=0; i<O.size(); ++i)
              if (O[i].wt>0) {Objs.push_back(O[i]);}
          } catch(std::exception& e) {myexception(e);}
        }
      }
    }
    // Now write the results to file, labeled by octant.
    std::ostringstream ss;  ss<<argv[8]<<"_oct"<<ic;
    write_objs(Objs,ss.str().c_str());
  }

  std::cout<<"# Finished."<<std::endl;std::cout.flush();
  return(0);
}
