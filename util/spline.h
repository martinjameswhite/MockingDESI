#ifndef	__SPLINE_H_
#define	__SPLINE_H_

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


/*

A simple cubic spline class.

Author:		Martin White	(UCB,LBNL)

*/




extern void	myexit(const int flag);
extern void	myexception(std::exception& e);





template<class Ttype> class Spline {
private:
  std::vector<Ttype>	xarr,yarr,yderiv2;
  int			nlen;
public:
  ~Spline() {};
  Spline(const std::vector<Ttype>& x, const std::vector<Ttype>& y) {
    if (x.size() != y.size()) {
      std::cout << "x and y must have the same dimensions." << std::endl;
      myexit(1);
    }
    nlen = x.size();
    if (x[0]>=x[nlen-1]) {
      std::cout << "x must be monotonically increasing." << std::endl;
      myexit(1);
    }
    std::vector<Ttype> u;
    try {
      xarr = x;
      yarr = y;
      yderiv2.resize(nlen);
      u.resize(nlen);
    } catch(std::exception& e) {myexception(e);}
    yderiv2[0]=u[0]=0.0;
    for (int i=1;i<nlen-1;++i) {
      Ttype ss,pp;
      ss=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      pp=ss*yderiv2[i-1]+2.0;
      yderiv2[i]=(ss-1.0)/pp;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-ss*u[i-1])/pp;
    }
    yderiv2[nlen-1]=0.0;
    for (int k=nlen-2;k>=0;k--)
      yderiv2[k]=yderiv2[k]*yderiv2[k+1]+u[k];
  }
  bool inrange(const Ttype x) const {
    if (x<xarr[0] || x>xarr[nlen-1])
      return(false);
    else
      return(true);
  }
  double operator()(const Ttype x, const Ttype eflag) const {
    // Returns the value of the spline, or "eflag" if out of range.
    if (x<xarr[0] || x>xarr[nlen-1]) {
      return(eflag);
    }
    else {
      int ii,ilo=0,ihi=nlen-1;
      while (ihi-ilo > 1) {       // Bisection search.
        ii=(ihi+ilo)/2;
        if (xarr[ii] > x)
          ihi=ii;
        else
          ilo=ii;
      }
      Ttype h=xarr[ihi]-xarr[ilo];
      if (h<=0.0) {std::cout<<"Error in bisection."<< std::endl;myexit(1);}
      Ttype a=(xarr[ihi]-x)/h;
      Ttype b=(x-xarr[ilo])/h;
      return( a*yarr[ilo]+b*yarr[ihi]+
              (a*(a*a-1)*yderiv2[ilo]+b*(b*b-1)*yderiv2[ihi])*(h*h)/6.0 );
    }
  }
};



#endif
