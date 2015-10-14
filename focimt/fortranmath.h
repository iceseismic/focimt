//-----------------------------------------------------------------------------
#ifndef TRILIB_FORTRANMATH_H
#define TRILIB_FORTRANMATH_H
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// fortranmath.h
//
// rev.
//  1.0.1 (2007.09.27) datan2(.) added.
//  1.0.0 (2007.08.06) Initial release.
//-----------------------------------------------------------------------------

namespace Taquart {
  double alog(double x);
  double alog10(double x);
  double amax1(double x, double y);
  double amin1(double x, double y);
  double sign(double x, double y);
  double datan2(double y, double x);
}

//-----------------------------------------------------------------------------
#endif
