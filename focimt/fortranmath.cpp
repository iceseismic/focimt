//-----------------------------------------------------------------------------
// fortranmath.cpp
//-----------------------------------------------------------------------------
#include <math.h>

//-----------------------------------------------------------------------------
#include "fortranmath.h"

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
double Taquart::alog(double Value) {
  return log(Value);
}

//-----------------------------------------------------------------------------
double Taquart::alog10(double Value) {
  return log10(Value);
}

//-----------------------------------------------------------------------------
double Taquart::amax1(double value1, double value2) {
  return value1 > value2 ? value1 : value2;
}

//-----------------------------------------------------------------------------
double Taquart::amin1(double value1, double value2) {
  return value1 > value2 ? value2 : value1;
}

//-----------------------------------------------------------------------------
double Taquart::sign(double value1, double value2) {
  if (value2 >= 0.0)
    return fabs(value1);
  else
    return -fabs(value1);
}

//-----------------------------------------------------------------------------
double Taquart::datan2(double y, double x) {
  double arctg;
  const double rdeg = 180.0 / M_PI;

  if (fabs(x) < 0.0001) {
    if (fabs(y) < 0.0001) {
      arctg = 0.0;
    }
    else
      arctg = y < 0.0 ? -90.0 : 90.0;
  }
  else if (x < 0.0)
    arctg = y < 0.0 ? atan(y / x) * rdeg - 180.0 : atan(y / x) * rdeg + 180.0;
  else
    arctg = atan(y / x) * rdeg;

  return (arctg);
}

