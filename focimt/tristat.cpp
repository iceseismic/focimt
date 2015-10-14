//-----------------------------------------------------------------------------
#include <math.h>
#include "tristat.h"
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------------
double Taquart::mean(double * X, unsigned int Size)
    throw (Taquart::TriENullPointer, Taquart::TriEOutOfRange) {
  if (X == 0L)
    throw Taquart::TriENullPointer("Trinity:mean(): Input pointer is NULL.");
  else if (Size == 0)
    throw Taquart::TriEOutOfRange(
        "Trinity:mean(): Size of input vector is zero.");
  else {
    double s = 0.0f;
    for (unsigned int i = 0; i < Size; i++)
      s += *(X + i);
    return s / double(Size);
  }
}

//---------------------------------------------------------------------------
double Taquart::std(double * X, unsigned int Size)
    throw (Taquart::TriENullPointer, Taquart::TriEOutOfRange) {
  if (X == 0L)
    throw Taquart::TriENullPointer("Trinity:std(): Input pointer is NULL.");
  else if (Size == 0)
    throw Taquart::TriEOutOfRange(
        "Trinity:std(): Size of input vector is zero.");
  else if (Size == 1)
    return 0.0f;
  else {
    const double MeanValue = mean(X, Size);
    double s = 0.0f;
    for (unsigned int i = 0; i < Size; i++)
      s += ((*(X + i) - MeanValue) * (*(X + i) - MeanValue));
    return sqrt(s / double(Size - 1));
  }
}

//---------------------------------------------------------------------------

