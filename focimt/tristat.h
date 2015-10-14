//-----------------------------------------------------------------------------
#ifndef TRILIB_TRISTAT_H
#define TRILIB_TRISTAT_H
//-----------------------------------------------------------------------------

#include "exceptions.h"

//-----------------------------------------------------------------------------
namespace Taquart {
  //! Returns the sample standard deviation for elements in an array.
  /*! Trinity::std calculates the sample standard deviation (the square root
   *  of the sample variance) of all values in the \p X array parameter.
   *  \p Size indicates the number of elements in the array.
   *  \param X Pointer to the array.
   *  \param Size Array size.
   *  \return Standard deviation value.
   *  \ingroup trilib
   */
  double std(double * X, unsigned int Size) throw (Taquart::TriENullPointer,
      Taquart::TriEOutOfRange);

  //! Returns the average of all values in an array.
  /*! Trinity::mean calculates the arithmetic average of all the values in
   *  the \p X array parameter. The \p Size parameter gives the number
   *  of array elements.
   *  \param X Pointer to the array.
   *  \param Size Array size.
   *  \return Mean value.
   *  \ingroup trilib
   */
  double mean(double * X, unsigned int Size) throw (Taquart::TriENullPointer,
      Taquart::TriEOutOfRange);
}

//-----------------------------------------------------------------------------
#endif
