//-----------------------------------------------------------------------------
// Source: timedist.h
// Module: focimt
// Class for storing time-distance pairs.
//
// Copyright (c) 2003-2015, Grzegorz Kwiatek.
//
// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
//
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//-----------------------------------------------------------------------------
#ifndef timedistH
#define timedistH
//---------------------------------------------------------------------------

namespace Taquart {
  //! Time-distance variable pair structure.
  /*! This structure contains time and distance data pair.
   *  \ingroup Foci
   */
  class TimeDist {
    public:
      double Time;
      double Distance;
      TimeDist(double ATime, double ADistance);
      TimeDist(void);
      TimeDist(const TimeDist &Src);
      const TimeDist& operator=(const TimeDist &Src);
  };

  //! Time-distance comparison structure.
  /*! Structure used to compare time and distance values.
   *  \ingroup Foci
   */
  struct TimeDistComparator {
      bool operator()(const TimeDist&a, const TimeDist& b) {
        return a.Distance < b.Distance;
      }
  };
}

bool operator<(const Taquart::TimeDist& X, const Taquart::TimeDist& Y);
bool operator==(const Taquart::TimeDist& X, const Taquart::TimeDist& Y);

#endif
