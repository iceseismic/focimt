//-----------------------------------------------------------------------------
// Source: traveltime.cpp
// Module: focimt
// 1D velocity model ray-tracing routines taken from HypoDD. Original code
// by Felix Waldhauser.
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
#ifndef TRAVELTIME_H_
#define TRAVELTIME_H_
//-----------------------------------------------------------------------------
#include <vector>

#define TT1D_RAYTRACE_MAXLAY (20)
#define TT1D_MAXT (100000.0)

//-----------------------------------------------------------------------------
void CalcTravelTime1D(double sta_elev, double depth, double delta,
    std::vector<double> Top, std::vector<double> Velocity, double &traveltime,
    double &takeoff, bool &directphase, double &aoi, int &kk, double &ray_dist);
void refract(const int nl, const double v[], const double vsq[],
    const double thk[], const int jl, const double tkj, const double delta,
    int &kk, double &tref, double &xovmax, double &ray_refract);
void tiddid(const int &jl, const int &nl, const double v[], const double vsq[],
    const double thk[], double tid[], double did[], double rid[]);
void ttime(const double &delta, const double &depth, const int &nl,
    const double v[], const double top[], double &t, double &ain,
    bool &directphase, double &aoi, int &kk, double &ray_dist);
void vmodel(const int& nl, const double v[], const double top[],
    const double &depth, double vsq[], double thk[], int &jl, double &tkj);
void direct1(const int &nl, const double v[], const double vsq[],
    const double thk[], const int& jl, const double& tkj, const double& delta,
    const double& depth, double& tdir, double& u, double& x,
    double &ray_direct);

//---------------------------------------------------------------------------
#endif /* TRAVELTIME_H_ */
