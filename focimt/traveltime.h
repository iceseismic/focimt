//---------------------------------------------------------------------------
#ifndef TRAVELTIME_H_
#define TRAVELTIME_H_
//---------------------------------------------------------------------------
#include <vector>

#define TT1D_RAYTRACE_MAXLAY (20)
#define TT1D_MAXT (100000.0)

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
