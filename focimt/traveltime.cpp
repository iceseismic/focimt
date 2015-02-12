//-----------------------------------------------------------------------------
#include <math.h>
//-----------------------------------------------------------------------------
#include "traveltime.h"

//-----------------------------------------------------------------------------
void CalcTravelTime1D(double sta_elev, double depth, double delta,
    std::vector<double> Top, std::vector<double> Velocity, double &traveltime,
    double &takeoff, bool &directphase, double &aoi, int &kk,
    double &ray_dist) {

  if (depth < 0.0) {
    depth = 0.0;
  }

  // The source depth cannot be on layer boundary.
  for (unsigned int i = 0; i < Top.size(); i++) {
    if (fabs(depth - Top[i]) < 0.0001) depth = depth - 0.001;
  }

  // Shift original velocity model according to station elevation.
  std::vector<double> TopAdj;
  TopAdj.clear();
  TopAdj.push_back(0.0);
  for (unsigned int k = 1; k < Top.size(); k++) {
    TopAdj.push_back(Top[k] + sta_elev);
  }

  // Fill in velocity model for travel time calculation.
  double VELOCITY[TT1D_RAYTRACE_MAXLAY];
  double TOP[TT1D_RAYTRACE_MAXLAY];
  VELOCITY[0] = 0.0;
  TOP[0] = 0.0;
  for (unsigned int i = 0; i < TopAdj.size(); i++) {
    VELOCITY[i + 1] = Velocity[i];
    TOP[i + 1] = TopAdj[i];
  }

  // Number of layers.
  const int nl = TopAdj.size();
  ttime(delta, depth + sta_elev, nl, VELOCITY, TOP, traveltime, takeoff,
      directphase, aoi, kk, ray_dist);
}

//-----------------------------------------------------------------------------
void ttime(const double &delta, const double &depth, const int &nl,
    const double v[], const double top[], double &t, double &takeoff,
    bool &directphase, double &aoi, int &kk, double &ray_dist) {
  int jl = 0;
  double tdir = 0.0;
  double thk[TT1D_RAYTRACE_MAXLAY];
  double tkj = 0.0;
  double tref = 0.0;
  double u = 0.0;
  double vsq[TT1D_RAYTRACE_MAXLAY];
  double x = 0.0;
  double xovmax = 0.0;
  double ray_dist_temp = 0.0; // Reset ray distance to 0.0f.

  kk = 0;
  directphase = false;

  vmodel(nl, v, top, depth, vsq, thk, jl, tkj);
  refract(nl, v, vsq, thk, jl, tkj, delta, kk, tref, xovmax, ray_dist_temp);
  t = tref;
  ray_dist = ray_dist_temp;
  if (kk > 0) {
    u = v[jl] / v[kk];
    takeoff = asin(u) * 180.0 / M_PI;
    aoi = asin(v[1] / v[kk]) * 180.0 / M_PI;
  }
  if (delta <= xovmax) {
    double ray_dist_temp2 = 0.0;
    direct1(nl, v, vsq, thk, jl, tkj, delta, depth, tdir, u, x, ray_dist_temp2);
    if (tref > tdir) {
      t = tdir;
      ray_dist = ray_dist_temp2;
      takeoff = 180.0 - asin(u) * 180.0 / M_PI;
      aoi = asin(sin(takeoff * M_PI / 180.0) * v[1] / v[jl]) * 180.0 / M_PI;
      directphase = true;
    }
  }
}

//-----------------------------------------------------------------------------
void vmodel(const int& nl, const double v[], const double top[],
    const double &depth, double vsq[], double thk[], int &jl, double &tkj) {
  int i = 0;
  for (i = 1; i <= nl; i++)
    vsq[i] = v[i] * v[i];
  jl = nl;
  for (i = 1; i <= nl; i++) {
    if (depth <= top[i]) {
      jl = i - 1;
      break;
    }
  }
  for (i = 1; i <= nl - 1; i++) {
    double a = top[i + 1] - top[i];
    thk[i] = a;
  }
  tkj = depth - top[jl];
}

//-----------------------------------------------------------------------------
void direct1(const int &nl, const double v[], const double vsq[],
    const double thk[], const int& jl, const double& tkj, const double& delta,
    const double& depth, double& tdir, double& u, double& x,
    double &ray_direct) {
  double r = 0.0, tklmax = 0.0, usq = 0.0;
  int lmax = 0, j1 = 0, l = 0;
  double vlmax = 0, del = 0.0, xtest = 0.0;
  double ua = 0.0, ub = 0.0, uasq = 0.0, ubsq = 0.0, xa = 0.0, xb = 0.0, dela =
      0.0, delb = 0.0, ubdiv = 0.0;
  ray_direct = 0.0;
  if (jl == 1) {
    r = hypot(depth, delta);
    tdir = r / v[1];
    u = delta / r;
    x = delta;
    ray_direct = r;
    return;
  }
  lmax = jl;
  tklmax = tkj;
  vlmax = v[jl];
  j1 = jl - 1;
  for (l = 1; l <= j1; l++) {
    if (v[l] > vlmax) {
      lmax = l;
      tklmax = thk[l];
      vlmax = v[l];
    }
  }
  if (tklmax <= 0.05) tklmax = 0.05;

  ua = (v[jl] / vlmax) * delta / sqrt(delta * delta + depth * depth);
  ub = (v[jl] / vlmax) * delta / sqrt(delta * delta + tklmax * tklmax);
  uasq = ua * ua;
  ubsq = ub * ub;
  if (ubsq >= 1.0f) ubsq = 0.99999;
  if (uasq >= 1.0f) uasq = 0.99999;
  xa = tkj * ua / sqrt(1.0 - uasq);
  if (lmax == jl) {
    xb = delta;
  }
  else {
    xb = tkj * ub / sqrt(1.0 - ubsq);
  }
  dela = xa;
  delb = xb;
  for (l = 1; l <= j1; l++) {
    dela = dela + thk[l] * ua / sqrt(vsq[jl] / vsq[l] - uasq);
    ubdiv = sqrt(vsq[jl] / vsq[l] - ubsq);
    if (ubdiv <= 1.0e-20) {
      ubdiv = 1.0e-20;
    }
    delb = delb + thk[l] * ub / sqrt(vsq[jl] / vsq[l] - ubsq);
  }
  for (int kount = 1; kount <= 25; kount++) {
    if (((delb - dela) < 0.02)) {
      x = 0.5 * (xa + xb);
      u = x / sqrt(x * x + tkj * tkj);
      usq = u * u;
      goto L23193;
    }
    x = xa + (delta - dela) * (xb - xa) / (delb - dela);
    u = x / sqrt(x * x + tkj * tkj);
    usq = u * u;
    del = x;
    for (int l = 1; l <= j1; l++) {
      del = del + thk[l] * u / sqrt(vsq[jl] / vsq[l] - usq);
    }
    xtest = del - delta;
    if (fabs(xtest) < 0.02) goto L23193;
    if (xtest < 0.0) {
      xa = x;
      dela = del;
    }
    else {
      xb = x;
      delb = del;
    }
  } // 23192
  L23193: tdir = sqrt(x * x + tkj * tkj) / v[jl];
  ray_direct = sqrt(x * x + tkj * tkj);

  for (int l = 1; l <= j1; l++) {
    tdir = tdir + thk[l] * v[jl] / (vsq[l] * sqrt(vsq[jl] / vsq[l] - usq));
    ray_direct = ray_direct
        + thk[l] * v[jl] / (v[l] * sqrt(vsq[jl] / vsq[l] - usq));
  }

  // PM&GK: It is not clear what's the purpose of the following lines!
  tdir = tdir - (u / v[jl]) * (del - delta);
  ray_direct = ray_direct - u * (del - delta);
}

//-----------------------------------------------------------------------------
void refract(const int nl, const double v[], const double vsq[],
    const double thk[], const int jl, const double tkj, const double delta,
    int &kk, double &tref, double &xovmax, double &ray_refract) {

  // tkj - depth of event in certain layer.
  // delta - horizontal distance between event and receiver
  double tid[TT1D_RAYTRACE_MAXLAY];
  double did[TT1D_RAYTRACE_MAXLAY]; // Terms in critical distance independent of tkj
  double rid[TT1D_RAYTRACE_MAXLAY];
  double tinj[TT1D_RAYTRACE_MAXLAY]; // Travel time intercept
  double didj[TT1D_RAYTRACE_MAXLAY]; // Critical distance
  double rinj[TT1D_RAYTRACE_MAXLAY];
  double tr[TT1D_RAYTRACE_MAXLAY]; // Travel time for refraction in layer X
  double rr[TT1D_RAYTRACE_MAXLAY];
  int j1 = 0;
  int jx = 0;
  int l = 0;
  int lx = 0;
  int m = 0;
  int m1 = 0;
  double sqt = 0.0;
  double tim = 0.0;

  tiddid(jl, nl, v, vsq, thk, tid, did, rid);
  tref = TT1D_MAXT;
  ray_refract = TT1D_MAXT;
  j1 = jl + 1;
  for (m = j1; m <= nl; m++) {
    if (tid[m] == TT1D_MAXT) {
      tr[m] = TT1D_MAXT;
      rr[m] = TT1D_MAXT;
    }
    else {
      sqt = sqrt(vsq[m] - vsq[jl]);
      tinj[m] = tid[m] - tkj * sqt / (v[m] * v[jl]);
      didj[m] = did[m] - tkj * v[jl] / sqt;
      rinj[m] = rid[m] - tkj * v[m] / sqt;
      tr[m] = tinj[m] + delta / v[m];
      rr[m] = rinj[m] + delta - didj[m];
      if (didj[m] > delta) {
        tr[m] = TT1D_MAXT;
        rr[m] = TT1D_MAXT;
      }
    }
    if (tr[m] < tref) {
      tref = tr[m];
      ray_refract = rr[m];
      kk = m;
    } //23157 continue
  } // 23151

  if (tref == TT1D_MAXT) {
    xovmax = TT1D_MAXT;
    kk = 0;
    return;
  } //23159 continue
  m = jl + 1;
  while (1) {
    if (!(tid[m] == TT1D_MAXT)) break;
    m = m + 1;
  }
  lx = m;

  if (jl == 1) {
    xovmax = tinj[lx] * v[lx] * v[1] / (v[lx] - v[1]);
    return;
  }
  m = jl;
  l23165: tid[m] = 0.0f;
  m1 = m - 1;
  for (l = 1; l <= m1; l++) {
    if (vsq[m] <= vsq[l]) {
      tid[m] = TT1D_MAXT;
    }
    else {
      sqt = sqrt(vsq[m] - vsq[l]);
      tim = thk[l] * sqt / (v[l] * v[m]);
      tid[m] = tid[m] + tim;
    }
  }
  m = m - 1;

  if (!(tid[m + 1] < TT1D_MAXT || m == 1)) goto l23165;

  if (tid[m + 1] < TT1D_MAXT) {
    jx = m + 1;
    xovmax = (tinj[lx] - tid[jx]) * v[lx] * v[jx] / (v[lx] - v[jx]);
  }
  else {
    xovmax = tinj[lx] * v[lx] * v[1] / (v[lx] - v[1]);
  }
}

//-----------------------------------------------------------------------------
void tiddid(const int &jl, const int &nl, const double v[], const double vsq[],
    const double thk[], double tid[], double did[], double rid[]) {
  double did1 = 0.0, did2 = 0.0, dimm = 0.0, rid1 = 0.0, rid2 = 0.0;
  int j1 = 0, l = 0, m = 0, m1 = 0;
  double sqt = 0.0, tid1 = 0.0, tid2 = 0.0, tim = 0.0, rim = 0.0;
  j1 = jl + 1;
  for (m = j1; m <= nl; m++) {
    tid[m] = 0.0;
    did[m] = 0.0;
    tid1 = 0.0;
    tid2 = 0.0;
    did1 = 0.0;
    did2 = 0.0;
    rid1 = 0.0;
    rid2 = 0.0;
    m1 = m - 1;
    for (l = 1; l <= m1; l++) {
      if ((vsq[m] <= vsq[l])) {
        tid[m] = 100000.0;
        did[m] = 100000.0;
        rid[m] = 100000.0;
        continue;
      } //23178 continue
      sqt = sqrt(vsq[m] - vsq[l]);
      tim = thk[l] * sqt / (v[l] * v[m]);
      dimm = thk[l] * v[l] / sqt;
      rim = thk[l] * v[m] / sqt;
      if ((l < jl)) {
        tid1 = tid1 + tim;
        did1 = did1 + dimm;
        rid1 = rid1 + rim;
      }
      else {
        tid2 = tid2 + tim;
        did2 = did2 + dimm;
        rid2 = rid2 + rim;
      }
    } //23176 continue

    if ((tid[m] != 100000.0)) {
      tid[m] = tid1 + 2 * tid2;
      did[m] = did1 + 2 * did2;
      rid[m] = rid1 + 2 * rid2;
    }
  } //23174 continue
}

//-----------------------------------------------------------------------------



