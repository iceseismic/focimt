//-----------------------------------------------------------------------------
#include <math.h>
//-----------------------------------------------------------------------------
#include "georoutines.h"
#include "fortranmath.h"

//-----------------------------------------------------------------------------
void Taquart::StrikeDipRake2MT(double S, double D, double R, double &M11,
    double &M22, double &M33, double &M12, double &M13, double &M23) {
  M11 = -0.0001
      + -1.0
          * (sin(D) * cos(R) * sin(2.0 * S)
              + sin(2.0 * D) * sin(R) * sin(S) * sin(S));
  M22 = -0.0001
      + (sin(D) * cos(R) * sin(2.0 * S)
          - sin(2.0 * D) * sin(R) * cos(S) * cos(S));
  M33 = 0.0002 + -1.0 * (M11 + M22);
  M12 = -0.0001
      + (sin(D) * cos(R) * cos(2 * S)
          + 0.5 * sin(2.0 * D) * sin(R) * sin(2.0 * S)); /* Mxy */
  M13 = 0.0001
      + -1.0 * (cos(D) * cos(R) * cos(S) + cos(2 * D) * sin(R) * sin(S));
  M23 = -0.0001
      + -1.0 * (cos(D) * cos(R) * sin(S) - cos(2 * D) * sin(R) * cos(S));
}

//-----------------------------------------------------------------------------
double Taquart::zero_360(double str) {
  if (str >= 360.0)
    str -= 360.0;
  else if (str < 0.0)
    str += 360.0;
  return (str);
}

//-----------------------------------------------------------------------------
void Taquart::sincosd(double i, double *s, double *c) {
  *s = sin(i * DEG2RAD);
  *c = cos(i * DEG2RAD);
}

//-------------------------------------------------------------------------------------------------
// TRICAIRO
void Taquart::sincos(double i, double *s, double *c) {
  *s = sin(i);
  *c = cos(i);
}

//-------------------------------------------------------------------------------------------------
double Taquart::computed_dip1(nodal_plane NP1) {
  /* Genevieve Patau */
  /* Compute second nodal plane dip when are given strike, dip and rake for
   the first nodal plane with AKI & RICHARD's convention.
   Angles are in degrees. */

  double am = NP1.rake == 0.0 ? 1.0 : NP1.rake / fabs(NP1.rake);
  return acos(am * sind(NP1.rake) * sind(NP1.dip)) / DEG2RAD;

  // Last tested: $Id: utilmeca.c 12822 2014-01-31 23:39:56Z remko $ OK
}

/* GMT 5.1.2
 double computed_dip1 (struct nodal_plane NP1)
 {
 Compute second nodal plane dip when are given strike,
 dip and rake for the first nodal plane with AKI & RICHARD's
 convention.  Angles are in degrees.
 Genevieve Patau

 double am = (GMT_IS_ZERO (NP1.rake) ? 1.0 : NP1.rake / fabs (NP1.rake));
 double dip2;

 dip2 = acosd (am * sind (NP1.rake) * sind (NP1.dip));

 return (dip2);
 }
 */

//-------------------------------------------------------------------------------------------------
double Taquart::computed_rake1(struct nodal_plane NP1) {
  /* Compute rake in the second nodal plane when strike ,dip
   and rake are given for the first nodal plane with AKI &
   RICHARD's convention.
   Angles are in degrees. */
  /* Genevieve Patau */

  double sinrake2;
  double str2 = Taquart::computed_strike1(NP1);
  double dip2 = Taquart::computed_dip1(NP1);
  double am = NP1.rake == 0.0 ? 1.0 : NP1.rake / fabs(NP1.rake);
  double sd, cd, ss, cs;
  double sd2, cd2;
  Taquart::sincos(NP1.dip * DEG2RAD, &sd, &cd);
  Taquart::sincos(dip2 * DEG2RAD, &sd2, &cd2); /*  Onur 25Sep02,  */
  Taquart::sincos((NP1.str - str2) * DEG2RAD, &ss, &cs);

  if (fabs(dip2 - 90.0) < EPSIL)
    sinrake2 = am * cd;
  else
    sinrake2 = -am * sd * cs / cd2; /* ONUR: cd2 [cos(DIP2)] must be used not cd [cos(DIP1)] */
  // This part is different in GMT 5.1, as there cd is used instead of cd2

  return Taquart::datan2(sinrake2, -am * sd * ss);
}

/* GMT 5.1.2
 double computed_rake1 (struct nodal_plane NP1)
 {

 Compute rake in the second nodal plane when strike ,dip
 and rake are given for the first nodal plane with AKI &
 RICHARD's convention.  Angles are in degrees.

 Genevieve Patau


 double computed_strike1(), computed_dip1();
 double rake2, sinrake2;
 double str2 = computed_strike1(NP1);
 double dip2 = computed_dip1(NP1);
 double am = (GMT_IS_ZERO (NP1.rake) ? 1.0 : NP1.rake / fabs (NP1.rake));
 double sd, cd, ss, cs;
 sincosd (NP1.dip, &sd, &cd);
 sincosd (NP1.str - str2, &ss, &cs);

 if (fabs (dip2 - 90.0) < EPSIL)
 sinrake2 = am * cd;
 else
 sinrake2 = -am * sd * cs / cd;

 rake2 = d_atan2d (sinrake2, -am * sd * ss);

 return (rake2);
 }
 */

//-------------------------------------------------------------------------------------------------
double Taquart::computed_strike1(struct nodal_plane NP1) {
  /* Compute the strike of the decond nodal plane when are given strike, dip and rake for the first 
   nodal plane with AKI & RICHARD's convention. Angles are in degrees. */
  /* Genevieve Patau */

  double str2;
  double cd1 = cosd(NP1.dip);
  double temp;
  double cp2, sp2;
  double am = NP1.rake == 0.0 ? 1.0 : NP1.rake / fabs(NP1.rake);
  double ss, cs, sr, cr;

  Taquart::sincos(NP1.rake * DEG2RAD, &sr, &cr);
  Taquart::sincos(NP1.str * DEG2RAD, &ss, &cs);
  if (cd1 < EPSIL && fabs(cr) < EPSIL) {
    str2 = NP1.str + 180.0;
  } else {
    temp = cr * cs;
    temp += sr * ss * cd1;
    sp2 = -am * temp;
    temp = ss * cr;
    temp -= sr * cs * cd1;
    cp2 = am * temp;
    str2 = datan2(sp2, cp2);
    str2 = Taquart::zero_360(str2);
  }
  return (str2);

}
/* GMT 5.1.2
 double computed_strike1 (struct nodal_plane NP1)
 {

 Compute the strike of the decond nodal plane when are given
 strike, dip and rake for the first nodal plane with AKI & RICHARD's
 convention.  Angles are in degrees.
 Genevieve Patau


 double str2, temp, cp2, sp2, ss, cs, sr, cr;
 double cd1 = cosd (NP1.dip);
 double am = (GMT_IS_ZERO (NP1.rake) ? 1. : NP1.rake /fabs (NP1.rake));

 sincosd (NP1.rake, &sr, &cr);
 sincosd (NP1.str, &ss, &cs);
 if (cd1 < EPSIL && fabs (cr) < EPSIL) {
 #if 0
 GMT_Report (API, GMT_MSG_DEBUG, "\nThe second plane is horizontal;");
 GMT_Report (API, GMT_MSG_DEBUG, "\nStrike is undetermined.");
 GMT_Report (API, GMT_MSG_DEBUG, "\nstr2 = NP1.str + 180. is taken to define");
 GMT_Report (API, GMT_MSG_DEBUG, "\nrake in the second plane.\n");
 #endif
 str2 = NP1.str + 180.0;
 }
 else {
 temp = cr * cs;
 temp += sr * ss * cd1;
 sp2 = -am * temp;
 temp = ss * cr;
 temp -= sr *  cs * cd1;
 cp2 = am * temp;
 str2 = d_atan2d(sp2, cp2);
 str2 = zero_360(str2);
 }
 return (str2);
 }
 */

//-------------------------------------------------------------------------------------------------
// ONUR
void Taquart::define_second_plane(struct nodal_plane NP1,
    struct nodal_plane *NP2) {
  /* Compute strike, dip, slip for the second nodal plane
   when are given strike, dip and rake for the first one. */
  /*  Genevieve Patau */

  NP2->str = Taquart::computed_strike1(NP1);
  NP2->dip = Taquart::computed_dip1(NP1);
  NP2->rake = Taquart::computed_rake1(NP1);
}

//-------------------------------------------------------------------------------------------------
// ONUR
double Taquart::null_axis_dip(double str1, double dip1, double str2,
    double dip2) {
  /* Compute null axis dip when strike and dip are given for each 
   nodal plane. Angles are in degrees. */
  /* Genevieve Patau */

  double den = asin(sind(dip1) * sind(dip2) * sind(str1 - str2)) / DEG2RAD;
  if (den < 0.0)
    den = -den;
  return (den);
}

/* GMT 5.1.2
 double null_axis_dip (double str1, double dip1, double str2, double dip2)
 {
 compute null axis dip when strike and dip are given
 for each nodal plane.  Angles are in degrees.

 Genevieve Patau

 double den;

 den = asind (sind (dip1) * sind (dip2) * sind (str1 - str2));
 if (den < 0.) den = -den;
 return (den);
 }
 */

//-------------------------------------------------------------------------------------------------
double Taquart::null_axis_strike(double str1, double dip1, double str2,
    double dip2) {
  /* Compute null axis strike when strike and dip are given for each 
   nodal plane. Angles are in degrees. */

  /* Genevieve Patau */
  double phn, cosphn, sinphn;
  double sd1, cd1, sd2, cd2, ss1, cs1, ss2, cs2;

  Taquart::sincos(dip1 * DEG2RAD, &sd1, &cd1);
  Taquart::sincos(dip2 * DEG2RAD, &sd2, &cd2);
  Taquart::sincos(str1 * DEG2RAD, &ss1, &cs1);
  Taquart::sincos(str2 * DEG2RAD, &ss2, &cs2);

  cosphn = sd1 * cs1 * cd2 - sd2 * cs2 * cd1;
  sinphn = sd1 * ss1 * cd2 - sd2 * ss2 * cd1;
  if (sind(str1 - str2) < 0.0) {
    cosphn = -cosphn;
    sinphn = -sinphn;
  }
  phn = datan2(sinphn, cosphn);
  if (phn < 0.0)
    phn += 360.0;
  return (phn);
}

/* GMT 5.1.2
 double null_axis_strike (double str1, double dip1, double str2, double dip2)
 {

 Compute null axis strike when strike and dip are given
 for each nodal plane.   Angles are in degrees.

 Genevieve Patau


 double phn, cosphn, sinphn, sd1, cd1, sd2, cd2, ss1, cs1, ss2, cs2;

 sincosd (dip1, &sd1, &cd1);
 sincosd (dip2, &sd2, &cd2);
 sincosd (str1, &ss1, &cs1);
 sincosd (str2, &ss2, &cs2);

 cosphn = sd1 * cs1 * cd2 - sd2 * cs2 * cd1;
 sinphn = sd1 * ss1 * cd2 - sd2 * ss2 * cd1;
 if (sind(str1 - str2) < 0.0) {
 cosphn = -cosphn;
 sinphn = -sinphn;
 }
 phn = d_atan2d(sinphn, cosphn);
 if (phn < 0.0) phn += 360.0;
 return (phn);
 }
 */

//-------------------------------------------------------------------------------------------------
double Taquart::computed_rake2(double str1, double dip1, double str2,
    double dip2, double fault) {
  /* Genevieve Patau */
  /* Compute rake in the second nodal plane when strikes and dips for the first
   and the second nodal plane are given with an additional variabled
   characterizing the fault: +1.0 for the inverse fault, -1.0 for the normal
   fault. Angles are in degrees. */

  double sinrake2;
  double sd, cd2, ss, cs;
  double cd;

  Taquart::sincos((str1 - str2) * DEG2RAD, &ss, &cs);
  sd = sind(dip1);
  cd = cosd(dip1); /* cos(dip1) not dip2   Onur Tan*/
  cd2 = cosd(dip2);

  if (fabs(dip2 - 90.0) < EPSIL)
    sinrake2 = fault * cd; // This part is different in GMT as cos(dip2) is used instead of cos(dip1)
  else
    sinrake2 = -fault * sd * cs / cd2; // This is the same as in GMT 5.1

  return Taquart::datan2(sinrake2, -fault * sd * ss);

}

/* GMT 5.1.2
 double computed_rake2 (double str1, double dip1, double str2, double dip2, double fault)
 {
 Compute rake in the second nodal plane when strike and dip
 for first and second nodal plane are given with a double
 characterizing the fault :
 +1. inverse fault
 -1. normal fault.
 Angles are in degrees.

 Genevieve Patau

 double rake2, sinrake2, sd, cd, ss, cs;

 sincosd (str1 - str2, &ss, &cs);

 sd = sind(dip1);        cd = cosd(dip2);
 if (fabs (dip2 - 90.0) < EPSIL)
 sinrake2 = fault * cd;
 else
 sinrake2 = -fault * sd * cs / cd;

 rake2 = d_atan2d (sinrake2, - fault * sd * ss);

 return (rake2);
 }
 */

//-------------------------------------------------------------------------------------------------
/* GMT 5.1.2
 void dc2axe (st_me meca, struct AXIS *T, struct AXIS *N, struct AXIS *P)
 {

 From FORTRAN routines of Anne Deschamps :
 compute azimuth and plungement of P-T axis
 from nodal plane strikes, dips and rakes.


 double cd1, sd1, cd2, sd2, cp1, sp1, cp2, sp2;
 double amz, amx, amy, dx, px, dy, py;

 cd1 = cosd (meca.NP1.dip) * M_SQRT2;
 sd1 = sind (meca.NP1.dip) * M_SQRT2;
 cd2 = cosd (meca.NP2.dip) * M_SQRT2;
 sd2 = sind (meca.NP2.dip) * M_SQRT2;
 cp1 = - cosd (meca.NP1.str) * sd1;
 sp1 = sind (meca.NP1.str) * sd1;
 cp2 = - cosd (meca.NP2.str) * sd2;
 sp2 = sind (meca.NP2.str) * sd2;

 amz = - (cd1 + cd2);
 amx = - (sp1 + sp2);
 amy = cp1 + cp2;
 dx = atan2d (hypot(amx, amy), amz) - 90.0;
 px = atan2d (amy, -amx);
 if (px < 0.0) px += 360.0;
 if (dx < EPSIL) {
 if (px > 90.0 && px < 180.0) px += 180.0;
 if (px >= 180.0 && px < 270.0) px -= 180.0;
 }

 amz = cd1 - cd2;
 amx = sp1 - sp2;
 amy = - cp1 + cp2;
 dy = atan2d (hypot(amx, amy), -fabs(amz)) - 90.0;
 py = atan2d (amy, -amx);
 if (amz > 0.0) py -= 180.0;
 if (py < 0.0) py += 360.0;
 if (dy < EPSIL) {
 if (py > 90.0 && py < 180.0) py += 180.0;
 if (py >= 180.0 && py < 270.0) py -= 180.0;
 }

 if (meca.NP1.rake > 0.0) {
 P->dip = dy; P->str = py;
 T->dip = dx; T->str = px;
 }
 else {
 P->dip = dx; P->str = px;
 T->dip = dy; T->str = py;
 }

 N->str = null_axis_strike (T->str, T->dip, P->str, P->dip);
 N->dip = null_axis_dip (T->str, T->dip, P->str, P->dip);
 }
 */

/*
 void Trinity::dc_to_axe(st_me meca,struct AXIS *T,struct AXIS *N,struct AXIS *P)
 {
 // From FORTRAN routines of Anne Deschamps : compute azimuth and
 //   plungement of P-T axis from nodal plane strikes, dips and rakes.

 int im;
 int pure_strike_slip = 0;
 double cd1, sd1, cd2, sd2;
 double cp1, sp1, cp2, sp2;
 double amz, amx, amy, dx, px, dy, py;

 if(fabs(sind(meca.NP1.rake)) > EPSIL) im = (int) (meca.NP1.rake / fabs(meca.NP1.rake));
 else if(fabs(sind(meca.NP2.rake)) > EPSIL) im = (int) (meca.NP2.rake / fabs(meca.NP2.rake));
 else pure_strike_slip = 1;

 if(pure_strike_slip)
 {
 if(cosd(meca.NP1.rake) < 0.)
 {
 P->str = zero_360(meca.NP1.str + 45.);
 T->str = zero_360(meca.NP1.str - 45.);
 }
 else
 {
 P->str = zero_360(meca.NP1.str - 45.);
 T ->str = zero_360(meca.NP1.str + 45.);
 }
 P->dip = 0.;
 T->dip = 0.;
 }
 else
 {
 cd1 = cosd(meca.NP1.dip) *  M_SQRT2;
 sd1 = sind(meca.NP1.dip) *  M_SQRT2;
 cd2 = cosd(meca.NP2.dip) *  M_SQRT2;
 sd2 = sind(meca.NP2.dip) *  M_SQRT2;
 cp1 = - cosd(meca.NP1.str) * sd1;
 sp1 = sind(meca.NP1.str) * sd1;
 cp2 = - cosd(meca.NP2.str) * sd2;
 sp2 = sind(meca.NP2.str) * sd2;

 amz = - (cd1 + cd2);
 amx = - (sp1 + sp2);
 amy = cp1 + cp2;
 dx = atan2(sqrt(amx * amx + amy * amy), amz) - M_PI_2;
 px = atan2(amy, - amx);
 if(px < 0.)
 px += TWO_PI;

 amz = cd1 - cd2;
 amx = sp1 - sp2;
 amy = - cp1 + cp2;
 dy = atan2(sqrt(amx * amx + amy * amy), - fabs(amz)) - M_PI_2;
 py = atan2(amy, - amx);
 if(amz > 0.)
 py -= M_PI;

 if(py < 0.)
 py += TWO_PI;

 if(im == 1)
 {
 P->dip = dy;
 P->str = py;
 T->dip = dx;
 T->str = px;
 }
 else
 {
 P->dip = dx;
 P->str = px;
 T->dip = dy;
 T->str = py;
 }
 }

 T->str /= D2R;
 T->dip /= D2R;
 P->str /= D2R;
 P->dip /= D2R;

 N->str =  null_axis_strike(T->str, T->dip, P->str, P->dip);
 N->dip =  null_axis_dip(T->str, T->dip, P->str, P->dip);
 }


 //-------------------------------------------------------------------------------------------------
 */
 void Taquart::axe2dc(AXIS T, AXIS P, nodal_plane *NP1,
 nodal_plane *NP2)
 {
 double pp, dp, pt, dt;
 double p1, d1, p2, d2;
 double PII = M_PI * 2.;
 double cdp, sdp, cdt, sdt;
 double cpt, spt, cpp, spp;
 double amz, amy, amx;
 double im;

 pp = P.str * DEG2RAD; dp = P.dip * DEG2RAD;
 pt = T.str * DEG2RAD; dt = T.dip * DEG2RAD;

 sincos (dp, &sdp, &cdp);
 sincos (dt, &sdt, &cdt);
 sincos (pt, &spt, &cpt);
 sincos (pp, &spp, &cpp);

 cpt *= cdt; spt *= cdt;
 cpp *= cdp; spp *= cdp;

 amz = sdt + sdp; amx = spt + spp; amy = cpt + cpp;
 d1 = atan2(sqrt(amx*amx + amy*amy), amz);
 p1 = atan2(amy, -amx);
 if(d1 > M_PI_2) {
 d1 = M_PI - d1;
 p1 += M_PI;
 if(p1 > PII) p1 -= PII;
 }
 if(p1 < 0.) p1 += PII;

 amz = sdt - sdp; amx = spt - spp; amy = cpt - cpp;
 d2 = atan2(sqrt(amx*amx + amy*amy), amz);
 p2 = atan2(amy, -amx);
 if(d2 > M_PI_2) {
 d2 = M_PI - d2;
 p2 += M_PI;
 if(p2 > PII) p2 -= PII;
 }
 if(p2 < 0.) p2 += PII;

 NP1->dip = d1 / DEG2RAD; NP1->str = p1 / DEG2RAD;
 NP2->dip = d2 / DEG2RAD; NP2->str = p2 / DEG2RAD;

 im = 1;
 if(dp > dt) im = -1;
 NP1->rake = computed_rake2(NP2->str,NP2->dip,NP1->str,NP1->dip,im);
 NP2->rake = computed_rake2(NP1->str,NP1->dip,NP2->str,NP2->dip,im);
 }
 /*

/* GMT 5.1.2
 void axe2dc (struct AXIS T, struct AXIS P, struct nodal_plane *NP1, struct nodal_plane *NP2)
 {

 Calculate double couple from principal axes.
 Angles are in degrees.

 Genevieve Patau, 16 juin 1997


 double p1, d1, p2, d2;
 double cdp, sdp, cdt, sdt, cpt, spt, cpp, spp;
 double amz, amy, amx, im;

 sincosd (P.dip, &sdp, &cdp);
 sincosd (P.str, &spp, &cpp);
 sincosd (T.dip, &sdt, &cdt);
 sincosd (T.str, &spt, &cpt);

 cpt *= cdt; spt *= cdt;
 cpp *= cdp; spp *= cdp;

 amz = sdt + sdp; amx = spt + spp; amy = cpt + cpp;
 d1 = atan2d (hypot(amx, amy), amz);
 p1 = atan2d (amy, -amx);
 if (d1 > 90.0) {
 d1 = 180.0 - d1;
 p1 -= 180.0;
 }
 if (p1 < 0.0) p1 += 360.0;

 amz = sdt - sdp; amx = spt - spp; amy = cpt - cpp;
 d2 = atan2d (hypot(amx, amy), amz);
 p2 = atan2d (amy, -amx);
 if (d2 > 90.0) {
 d2 = 180.0 - d2;
 p2 -= 180.0;
 }
 if (p2 < 0.0) p2 += 360.0;

 NP1->dip = d1; NP1->str = p1;
 NP2->dip = d2; NP2->str = p2;

 im = 1;
 if (P.dip > T.dip) im = -1;
 NP1->rake = computed_rake2 (NP2->str, NP2->dip, NP1->str, NP1->dip, im);
 NP2->rake = computed_rake2 (NP1->str, NP1->dip, NP2->str, NP2->dip, im);
 }
 */

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
