//-----------------------------------------------------------------------------
#ifndef TRILIB_GEOROUTINES_H
#define TRILIB_GEOROUTINES_H
//-----------------------------------------------------------------------------

#define EPSIL 0.0001
#define DEG2RAD (M_PI / 180.0)
#define sind(x) sin ((x) * DEG2RAD)
#define cosd(x) cos ((x) * DEG2RAD)
#define tand(x) tan ((x) * DEG2RAD)

namespace Taquart {

  typedef struct nodal_plane {
      double str;
      double dip;
      double rake;
  } TriCairo_NodalPlane;

  typedef struct TriCairo_Axis {
      double str;
      double dip;
      double val;
      TriCairo_Axis(void) {
        str = 0.0;
        dip = 0.0;
        val = 0.0;
      }
  } AXIS;

  struct MOMENT {
      double mant;
      int exponent;
  };

  struct MECHANISM {
      struct nodal_plane NP1;
      struct nodal_plane NP2;
      struct MOMENT moment;
      double magms;
  };

  typedef struct MECHANISM st_me;

  void axe2dc(AXIS T, AXIS P, nodal_plane *NP1,
   nodal_plane *NP2);
  void StrikeDipRake2MT(double S, double D, double R, double &M11, double &M22,
      double &M33, double &M12, double &M13, double &M23);
  double zero_360(double str);
  void sincosd(double i, double *s, double *c);
  void sincos(double i, double *s, double *c);
  double computed_dip1(nodal_plane NP1);
  double computed_rake1(struct nodal_plane NP1);
  double computed_strike1(struct nodal_plane NP1);
  void define_second_plane(struct nodal_plane NP1, struct nodal_plane *NP2);
  double null_axis_dip(double str1, double dip1, double str2, double dip2);
  double null_axis_strike(double str1, double dip1, double str2, double dip2);
  double computed_rake2(double str1, double dip1, double str2, double dip2,
      double fault);
}

//-----------------------------------------------------------------------------
#endif
