//---------------------------------------------------------------------------
#ifndef tricairo_mecaH
#define tricairo_mecaH
//---------------------------------------------------------------------------
#include "fortranmath.h"
#include "georoutines.h"
#include "tricairo_main.h"

namespace Taquart {
  typedef int GMT_LONG;

  //! Hemisphere projection.
  /*! \ingroup tricairo
   */
  enum TriCairo_Hemisphere {
    heLower = 0, heUpper = 1
  };

  //! Type of station marker.
  /*! \ingroup tricairo
   */
  enum TriCairo_MarkerType {
    mtCircle = 0, mtSquare = 1, mtPlusMinus = 2, mtBWCircle = 3
  };

  //! Type of network projection.
  /*! \ingroup tricairo
   */
  enum TriCairo_Projection {
    prWulff = 0, prSchmidt = 1
  };

  //! Structure stores information about plunge and trend of axis.
  /*! \ingroup tricairo
   */

  //! Structure stores information about strike, dip and rake of a fault.
  /*! \ingroup tricairo
   */

  /*
   typedef struct DLL_EXP TriCairo_NodalPlane
   {
   double str;
   double dip;
   double rake;
   } nodal_plane;
   */
  //! Structure stores information moment tensor.
  /*! The corresponding matrix elements correspond to the moment tensor
   *  components according to CMT convention.
   *  \ingroup tricairo
   */
  typedef struct TriCairo_MomentTensor {
      double f[6]; /* mrr mtt mff mrt mrf mtf in 10**expo dynes-cm */
      TriCairo_MomentTensor(void) {
        for (int i = 0; i < 6; i++)
          f[i] = 0;
      }
  } M_TENSOR;

  //! Class for producing the graphical representation of moment tensor component.
  /*! This class is capable to produce a graphical representation of the seismic
   *   moment tensor (so called beach balls).
   *  \ingroup tricairo
   */
  class TriCairo_Meca: public TriCairo {
    public:
      // Constructor.
      TriCairo_Meca(unsigned int width, unsigned int height,
          TriCairo_CanvasType type, Taquart::String filename = "");
      // Destructor.
      virtual ~TriCairo_Meca(void);
      //void Draw(double M11, double M12, double M13, double M22, double M23, double M33,
      //  double s1, double d1, double r1, double s2, double d2, double r2);

      // Public variables.
      bool DrawCross;
      bool DrawDC;

      TriCairo_Hemisphere Hemisphere;
      unsigned int Margin; // Margin size.
      unsigned int BWidth;
      unsigned int BHeight;
      unsigned int BRadius;
      unsigned int BXo;
      unsigned int BYo;
      double BOutlineWidth;
      TCColor BOutlineColor;
      TCColor BPlusColor; // For compressional part.
      TCColor BMinusColor; // For dilatational part.
      TCColor BTensorOutline;
      TCColor BDCColor;
      double BDCWidth;

      bool DrawAxis;
      double AxisFontSize;
      Taquart::String AxisFontFace;

      bool DrawStations;
      bool DrawStationName;
      double StationMarkerSize;
      double StationFontSize;
      TCColor StationPlusColor;
      TCColor StationMinusColor;
      TriCairo_MarkerType StationMarkerType;
      Taquart::String StationFontFace;
      TCColor StationTextColor;

      TriCairo_Projection Projection;

      // Main drawing routines.
      void Station(double Azimuth, double Takeoff, double Disp,
          Taquart::String Label, double &mx, double &my);
      void Tensor(AXIS T, AXIS N, AXIS P);
      void Axis(AXIS A, Taquart::String Text);
      void DoubleCouple(double Strike, double Dip);
      void DrawStationMarker(double x, double y, double disp,
          Taquart::String Label);
      void CenterCross(void);
      void GMT_momten2axe(M_TENSOR mt, AXIS *T, AXIS *N, AXIS *P);
      void Station(double GA[], double Disp, Taquart::String Label, double &mx,
          double &my);

    protected:

      // Upper or lower hemisphere projection routines.
      void Project(double &X, double &Y);
      void Project(double * X, double * Y, unsigned int npoints);

    private:
      // Additional routines.
      double squared(double v);
      void axe2dc(AXIS T, AXIS P, nodal_plane *NP1, nodal_plane *NP2);
      double proj_radius2(double str1, double dip1, double str);
      void ps_circle(double x0, double y0, double radius_size, TCColor fc);
      void Polygon(double xp1[], double yp1[], int npoints, TCColor oc,
          bool fill, TCColor fc = TCColor(), double OutlineWidth = 1.0);

      GMT_LONG GMT_jacobi(double *a, GMT_LONG *n, GMT_LONG *m, double *d,
          double *v, double *b, double *z, GMT_LONG *nrots);

      void * GMT_memory(void *prev_addr, GMT_LONG nelem, size_t size);
      void GMT_free(void *addr);

  };
// class TriCairo_Meca
}
#endif
