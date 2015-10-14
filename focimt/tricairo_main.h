//---------------------------------------------------------------------------
#ifndef tricairo_mainH
#define tricairo_mainH
//---------------------------------------------------------------------------
#include <cairo/cairo.h>
//#include <cairo/cairo-win32.h>
#include <cairo/cairo-svg.h>
#include <cairo/cairo-pdf.h>
#include <cairo/cairo-ps.h>
#include "tricairo_color.h"
#include <math.h>
#include "tstring.h"

//! \defgroup tricairo Interface to CAIRO library.

namespace Taquart {

  //! Type of line joint
  /*! \ingroup tricairo
   */
  enum TriCairo_LineJoin {
    ljMiter, ljBevel, ljRound
  };

  //! Type of line cap.
  /*! \ingroup tricairo
   */
  enum TriCairo_LineCap {
    lcButt, lcRound, lcSquare
  };

  //! Type of output canvas.
  /*! \ingroup tricairo
   */
  enum TriCairo_CanvasType {
    ctBitmap, ctSurface, ctSVG, ctPDF, ctPS
  };

  //! Font style.
  /*! \ingroup tricairo
   */
  enum TriCairo_FontStyle {
    fsNormal, fsItalic, fsOblique, fsNormalBold, fsItalicBold, fsObliqueBold
  };

  //! Horizontal alignment of text.
  /*! \ingroup tricairo
   */
  enum TriCairo_HorizontalAlignment {
    haLeft, haCenter, haRight
  };

  //! Vertical alignment of text.
  /*! \ingroup tricairo
   */
  enum TriCairo_VerticalAlignment {
    vaTop, vaMiddle, vaBottom
  };

  //! Base class wrapping the interface between BDS2006 and Cairo library.
  /*! \ingroup tricairo
   */
  class TriCairo {
    public:
      TriCairo(unsigned int width, unsigned int height,
          Taquart::TriCairo_CanvasType canvastype,
          Taquart::String Filename = "");
      virtual ~TriCairo(void);
      virtual void Save(Taquart::String filename);
      //Graphics::TBitmap * GetBitmap(void);
      //Graphics::TBitmap * CreateBitmap(void);
      //void DrawCanvas(TCanvas * Canvas, int Left = 0, int Top = 0);

      // Drawing functions.
      void Color(Taquart::TriCairo_Color Color);
      void ColorD(double r, double g, double b, double a = 1.0);
      void ColorB(unsigned int r, unsigned int g, unsigned int b,
          unsigned int a = 255);
      void Clear(double r, double g, double b, double a = 1.0);
      void Clear(Taquart::TriCairo_Color Color);
      void MoveTo(double x, double y);
      void LineTo(double x, double y);
      void LineToRel(double x, double y);
      void LineWidth(double w);
      void LineCap(TriCairo_LineCap lc);
      void LineJoin(TriCairo_LineJoin lj);
      void Font(Taquart::String Name, double Size, TriCairo_FontStyle fs);
      void Text(double x, double y, Taquart::String Text,
          TriCairo_HorizontalAlignment ha = haLeft,
          TriCairo_VerticalAlignment va = vaTop);
      void Stroke(void);
      void StrokePreserve(void);
      void Fill(void);
      void FillPreserve(void);
      void Arc(double x, double y, double r, double start = 0.0,
          double end = 2 * M_PI);
      void Rectangle(double x, double y, double w, double h);
      void Circle(double x, double y, double r);
      void ClosePath(void);

      Taquart::String Filename;
    private:
      //Graphics::TBitmap * Bitmap;

    protected:
      cairo_surface_t * CreateSurface(TriCairo_CanvasType CanvasType);

      // Additional functions.
      TriCairo_CanvasType CanvasType;
      const unsigned int Width;
      const unsigned int Height;
      cairo_surface_t * surface;
      cairo_t * cr;
  };

}

#endif
