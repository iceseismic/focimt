//---------------------------------------------------------------------------
#include "tricairo_main.h"

using namespace Taquart;

//---------------------------------------------------------------------------
TriCairo::TriCairo(unsigned int width, unsigned int height,
    TriCairo_CanvasType canvastype, String filename) :
    Width(width), Height(height) {
  Filename = filename;

  // Default constructor.
  CreateSurface(canvastype);
  cr = cairo_create(surface);

  // Clean up the surface.
  Clear(1.0, 1.0, 1.0);
  ColorD(0.0, 0.0, 0.0);

  //Bitmap = NULL;
}

//---------------------------------------------------------------------------
cairo_surface_t * TriCairo::CreateSurface(TriCairo_CanvasType canvastype) {
  CanvasType = canvastype;
  switch (CanvasType) {
    case ctBitmap:
      //surface = cairo_win32_surface_create_with_dib(CAIRO_FORMAT_ARGB32, Width,
      //    Height);
      break;
    case ctSurface:
      surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, Width, Height);
      break;
    case ctSVG:
      surface = cairo_svg_surface_create(Filename.c_str(), Width, Height);
      break;
    case ctPDF:
      surface = cairo_pdf_surface_create(Filename.c_str(), Width, Height);
      break;
    case ctPS:
      surface = cairo_ps_surface_create(Filename.c_str(), Width, Height);
      break;

  };
  return surface;
}

//---------------------------------------------------------------------------
TriCairo::~TriCairo(void) {
  /*
   if(Bitmap)
   {
   delete Bitmap;
   Bitmap = NULL;
   }
   */
  // Default destructor.
  if (cr) {
    cairo_destroy(cr);
    cr = NULL;
  }
  if (surface) {
    cairo_surface_destroy(surface);
    surface = NULL;
  }
}

//---------------------------------------------------------------------------
void TriCairo::Text(double x, double y, String Text,
    TriCairo_HorizontalAlignment ha, TriCairo_VerticalAlignment va) {
  double xo = 0.0;
  double yo = 0.0;

  cairo_text_extents_t extents;
  cairo_text_extents(cr, Text.c_str(), &extents);

  switch (ha) {
    case haLeft:
      xo = 0;
      break;
    case haCenter:
      xo = -extents.width / 2.0;
      break;
    case haRight:
      xo = -extents.width;
      break;
  }

  switch (va) {
    case vaBottom:
      yo = 0;
      break;
    case vaMiddle:
      yo = extents.height / 2.0;
      break;
    case vaTop:
      yo = extents.height;
      break;
  }

  cairo_move_to(cr, x + xo, y + yo);
  cairo_show_text(cr, Text.c_str());
}

//---------------------------------------------------------------------------
void TriCairo::Font(String Name, double Size, TriCairo_FontStyle Style) {
  switch (Style) {
    case Taquart::fsNormal:
      cairo_select_font_face(cr, Name.c_str(), CAIRO_FONT_SLANT_NORMAL,
          CAIRO_FONT_WEIGHT_NORMAL);
      break;
    case Taquart::fsItalic:
      cairo_select_font_face(cr, Name.c_str(), CAIRO_FONT_SLANT_ITALIC,
          CAIRO_FONT_WEIGHT_NORMAL);
      break;
    case Taquart::fsOblique:
      cairo_select_font_face(cr, Name.c_str(), CAIRO_FONT_SLANT_OBLIQUE,
          CAIRO_FONT_WEIGHT_NORMAL);
      break;
    case Taquart::fsNormalBold:
      cairo_select_font_face(cr, Name.c_str(), CAIRO_FONT_SLANT_NORMAL,
          CAIRO_FONT_WEIGHT_BOLD);
      break;
    case Taquart::fsItalicBold:
      cairo_select_font_face(cr, Name.c_str(), CAIRO_FONT_SLANT_ITALIC,
          CAIRO_FONT_WEIGHT_BOLD);
      break;
    case Taquart::fsObliqueBold:
      cairo_select_font_face(cr, Name.c_str(), CAIRO_FONT_SLANT_OBLIQUE,
          CAIRO_FONT_WEIGHT_BOLD);
      break;
  }

  cairo_set_font_size(cr, Size);
}

//---------------------------------------------------------------------------
void TriCairo::ColorD(double r, double g, double b, double a) {
  cairo_set_source_rgba(cr, r, g, b, a);
}

//---------------------------------------------------------------------------
void TriCairo::ColorB(unsigned int r, unsigned int g, unsigned int b,
    unsigned int a) {
  cairo_set_source_rgba(cr, double(r) / 255.0f, double(g) / 255.0f,
      double(b) / 255.0f, double(a) / 255.0f);
}

//---------------------------------------------------------------------------
void TriCairo::Color(TriCairo_Color Color) {
  double r = 0.0, g = 0.0, b = 0.0, a = 0.0;
  Color.Dispatch(r, g, b, a);
  cairo_set_source_rgba(cr, r, g, b, a);
}

//---------------------------------------------------------------------------
void TriCairo::Clear(double r, double g, double b, double a) {
  cairo_save(cr);
  if (a == 1.0f) {
    // Fully opaque.
    cairo_set_source_rgb(cr, r, g, b);
  }
  else {
    // Transparent.
    cairo_set_source_rgba(cr, r, g, b, a);
    cairo_set_operator(cr, CAIRO_OPERATOR_SOURCE);
  }
  cairo_paint(cr);
  cairo_restore(cr);
}

//---------------------------------------------------------------------------
void TriCairo::Clear(TriCairo_Color Color) {
  double r = 0.0, g = 0.0, b = 0.0, a = 0.0;
  Color.Dispatch(r, g, b, a);
  Clear(r, g, b, a);
}

//---------------------------------------------------------------------------
void TriCairo::MoveTo(double x, double y) {
  cairo_move_to(cr, x, y);
}

//---------------------------------------------------------------------------
void TriCairo::LineTo(double x, double y) {
  cairo_line_to(cr, x, y);
}

//---------------------------------------------------------------------------
void TriCairo::LineToRel(double x, double y) {
  cairo_rel_line_to(cr, x, y);
}

//---------------------------------------------------------------------------
void TriCairo::Arc(double x, double y, double r, double start, double end) {
  cairo_arc(cr, x, y, r, start, end);
}

void TriCairo::ClosePath(void) {
  cairo_close_path(cr);
}

//---------------------------------------------------------------------------
void TriCairo::Circle(double x, double y, double r) {
  Arc(x, y, r);
}

void TriCairo::Rectangle(double x, double y, double w, double h) {
  cairo_rectangle(cr, x, y, w, h);
}

//---------------------------------------------------------------------------
void TriCairo::LineWidth(double w) {
  cairo_set_line_width(cr, w);
}

//---------------------------------------------------------------------------
void TriCairo::Stroke(void) {
  cairo_stroke(cr);
}

//---------------------------------------------------------------------------
void TriCairo::StrokePreserve(void) {
  cairo_stroke_preserve(cr);
}

//---------------------------------------------------------------------------
void TriCairo::Fill(void) {
  cairo_fill(cr);
}

//---------------------------------------------------------------------------
void TriCairo::FillPreserve(void) {
  cairo_fill_preserve(cr);
}

//---------------------------------------------------------------------------
void TriCairo::LineCap(TriCairo_LineCap lc) {
  switch (lc) {
    case lcButt:
      cairo_set_line_cap(cr, CAIRO_LINE_CAP_BUTT); /* default */
      break;
    case lcRound:
      cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
      break;
    case lcSquare:
      cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
      break;
  }
}

//---------------------------------------------------------------------------
void TriCairo::LineJoin(TriCairo_LineJoin lj) {
  switch (lj) {
    case ljMiter:
      cairo_set_line_join(cr, CAIRO_LINE_JOIN_MITER); /* default */
      break;
    case ljBevel:
      cairo_set_line_join(cr, CAIRO_LINE_JOIN_BEVEL);
      break;
    case ljRound:
      cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
      break;
  }
}

//---------------------------------------------------------------------------
void TriCairo::Save(String filename) {
  if (CanvasType == ctBitmap) {
    /*
     Graphics::TBitmap * Bitmap = new Graphics::TBitmap;
     Bitmap -> Width = Width;
     Bitmap -> Height = Height;
     BitBlt(Bitmap -> Canvas -> Handle, 0, 0, Width, Height,
     cairo_win32_surface_get_dc(surface), 0, 0, SRCCOPY);

     // Save result to output file.
     Bitmap -> SaveToFile(filename);
     delete Bitmap;
     */
  }
  else if (CanvasType == ctSurface) {
    cairo_surface_write_to_png(surface, filename.c_str());
  }
}

//---------------------------------------------------------------------------
/*
 Graphics::TBitmap * TriCairo::GetBitmap(void)
 {
 if(CanvasType == ctBitmap)
 {
 if(Bitmap == NULL)
 {
 Bitmap = new Graphics::TBitmap;
 Bitmap -> Width = Width;
 Bitmap -> Height = Height;
 }
 BitBlt(Bitmap -> Canvas -> Handle, 0, 0, Width, Height,
 cairo_win32_surface_get_dc(surface), 0, 0, SRCCOPY);
 return Bitmap;
 }
 else return NULL;
 }
 */

//---------------------------------------------------------------------------
/*
 void TriCairo::DrawCanvas(TCanvas * Canvas, int Left, int Top)
 {
 if(CanvasType == ctBitmap)
 {
 BitBlt(Canvas -> Handle, Left, Top, Width, Height,
 cairo_win32_surface_get_dc(surface), 0, 0, SRCCOPY);

 }
 }
 */

//---------------------------------------------------------------------------
/*
 Graphics::TBitmap * TriCairo::CreateBitmap(void)
 {
 if(CanvasType == ctBitmap)
 {
 Graphics::TBitmap *const Bitmap = new Graphics::TBitmap;
 Bitmap -> Width = Width;
 Bitmap -> Height = Height;
 BitBlt(Bitmap -> Canvas -> Handle, 0, 0, Width, Height,
 cairo_win32_surface_get_dc(surface), 0, 0, SRCCOPY);
 return Bitmap;
 }
 else return NULL;
 }
 */
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
