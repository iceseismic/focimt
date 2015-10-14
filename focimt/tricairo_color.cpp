//---------------------------------------------------------------------------
#include "tricairo_color.h"

//---------------------------------------------------------------------------
using namespace Taquart;

//---------------------------------------------------------------------------
TriCairo_Color::TriCairo_Color(void) {
  R = 0.0;
  G = 0.0;
  B = 0.0;
  A = 1.0;
}

//---------------------------------------------------------------------------
TriCairo_Color::TriCairo_Color(double r, double g, double b, double a) {
  R = r;
  G = g;
  B = b;
  A = a;
}

//---------------------------------------------------------------------------
TriCairo_Color::TriCairo_Color(unsigned char r, unsigned char g,
    unsigned char b, unsigned char a) {
  R = (double) r / 255.0;
  G = (double) g / 255.0;
  B = (double) b / 255.0;
  A = (double) a / 255.0;
}

//---------------------------------------------------------------------------
/*
 TriCairo_Color::TriCairo_Color(TColor t, double a)
 {
 R = (double)(t & 0x000000ff) / 255.0;
 G = (double)( (t & 0x0000ff00) >> 8) / 255.0;
 B = (double)( (t & 0x00ff0000) >> 16) / 255.0;
 A = a;
 }
 */

//---------------------------------------------------------------------------
void TriCairo_Color::Dispatch(double &r, double &g, double &b, double &a) {
  r = R;
  g = G;
  b = B;
  a = A;
}

