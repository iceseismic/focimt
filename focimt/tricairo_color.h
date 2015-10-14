//---------------------------------------------------------------------------
#ifndef tricairo_colorH
#define tricairo_colorH
//---------------------------------------------------------------------------

namespace Taquart {
  //! Color description class.
  /*! \ingroup tricairo
   */
  typedef class TriCairo_Color {
    public:
      double R;
      double G;
      double B;
      double A;
      TriCairo_Color(void);
      TriCairo_Color(double R, double G, double B, double A = 1.0);
      TriCairo_Color(unsigned char R, unsigned char G, unsigned char B,
          unsigned char A = 255);
      //TriCairo_Color(TColor t, double A=1.0);
      void Dispatch(double &r, double &g, double &b, double &a);
    private:
    protected:
  } TCColor;
}
#endif
