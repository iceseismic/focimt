//-----------------------------------------------------------------------------
#ifndef TSTRING_H_
#define TSTRING_H_
//-----------------------------------------------------------------------------

#include <string>

namespace Taquart {

  class String {
    public:
      String(void);
      String(const std::string &AText);
      String(const char *AText);
      String(const Taquart::String &Source);
      String(const double Source);
      Taquart::String &operator= (const Taquart::String &Source);
      Taquart::String operator+ (const Taquart::String &Source);
      char &operator[] (int i);
      void Assign(const Taquart::String &Source);
      virtual ~String(void);
      const char * c_str(void);
      std::string c_txt(void);
      int Pos(const std::string &Needle);
      int Pos(const char *Needle);
      int Pos(Taquart::String &Needle);
      //int Pos(char Needle);
      int Length(void);
      Taquart::String SubString(int Pos, int Length);
      Taquart::String Trim(void);
      Taquart::String TrimLeft(void);
      Taquart::String TrimRight(void);
      Taquart::String UpperCase(void);
      Taquart::String LowerCase(void);
      double ToDouble(void);
      int ToInt(void);
      //friend bool operator!= (Taquart::String &Object1, Taquart::String &Object2);
      friend bool operator== (const Taquart::String &Object1, const Taquart::String &Object2);
    private:
      std::string Text;
    protected:
  };

  Taquart::String ExtractFileName(Taquart::String Filename);

  Taquart::String FormatFloat(Taquart::String Format, double Value);
  //Taquart::String FormatFloat(Taquart::String Format, int Value);

  bool operator== (const Taquart::String &Object1, const Taquart::String &Object2);
} /* namespace Taquart */

#endif /* TSTRING_H_ */
