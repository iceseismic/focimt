//-----------------------------------------------------------------------------
#include <algorithm>
#include <stdio.h>
//-----------------------------------------------------------------------------
#include "tstring.h"

Taquart::String Taquart::FormatFloat(Taquart::String Format, double Value) {
  char buffer[50];
  snprintf(buffer, 49, Format.c_str(), Value);
  return Taquart::String(buffer);
}

//-----------------------------------------------------------------------------
Taquart::String Taquart::ExtractFileName(Taquart::String Filename) {
  // Remove directory if present.
  // Do this before extension removal incase directory has a period character.
  std::string filename(Filename.c_str());
  const size_t last_slash_idx = filename.find_last_of("\\/");
  if (std::string::npos != last_slash_idx) {
    filename.erase(0, last_slash_idx + 1);
  }

  // Remove extension if present.
  const size_t period_idx = filename.rfind('.');
  if (std::string::npos != period_idx) {
    filename.erase(period_idx);
  }
  return Taquart::String(filename);
}

//-----------------------------------------------------------------------------
Taquart::String::String(void) {
  // Default constructor.
}

Taquart::String Taquart::String::UpperCase(void) {
  std::string Out(Text);
  std::transform(Out.begin(), Out.end(), Out.begin(), ::toupper);
  return Out;
}

Taquart::String Taquart::String::LowerCase(void) {
  std::string Out(Text);
  std::transform(Out.begin(), Out.end(), Out.begin(), ::tolower);
  return Out;
}

Taquart::String::String(const std::string &AText) {
  Text = AText;
}

Taquart::String::String(const char *AText) {
  Text = std::string(AText);
}

std::string Taquart::String::c_txt(void) {
  return Text;
}

Taquart::String::String(const double Source) {
  Text = FormatFloat("%f", Source).c_txt();
}

const char * Taquart::String::c_str(void) {
  return Text.c_str();
}

int Taquart::String::Length(void) {
  return Text.length();
}

int Taquart::String::Pos(const char *Needle) {
  std::string a(Needle);
  return Pos(a);
}

/*
 int Taquart::String::Pos(char character) {
 std::string a(&character);
 return Pos(a);
 }
 */

double Taquart::String::ToDouble(void) {
  return atof(Text.c_str());
}

int Taquart::String::ToInt(void) {
  return atoi(Text.c_str());
}

int Taquart::String::Pos(Taquart::String &Needle) {
  return Pos(Needle.c_str());
}

int Taquart::String::Pos(const std::string &Needle) {
  std::size_t position = Text.find(Needle);
  if (position != std::string::npos)
    return position + 1;
  else
    return 0;
}

Taquart::String::String(const Taquart::String &Source) {
  Assign(Source);
}

Taquart::String& Taquart::String::operator=(const Taquart::String &Source) {
  Assign(Source);
  return *this;
}

Taquart::String Taquart::String::operator+(const Taquart::String &Source) {
  return Taquart::String(Text + Source.Text);
}

void Taquart::String::Assign(const Taquart::String &Source) {
  Text = Source.Text;
}

Taquart::String Taquart::String::SubString(int pos, int len) {
  return Taquart::String(Text.substr(pos - 1, len).c_str());
}

Taquart::String::~String(void) {
  // Default destructor.
}

Taquart::String Taquart::String::Trim(void) {
  return TrimRight().TrimLeft();
}

Taquart::String Taquart::String::TrimLeft(void) {
  size_t startpos = Text.find_first_not_of(" \n\r\t");
  if (startpos == std::string::npos)
    return Taquart::String();
  else
    return Taquart::String(Text.substr(startpos));
}

Taquart::String Taquart::String::TrimRight(void) {
  size_t endpos = Text.find_last_not_of(" \n\r\t");
  return
      (endpos == std::string::npos) ?
          Taquart::String() : Taquart::String(Text.substr(0, endpos + 1));
}

char & Taquart::String::operator[](int i) {
  return Text[i - 1];
}

bool Taquart::operator==(const Taquart::String &Object1,
    const Taquart::String &Object2) {
  return Object1.Text == Object2.Text;
}

/*
 bool operator!= (Taquart::String &Object1, Taquart::String &Object2) {
 return !(Object1 == Object2);
 }

 */
