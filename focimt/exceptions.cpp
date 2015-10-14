//---------------------------------------------------------------------------
// EXCEPTIONS
//  Exception classes
//---------------------------------------------------------------------------

#include "exceptions.h"

using namespace Taquart;

//---------------------------------------------------------------------------
TriException::TriException(void) {
  Message = 0;
}

//---------------------------------------------------------------------------
TriException::TriException(const char * AMessage) {
  Message = 0;
  fill(AMessage);
}

//---------------------------------------------------------------------------
TriException::TriException(const TriException& ASource) {
  Message = 0;
  fill(ASource.Message);
}

//---------------------------------------------------------------------------
TriException& TriException::operator=(const TriException& ASource) {
  if (this != &ASource) {
    fill(ASource.Message);
  }
  return *this;
}

//---------------------------------------------------------------------------
TriException::~TriException(void) {
  clear();
}

//---------------------------------------------------------------------------
unsigned int TriException::length(const char * Text) {
  unsigned int i = 0;
  while (*(Text + i) != 0)
    i++;
  return i;
}

//---------------------------------------------------------------------------
void TriException::clear(void) {
  if (Message) {
    delete Message;
    Message = 0;
  }
}

//---------------------------------------------------------------------------
void TriException::fill(const char * AMessage) {
  clear();
  int Len = length(AMessage);
  Message = new char[Len + 1];
  for (int i = 0; i < Len; i++) {
    Message[i] = AMessage[i];
  }
  Message[Len] = 0;
}

//---------------------------------------------------------------------------
TriEWrongOperation::TriEWrongOperation(const char * AMessage) :
    TriException(AMessage) {
}

//---------------------------------------------------------------------------
TriEOutOfRange::TriEOutOfRange(const char * AMessage) :
    TriException(AMessage) {
}

//---------------------------------------------------------------------------
TriEOperationAborted::TriEOperationAborted(const char * AMessage) :
    TriException(AMessage) {
}

//---------------------------------------------------------------------------
TriEIOError::TriEIOError(const char * AMessage) :
    TriException(AMessage) {
}

//---------------------------------------------------------------------------
TriEConversionError::TriEConversionError(const char * AMessage) :
    TriException(AMessage) {
}

//---------------------------------------------------------------------------
TriEEmptyString::TriEEmptyString(const char * AMessage) :
    TriException(AMessage) {
}

//---------------------------------------------------------------------------
TriENullPointer::TriENullPointer(const char * AMessage) :
    TriException(AMessage) {
}

