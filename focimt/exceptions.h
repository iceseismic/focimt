//---------------------------------------------------------------------------
#ifndef trinityexceptionsH
#define trinityexceptionsH
//---------------------------------------------------------------------------

namespace Taquart {
  /*! \defgroup stdexception STD-compliant exception classes (libexceptions.a)
   *  \ingroup trilib
   */

//! Parent class for simple exception classes.
  /*! TriException class is a parent class for various exception classes
   *  in Trinity namespace.
   *
   *  \version 1.4.0 [2006.10.01] Removed <stdio.h> header file, TriException
   *   class was totally rebuilt (copy constructor and assignment operator
   *   added).
   *  \version 1.2.0 [2006.09.16] Removed <string> header file.
   *  \version 1.1.2 [2006.06.13] Corrected header file.
   *  \version 1.1.1 [2005.11.06] Fixed a few errors in the code, full
   *   documentation included.
   *  \version 1.1.0 [2005.11.05] Re-coded from scratch
   *  \version 1.0.0 [2004.11.10] First version released
   *
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriException {
    public:
      //! Pointer to the exception message.
      char * Message;

      //! Default constructor.
      TriException(void);

      //! Constructor.
      /*! \param AMessage Custom exception message (will be merged with
       *  a default error message for parent class.)
       */
      TriException(const char * AMessage);

      //! Copy constructor.
      /*! \param ASource Reference to the source object.
       */
      TriException(const TriException& ASource);

      //! Assignment operator.
      /*! \param ASource Reference to the source object.
       *  \return Reference to the current object.
       */
      TriException& operator=(const TriException& ASource);

      //! Default destructor.
      virtual ~TriException(void);

    private:

      unsigned int length(const char * Text);
      void clear(void);
      void fill(const char * AMessage);

    protected:
  };

//! "Wrong operation" exception class.
  /*! Trinity::TriEWrongOperation is thrown when forbidden operation occurs.
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriEWrongOperation: public TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriEWrongOperation(const char * AMessage = "Wrong operation");

    private:
    protected:
  };

//! "Out of range" exception class.
  /*! Trinity::TriEOutOfRange is thrown when the number is out of range.
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriEOutOfRange: public Taquart::TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriEOutOfRange(const char * AMessage = "Out of range");

    private:
    protected:
  };

//! "Operation aborted" exception class.
  /*! Trinity::TriEOperationAborted is thrown when current operation has
   *  been aborted.
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriEOperationAborted: public TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriEOperationAborted(const char * AMessage = "Operation aborted");

    private:
    protected:
  };

//! "Input/output" exception class.
  /*! Trinity::TriEIOError is thrown when a general input/output error occurs.
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriEIOError: public TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriEIOError(const char * AMessage = "Input / Output error");

    private:
    protected:
  };

//! "Conversion failed" exception class.
  /*! Trinity::TriEConversionError is thrown when conversion from a one type to
   *  another one is not available (e.g. while converting string to
   *  floating-point number is impossible).
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriEConversionError: public TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriEConversionError(const char * AMessage = "Conversion error");

    private:
    protected:
  };

//! "Empty string" exception class.
  /*! Trinity::TriEEmptyString is thrown when some string is empty (e.g. its
   *  length is equal to zero). To signal that the pointer is NULL, use
   *  Trinity::TriENullPointer instead.
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriEEmptyString: public TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriEEmptyString(const char * AMessage = "Empty string");

    private:
    protected:
  };

//! "NULL pointer" exception class.
  /*! Trinity::TriENullPointer is thrown when parameter specified in e.g.
   *  function call is NULL. To signal that the string is empty, use
   *  Trinity::TriENullPointer instead.
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriENullPointer: public TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriENullPointer(const char * AMessage = "NULL Pointer");

    private:
    protected:
  };
} // namespace Trinity

#endif
