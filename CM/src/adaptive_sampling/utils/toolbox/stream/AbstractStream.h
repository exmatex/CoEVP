//
// File:	AbstractStream.h
// Package:	toolbox
// 
// 
// 
// Description:	Abstract base class for data stream abstraction
//

#ifndef included_toolbox_AbstractStream
#define included_toolbox_AbstractStream

#if 0
#ifndef included_toolbox_Complex
#include "toolbox/Complex.h"
#endif
#endif

namespace toolbox {


/*!
 * @brief AbstractStream is an abstract virtual base class that
 * provides a stream abstraction for packing and unpacking data.  
 *
 * Data sizes are defined to be compliant with the XDR specs (the data buffer
 * is always four-byte aligned).  This alignment may make buffers a bit
 * larger than necessary, but allows for XDR translation.
 *
 * The sizeofXXX() functions should be used to calculate the amount of
 * buffer space needed to pack primitive data type XXX, where XXX is one
 * of bool, char, double, float, or int.
 */

class AbstractStream 
{
public:
   /*!
    * Default constructor for an abstract stream.
    */
   AbstractStream();

   /*!
    * Virtual destructor for an abstract stream.
    */
   virtual ~AbstractStream();

   /*!
    * @name Stream Space Calculation Primitives
    * Calculate the stream space needed for various data types.
    */
   //@{
   /// Calculate the stream space needed for n bools.
   static int sizeofBool(const int n = 1);
   /// Calculate the stream space needed for n chars.
   static int sizeofChar(const int n = 1);
#if 0
   /// Calculate the stream space needed for n double complex.
   static int sizeofDoubleComplex(const int n = 1);
#endif
   /// Calculate the stream space needed for n doubles.
   static int sizeofDouble(const int n = 1);
   /// Calculate the stream space needed for n floats.
   static int sizeofFloat(const int n = 1);
   /// Calculate the stream space needed for n ints.
   static int sizeofInt(const int n = 1);
   //@}

   /*!
    * @name Boolean Stream Primitives
    * Pack and unpack booleans into and out of the abstract stream.
    */
   //@{
   /// Pack a single bool into the abstract stream.
   virtual AbstractStream& operator<<(const bool& data) = 0;
   /// Remove a single bool from the abstract stream.
   virtual AbstractStream& operator>>(bool& data) = 0;
   /// Pack an array of bools into the abstract stream.
   virtual void pack(const bool* data, const int n = 1) = 0;
   /// Remove an array of bools from the abstract stream.
   virtual void unpack(bool* data, const int n = 1) = 0;
   //@}

   /*!
    * @name Character Stream Primitives
    * Pack and unpack chars into and out of the abstract stream.
    */
   //@{
   /// Pack a single char into the abstract stream.
   virtual AbstractStream& operator<<(const char& data) = 0;
   /// Remove a single char from the abstract stream.
   virtual AbstractStream& operator>>(char& data) = 0;
   /// Pack an array of chars into the abstract stream.
   virtual void pack(const char* data, const int n = 1) = 0;
   /// Remove an array of chars from the abstract stream.
   virtual void unpack(char* data, const int n = 1) = 0;
   //@}

#if 0
   /*!
    * @name Double Complex Stream Primitives
    * Pack and unpack double complex into and out of the abstract stream.
    */
   //@{
   /// Pack a single double complex into the abstract stream.
   virtual AbstractStream& operator<<(const dcomplex& data) = 0;
   /// Remove a single double complex from the abstract stream.
   virtual AbstractStream& operator>>(dcomplex& data) = 0;
   /// Pack an array of double complex into the abstract stream.
   virtual void pack(const dcomplex* data, const int n = 1) = 0;
   /// Remove an array of double complex from the abstract stream.
   virtual void unpack(dcomplex* data, const int n = 1) = 0;
   //@}
#endif

   /*!
    * @name Double Stream Primitives
    * Pack and unpack doubles into and out of the abstract stream.
    */
   //@{
   /// Pack a single double into the abstract stream.
   virtual AbstractStream& operator<<(const double& data) = 0;
   /// Remove a single double from the abstract stream.
   virtual AbstractStream& operator>>(double& data) = 0;
   /// Pack an array of doubles into the abstract stream.
   virtual void pack(const double* data, const int n = 1) = 0;
   /// Remove an array of doubles from the abstract stream.
   virtual void unpack(double* data, const int n = 1) = 0;
   //@}

   /*!
    * @name Float Stream Primitives
    * Pack and unpack floats into and out of the abstract stream.
    */
   //@{
   /// Pack a single float into the abstract stream.
   virtual AbstractStream& operator<<(const float& data) = 0;
   /// Remove a single float from the abstract stream.
   virtual AbstractStream& operator>>(float& data) = 0;
   /// Pack an array of floats into the abstract stream.
   virtual void pack(const float* data, const int n = 1) = 0;
   /// Remove an array of floats from the abstract stream.
   virtual void unpack(float* data, const int n = 1) = 0;
   //@}

   /*!
    * @name Integer Stream Primitives
    * Pack and unpack integers into and out of the abstract stream.
    */
   //@{
   /// Pack a single integer into the abstract stream.
   virtual AbstractStream& operator<<(const int& data) = 0;
   /// Remove a single integer from the abstract stream.
   virtual AbstractStream& operator>>(int& data) = 0;
   /// Pack an array of integers into the abstract stream.
   virtual void pack(const int* data, const int n = 1) = 0;
   /// Remove an array of integers from the abstract stream.
   virtual void unpack(int* data, const int n = 1) = 0;
   //@}

private:
   static int roundXDR(const int n);
};

}

#ifndef DEBUG_NO_INLINE
#include "toolbox/stream/AbstractStream.I"
#endif
#endif
