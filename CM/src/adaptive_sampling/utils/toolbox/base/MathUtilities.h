//
// File:	MathUtilities.h
// Package:	toolbox
// 
// 
// 
// Description:	Routines to get signaling NaNs and POSIX contants
//

#ifndef included_toolbox_MathUtilities
#define included_toolbox_MathUtilities

#if 0
#ifndef included_toolbox_Complex
#include "Complex.h"
#endif
#endif

namespace toolbox {

/*!
 * @brief MathUtilities is a utility for managing basic math functions
 * and the initialization of data to signaling NaNs, and for managing
 * Posix constants like INT_MAX, FLT_MAX, DBL_MAX, etc.  Signaling
 * NaNs force a trap if they are used in a numerical operation, so
 * they are a useful way to track uninitialized floating point data.
 * Setting integer values to INT_MAX is a useful way to track
 * uninitialized integer values.
 *
 * The implementation of this class depends heavily on the particular
 * computer architecture and how it implements floating point arithmetic
 * and hardware traps.  
 */

template<class TYPE> class MathUtilities
{
  public:

   /*!
    * @brief Get the value 0 for TYPE.
    *
    */
   static TYPE getZero();

   /*!
    * @brief Get the value 1 for TYPE.
    *
    */
   static TYPE getOne();


   /*!
    * @brief Get the IEEE signaling NaN for TYPE on architectures that
    * support it.  
    *
    * Using this value in a numerical expression will
    * cause a program abort.
    *
    * Valid for float/double only.  For non float/double will return 
    * false.
    */
   static TYPE getSignalingNaN();

   /*!
    * @brief Indicates whether the supplied value is NaN.
    *
    * Valid for float/double only.  For non float/double will return 
    * false.
    *
    * @param value Value to test
    */
   static bool isNaN(const TYPE& value);

   /*!
    * @brief Return true if values have a relative difference
    * smaller than  sqrt(mach_eps).
    *
    * Valid for float/double/dcomplex only.  For dcomplex will
    * return true if both real and imaginary parts numbers
    * have a relative difference smaller than sqrt(mach_eps).
    * For non float/double/dcomplex will just return ( a == b ).
    *
    * @param a
    * @param b
    */
   static bool equalEps(const TYPE& a, const TYPE& b ); 

   /*!
    * @brief Get max for the templated type.
    */
   static TYPE getMax();

   /*!
    * @brief Get min for the templated type.
    */
   static TYPE getMin();

   /**
    * @brief Get epsilon for the templated type.
    */
   static TYPE getEpsilon();

   /*!
    * @brief Get value to set for undefined data.  
    *
    * Typicaly used for initialziation during debugging.
    */
   static TYPE getUndefined();

   /*!
    * @brief Compute the minimum of a and b.
    *
    * @param a
    * @param b
    */
   static TYPE Min(TYPE a, TYPE b);

   /*!
    * @brief Compute the maximum of a and b.
    *
    * @param a
    * @param b
    */
   static TYPE Max(TYPE a, TYPE b);

   /*!
    * @brief Return absolute value of a.
    *
    * Valid for float/double/int only.
    * For non float/double/int will return input value.
    *
    * @param a
    */
   static TYPE Abs(TYPE a);

   /*!
    * @brief Generate a random value from low to low+width.
    *
    * @param low   Starting value for range
    * @param width Width of the range.
    */
   static TYPE Rand(const TYPE& low, const TYPE& width);

private:
   static TYPE  s_zero;
   static TYPE  s_one;
   static TYPE  s_signaling_nan;
   static TYPE  s_max;
   static TYPE  s_min;
   static TYPE  s_epsilon;
   static TYPE  s_undefined;

};

/*
 * Template specializations.
 */

template<> bool MathUtilities<float>::isNaN(const float& value);
template<> bool MathUtilities<double>::isNaN(const double& value);

template<> bool MathUtilities<float>::equalEps(const float& a, const float& b);
template<> bool MathUtilities<double>::equalEps(const double& a, const double& b);

template<> bool     MathUtilities<bool>::Rand(const bool& low, const bool& width);
template<> char     MathUtilities<char>::Rand(const char& low, const char& width);
template<> int      MathUtilities<int>::Rand(const int& low, const int& width);
template<> float    MathUtilities<float>::Rand(const float& low, const float& width);
template<> double   MathUtilities<double>::Rand(const double& low, const double& width);

template<> int      MathUtilities<int>::Abs(int a);
template<> float    MathUtilities<float>::Abs(float a);
template<> double   MathUtilities<double>::Abs(double a);

#if 0
template<> bool MathUtilities<double>::isNaN(const double& value);
template<> bool MathUtilities<dcomplex>::equalEps(const dcomplex& a, const dcomplex& b);
template<> dcomplex MathUtilities<dcomplex>::Rand(const dcomplex& low, const dcomplex& width);

template<> dcomplex MathUtilities<dcomplex>::Min(dcomplex a, dcomplex b);
template<> dcomplex MathUtilities<dcomplex>::Max(dcomplex a, dcomplex b);
#endif

}

#ifndef DEBUG_NO_INLINE
#include "MathUtilities.I"
#endif

#endif
