// DO-NOT-DELETE revisionify.begin() 
/*

                            Copyright (c) 2014.
               Lawrence Livermore National Security, LLC.
         Produced at the Lawrence Livermore National Laboratory
                             LLNL-CODE-656392.
                           All rights reserved.

This file is part of CoEVP, Version 1.0. Please also read this link -- http://www.opensource.org/licenses/index.php

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the disclaimer below.

   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the disclaimer (as noted below)
     in the documentation and/or other materials provided with the
     distribution.

   * Neither the name of the LLNS/LLNL nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
   Department of Energy (DOE). This work was produced at Lawrence Livermore
   National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
   Security, LLC nor any of their employees, makes any warranty, express
   or implied, or assumes any liability or responsibility for the accuracy,
   completeness, or usefulness of any information, apparatus, product, or
   process disclosed, or represents that its use would not infringe
   privately-owned rights.

3. Also, reference herein to any specific commercial products, process, or
   services by trade name, trademark, manufacturer or otherwise does not
   necessarily constitute or imply its endorsement, recommendation, or
   favoring by the United States Government or Lawrence Livermore National
   Security, LLC. The views and opinions of authors expressed herein do not
   necessarily state or reflect those of the United States Government or
   Lawrence Livermore National Security, LLC, and shall not be used for
   advertising or product endorsement purposes.

*/
// DO-NOT-DELETE revisionify.end() 
//
// File:	MathUtilitiesSpecial.cc
// Package:	toolbox
// 
// 
// 
// Description:	MathUtilities routines to set up handlers and get signaling NaNs
//

#include "MathUtilities.h"

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <limits.h>

#if 0
#include "toolbox/Complex.h"
#endif

/*
 * Floating point exception handling.  
 * The following lines setup exception handling header files on 
 * systems other than solaris.
 */
#if defined(HAVE_EXCEPTION_HANDLING) 
#include <stdio.h>
#include <fpu_control.h>
#include <signal.h>
#endif

/*
 * The following lines setup exception handling headers on the Sun.  If we
 * use Sun's native compiler, just pull in the <sunmath.h> include file.
 * If we are under solaris but use a different compiler (e.g. KCC, g++)
 * we have to explicitly define the functions that <sunmath.h> defines,
 * since we don't have access to this file.
 */
#ifdef __SUNPRO_CC
#include <sunmath.h>
#endif

namespace toolbox {

/*
 *  Settings for the various signaling NaNs on different systems
 */
#if !defined(FLT_SNAN_IS_BROKEN)  
#define MPTCOUPLER_FLT_SNAN   FLT_SNAN
#elif !defined(FLT_MAX_IS_BROKEN)
#define MPTCOUPLER_FLT_SNAN   FLT_MAX
#else
#define MPTCOUPLER_FLT_SNAN   NAN
#endif

#if !defined(DBL_SNAN_IS_BROKEN)
#define MPTCOUPLER_DBL_SNAN   DBL_SNAN
#elif !defined(DBL_MAX_IS_BROKEN)
#define MPTCOUPLER_DBL_SNAN   DBL_MAX
#else
#define MPTCOUPLER_DBL_SNAN   NAN
#endif

template<> bool   MathUtilities<bool>::s_zero           = false;
template<> bool   MathUtilities<bool>::s_one            = true;
template<> bool   MathUtilities<bool>::s_signaling_nan  = false;
template<> bool   MathUtilities<bool>::s_max            = true;
template<> bool   MathUtilities<bool>::s_min            = false;
template<> bool   MathUtilities<bool>::s_epsilon        = true;
template<> bool   MathUtilities<bool>::s_undefined      = false;

template<> char   MathUtilities<char>::s_zero           = 0;
template<> char   MathUtilities<char>::s_one            = 1;
template<> char   MathUtilities<char>::s_signaling_nan  = CHAR_MAX;
template<> char   MathUtilities<char>::s_max            = CHAR_MAX;
template<> char   MathUtilities<char>::s_min            = CHAR_MIN;
template<> char   MathUtilities<char>::s_epsilon        = 1;
template<> char   MathUtilities<char>::s_undefined      = (char) 0xff;

template<> int    MathUtilities<int>::s_zero           = 0;
template<> int    MathUtilities<int>::s_one            = 1;
template<> int    MathUtilities<int>::s_signaling_nan   = INT_MAX;
template<> int    MathUtilities<int>::s_max             = INT_MAX;
template<> int    MathUtilities<int>::s_min             = INT_MIN;
template<> int    MathUtilities<int>::s_epsilon         = 1;
template<> int    MathUtilities<int>::s_undefined       = INT_MAX;

template<> float  MathUtilities<float>::s_zero          = 0.0;
template<> float  MathUtilities<float>::s_one           = 1.0;
template<> float  MathUtilities<float>::s_signaling_nan = MPTCOUPLER_FLT_SNAN;
template<> float  MathUtilities<float>::s_max           = FLT_MAX;
template<> float  MathUtilities<float>::s_min           = FLT_MIN;
template<> float  MathUtilities<float>::s_epsilon       = FLT_EPSILON;
template<> float  MathUtilities<float>::s_undefined     = MPTCOUPLER_FLT_SNAN;

template<> double MathUtilities<double>::s_zero          = 0.0;
template<> double MathUtilities<double>::s_one           = 1.0;
template<> double MathUtilities<double>::s_signaling_nan = MPTCOUPLER_DBL_SNAN;
template<> double MathUtilities<double>::s_max           = DBL_MAX;
template<> double MathUtilities<double>::s_min           = DBL_MIN;
template<> double MathUtilities<double>::s_epsilon       = DBL_EPSILON;
template<> double MathUtilities<double>::s_undefined     = MPTCOUPLER_DBL_SNAN;

#if 0
template<> dcomplex   MathUtilities<dcomplex>::s_zero             = dcomplex(0.0,0.0);
template<> dcomplex   MathUtilities<dcomplex>::s_one              = dcomplex(1.0,0.0);
template<> dcomplex   MathUtilities<dcomplex>::s_signaling_nan  = dcomplex(MPTCOUPLER_DBL_SNAN,MPTCOUPLER_DBL_SNAN);
template<> dcomplex   MathUtilities<dcomplex>::s_max            = dcomplex(DBL_MAX,DBL_MAX);
template<> dcomplex   MathUtilities<dcomplex>::s_min            = dcomplex(DBL_MIN,DBL_MIN);
template<> dcomplex   MathUtilities<dcomplex>::s_epsilon        = dcomplex(DBL_MIN,0.0);
template<> dcomplex   MathUtilities<dcomplex>::s_undefined      = dcomplex(MPTCOUPLER_DBL_SNAN,MPTCOUPLER_DBL_SNAN);
#endif


template<>
bool MathUtilities<float>::isNaN(const float& value)
{
   int i = isnan(value);
   if (i != 0) {
     return(true);
   } else {
     return(false);
   }
}

template<> 
bool MathUtilities<double>::isNaN(const double& value)
{
   int i = isnan(value);
   if (i != 0) {
     return(true);
   } else {
     return(false);
   }
}

template<>
bool MathUtilities<float>::equalEps(const float& a, const float& b)
{
   float absmax = MathUtilities<float>::Max(
                     MathUtilities<float>::Abs(a),
                     MathUtilities<float>::Abs(b) );
   float numerator = MathUtilities<float>::Abs(a-b);
   float denomenator =
      MathUtilities<float>::Max(absmax,
           MathUtilities<float>::s_epsilon);

   return( numerator/denomenator < sqrt(MathUtilities<float>::s_epsilon) );
}

template<>
bool MathUtilities<double>::equalEps(const double& a, const double& b)
{
   double absmax = MathUtilities<double>::Max(
                      MathUtilities<double>::Abs(a),
                      MathUtilities<double>::Abs(b) );
   double numerator = MathUtilities<double>::Abs(a-b);
   double denomenator =
      MathUtilities<double>::Max(absmax,
           MathUtilities<double>::s_epsilon);

   return( numerator/denomenator < sqrt(MathUtilities<double>::s_epsilon) );
}

#if 0
template<>
bool MathUtilities<dcomplex>::equalEps(const dcomplex& a, const dcomplex& b)
{
   double a_re = real(a);
   double a_im = imag(a);
   double b_re = real(b);
   double b_im = imag(b);

   return( equalEps(a_re,b_re) && equalEps(a_im,b_im) );
}

template<> 
dcomplex MathUtilities<dcomplex>::Min(dcomplex a, dcomplex b)
{
   return(norm(a) < norm(b) ? a : b);
}

template<> 
dcomplex MathUtilities<dcomplex>::Max(dcomplex a, dcomplex b)
{
   return(norm(a) > norm(b) ? a : b);
}
#endif

template<> 
int MathUtilities<int>::Abs(int a)
{
   return(a > 0 ? a : -a); 
}

template<> 
float MathUtilities<float>::Abs(float a)
{
   return(a > 0.0 ? a : -a); 
}

template<> 
double MathUtilities<double>::Abs(double a)
{
   return(a > 0.0 ? a : -a); 
}

template<> 
bool MathUtilities<bool>::Rand(const bool& low, const bool& width)
{
   return mrand48() > 0 ? true : false;
}

template<> 
char MathUtilities<char>::Rand(const char& low, const char& width)
{
   return (char)((double)width * drand48()) + low;
}

template<> 
int MathUtilities<int>::Rand(const int& low, const int& width)
{
   return (int)((double)width * drand48()) + low;
}

template<> 
float MathUtilities<float>::Rand(const float& low, const float& width)
{
   return width * drand48() + low;
}

template<> 
double MathUtilities<double>::Rand(const double& low, const double& width)
{
   return width * drand48() + low;
}

#if 0
template<> 
dcomplex MathUtilities<dcomplex>::Rand(const dcomplex& low, const dcomplex& width)
{
   double real_part = real(width) * drand48() + real(low);
   double imag_part = imag(width) * drand48() + imag(low);
   return dcomplex(real_part, imag_part);
}
#endif

}




