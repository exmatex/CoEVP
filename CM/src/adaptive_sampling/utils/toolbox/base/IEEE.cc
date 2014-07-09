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
// File:	IEEE.cc
// Package:	toolbox
// 
// 
// 
// Description:	IEEE routines to set up handlers and get signaling NaNs
//

#include "toolbox/base/IEEE.h"
#include "toolbox/parallel/MPI.h"
#include <float.h>
#include <math.h>
#include <limits.h>

/*
 * Floating point exception handling.  
 * 
 * The following lines setup exception handling header files.
 */
#if defined(HAVE_EXCEPTION_HANDLING)
#include <stdlib.h>
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

#ifdef DEBUG_NO_INLINE
#include "toolbox/base/IEEE.I"
#endif

namespace toolbox {


/*
 * Create the function invoked when an exception is tripped. 
 */
#if defined(HAVE_EXCEPTION_HANDLING)
static void error_action(int error) 
{
   fprintf(stderr, "floating point exception -- program abort!\n");
   abort();
   MPI::abort();
}
#endif

/*
 *  Settings for the various signaling NaNs on different systems
 */

#if !defined(FLT_SNAN_IS_BROKEN)  
float  IEEE::s_signaling_nan_float  = FLT_SNAN;
#elif !defined(FLT_MAX_IS_BROKEN)
float  IEEE::s_signaling_nan_float  = FLT_MAX;
#else
float  IEEE::s_signaling_nan_float  = NAN;
#endif

#if !defined(DBL_SNAN_IS_BROKEN)
double  IEEE::s_signaling_nan_double  = DBL_SNAN;
#elif !defined(DBL_MAX_IS_BROKEN)
double  IEEE::s_signaling_nan_double  = DBL_MAX;
#else
double  IEEE::s_signaling_nan_double  = NAN;
#endif

int    IEEE::s_int_max = INT_MAX;
int    IEEE::s_int_min = INT_MIN;
float  IEEE::s_flt_max = FLT_MAX;
float  IEEE::s_flt_min = FLT_MIN;
float  IEEE::s_flt_epsilon = FLT_EPSILON;
double IEEE::s_dbl_max = DBL_MAX;
double IEEE::s_dbl_min = DBL_MIN;
double IEEE::s_dbl_epsilon = DBL_EPSILON;
   


/*
*************************************************************************
*									*
* Set up the IEEE exception handlers so that normal IEEE exceptions	*
* will cause a program abort.  How this is done varies wildly from	*
* architecture to architecture. 					*
*************************************************************************
*/

void IEEE::setupExceptionHandlers()
{

#if defined(HAVE_EXCEPTION_HANDLING)
   unsigned short fpu_flags = _FPU_DEFAULT;          
   fpu_flags &= ~_FPU_MASK_IM;  /* Execption on Invalid operation */
   fpu_flags &= ~_FPU_MASK_ZM;  /* Execption on Division by zero  */
   fpu_flags &= ~_FPU_MASK_OM;  /* Execption on Overflow */
   _FPU_SETCW(fpu_flags);
   signal(SIGFPE, error_action);
#endif
}

/*
*************************************************************************
*									*
* Initialize float and double values to the signaling nan.              *
* Initialize int to INT_MAX.                                            *
*									*
*************************************************************************
*/

void IEEE::setNaN(float &f)
{  
   f = s_signaling_nan_float;
}

void IEEE::setNaN(double &d)
{  
   d = s_signaling_nan_double;
}

/*
*************************************************************************
*									*
* Initialize float and double arrays to signaling NaNs.			*
* Initialize int array to INT_MAX.                                      *
*									*
*************************************************************************
*/

void IEEE::initializeArrayToSignalingNaN(float *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_signaling_nan_float;
   }
}

void IEEE::initializeArrayToSignalingNaN(double *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_signaling_nan_double;
   }
}

void IEEE::initializeArrayToINT_MAX(int *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_int_max;
   }
}

void IEEE::initializeArrayToINT_MIN(int *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_int_min;
   }
}

void IEEE::initializeArrayToFLT_MAX(float *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_flt_max;
   }
}

void IEEE::initializeArrayToFLT_MIN(float *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_flt_min;
   }
}

void IEEE::initializeArrayToDBL_MAX(double *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_dbl_max;
   }
}

void IEEE::initializeArrayToDBL_MIN(double *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_dbl_min;
   }
}

/*
*************************************************************************
*									*
* Return whether or not the value is a NaN.     	                *
*									*
*************************************************************************
*/

bool IEEE::isNaN(const float &f) 
{
   int i = isnan(f);
   if (i != 0) {
     return(true);
   } else {
     return(false);
   }
}

bool IEEE::isNaN(const double &d) 
{
   int i = isnan(d);
   if (i != 0) {
     return(true);
   } else {
     return(false);
   }
}


}



