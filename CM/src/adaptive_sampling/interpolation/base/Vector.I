/* DO-NOT-DELETE revisionify.begin() */
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
/* DO-NOT-DELETE revisionify.end() */
//
// File:        Vector.I
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing MTL-aware vector.
//
// $Id$
//
// $Log$
//

#include <algorithm>

#ifndef included_mtl_mtl_h
#define included_mtl_mtl_h
#include <mtl/mtl.h>
#endif

#ifndef included_mtl_blas_interface_h
#define included_mtl_blas_interface_h
#include <mtl/blas_interface.h>
#endif 

#ifdef DEBUG_NO_INLINE
#define inline
#endif // DEBUG_NO_INLINE

namespace krigalg {

    //
    // construction
    //

    inline
    Vector::Vector()
    {

      return;

    }

    inline
    Vector::Vector(int n)
      : mtl::dense1D<double>(n)
    {

      return;

    }
  
    inline
    Vector::Vector(int n, double init)
      : mtl::dense1D<double>(n, init)
    {

      return;

    }

    inline
    Vector::Vector(int n, const double * value)
      : mtl::dense1D<double>(n) 
    {

      std::copy(&(value[0]),
		&(value[n]),
		&((*this)[0]));

      return;

    }

    //
    // addition
    //

    inline Vector & 
    operator+=(Vector       & vector1,
	       const Vector & vector2)
    {

      mtl::add(vector1,
	       vector2,
	       vector1);

      return vector1;

    }
  

    inline Vector 
    operator+(const Vector & vector1,
	      const Vector & vector2)
    {

      Vector result(vector1.size());
      mtl::copy(vector1, result);
      return result += vector2;
	     
    }

    //
    // subtraction
    //

    inline Vector &
    operator-=(Vector       & vector1,
	       const Vector & vector2)
    {

      mtl::add(vector1,
	       mtl::scaled(vector2, -1.0),
	       vector1);
    
      return vector1;

    }

    inline Vector
    operator-(const Vector & vector1,
	      const Vector & vector2)
    {

      Vector result(vector1.size());
      mtl::copy(vector1, result);
      return result -= vector2;

    }
     
     //
     // scalar product
     //

     inline Vector &
     operator*=(Vector & vector1,
		double alpha)
     {
	
	mtl::scale(vector1, alpha);
	return vector1;

     }

     inline Vector
     operator*(const Vector & vector1,
	       double alpha)
     {
	
	Vector result(vector1.size());
	mtl::copy(vector1, result);
	mtl::scale(result, alpha);
	return result;

     }

    //
    // compute dot product of two vectors
    //

    inline double
    dot(const Vector & vector1,
	const Vector & vector2)
    {

      //
      // firewalls
      //

      assert(vector1.size() == vector2.size());
      
      //
      // setup parameters
      //

      int n = static_cast<int>(vector1.size());
      double * x = const_cast<double *>(&(vector1[0]));
      int incx = 1;
      double * y = const_cast<double *>(&(vector2[0]));
      int incy = 1;

      //
      //
      //

      return ddot_(&n,
		   x,
		   &incx,
		   y,
		   &incy);
      
    }

}


