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

#include <assert.h>
#include "Plasticity.h"


void
Plasticity::evaluateNative( const Tensor2Sym& in,
                            Tensor2Sym&       out_value,
                            Tensor4LSym&      out_derivative ) const
{
   out_value = tensorFunction(in);

   out_derivative = tensorFunctionDerivative(in);
}


void
Plasticity::evaluate( const vector<double>& point,
                      vector<double>&       value ) const
{
   Tensor2Sym in;
   unpackInputVector(point, in);

   Tensor2Sym out_value;
   Tensor4LSym out_derivative;
   evaluateNative(in, out_value, out_derivative);

   value.resize(valueAndDerivativeDimension());
   packOutputVector(out_value, out_derivative, value);
}



Tensor4LSym
Plasticity::tensorFunctionDerivative( const Tensor2Sym& in ) const
{
   Tensor2Sym tensor_function_in(tensorFunction(in)) ;

   double eps = 1.e-3;
   double delta = eps * norm(in);
   if (delta == 0.) delta = eps;  // Use an absolute delta since there's no way
                                  // to determine the scale of the argument    

   Tensor4LSym out;

   for (int k=1; k<=3; ++k) {
      for (int l=1; l<=k; ++l) {
         Tensor2Sym increment(0);
         increment(k,l) = delta; 

         // This factor accounts for the fact that we are differentiating a
         // symmetric tensor with respect to a symmetric argument
         double symmetry_factor = k==l? 1.: 0.5;

         Tensor2Sym difference = (tensorFunction(in + increment) - tensor_function_in) / delta;

         for (int i=1; i<=3; ++i) {
            for (int j=1; j<=i; ++j) {
               out(i,j,k,l) = symmetry_factor * difference(i,j);
            }
         }
      }
   }
   
   //   cout << setiosflags(ios::fixed | ios::showpoint) << setw(12) << out << endl;

   return out;
}



void
Plasticity::packInputVector( const Tensor2Sym& point,
                             vector<double>&   in ) const
{
   int index = 0;

   for (int i=1; i<=3; ++i) {
      for (int j=1; j<=i; ++j) {
         in[index++] = point(i,j);
      }
   }

   assert(index == pointDimension());
}


void
Plasticity::unpackInputVector( const vector<double>& in,
                               Tensor2Sym&           point ) const
{
   int index = 0;

   for (int i=1; i<=3; ++i) {
      for (int j=1; j<=i; ++j) {
         point(i,j) = in[index++];
      }
   }

   assert(index == pointDimension());
}

#if 0
void
Plasticity::packOutputVector( Tensor2Sym&     value,
                              Tensor4LSym&    derivative,
                              vector<double>& out ) const
{
   int index = 0;

   for (int i=1; i<=3; ++i) {
      for (int j=1; j<=i; ++j) {
          out[index++] = value(i,j);
      }
   }

   double sqrt2 = sqrt(2.);

   for (int k=1; k<=3; ++k) {
      for (int l=1; l<=k; ++l) {
         for (int i=1; i<=3; ++i) {
            for (int j=1; j<=i; ++j) {
               double fac = 1.;
               if ( j<i ) fac *= sqrt2;
               if ( l<k ) fac *= sqrt2;
               out[index++] = fac * derivative(i,j,k,l);
            }
         }
      }
   }

   assert(index == valueAndDerivativeDimension());
}
#else
void
Plasticity::packOutputVector( Tensor2Sym&     value,
                              Tensor4LSym&    derivative,
                              vector<double>& out ) const
{
   int index = 0;

   for (int i=1; i<=3; ++i) {
      for (int j=1; j<=i; ++j) {
          out[index++] = value(i,j);
      }
   }

   for (int k=1; k<=3; ++k) {
      for (int l=1; l<=k; ++l) {

         double symmetry_factor = (k==l)? 1.: 2.;

         for (int i=1; i<=3; ++i) {
            for (int j=1; j<=i; ++j) {
               out[index++] = symmetry_factor * derivative(i,j,k,l);
            }
         }
      }
   }

   assert(index == valueAndDerivativeDimension());
}
#endif

#if 0
void
Plasticity::unpackOutputVector( const vector<double>& out,
                                Tensor2Sym&           value,
                                Tensor4LSym&          derivative ) const
{
   int index = 0;

   for (int i=1; i<=3; ++i) {
      for (int j=1; j<=i; ++j) {
         value(i,j) = out[index++];
      }
   }

   double sqrt2inv = 1./sqrt(2.);

   for (int k=1; k<=3; ++k) {
      for (int l=1; l<=k; ++l) {
         for (int i=1; i<=3; ++i) {
            for (int j=1; j<=i; ++j) {
               double fac = 1.;
               if ( j<i ) fac *= sqrt2inv;
               if ( l<k ) fac *= sqrt2inv;
               derivative(i,j,k,l) = fac * out[index++];
            }
         }
      }
   }

   assert(index == valueAndDerivativeDimension());
}
#else
void
Plasticity::unpackOutputVector( const vector<double>& out,
                                Tensor2Sym&           value,
                                Tensor4LSym&          derivative ) const
{
   int index = 0;

   for (int i=1; i<=3; ++i) {
      for (int j=1; j<=i; ++j) {
         value(i,j) = out[index++];
      }
   }

   for (int k=1; k<=3; ++k) {
      for (int l=1; l<=k; ++l) {

         double symmetry_factor = (k==l)? 1: 0.5;
         for (int i=1; i<=3; ++i) {
            for (int j=1; j<=i; ++j) {
               derivative(i,j,k,l) = symmetry_factor * out[index++];
            }
         }
      }
   }

   assert(index == valueAndDerivativeDimension());
}
#endif
