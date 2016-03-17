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

#include "Taylor.h"
#include <stdlib.h>

#define USE_FORT 1

extern "C"
{
   void taylorf_(int*, double*, double*, const double*, const double*, const double*);
   void printtensor2_(double* in);
}

Tensor2Sym
Taylor::tensorFunction( const Tensor2Sym& tau_prime ) const
{
   Tensor2Sym val;
   Tensor2Sym tau_dev;

   // project onto deviatoric 
   tau_dev = dev(tau_prime);

#if USE_FORT
   int size = 6;
   double in[6];

   double out[6];

   // copy tensor values to flat array
   for (int i = 0; i < 6; i++) { 
      in[i] = tau_dev.a[i];
   }

   // printf("Calling Taylor fortran\n");

   taylorf_(&size, in, out, &m_m, &m_D_0, &m_g);
   for (int i = 0; i < 6; i++) { 
      val.a[i] = out[i];
   }

#else

   if ( m_m == 1. ) {
      val = m_D_0 * tau_dev / m_g;
   }
   else {

      double norm_tau_dev = norm(tau_dev);

//      cout << norm_tau_dev << endl;
      if (norm_tau_dev > 0.) {
         val = m_D_0 * tau_dev * ( pow(norm_tau_dev/m_g, 1./m_m) / norm_tau_dev );
      }
      else {
         val = 0.;
      }
   }

#endif
   return val;
}

#if 0
Tensor4LSym
Taylor::tensorFunctionDerivative( const Tensor2Sym& in ) const
{
   Tensor4LSym val;

   if ( m_m == 1. ) {
      val = m_D_0 * Tensor4LSym(1) / m_g;
   }
   else {
      double in_norm = norm(in);

      if (in_norm > 0.) {

         for (int i=1; i<=3; ++i) {
            for (int j=1; j<=i; ++j) {
               for (int k=1; k<=3; ++k) {
                  for (int l=1; l<=k; ++l) {
                     val(i,j,k,l) = in(i,j) * in(k,l);
                  }
               }
            }
         }

         val *= (1. - m_m) / (m_m * in_norm * in_norm);

         val += Tensor4LSym(1);

         val *= pow(in_norm, (1. - m_m)/m_m ) * m_D_0 * pow(m_g, -1./m_m);
      }
      else {
         val = 0.;
      }
   }

   return val;
}
#endif



void
Taylor::getScalingsForSampling( vector<double>& input_scaling,
                                vector<double>& output_scaling ) const
{
   assert(input_scaling.size() == m_pointDimension);
   assert(output_scaling.size() == m_valueDimension);

   if ( m_m == 1. ) {
      for (int i=0; i<m_pointDimension; ++i) {
         input_scaling[i] = 1.e-1;
      }
      for (int i=0; i<m_valueDimension; ++i) {
         output_scaling[i] = 1.e3;
      }
   }
   else if ( m_m == 1./2. ) {
      for (int i=0; i<m_pointDimension; ++i) {
         //         input_scaling[i] = 1.e-2;
         input_scaling[i] = 1.e-4;
      }
      for (int i=0; i<m_valueDimension; ++i) {
         output_scaling[i] = 1.e1;
      }
   }
   else {
      for (int i=0; i<m_pointDimension; ++i) {
         input_scaling[i] = m_g;
      }
      for (int i=0; i<m_valueDimension; ++i) {
         output_scaling[i] = m_D_0;
      }
   }
}


void
Taylor::advance( const double delta_t, void* state)
{
   // This doesn't actually do anything but copy three doubles
   // back and forth, but is here as an example of how a more
   // sophisticated plasticity model might work

   // Get the state
   double* local_state = (double*)state;

   m_D_0 = *local_state; local_state++;
   m_m   = *local_state; local_state++;
   m_g   = *local_state;

   // Advance the fine-scale model (trivially)

   // Update the state
   local_state = (double*)state;

   *local_state = m_D_0; local_state++;
   *local_state = m_m;   local_state++;
   *local_state = m_g;
}


void
Taylor::getState( void* state ) const
{
   double* local_state = (double*)state;

   *local_state = m_D_0; local_state++;
   *local_state = m_m;   local_state++;
   *local_state = m_g;
}


