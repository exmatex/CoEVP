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

#include <stdlib.h>
#include "ElastoViscoPlasticity.h"
#include "xtensor.h"


ElastoViscoPlasticity::ElastoViscoPlasticity( const Tensor2Gen&      L,
                                              const double           bulk_modulus,
                                              const double           shear_modulus,
                                              const EOS*             eos_model,
                                              const Plasticity*      plasticity_model,
                                              const bool             use_adaptive_sampling )
   : m_D_old(sym(L)),
     m_W_old(skew(L)),
     m_R(1),
     m_J(1.),
     m_Vbar_prime(0),
     m_Vbar_prime_dot(0),
     m_Dbar_prime(0),
     m_Wbar(0),
     m_K(bulk_modulus),
     m_G(shear_modulus),
     m_eos_model(eos_model),
     m_plasticity_model(plasticity_model)
{
   assert(eos_model != NULL);
   assert(plasticity_model != NULL);

   if ( use_adaptive_sampling ) {

      int pointDimension           = m_plasticity_model->pointDimension();
      int valueDimension           = m_plasticity_model->valueDimension();
      int maxKrigingModelSize      = 6;
      int maxNumberSearchModels    = 1;
      double theta                 = 5.e3;
      double meanErrorFactor       = 1.;
      double tolerance             = 3.e-3;
      double maxQueryPointModelDistance = 0.2;

      std::vector<double> pointScaling( pointDimension );
      std::vector<double> valueScaling( valueDimension );
      plasticity_model->getScalingsForSampling( pointScaling, valueScaling );

      enableAdaptiveSampling( pointDimension, valueDimension, pointScaling, valueScaling,
                              maxKrigingModelSize, maxNumberSearchModels, theta, meanErrorFactor,
                              tolerance, maxQueryPointModelDistance );
   }
}


ElastoViscoPlasticity::~ElastoViscoPlasticity()
{
}


Tensor2Sym
ElastoViscoPlasticity::stress( const double compression,
                               const double e,
                               const double q ) const
{
   if ( m_eos_model != NULL ) {
      double p = m_eos_model->evaluate( compression, e );

      return -(p + q)*Tensor2Sym(1) + stressDeviator();
   }
   else {
      cout << "ElastoViscoPlasticity::stress(): Can't compute full stress tensor. "
           << "No eos model was specified in the constructor." << endl;
      exit(1);
   }
}


Tensor2Sym
ElastoViscoPlasticity::stressDeviator() const
{
   Tensor2Sym sigma_bar_prime( tauBarPrime(a(m_J), m_Vbar_prime) );
   sigma_bar_prime /= m_J;

   Tensor2Sym sigma_prime;
   convertToCoarse( sigma_bar_prime, m_R, sigma_prime );

   return sigma_prime;
}


void
ElastoViscoPlasticity::setNewVelocityGradient( const Tensor2Gen& L_new )
{
   m_D_new = sym(L_new);
   m_W_new = skew(L_new);
}


double
ElastoViscoPlasticity::soundSpeedSquared( const double reference_density,
                                          const double relativeVolume,
                                          const double energy ) const
{
   double mu = 1./relativeVolume - 1.;

   double K;
#if 0   
   // Estimate the bulk modulus
   K = m_eos_model->evaluate_dpdrho( mu, reference_density, energy );
#else
   K = m_K;
#endif

   // p-wave modulus
   return (K + 4.*m_G/3.) / reference_density;
}


void
ElastoViscoPlasticity::advance( const double delta_t )
{
   // Get the deviatoric strain rate at the new time
   Tensor2Sym Dprime_new( dev(m_D_new) );

   // Evaluate W^R per equation (29)
   Tensor2Gen WR;
   evaluateWR( m_Vbar_prime, m_Vbar_prime_dot, m_Dbar_prime, m_W_old,
               m_Wbar, m_R, a(m_J), WR );

   // Update the rotation per equation (32)
   Tensor2Gen R;
   updateR( m_R, WR, delta_t, R );

   // Update the stretch determinant per equation (31)
   double J;
   updateJ( m_J, m_D_old, delta_t, J);

   // Use the current value of the stretch deviator and its
   // approximate derivative to estimate the new time value
   Tensor2Sym Vbar_prime;
   Vbar_prime = delta_t * m_Vbar_prime_dot;
   Vbar_prime += m_Vbar_prime;

   // Update the stretch deviator per equation (30)
   Tensor2Sym Dbar_prime;
   Tensor2Gen Wbar;
   updateVbar_prime( m_Vbar_prime, Dprime_new, R, a(J), delta_t, Vbar_prime, Dbar_prime, Wbar );

   // Update the internal state in preparation for the next call

   m_Vbar_prime_dot = Vbar_prime - m_Vbar_prime;
   m_Vbar_prime_dot /= delta_t;

   m_Vbar_prime = Vbar_prime;
   m_J = J;
   m_R = R;
   m_Dbar_prime = Dbar_prime;
   m_Wbar = Wbar;
   m_D_old = m_D_new;
   m_W_old = m_W_new;
}


Tensor2Sym
ElastoViscoPlasticity::tauBarPrime( const double      a,
                                    const Tensor2Sym& Vbar_prime ) const
{
   return (2.*m_G/a)*Vbar_prime;
}


void
ElastoViscoPlasticity::convertToFine( const Tensor2Sym& in,
                                      const Tensor2Gen& R,
                                      Tensor2Sym&       out ) const

{
   //   out = R^T * in * R 

   out = sym( R.transpose() * in * R );
}


void
ElastoViscoPlasticity::convertToCoarse( const Tensor2Gen& in,
                                        const Tensor2Gen& R,
                                        Tensor2Gen&       out ) const

{
   //   out = R * in * R^T 

   out = R * in * R.transpose();
}


void
ElastoViscoPlasticity::convertToCoarse( const Tensor2Sym& in,
                                        const Tensor2Gen& R,
                                        Tensor2Sym&       out ) const

{
   //   out = R * in * R^T 

   out = sym( R * in * R.transpose() );
}


void
ElastoViscoPlasticity::evaluateStretchRHS( const Tensor2Sym& D_prime,
                                           const Tensor2Sym& Dbar_prime,
                                           const Tensor2Gen& R,
                                           Tensor2Sym&       rhs ) const
{
   /*
     Evaluates

       rhs = R^T * D_prime * R - Dbar_prime

     per equation (27) of the specification.

   */

   convertToFine( D_prime, R, rhs );

   rhs -= Dbar_prime;
}


void
ElastoViscoPlasticity::updateVbar_prime( const Tensor2Sym& Vbar_prime_old,
                                         const Tensor2Sym& Dprime_new,
                                         const Tensor2Gen& R_new,
                                         const double      a_new,
                                         const double      delta_t,
                                         Tensor2Sym&       Vbar_prime_new,
                                         Tensor2Sym&       Dbar_prime_new,
                                         Tensor2Gen&       Wbar )
{
   /*
     Uses Newton iteration to update Vbar_prime_new per equations (33) through (37)
     of the specification.  It is assumed that an initial guess is being passed in.
   */

   double tol = 1.e-4;        // Read from input if we find some reason to later
   int max_iter = 15;
   double saved_residual[max_iter];

   Tensor2Sym residual;

   Tensor4LSym Dbar_prime_deriv;
   evaluateFineScaleModel( tauBarPrime(a_new, Vbar_prime_new), Dbar_prime_new, Dbar_prime_deriv );
   Wbar = 0.;

   computeResidual( Vbar_prime_new, Vbar_prime_old, Dprime_new, Dbar_prime_new,
                    R_new, a_new * delta_t, residual );

   double relative_residual;
   double residual_scale = 0.5 * (norm(Dprime_new) + norm(Dbar_prime_new));

   if (residual_scale > 0.) {
      relative_residual = norm(residual) / residual_scale;
   }
   else {
      relative_residual = 0.;
   }

   bool converged = relative_residual < tol;
   m_num_iters = 0;

   double theta = 0.1;
   double t = 1.e-4;
   double old_residual_norm = norm(residual);

   while ( !converged && m_num_iters < max_iter) {

      // Compute the Jacobian to get the update
      Tensor4LSym jacobian;
      computeJacobian( Dbar_prime_deriv, a_new, delta_t, jacobian );

      Tensor2Sym delta;
      solveJacobianSystem( jacobian, residual, delta );

      double eta = 0.;
      Tensor2Sym proposed_solution = Vbar_prime_new + delta;

      evaluateFineScaleModel( tauBarPrime(a_new, proposed_solution), Dbar_prime_new, Dbar_prime_deriv );
      Wbar = 0.;

      computeResidual( proposed_solution, Vbar_prime_old, Dprime_new, Dbar_prime_new,
                       R_new, a_new * delta_t, residual );

#if 0
      while ( norm(residual) > (1. - t*(1. - eta)) * old_residual_norm ) {

         delta *= theta;
         proposed_solution = Vbar_prime_new + delta;
         eta = 1. - theta * (1. - eta);

         evaluateFineScaleModel( tauBarPrime(a_new, proposed_solution), Dbar_prime_new, Dbar_prime_deriv );
         Wbar = 0.;

         computeResidual(proposed_solution, Vbar_prime_old, Dprime_new, Dbar_prime_new,
                         R_new, a_new * delta_t, residual );
      }
#endif

      // Update solution
      Vbar_prime_new = proposed_solution;


      double new_residual_norm ;

      relative_residual = (new_residual_norm = norm(residual)) * 2. /
                          (norm(Dprime_new) + norm(Dbar_prime_new));

      saved_residual[m_num_iters] = relative_residual ;

      // Test convergence
      converged = relative_residual <= tol;

      m_num_iters++;
      old_residual_norm = new_residual_norm;
   }

#if 0
   if ( !converged && m_num_iters >= max_iter ) {
      std::cout << "Newton solve did not converge, num_iters = " << m_num_iters << ", relative residual = " << relative_residual << std::endl;
      std::cout << std::endl;

      for (int i=0; i<m_num_iters; ++i) {
         cout << i << " " << saved_residual[i] << endl;
      }
      exit(1);

   }
   else {
      //      cout << "   Number of Newton iterations = " << m_num_iters << ", relative residual = " << relative_residual << endl;
   }
#endif
}


void
ElastoViscoPlasticity::computeResidual( const Tensor2Sym& Vbar_prime_new,
                                        const Tensor2Sym& Vbar_prime_old,
                                        const Tensor2Sym& Dprime_new,
                                        const Tensor2Sym& Dbar_prime_new,
                                        const Tensor2Gen& R_new,
                                        const double      a_delta_t,
                                        Tensor2Sym&       residual ) const
{
   /*
     Computes the Newton residual

       F = (1/a_delta_t)(Vbar_prime_new - Vbar_prime_old) - R^T * Dprime_new * R + Dbar_prime_new

     per equation (34) of the specification.

   */

   residual = (Vbar_prime_new - Vbar_prime_old) / a_delta_t;

   Tensor2Sym RT_Dprime_new_R;
   convertToFine( Dprime_new, R_new, RT_Dprime_new_R );

   residual += Dbar_prime_new - RT_Dprime_new_R;
}


void
ElastoViscoPlasticity::computeJacobian( const Tensor4LSym& Dbar_deriv,
                                        const double       a,
                                        const double       delta_t,
                                        Tensor4LSym&       jacobian ) const
{
   /*
     Computes the Jacobian

       jacobian = (1/a*delta_t)I + (2G/a)*Dbar_deriv

     per equation (37) of the specification.

   */

   jacobian = Tensor4LSym(1) / (a * delta_t) + (2. * m_G/ a) * Dbar_deriv;
}


void
ElastoViscoPlasticity::solveJacobianSystem( const Tensor4LSym& jacobian,
                                            const Tensor2Sym&  residual,
                                            Tensor2Sym&        delta ) const
{
   /*
     Solves the Jacobian system

       jacobian * delta = - residual

     per equation (36) of the specification.

   */

   delta = solve( jacobian, -residual );
}


void
ElastoViscoPlasticity::updateJ( const double      J_old,
                                const Tensor2Sym& D,
                                const double      delta_t,
                                double&           J_new ) const
{
#if 0
   /*
     Updates J per equation (31) of the specification.
   */

   J_new = exp( trace(D) * delta_t ) * J_old;
#else
   J_new = m_volume_change * J_old;
#endif
}


void
ElastoViscoPlasticity::evaluateWR( const Tensor2Gen& Vbar_prime,
                                   const Tensor2Gen& Vbar_prime_dot,
                                   const Tensor2Gen& Dbar_prime,
                                   const Tensor2Gen& W,
                                   const Tensor2Gen& Wbar,
                                   const Tensor2Gen& R,
                                   const double      a,
                                   Tensor2Gen&       WR ) const
{
   /*
     Evaluates

       WR = W - R * 

          { Wbar 
 
              - (1/a) R * [  Vbar_prime * (Dbar_prime + (1/2a)Vbar_prime_dot) 

                         +  (Dbar_prime + (1/2a)Vbar_prime_dot) * Vbar_prime ]

          } * R^T

     per equation (29) of the specification.

   */

   Tensor2Gen temp( Dbar_prime + Vbar_prime_dot/(2.*a) );

   Tensor2Gen term1(Vbar_prime * temp);
   Tensor2Gen term2(temp * Vbar_prime);

   temp = (term2 - term1)/a - Wbar;

   convertToCoarse( temp, R, WR );
   WR += W;
}


void
ElastoViscoPlasticity::updateR( const Tensor2Gen& R_old,
                                const Tensor2Gen& WR,
                                const double      delta_t,
                                Tensor2Gen&       R_new ) const
{
   /*
     Updates R per equation (32) of the specification.
   */

   R_new = expW( WR*delta_t ) * R_old;
}


void
ElastoViscoPlasticity::evaluateFineScaleModel( const Tensor2Sym& tau_bar_prime,
                                               Tensor2Sym&       Dbar_prime,
                                               Tensor4LSym&      Dbar_prime_deriv ) const
{
   if ( adaptiveSamplingEnabled() ) {

      std::vector<double> point(m_plasticity_model->pointDimension());
      std::vector<double> value(m_plasticity_model->valueAndDerivativeDimension());

      m_plasticity_model->packInputVector( tau_bar_prime, point );

      sample( *m_plasticity_model, point, value );

      m_plasticity_model->unpackOutputVector( value, Dbar_prime, Dbar_prime_deriv );

   }
   else {

      m_plasticity_model->evaluateNative( tau_bar_prime, Dbar_prime, Dbar_prime_deriv );

   }
}

