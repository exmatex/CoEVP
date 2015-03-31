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
     m_plasticity_model(plasticity_model),
     m_Delta_max(1.e-2),
     m_Delta(0.99*m_Delta_max)
{
   assert(eos_model != NULL);
   assert(plasticity_model != NULL);

   if ( use_adaptive_sampling ) {

      int pointDimension           = m_plasticity_model->pointDimension();
      int valueDimension           = m_plasticity_model->valueDimension();
      int maxKrigingModelSize      = 6;
      int maxNumberSearchModels    = 4;
      double theta                 = 5.e2;
      double meanErrorFactor       = 1.;
      double tolerance             = 5.e-4;
      double maxQueryPointModelDistance = 5.e-2;

      std::vector<double> pointScaling( pointDimension );
      std::vector<double> valueScaling( valueDimension );
      plasticity_model->getScalingsForSampling( pointScaling, valueScaling );

      enableAdaptiveSampling( pointDimension, valueDimension, pointScaling, valueScaling,
                              maxKrigingModelSize, maxNumberSearchModels, theta, meanErrorFactor,
                              tolerance, maxQueryPointModelDistance );

      m_tol = tolerance;
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

   Tensor2Sym Vbar_prime;
   // Use the current value of the stretch deviator to 
   // estimate the new time value
   Vbar_prime = m_Vbar_prime;

   // Update the stretch deviator per equation (30)
   Tensor2Sym Dbar_prime;
   Tensor2Gen Wbar_new;
   updateVbar_prime( m_Vbar_prime, Dprime_new, R, a(J), delta_t, Vbar_prime, Dbar_prime, Wbar_new );

   // Update the internal state in preparation for the next call

   m_Vbar_prime_dot = Vbar_prime - m_Vbar_prime;
   m_Vbar_prime_dot /= delta_t;

   m_Vbar_prime = Vbar_prime;
   m_J = J;
   m_R = R;
   m_Dbar_prime = Dbar_prime;
   m_Wbar = Wbar_new;
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
ElastoViscoPlasticity::getDoglegStep( const Tensor2Sym&  r,
                                      const Tensor4LSym& J,
                                      const double       Delta,
                                      Tensor2Sym&        p ) const
{
   // Compute the Cauchy point

   Tensor2Sym Jtranspose_r;

   for (int k=1; k<=3; ++k) {
      for (int l=1; l<=k; ++l) {

         double symmetry_factor = k==l? 1.: 2.;

         double sum = 0.;
         for (int i=1; i<=3; ++i) {
            for (int j=1; j<=i; ++j) {
               sum += symmetry_factor * J(i,j,k,l) * r(i,j);
            }
         }
         Jtranspose_r(k,l) = sum;
      }
   }

   double norm_Jtranspose_r = sqrt(norm2(Jtranspose_r));

   Tensor2Sym J_Jtranspose_r = J * Jtranspose_r;

   double tau = pow(norm_Jtranspose_r,3) / (Delta * norm2(J_Jtranspose_r));
   if ( tau > 1. ) tau = 1.;

   Tensor2Sym Cauchy_p = - tau * (Delta / norm_Jtranspose_r) * Jtranspose_r;

   // If the Cauchy point is within the trust region, compute a Newton step
   // and interpolate the two.  Otherwise, return the Cauchy point.

   if ( sqrt(norm2(Cauchy_p)) >= Delta ) {  // Steepest descent step reaches trust region boundary
      p = Cauchy_p;
   }
   else {

      Tensor2Sym Newton_p;
      solveJacobianSystem( J, r, Newton_p );

      Tensor2Sym p_diff = Newton_p - Cauchy_p;

      double a = norm2(p_diff);

      if (a == 0.) {  // Newton step and Cauchy point are the same
         p = Newton_p;
      }
      else {
         // Solve a quadratic equation to find the largest tau in the unit interval
         // such that norm(Cauchy_p + tau * (Newton_p - Cauchy_p)) <= Delta

         double b = 2. * dot(Cauchy_p,p_diff);
         double c = norm2(Cauchy_p);

         tau = 0.5 * (-b + sqrt(b*b - 4.*a*(c - Delta*Delta))) / a;
         if (tau > 1.) tau = 1.;

         p = Cauchy_p + tau * p_diff;
      }
   }
}


void
ElastoViscoPlasticity::updateVbar_prime( const Tensor2Sym& Vbar_prime_old,
                                         const Tensor2Sym& Dprime_new,
                                         const Tensor2Gen& R_new,
                                         const double      a_new,
                                         const double      delta_t,
                                         Tensor2Sym&       Vbar_prime_new,
                                         Tensor2Sym&       Dbar_prime_new,
                                         Tensor2Gen&       Wbar_new )
{
   /*
     This function solves the Vbar_prime equation using the Dogleg Newton algorithm
     described in Jorge Nocedal and Stephen J. Wright, "Numerical Optimzation", Second Edition,
     Springer Series in Operation Research and Financial Engineering, Springer 2006.
   */ 

   double eta = 0.24;  // Tolerance for the actual merit function decrease relative
                       // to the prediction given by the local model.  A larger value
                       // corresponds to a tighter tolerance, but the maximum value
                       // allowed is 0.25;
   int max_iter = 40;  // Maximum number of iterations allowed, including those due only
                       // to trust region modifications
   double tol = 1.e-4; // Requested merit function tolerance.  The merit function is
                       // one-half the square of the residual norm

   // end of options (need to move this elsewhere)

   // For now, we assume the fine scale model spin is always zero, but we nevertheless
   // include an "update" stub here in case we want to add a non-zero spin in the future.
   Wbar_new = Tensor2Gen(0);

   assert(eta < 0.25);

   // Quantities saved for convergence failure post mortem
   double* diagnostic_merit_value = new double[max_iter+1];
   double* diagnostic_Delta = new double[max_iter+1];
   double* diagnostic_rho = new double[max_iter+1];
   double* diagnostic_p_norm = new double[max_iter+1];
   double* diagnostic_error = new double[max_iter+1];
   int*    diagnostic_model = new int[max_iter+1];

   double tol2 = 2. * tol;

   Tensor4LSym Dbar_prime_deriv;
   evaluateFineScaleModel( tauBarPrime(a_new, Vbar_prime_new), Dbar_prime_new, Dbar_prime_deriv );

   double saved_estimate = m_error_estimate;
   int current_model = m_hint;
   diagnostic_model[0] = m_hint;
   diagnostic_rho[0] = 0.;
   diagnostic_error[0] = m_error_estimate;

   Tensor2Sym r;
   computeResidual( Vbar_prime_new, Vbar_prime_old, Dprime_new, Dbar_prime_new,
                    R_new, a_new * delta_t, r );
   double rnorm2 = norm2(r);

   diagnostic_merit_value[0] = 0.5 * rnorm2;
   diagnostic_Delta[0] = m_Delta;

   Tensor4LSym jacobian;

   m_num_iters = 0;
   bool converged = rnorm2 < tol2;
   bool new_local_model = true;
   
   double Delta_saved = m_Delta;

   while ( !converged && m_num_iters < max_iter ) {

      if ( new_local_model ) {
         computeJacobian( Dbar_prime_deriv, a_new, delta_t, jacobian );
      }

      // Calculate dogleg step
      Tensor2Sym p;
      getDoglegStep(r, jacobian, m_Delta, p);

      Tensor2Sym proposed_solution = Vbar_prime_new + p;

      diagnostic_p_norm[m_num_iters] = norm(p);

      evaluateFineScaleModel( tauBarPrime(a_new, proposed_solution), Dbar_prime_new, Dbar_prime_deriv );

      diagnostic_model[m_num_iters+1] = m_hint;
      diagnostic_error[m_num_iters+1] = m_error_estimate;

#if 0
      if ( m_hint == -1 && current_model != -1 ) {

         // A fine-scale model evaluation was just performed and the result was either added to an existing
         // kriging model in the database or a new model was created.  If the modified model happens to be the
         // current_model (we need some way to determine this for certain) we need to re-evaluate the current
         // residual.  Otherwise, the backtracking resulting from reductions of the trust region size may
         // not succeed.

         Tensor2Sym Dbar_prime_new_tmp;
         Tensor4LSym Dbar_prime_deriv_tmp;
         evaluateSpecificModel( current_model, tauBarPrime(a_new, Vbar_prime_new), Dbar_prime_new_tmp, Dbar_prime_deriv_tmp );

         Tensor2Sym r_new;
         computeResidual( Vbar_prime_new, Vbar_prime_old, Dprime_new, Dbar_prime_new_tmp,
                          R_new, a_new * delta_t, r_new );

         if ( norm(r-r_new) != 0. ) {  // The current residual has changed, presumably as the result of adding
                                       // the new fine-scale model evaluation to it            

            //            cout << "residual changed by " << norm(r-r_new) << ", restarting at iteration " << m_num_iters;

            if ( m_error_estimate >= 5.e-4 ) {

               //               cout << ", old estimate = " << saved_estimate << ", new estimate = " << m_error_estimate;

               // For some reason, re-evaluation of the current_model has produced an error estimate that now
               // exceeds the tolerance, whereas it presumably had originally. In this case, we start completely over.

               evaluateFineScaleModel( tauBarPrime(a_new, Vbar_prime_new), Dbar_prime_new_tmp, Dbar_prime_deriv_tmp );

               computeResidual( Vbar_prime_new, Vbar_prime_old, Dprime_new, Dbar_prime_new_tmp,
                                R_new, a_new * delta_t, r_new );

               current_model = m_hint;
            }

            //            cout << endl;

            m_Delta = Delta_saved;
            m_num_iters = 0;
            Dbar_prime_deriv = Dbar_prime_deriv_tmp;
            r = r_new;
            rnorm2 = norm2(r);
            new_local_model = true;
            m_hint = current_model;
            saved_estimate = m_error_estimate;
            continue;
         }
   }
#endif

      Tensor2Sym proposed_r;
      computeResidual( proposed_solution, Vbar_prime_old, Dprime_new, Dbar_prime_new,
                       R_new, a_new * delta_t, proposed_r );
      double rnorm2_proposed = norm2(proposed_r);

      // Evaluate the local linear model
      Tensor2Sym lm_proposed = r + jacobian * p;
      double lm_proposed_norm2 = norm2(lm_proposed);

      // Compute the ratio of the actual merit function reduction to the predicted reduction
      // to estimate the accuracy of the local model
      double rho_numerator = rnorm2 - rnorm2_proposed;
      double rho_denominator = rnorm2 - lm_proposed_norm2;

      if( rho_denominator < 0. ) {
         cout << "Negative rho denominator" << endl;
         exit(1);
      }

      double rho = rho_numerator / rho_denominator;

      // Evaluate the effectiveness of the current trust region
      if ( rho < 0.25 ) {  // Local model is not sufficiently predictive; shrink the trust region
         m_Delta *= 0.25;
      }
      else {
         if ( rho > 0.75 ) { // Local model works well; try increasing the trust region
            m_Delta *= 2.;
            if ( m_Delta > m_Delta_max ) m_Delta = m_Delta_max;
         }
         else {
            // Keep the current trust region
         }
      }

      if ( rho > eta ) {  // Norm decrease was satisfactory; accept the step
         Vbar_prime_new = proposed_solution;
         r = proposed_r;
         rnorm2 = rnorm2_proposed;
            
         converged = rnorm2 < tol2;
         new_local_model = true;
         current_model = m_hint;

         saved_estimate = m_error_estimate;
         Delta_saved = m_Delta;
      }
      else {  // Step was rejected; keep the current solution, residual and Jacobian
         new_local_model = false;
      }

      m_num_iters++;

      diagnostic_rho[m_num_iters] = rho;
      diagnostic_Delta[m_num_iters] = m_Delta;
      diagnostic_merit_value[m_num_iters] = 0.5 * rnorm2;
   }

   if ( m_num_iters >= max_iter && !converged ) {
      cout << "ElastoViscoPlasticity::updateVbar_primeTR() failed to converge" << endl;
      
      for (int n=0; n<m_num_iters; ++n) {
         cout << "merit[" << n << "] = " << diagnostic_merit_value[n] << ", Delta = "
              << diagnostic_Delta[n] << ", model = " << diagnostic_model[n]
              << ", rho = " << diagnostic_rho[n] << ", p_norm = " << diagnostic_p_norm[n]
              << ", error estimate = " << diagnostic_error[n] << endl;
      }

      delete [] diagnostic_error;
      delete [] diagnostic_p_norm;
      delete [] diagnostic_rho;
      delete [] diagnostic_model;
      delete [] diagnostic_Delta;
      delete [] diagnostic_merit_value;
      exit(1);
   }

   delete [] diagnostic_error;
   delete [] diagnostic_p_norm;
   delete [] diagnostic_rho;
   delete [] diagnostic_model;
   delete [] diagnostic_Delta;
   delete [] diagnostic_merit_value;
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


void
ElastoViscoPlasticity::evaluateSpecificModel( const int model,
                                              const Tensor2Sym& tau_bar_prime,
                                              Tensor2Sym& Dbar_prime,
                                              Tensor4LSym& Dbar_prime_deriv ) const
{
   if ( adaptiveSamplingEnabled() ) {

      std::vector<double> point(m_plasticity_model->pointDimension());
      std::vector<double> value(m_plasticity_model->valueAndDerivativeDimension());

      m_plasticity_model->packInputVector( tau_bar_prime, point );

      Constitutive::evaluateSpecificModel( model, *m_plasticity_model, point, value );
            
      m_plasticity_model->unpackOutputVector( value, Dbar_prime, Dbar_prime_deriv );
   }
   else {
      m_plasticity_model->evaluateNative( tau_bar_prime, Dbar_prime, Dbar_prime_deriv );
   }
}


double
ElastoViscoPlasticity::dot( const Tensor2Sym& tensor1,
                            const Tensor2Sym& tensor2 ) const
{
   double dotprod = 0.;

   for (int i=1; i<=3; i++) {
      for (int j=1; j<=i; ++j) {
         dotprod += tensor1(i,j) * tensor2(i,j);
      }
   }

   return dotprod;
}


void
ElastoViscoPlasticity::printTensor2Sym( const Tensor2Sym& tensor ) const
{
   for (int i=1; i<=3; i++) {
      for (int j=1; j<=i; ++j) {
         cout << "    (" << i << ", " << j << ")  " << tensor(i,j) << endl;
      }
   }
}


void
ElastoViscoPlasticity::printTensor4LSym( const Tensor4LSym& tensor ) const
{
   for (int i=1; i<=3; i++) {
      for (int j=1; j<=i; ++j) {
         for (int k=1; k<=3; ++k) {
            for (int l=1; l<=k; ++l) {
               cout << "    (" << i << ", " << j << ", " << k << ", " << l << ")  " << tensor(i,j,k,l) << endl;
            }
         }
      }
   }
   cout.flush();
}


