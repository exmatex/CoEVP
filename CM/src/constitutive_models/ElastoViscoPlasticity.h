#ifndef _ELASTO_VISCO_PLASTICITY_
#define _ELASTO_VISCO_PLASTICITY_

#include <fstream>
#include "Constitutive.h"
#include "Plasticity.h"
#include "EOS.h"


class ElastoViscoPlasticity
   : public Constitutive
{
   public:

      ElastoViscoPlasticity( const Tensor2Gen& L,
                             const double      bulk_modulus,
                             const double      shear_modulus,
                             const EOS*        eos,
                             const Plasticity* fine_scale_model,
                             const bool        use_adaptive_sampling = false );

      ~ElastoViscoPlasticity();

      virtual Tensor2Sym stress( const double compression,
                                 const double e,
                                 const double q ) const;

      virtual Tensor2Sym stressDeviator() const;
      
      virtual void advance( const double delta_t );

      virtual void setNewVelocityGradient( const Tensor2Gen& L_new );

      virtual void setVolumeChange( const double volume_change ) {m_volume_change = volume_change;}

      virtual double pressure( const double compression,
                               const double internal_energy ) const
            {return m_eos_model->evaluate( compression, internal_energy );}

      virtual double soundSpeedSquared( const double reference_density,
                                        const double relativeVolume,
                                        const double energy ) const;
   
      int numNewtonIterations() const {return m_num_iters;}

   private:

      inline Tensor2Sym tauBarPrime( const double      a,
                                     const Tensor2Sym& Vbar_prime ) const;

      inline double a(double J) const {return pow(J,1./3.);}

      inline double J(const Tensor2Sym& V) const {return det(V);}

      void updateVbar_prime( const Tensor2Sym& Vbar_prime_old,
                             const Tensor2Sym& Dprime_new,
                             const Tensor2Gen& R_new,
                             const double      a_new,
                             const double      delta_t,
                             Tensor2Sym&       Vbar_prime_new,
                             Tensor2Sym&       Dbar_prime,
                             Tensor2Gen&       Wbar );

      void getDoglegStep( const Tensor2Sym&  r,
                          const Tensor4LSym& J,
                          const double       Delta,
                          Tensor2Sym&        p ) const;

      void convertToFine( const Tensor2Sym& in,
                          const Tensor2Gen& R,
                          Tensor2Sym&       out ) const;

      void convertToCoarse( const Tensor2Gen& in,
                            const Tensor2Gen& R,
                            Tensor2Gen&       out ) const;

      void convertToCoarse( const Tensor2Sym& in,
                            const Tensor2Gen& R,
                            Tensor2Sym&       out ) const;

      void evaluateStretchRHS( const Tensor2Sym& D_prime,
                               const Tensor2Sym& Dbar_prime,
                               const Tensor2Gen& R,
                               Tensor2Sym&       rhs ) const;

      void computeResidual( const Tensor2Sym& Vbar_prime_new,
                            const Tensor2Sym& Vbar_prime_old,
                            const Tensor2Sym& Dprime_new,
                            const Tensor2Sym& Dbar_prime_new,
                            const Tensor2Gen& R_new,
                            const double      a_delta_t,
                            Tensor2Sym&       residual ) const;

      void computeJacobian( const Tensor4LSym& Dbar_deriv,
                            const double       a,
                            const double       delta_t,
                            Tensor4LSym&       jacobian ) const;

      void solveJacobianSystem( const Tensor4LSym& jacobian,
                                const Tensor2Sym&  rhs,
                                Tensor2Sym&        solution ) const;

      void updateJ( const double      J_old,
                    const Tensor2Sym& D,
                    const double      delta_t,
                    double&           J_new ) const;

      void evaluateWR( const Tensor2Gen& Vbar_prime,
                       const Tensor2Gen& Vbar_prime_dot,
                       const Tensor2Gen& Dbar_prime,
                       const Tensor2Gen& W,
                       const Tensor2Gen& Wbar,
                       const Tensor2Gen& R,
                       const double      a,
                       Tensor2Gen&       WR ) const;

      void updateR( const Tensor2Gen& R_old,
                    const Tensor2Gen& WR,
                    const double      delta_t,
                    Tensor2Gen&       R_new ) const;

      void evaluateFineScaleModel( const Tensor2Sym& tau_bar_prime,
                                   Tensor2Sym&       Dbar_prime,
                                   Tensor4LSym&      Dbar_prime_deriv ) const;

      void evaluateSpecificModel( const int         model,
                                  const Tensor2Sym& tau_bar_prime,
                                  Tensor2Sym&       Dbar_prime,
                                  Tensor4LSym&      Dbar_prime_deriv ) const;

      double dot( const Tensor2Sym& tensor1,
                  const Tensor2Sym& tensor2 ) const;
      double norm2( const Tensor2Sym& tensor ) const {return dot(tensor,tensor);}
      void printTensor2Sym( const Tensor2Sym& tensor ) const;
      void printTensor4LSym( const Tensor4LSym& tensor ) const;

      double m_Delta_max;
      double m_Delta;

      Tensor2Sym m_D_old;    // Velocity gradient symmetric part (coarse-scale strain rate) at old time
      Tensor2Gen m_W_old;    // Velocity gradient skew part (coarse-scale spin) at old time
      Tensor2Sym m_D_new;    // Velocity gradient symmetric part (coarse-scale strain rate) at new time
      Tensor2Gen m_W_new;    // Velocity gradient skew part (coarse-scale spin) at new time
      Tensor2Gen m_R;    // Rotation from fine-scale to coarse-scale models 
      double m_J;        // Stretch determinant

      int m_num_iters;   // Number of Newton iterations performed in the most recent solve

      Tensor2Sym m_Vbar_prime;       // Stretch deviator in fine-scale frame
      Tensor2Sym m_Vbar_prime_dot;   // Time derivative of stretch deviator in fine-scale frame
      Tensor2Sym m_Dbar_prime;       // Deviator of the fine-scale strain rate
      Tensor2Gen m_Wbar;             // Fine-scale spin
      double m_volume_change;        // Ratio of new to old volume;

      double m_tol;

      double m_K;         // Bulk modulus
      double m_G;         // Shear modulus

      const EOS* m_eos_model;   // Equation of state model

      const Plasticity* m_plasticity_model;  // Fine-scale model
};

#endif
