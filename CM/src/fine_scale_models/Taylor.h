#ifndef _TAYLOR_
#define _TAYLOR_

#include "Plasticity.h"


class Taylor
  : public Plasticity
{
   public:

        Taylor(const double D_0,
               const double m,
               const double g)
           : m_D_0(D_0), m_m(m), m_g(g) {};

      ~Taylor() {;}
   
      virtual Tensor2Sym tensorFunction( const Tensor2Sym& in ) const;

      //      virtual Tensor4LSym tensorFunctionDerivative( const Tensor2Sym& in ) const;

      virtual void getScalingsForSampling( vector<double>& input_scaling,
                                           vector<double>& output_scaling ) const;

      virtual void advance( const double delta_t, void* state) const {};

      virtual size_t getStateSize() const {return 0;}

   private:

      double m_D_0;      // Reference strain rate

      double m_m;        // Sensitivity

      double m_g;        // Hardness
};

#endif
