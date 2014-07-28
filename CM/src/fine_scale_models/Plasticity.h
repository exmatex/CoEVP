#ifndef _PLASTICITY_
#define _PLASTICITY_


#include "FineScale.h"
#include "tensor.h"

class Plasticity
   : public FineScale
{
   public:

      Plasticity()
         : FineScale(6,6) {};

   ~Plasticity() {};
   
      virtual Tensor2Sym tensorFunction( const Tensor2Sym& in ) const = 0;

      void evaluateNative( const Tensor2Sym& in,
                           Tensor2Sym&       out_value,
                           Tensor4LSym&      out_derivative ) const;

      void evaluate( const vector<double>& a_point,
                     vector<double>&       a_value ) const;

      virtual Tensor4LSym tensorFunctionDerivative( const Tensor2Sym& in ) const;

      virtual void getScalingsForSampling( vector<double>& input_scaling,
                                           vector<double>& output_scaling ) const = 0;

      void packInputVector( const Tensor2Sym& point,
                            vector<double>&   in ) const;

      void unpackInputVector( const vector<double>& in,
                              Tensor2Sym&           point ) const;

      void packOutputVector( Tensor2Sym&     value,
                             Tensor4LSym&    derivative,
                             vector<double>& out ) const;

      void unpackOutputVector( const vector<double>& out,
                               Tensor2Sym&           value,
                               Tensor4LSym&          derivative ) const;
};

#endif
