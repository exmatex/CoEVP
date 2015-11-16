#ifndef _IDEAL_GAS_
#define _IDEAL_GAS_

#include "Constitutive.h"
#include "GammaLawGas.h"
#include "AddPair.h"

class IdealGas
  : public Constitutive
{
 public:

   IdealGas( ConstitutiveGlobal&     global,
             ApproxNearestNeighbors* ann,
             const double            gamma,
             const bool              use_adaptive_sampling = false );

   virtual ConstitutiveData advance( const double delta_t, const Tensor2Gen& L_new, const double, void* state );

   virtual Tensor2Sym stress( const double compression,
                              const double e,
                              const double q ) const;

   virtual Tensor2Sym stressDeviator() const;

   virtual double pressure( const double compression,
                            const double internal_energy ) const
                {return m_eos.evaluate( compression, internal_energy );}

   virtual double soundSpeedSquared( const double reference_density,
                                     const double relativeVolume,
                                     const double energy ) const;

   virtual size_t getStateSize() const {return 0;}

   virtual void getState( void* buffer ) const {};

   virtual void setState( void* buffer ) {};

 private:
   
   const GammaLawGas m_eos;

   const AddPair m_fine_scale_model;

};

#endif
