#ifndef _GAMMA_LAW_GAS_
#define _GAMMA_LAW_GAS_

#include "EOS.h"

class GammaLawGas
   : public EOS
{
 public:

   GammaLawGas( const double gamma )
      : m_gamma(gamma) {;}

   virtual double evaluate( const double compression,
                            const double internal_energy ) const;

   virtual double evaluate_dpdrho( const double compression,
                                   const double reference_density,
                                   const double internal_energy ) const;

   double gamma() const {return m_gamma;}

 private:

   double m_gamma;

};

#endif
