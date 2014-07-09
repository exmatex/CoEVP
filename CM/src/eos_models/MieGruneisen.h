#ifndef _MIEGRUNEISEN_
#define _MIEGRUNEISEN_

#include "EOS.h"

class MieGruneisen
   : public EOS
{
   public:

      MieGruneisen( const double k1,
                    const double k2,
                    const double k3,
                    const double Gamma )
         : m_k1(k1),
           m_k2(k2),
           m_k3(k3),
           m_Gamma(Gamma) {;}

      ~MieGruneisen() {};

      virtual double evaluate( const double mu,
                               const double e ) const;

      virtual double evaluate_dpdrho( const double mu,
                                      const double rho0,
                                      const double e ) const;

   private:

      double m_k1, m_k2, m_k3, m_Gamma;
};

#endif
