#ifndef _BULK_PRESSURE_
#define _BULK_PRESSURE_

#include "EOS.h"

class BulkPressure
   : public EOS
{
 public:

    BulkPressure( const double K0,
                  const double a,
                  const double S,
                  const double Gamma0 )
       : m_K0(K0),
         m_a(a),
         m_S(S),
         m_Gamma0(Gamma0)
      {;}

   virtual double evaluate( const double mu,
                            const double e ) const;

   virtual double evaluate_dpdrho( const double mu,
                                   const double rho0,
                                   const double e ) const;

 private:

   double m_Gamma0, m_K0, m_a, m_S;

};

#endif
