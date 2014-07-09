#ifndef _EOS_
#define _EOS_

class EOS
{
 public:

   EOS() {;}

   ~EOS() {;}

   virtual double evaluate( const double mu,
                            const double e ) const = 0;

   virtual double evaluate_dpdrho( const double mu,
                                   const double rho0,
                                   const double e ) const = 0;
};

#endif
