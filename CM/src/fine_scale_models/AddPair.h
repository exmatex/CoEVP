#ifndef _ADD_PAIR_
#define _ADD_PAIR_

#include "FineScale.h"


class AddPair
  : public FineScale
{
   public:

      AddPair( const int pointDimension,
               const int valueDimension )
         : FineScale(pointDimension, valueDimension) {};

      ~AddPair() {};
   
      virtual void evaluate( const std::vector<double>& points,
                             std::vector<double>&       values ) const;

      virtual void advance( const double delta_t, void* state) {};

      virtual size_t getStateSize() const {return 0;}
};

#endif
