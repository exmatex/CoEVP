#ifndef _FINE_SCALE_
#define _FINE_SCALE_

using namespace std;

#include <vector>

class FineScale
{
   public:

      FineScale( const int pointDimension,
                 const int valueDimension )
         : m_pointDimension(pointDimension),
           m_valueDimension(valueDimension) {};

      ~FineScale() {};

      virtual void evaluate( const std::vector<double>& points,
                             vector<double>&            values ) const = 0;

      int pointDimension() const {return m_pointDimension;}

      int valueDimension() const {return m_valueDimension;}

      int valueAndDerivativeDimension() const {return m_valueDimension*(m_pointDimension + 1);}

   protected:

   int m_pointDimension;
   int m_valueDimension;
};

#endif
