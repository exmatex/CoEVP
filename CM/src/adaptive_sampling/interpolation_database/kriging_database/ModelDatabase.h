//
// File:        ModelDataBase.h
//
// Revision:    $Revision$
// Modified:    $Date$
// Description: Abstract class to define interface to (Kriging) Model Database
//

#ifndef included_ModelDatabase_h
#define included_ModelDatabase_h

#include <memory>
//#include "base/ResponsePoint.h"
//#include "base/InterpolationModelFactory.h"
//#include "base/InterpolationModel.h"

#ifndef uint128_t
#define uint128_t unsigned __int128
#endif

//Hacky forward declarations
///TODO: Find a better solution that isn't just exposing everything in include
#ifndef included_krigalg_InterpolationModel
namespace krigalg
{
  class InterpolationModel;
  typedef std::shared_ptr<InterpolationModel> InterpolationModelPtr;
  class InterpolationModelFactory;
  typedef std::shared_ptr<InterpolationModelFactory> InterpolationModelFactoryPointer;
}
#endif
#ifndef included_krigcpl_ResponsePoint
namespace krigcpl
{
  class ResponsePoint;
}
#endif

class ModelDatabase {

public:
    virtual void insert(uint128_t & model_key, krigalg::InterpolationModelPtr krigingModel, krigcpl::ResponsePoint * point) = 0;
    ///TODO: Consider making ResponsePoint optional as it is only needed for the SingletonDB version
    virtual krigalg::InterpolationModelPtr extract(uint128_t & model_key, krigalg::InterpolationModelFactoryPointer  * newFact=nullptr) = 0;
    ///TODO: Consider operator overloading to allow for ModelDatabase[key]
    virtual void erase(uint128_t & model_key) = 0;
};

#endif // included_ModelDatabase_h
