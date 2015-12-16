//
// File:        ModelDB_Hashmap.h
//
// Revision:    $Revision$
// Modified:    $Date$
// Description: Implementation of ModelDatabase interface using a Hashmap as backend
//

#ifndef included_ModelDB_HashMap_h
#define included_ModelDB_HashMap_h

#include "ModelDatabase.h"

#include <unordered_map>

namespace std {

   // Defining a hash function in order to use uint128_t as a key
   // for an std::unordered_map
   template <>
      struct hash<uint128_t>
      {
         std::size_t operator()(const uint128_t& in) const
            {
               return (size_t)in;
            }
      };


}

class ModelDB_HashMap : public ModelDatabase {
public:
    ModelDB_HashMap();
    virtual void insert(uint128_t & model_key, krigalg::InterpolationModelPtr krigingModel, krigcpl::ResponsePoint * point);
    virtual krigalg::InterpolationModelPtr extract(uint128_t & model_key, krigalg::InterpolationModelFactoryPointer  * newFact=nullptr);
    virtual void erase(uint128_t & model_key);
private:
    std::unordered_map<uint128_t, krigalg::InterpolationModelPtr> InterpolationModelDataBase;
};

#endif // included_ModelDB_HashMap_h
