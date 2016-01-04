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
#include "KeyHash.h"

#include <unordered_map>

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
