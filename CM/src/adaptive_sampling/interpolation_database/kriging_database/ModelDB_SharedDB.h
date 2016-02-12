//
// File:        ModelDB_SharedDB.h
//
// Revision:    $Revision$
// Modified:    $Date$
// Description: Implementation of ModelDatabase interface using a SharedDB as backend
//

#ifndef included_ModelDB_SharedDB_h
#define included_ModelDB_SharedDB_h


#include "ModelDatabase.h"
#include "KrigingDataBase.h"
#include "SingletonDB_Backend.h"
#include "SingletonDB.h"


class ModelDB_SharedDB : public ModelDatabase {
private:
    SingletonDB_Backend * dbRef;
public:
    ModelDB_SharedDB(SingletonDBBackendEnum backType, int nArgs = 0, ...);
    virtual void insert(uint128_t & model_key, krigalg::InterpolationModelPtr krigingModel, krigcpl::ResponsePoint * point);
    virtual krigalg::InterpolationModelPtr extract(uint128_t & model_key, krigalg::InterpolationModelFactoryPointer  * newFact=nullptr);
    virtual void erase(uint128_t & model_key);
};


#endif // included_ModelDB_SharedDB_h
