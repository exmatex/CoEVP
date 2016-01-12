//
// File:        ModelDB_SingletonDB.h
//
// Revision:    $Revision$
// Modified:    $Date$
// Description: Implementation of ModelDatabase interface using a SingletonDB as backend
//

#ifndef included_ModelDB_SingletonDB_h
#define included_ModelDB_SingletonDB_h


#include "ModelDatabase.h"
#include "KrigingDataBase.h"
#include "SingletonDB.h"


class ModelDB_SingletonDB : public ModelDatabase {
private:
    SingletonDB & dbRef;
public:
    ModelDB_SingletonDB();
    virtual void insert(uint128_t & model_key, krigalg::InterpolationModelPtr krigingModel, krigcpl::ResponsePoint * point);
    virtual krigalg::InterpolationModelPtr extract(uint128_t & model_key, krigalg::InterpolationModelFactoryPointer  * newFact=nullptr);
    virtual void erase(uint128_t & model_key);
};


#endif // included_ModelDB_SingletonDB_h
