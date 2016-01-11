#include "ModelDB_HashMap.h"
#include "KrigingDataBase.h"

ModelDB_HashMap::ModelDB_HashMap()
{
    //Nothing to do?
}

void ModelDB_HashMap::insert(uint128_t & model_key, krigalg::InterpolationModelPtr krigingModel, krigcpl::ResponsePoint * point)
{
    this->InterpolationModelDataBase[model_key] = krigingModel;
}

krigalg::InterpolationModelPtr ModelDB_HashMap::extract(uint128_t & model_key, krigalg::InterpolationModelFactoryPointer  * newFact)
{
    ///TODO: Error checking would probably be smart
    return this->InterpolationModelDataBase[model_key];
}

void ModelDB_HashMap::erase(uint128_t & model_key)
{
	///TODO: Do nothing because we still have the multiple NNs using the same DB issue
    //this->InterpolationModelDataBase.erase(model_key);
}
