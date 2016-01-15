#include "ModelDB_SingletonDB.h"

#include "base/ResponsePoint.h"
#include "base/InterpolationModelFactory.h"
#include "base/InterpolationModel.h"

#include <iostream>

ModelDB_SingletonDB::ModelDB_SingletonDB()
: dbRef(SingletonDB::getInstance())
{
    //Nothing in here as the SingletonDB ref needs to be in the initializer list
}

void ModelDB_SingletonDB::insert(uint128_t & model_key, krigalg::InterpolationModelPtr krigingModel, krigcpl::ResponsePoint * point)
{
	//double dKey = model_key;
	//std::cerr << "Key is " << dKey << std::endl;
    std::vector<double> packedContainer;
    krigingModel->pack(*point, packedContainer);
    dbRef.push(model_key, packedContainer, point->size());
}

krigalg::InterpolationModelPtr ModelDB_SingletonDB::extract(uint128_t & model_key,krigalg::InterpolationModelFactoryPointer  * newFact)
{
    ///TODO: Error checking would probably be smart
    std::vector<double> packedContainer = dbRef.pull(model_key);
	krigalg::InterpolationModelPtr retPtr;
	retPtr = (*newFact)->build();
    retPtr->unpack(packedContainer);
    return retPtr;

}

void ModelDB_SingletonDB::erase(uint128_t & model_key)
{
	///TODO: We should really erase the key, but this fixes our problems
    //this->dbRef.erase( model_key);
}

