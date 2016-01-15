#include "ModelDB_HashMap.h"
#include "KrigingDataBase.h"

#include "base/ResponsePoint.h"
#include "base/InterpolationModelFactory.h"
#include "base/InterpolationModel.h"

ModelDB_HashMap::ModelDB_HashMap()
{
    //Nothing to do?
}

void ModelDB_HashMap::insert(uint128_t & model_key, krigalg::InterpolationModelPtr krigingModel, krigcpl::ResponsePoint * point)
{
    std::vector<double> packedContainer;
    krigingModel->pack(*point, packedContainer); 
	this->InterpolationModelDataBase[model_key] = packedContainer;
}

krigalg::InterpolationModelPtr ModelDB_HashMap::extract(uint128_t & model_key, krigalg::InterpolationModelFactoryPointer  * newFact)
{
    ///TODO: Error checking would probably be smart
	std::unordered_map<uint128_t, std::vector<double>>::const_iterator iter = this->InterpolationModelDataBase.find(model_key);
	if(iter == this->InterpolationModelDataBase.end())
	{
		std::cerr << "Error: Searching for Non-Existing Key\n" << std::endl;
	}
    std::vector<double> packedContainer = iter->second;
	krigalg::InterpolationModelPtr retPtr;
	retPtr = (*newFact)->build();
    retPtr->unpack(packedContainer);
    return retPtr;
}

void ModelDB_HashMap::erase(uint128_t & model_key)
{
	///TODO: Do nothing because we still have the multiple NNs using the same DB issue
    //this->InterpolationModelDataBase.erase(model_key);
}
