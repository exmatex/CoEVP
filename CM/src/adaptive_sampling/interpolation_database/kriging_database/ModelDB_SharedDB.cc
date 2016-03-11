#include "ModelDB_SharedDB.h"

#include "base/ResponsePoint.h"
#include "base/InterpolationModelFactory.h"
#include "base/InterpolationModel.h"

#include <iostream>
#include <stdexcept>
#include <cstdarg>

#ifdef REDIS
#include "SingletonDB_Redis.h"
#else
#include "SingletonDB_Dummy.h"
typedef SingletonDB_Dummy SingletonDB_Redis;
#endif // REDIS

#include "SingletonDB_HashMap.h"


ModelDB_SharedDB::ModelDB_SharedDB(SingletonDBBackendEnum backType, int nArgs, ...)
{
	//Use enum like Singleton does
	if(backType == REDIS_DB)
	{
		if(nArgs != 0)
		{
			//If we have arguments, assume the user is capable of using them correctly. This is probably a bad idea
			va_list args;
			va_start(args, nArgs);
			this->dbRef = new SingletonDB_Redis(nArgs, args);
			va_end(args);
		}
		else
		{
			//No arguments, so we want undistributed with this one not spawning a server
			this->dbRef = new SingletonDB_Redis(2, false, true);
		}
	}
	else if(backType == HASHMAP_DB)
	{
		//For now, we haven't implemented multiple entry points to a shared hashmap, so throw up a warning that we will ignore
		std::cerr << "Warning: Non-singleton Hashmap is not Shared" << std::endl;
		this->dbRef = new SingletonDB_HashMap();
		//throw std::runtime_error("Unimplemented Backend (Hashmap) used in ModelDB_SharedDB.cc");
	}
	else if(backType == DIST_REDIS_DB)
	{
		if(nArgs != 0)
		{
			//If we have arguments, assume the user is capable of using them correctly. This is probably a bad idea
			va_list args;
			va_start(args, nArgs);
			this->dbRef = new SingletonDB_Redis(nArgs, args);
			va_end(args);
		}
		else
		{
			//No arguments, so we want distributed with this one not spawning a server
			this->dbRef = new SingletonDB_Redis(2, true, true);
		}
	}
	else
	{
		//Undefined backend specified
		throw std::runtime_error("Invalid DB Backend Used in ModelDB_SharedDB.cc");
	}

}

void ModelDB_SharedDB::insert(uint128_t & model_key, krigalg::InterpolationModelPtr krigingModel, krigcpl::ResponsePoint * point)
{
	//double dKey = model_key;
	//std::cerr << "Key is " << dKey << std::endl;
    std::vector<double> packedContainer;
    krigingModel->pack(*point, packedContainer);
    dbRef->push(model_key, packedContainer, point->size());
}

krigalg::InterpolationModelPtr ModelDB_SharedDB::extract(uint128_t & model_key,krigalg::InterpolationModelFactoryPointer  * newFact)
{
    ///TODO: Error checking would probably be smart
    std::vector<double> packedContainer = dbRef->pull(model_key);
	krigalg::InterpolationModelPtr retPtr;
	retPtr = (*newFact)->build();
    retPtr->unpack(packedContainer);
    return retPtr;

}

void ModelDB_SharedDB::erase(uint128_t & model_key)
{
	///TODO: We should really erase the key, but this fixes our problems
    //this->dbRef.erase( model_key);
}

