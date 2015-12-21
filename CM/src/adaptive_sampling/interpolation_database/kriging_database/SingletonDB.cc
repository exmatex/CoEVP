//
//  First attempt at using a global (singleton) redis connections.
//  This WILL evolve--it's mostly hack code now.  For instance,
//  we'll likeliy refactor to a base class so that we can suport
//  many types of databases.
//
// http://stackoverflow.com/questions/1008019/c-singleton-design-pattern

#ifdef REDIS

#include "SingletonDB.h"
#include "SingletonDB_Redis.h"
#include "SingletonDB_HashMap.h"

    
//  Will eventually be something like add_points
void  SingletonDB::push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length) {
  this->backend->push(key, buf, key_length);
}

void  SingletonDB::erase(const uint128_t &key){
  this->backend->erase(key);
}

std::vector<double> SingletonDB::pull(const uint128_t &key) {
  return this->backend->pull(key);
}

std::vector<double> SingletonDB::pull_key(const uint128_t &key) {
  return this->backend->pull_key(key);
}

SingletonDB::SingletonDB(SingletonDBBackendEnum backType) {
	if(backType == REDIS_DB)
	{
		this->backend = new SingletonDB_Redis();
	}
	else if(backType == HASHMAP_DB)
	{
		this->backend = new SingletonDB_HashMap();
	}
	else
	{
		///TODO: Throw error
	}
}

//  Shutdown redis databse and print some simple info about the accumulated database.
//
SingletonDB::~SingletonDB() {
  if(this->backend != nullptr)
  {
    delete this->backend;
  }
}

#endif
