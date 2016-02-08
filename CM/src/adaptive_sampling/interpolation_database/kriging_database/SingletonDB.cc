//
//  First attempt at using a global (singleton) redis connections.
//  This WILL evolve--it's mostly hack code now.  For instance,
//  we'll likeliy refactor to a base class so that we can suport
//  many types of databases.
//
// http://stackoverflow.com/questions/1008019/c-singleton-design-pattern


#include "SingletonDB.h"
#include "SingletonDB_HashMap.h"

#ifdef REDIS
#include "SingletonDB_Redis.h"
#else
#include "SingletonDB_Dummy.h"
typedef SingletonDB_Dummy SingletonDB_Redis;
#endif // REDIS

#include <iostream>
#include <fstream>
extern std::ofstream logFile;


std::string to_string_hack2(unsigned long long inVal)
{
	char buf[256];
	sprintf(buf, "%llu", inVal);
	std::string retString(buf);
	return buf;
}

static std::string uint128_to_string2(const uint128_t &in){
   uint64_t *in64 = (uint64_t *)&in; 
   return to_string_hack2((unsigned long long) *in64)+to_string_hack2((unsigned long long)*(in64+1));
}


    
//  Will eventually be something like add_points
void  SingletonDB::push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length) {
	logFile << "WRITE\t" << uint128_to_string2(key) << std::endl;
  this->backend->push(key, buf, key_length);
}

void  SingletonDB::erase(const uint128_t &key){
	logFile << "ERASE\t" << uint128_to_string2(key) << std::endl;
  this->backend->erase(key);
}

std::vector<double> SingletonDB::pull(const uint128_t &key) {
	logFile << "READ\t" << uint128_to_string2(key) << std::endl;
  return this->backend->pull(key);
}

std::vector<double> SingletonDB::pull_key(const uint128_t &key) {
  return this->backend->pull_key(key);
}

SingletonDB::SingletonDB(SingletonDBBackendEnum backType) {
	if(backType == REDIS_DB)
	{
		this->backend = new SingletonDB_Redis(false);
	}
	else if(backType == HASHMAP_DB)
	{
		this->backend = new SingletonDB_HashMap();
	}
	else if(backType == DIST_REDIS_DB)
	{
		this->backend = new SingletonDB_Redis(true);
	}
	else
	{
		std::cerr << "Invalid DB Backend Used in SingletonDB.cc" << std::endl;
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

