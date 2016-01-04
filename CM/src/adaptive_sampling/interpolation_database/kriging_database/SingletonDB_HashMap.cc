#ifdef REDIS

#include "SingletonDB_HashMap.h"

     
//  Will eventually be something like add_points
void  SingletonDB_HashMap::push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length) {
  this->DBMap[key] = buf;
}

void  SingletonDB_HashMap::erase(const uint128_t &key){
  this->DBMap.erase(key);
}

std::vector<double> SingletonDB_HashMap::pull(const uint128_t &key) {
  return this->DBMap[key];
}

std::vector<double> SingletonDB_HashMap::pull_key(const uint128_t &key) {
  ///TODO: Verify there isn't anything extra I have to do
  return this->DBMap[key];
}

SingletonDB_HashMap::SingletonDB_HashMap() {

}

//  Shutdown redis databse and print some simple info about the accumulated database.
//
SingletonDB_HashMap::~SingletonDB_HashMap() {

}

#endif
