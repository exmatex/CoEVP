#ifndef included_SingletonDB_Dummy_h
#define included_SingletonDB_Dummy_h

#include "SingletonDB_Backend.h"
#include "KeyHash.h"

#include <iostream>

class SingletonDB_Dummy : public SingletonDB_Backend{
 public:
 
  SingletonDB_Dummy(int dummyArg = 0, ...);
  ~SingletonDB_Dummy();

  virtual void  push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length);
  virtual void  erase(const uint128_t &key);
  virtual std::vector<double> pull(const uint128_t &key);
  virtual std::vector<double> pull_key(const uint128_t &key);

};

#endif // included_SingletonDB_Dummy_h
