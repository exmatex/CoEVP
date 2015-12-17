#ifndef included_SingletonDB_HashMap_h
#define included_SingletonDB_HashMap_h

#include "SingletonDB_Backend.h"
#include "KeyHash.h"

#include <vector>
#include <unordered_map>

class SingletonDB_HashMap : public SingletonDB_Backend{
 public:
 
  SingletonDB_HashMap();
  ~SingletonDB_HashMap();

  virtual void  push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length);
  virtual void  erase(const uint128_t &key);
  virtual std::vector<double> pull(const uint128_t &key);
  virtual std::vector<double> pull_key(const uint128_t &key);

private:
  std::unordered_map<uint128_t, std::vector<double> > DBMap;
};


#endif // included_SingletonDB_HashMap_h
