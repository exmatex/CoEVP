#ifndef included_SingletonDB_Backend_h
#define included_SingletonDB_Backend_h


#include <vector>
#ifdef REDIS
#include <hiredis.h>
#endif
#define uint128_t unsigned __int128

class SingletonDB_Backend {
 public:
  virtual void  push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length) = 0;
  virtual void  erase(const uint128_t &key) = 0;
  virtual std::vector<double> pull(const uint128_t &key) = 0;
  virtual std::vector<double> pull_key(const uint128_t &key) = 0;
};


#endif // included_SingletonDB_Backend_h
