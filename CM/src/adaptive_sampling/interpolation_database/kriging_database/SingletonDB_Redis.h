#ifndef included_SingletonDB_Redis_h
#define included_SingletonDB_Redis_h

#include "SingletonDB_Backend.h"

#include <vector>
#include <hiredis.h>
#define uint128_t unsigned __int128

class SingletonDB_Redis : public SingletonDB_Backend{
 public:
 
  SingletonDB_Redis();
  ~SingletonDB_Redis();

  virtual void  push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length);
  virtual void  erase(const uint128_t &key);
  virtual std::vector<double> pull(const uint128_t &key);
  virtual std::vector<double> pull_key(const uint128_t &key);

private:
  redisContext*   redis;
  FILE * redisServerHandle;

  redisReply *pull_data(const uint128_t &key);
};


#endif // included_SingletonDB_Redis_h
