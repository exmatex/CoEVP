#ifndef included_SingletonDB_POSIX_h
#define included_SingletonDB_POSIX_h

#include "SingletonDB_Backend.h"

#include <vector>
#define uint128_t unsigned __int128

class SingletonDB_POSIX : public SingletonDB_Backend{
 public:
 
  SingletonDB_POSIX();
  ~SingletonDB_POSIX();

  virtual void  push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length);
  virtual void  erase(const uint128_t &key);
  virtual std::vector<double> pull(const uint128_t &key);
  virtual std::vector<double> pull_key(const uint128_t &key);

private:
    string dbPath;
};


#endif // included_SingletonDB_POSIX_h
