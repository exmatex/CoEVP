//
//  First attempt at using a global (singleton) redis connections.
//  This WILL evolve--it's mostly hack code now.  For instance,
//  we'll likeliy refactor to a base class so that we can suport
//  many types of databases.
//
// http://stackoverflow.com/questions/1008019/c-singleton-design-pattern

#ifndef included_SingletonDB_h
#define included_SingletonDB_h


#include <vector>
#include <hiredis.h>
#define uint128_t unsigned __int128

class SingletonDB {
 public:
  
  //  Return the single instance that was initialized in the private constructor.
  static  SingletonDB&  getInstance() {
    static  SingletonDB   instance;
    return instance;
  }

  void  push(const uint128_t &key, const std::vector<double>& buf);
  std::vector<double> pull(const uint128_t &key);

private:
  redisContext*   redis;

  SingletonDB();
  ~SingletonDB();

  //  This technique requires C++11 (can do a C++03 version too)
  SingletonDB(SingletonDB const&)    = delete;
  void operator=(SingletonDB const&) = delete;
};


#endif // included_SingletonDB_h
