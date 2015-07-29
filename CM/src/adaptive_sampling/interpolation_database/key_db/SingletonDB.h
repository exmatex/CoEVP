//
//  First attempt at using a global (singleton) redis connections.
//  This WILL evolve--it's mostly hack code now.  For instance,
//  we'll likeliy refactor to a base class so that we can suport
//  many types of databases.
//
// http://stackoverflow.com/questions/1008019/c-singleton-design-pattern

#include <vector>
#include <hiredis.h>

class SingletonDB {
 public:
  
  //  Return the single instance that was initialized in the private constructor.
  static  SingletonDB&  getInstance() {
    static  SingletonDB   instance;
    return instance;
  }

  void  sadd_sb(const char *key, const std::vector<double>& buf);
  std::vector<double> smembers_s(const char *key);

private:
  redisContext*   redis;

  SingletonDB();
  ~SingletonDB();

  //  This technique requires C++11 (can do a C++03 version too)
  SingletonDB(SingletonDB const&)    = delete;
  void operator=(SingletonDB const&) = delete;
};

