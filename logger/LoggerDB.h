//
//  First attempt at using a global (singleton) redis connections.
//  This WILL evolve--it's mostly hack code now.  For instance,
//  we'll likeliy refactor to a base class so that we can suport
//  many types of databases.
//
// http://stackoverflow.com/questions/1008019/c-singleton-design-pattern

#ifndef included_LoggerDB_h
#define included_LoggerDB_h


#include <vector>
#include <hiredis.h>
#define uint128_t unsigned __int128

class LoggerDB {
 public:
  //  Return the single instance that was initialized in the private constructor.
  static  LoggerDB&  getInstance() {
    static  LoggerDB   instance;
    return instance;
  }

  void  push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length);
  std::vector<double> pull(const uint128_t &key);
  std::vector<double> pull_key(const uint128_t &key);

 private:
  redisContext*   redis;

  LoggerDB();
 ~LoggerDB();

  //  This technique requires C++11 (can do a C++03 version too)
  LoggerDB(LoggerDB const&)       = delete;
  void operator=(LoggerDB const&) = delete;

  redisReply *pull_data(const uint128_t &key);
};

#endif // included_LoggerDB_h
