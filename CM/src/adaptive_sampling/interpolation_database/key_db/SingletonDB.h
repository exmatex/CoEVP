//
//  First attempt at using a global (singleton) redis connections.
//  This WILL evolve--it's mostly hack code now.  For instance,
//  we'll likeliy refactor to a base class so that we can suport
//  many types of databases.
//
// Allen McPherson
// Los Alamos National Laboratory

// http://stackoverflow.com/questions/1008019/c-singleton-design-pattern

#define REDIS_PORT 6379
#define REDIS_HOST "localhost"

#include <iostream>
#include <string>
#include <vector>
#include <assert.h>
#include <hiredis.h>



class SingletonDB {
 public:
  static  SingletonDB&  getInstance() {
    static  SingletonDB   instance;
    return instance;
  }
  void  sadd_sb(const char *key, const double *buf, int sz) {
    redisReply* reply;
    reply = (redisReply *)redisCommand(redis, "SADD %s %b", key, buf, sz);
    if (reply==NULL) {
      std::cerr << "No connection to redis server" << std::endl;
      assert(reply==NULL);
    }
  }
  std::vector<double> smembers_s(const char *key, int sz) {
    redisReply* reply;
    reply = (redisReply *)redisCommand(redis, "SMEMBERS %s", key);
    if (reply==NULL) {
      std::cerr << "No connection to redis server" << std::endl;
      assert(reply==NULL);
    }
    assert(reply->type == REDIS_REPLY_ARRAY);
    //TODO not sure what to do with other replies
    assert(reply->elements == 1);
    double *raw=(double *)(reply->element)[0]->str;
    std::vector<double> packedContainer;
    packedContainer.assign(raw, raw + sz);
    freeReplyObject(reply);
    return packedContainer;
  }
private:
  redisContext*   redis;

  SingletonDB() {
    redis = redisConnect(REDIS_HOST, REDIS_PORT);
    std::cout << "Connecting to redis in SingletonDB..." << std::endl;
  }
  ~SingletonDB() {
    redisFree(redis);
    std::cout << "Closing redis in SingletonDB..." << std::endl;
  }
  //  This technique requires C++11 (can do a C++03 version too)
  SingletonDB(SingletonDB const&)    = delete;
  void operator=(SingletonDB const&) = delete;
};

