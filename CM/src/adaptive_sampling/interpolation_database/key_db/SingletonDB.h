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
  void  sadd_sb(const char *key, const std::vector<double>& buf) {
    redisReply* reply;
    reply = (redisReply *)redisCommand(redis, "SADD %s %lu", key, buf.size());
    if (!reply) {
      throw std::runtime_error("No connection to redis server, please start one on host'" + std::string(REDIS_HOST) + "' and port " + std::to_string(REDIS_PORT));
    }
    reply = (redisReply *)redisCommand(redis, "SADD %s %b", key, &buf[0],buf.size()*sizeof(double));
    if (!reply) {
      throw std::runtime_error("No connection to redis server, please start one on host'" + std::string(REDIS_HOST) + "' and port " + std::to_string(REDIS_PORT));
    }
    freeReplyObject(reply);
  }

  std::vector<double> smembers_s(const char *key) {
    redisReply* reply;
    //std::cout << "Connecting to redis in SingletonDB..." << std::endl;
    reply = (redisReply *)redisCommand(redis, "SMEMBERS %s", key);
    if (!reply) {
      throw std::runtime_error("No connection to redis server, please start one on host'" + std::string(REDIS_HOST) + "' and port " + std::to_string(REDIS_PORT));
    }
    if (reply->type != REDIS_REPLY_ARRAY){
      throw std::runtime_error("Wrong redis return type");
    }
    //TODO not sure what to do with other replies
    if(reply->elements < 2) {
      throw std::runtime_error("Number of redis reply elements wrong, got: " + std::to_string(reply->elements));
    }
    unsigned long sz=std::stoul((reply->element)[0]->str,nullptr,0);
    double *raw=(double *)(reply->element)[1]->str;
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

