//
//  First attempt at using a global (singleton) redis connections.
//  This WILL evolve--it's mostly hack code now.  For instance,
//  we'll likeliy refactor to a base class so that we can suport
//  many types of databases.
//
// http://stackoverflow.com/questions/1008019/c-singleton-design-pattern

#ifdef REDIS

#include "SingletonDB.h"

#define REDIS_PORT 6379
#define REDIS_HOST "localhost"

#include <iostream>
#include <string>
#include <assert.h>
#include <stdexcept>

//  Will eventually be something like add_points
void  SingletonDB::sadd_sb(const char *key, const std::vector<double>& buf) {
  redisReply* reply;
  unsigned long sz = buf.size();
  reply = (redisReply *)redisCommand(redis, "SADD %s %b%b",
                                     key, &sz, sizeof(unsigned long),
                                     &buf[0], buf.size()*sizeof(double));
  if (!reply) {
    throw std::runtime_error("No connection to redis server, please start one on host'"
                             + std::string(REDIS_HOST) + "' and port "
                             + std::to_string(REDIS_PORT));
  }
  freeReplyObject(reply);
}


//  Will eventually be something like get_points
std::vector<double> SingletonDB::smembers_s(const char *key) {
  redisReply* reply;
  reply = (redisReply *)redisCommand(redis, "SMEMBERS %s", key);
  if (!reply) {
    throw std::runtime_error("No connection to redis server, please start one on host'"
                             + std::string(REDIS_HOST) + "' and port "
                             + std::to_string(REDIS_PORT));
  }
  if (reply->type != REDIS_REPLY_ARRAY){
    throw std::runtime_error("Wrong redis return type");
  }
  //TODO not sure what to do with other replies
  if(reply->elements < 1) {
    throw std::runtime_error("Number of redis reply elements wrong, got: "
                             + std::to_string(reply->elements));
  }
  unsigned long *sz = (unsigned long *)(reply->element)[0]->str;
  double *raw=(double *)(sz+1);
  std::vector<double> packedContainer;
  packedContainer.assign(raw, raw + *sz);
  freeReplyObject(reply);
  return packedContainer;
}


SingletonDB::SingletonDB() {
  redis = redisConnect(REDIS_HOST, REDIS_PORT);
  std::cout << "Connecting to redis in SingletonDB..." << std::endl;
}


SingletonDB::~SingletonDB() {
  redisFree(redis);
  std::cout << "Closing redis in SingletonDB..." << std::endl;
}

#endif
