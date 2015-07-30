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

//  Open a connection to redis (and retain the connection).
//  If the databse is already populated, zero it out.  Existing entries will screw up
//  CoEVP.
//
SingletonDB::SingletonDB() {
  std::cout << "Connecting to redis in SingletonDB..." << std::endl;
  redis = redisConnect(REDIS_HOST, REDIS_PORT);
  if (redis != NULL && redis->err) {
    throw std::runtime_error("Error connecting to redis, please start one on host'"
                             + std::string(REDIS_HOST) + "' and port "
                             + std::to_string(REDIS_PORT));
  }    
  redisReply *reply = (redisReply *) redisCommand(redis, "DBSIZE");
  if (!reply) {
    throw std::runtime_error("No connection to redis server, please start one on host'"
                             + std::string(REDIS_HOST) + "' and port "
                             + std::to_string(REDIS_PORT));
  }
  if (reply->type != REDIS_REPLY_INTEGER){
    throw std::runtime_error("Wrong redis return type");
  }
  if(reply->integer > 0) {
    std::cout << "...existing database found (" << reply->integer << " keys), erasing it" << std::endl;
    redisCommand(redis, "FLUSHDB");    // Never supposed to fail :-)
  }
}

//  Shutdown redis databse and print some simple info about the accumulated database.
//
SingletonDB::~SingletonDB() {
  std::cout << "Closing redis in SingletonDB..." << std::endl;
  redisReply *reply = (redisReply *) redisCommand(redis, "INFO");
  if (!reply) {
    throw std::runtime_error("No connection to redis server, please start one on host'"
                             + std::string(REDIS_HOST) + "' and port "
                             + std::to_string(REDIS_PORT));
  }
  if (reply->type != REDIS_REPLY_STRING){
    throw std::runtime_error("Wrong redis return type");
  }
  std::cout << reply->str << std::endl;
  redisFree(redis);
}

#endif
