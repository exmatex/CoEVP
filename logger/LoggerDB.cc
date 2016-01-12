//
//  First attempt at using a global (singleton) redis connections.
//  This WILL evolve--it's mostly hack code now.  For instance,
//  we'll likeliy refactor to a base class so that we can suport
//  many types of databases.
//
// http://stackoverflow.com/questions/1008019/c-singleton-design-pattern

#ifdef REDIS

#include "LoggerDB.h"

#define REDIS_PORT 6379
#define REDIS_HOST "localhost"

#include <iostream>
#include <string>
#include <assert.h>
#include <stdexcept>

static std::string uint128_to_string(const uint128_t &in){
  uint64_t *in64 = (uint64_t *)&in;
  return std::to_string(*in64)+std::to_string(*(in64+1));
}

//  Will eventually be something like add_points
void  LoggerDB::push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length) {
  redisReply* reply;
  unsigned long sz = buf.size();
  std::string skey=uint128_to_string(key);
  reply = (redisReply *)redisCommand(redis, "SADD %s %b%b%b",
                                     skey.c_str(),
                                     &sz, sizeof(unsigned long),
                                     &key_length, sizeof(unsigned long),
                                     &buf[0], buf.size()*sizeof(double));
  if (!reply) {
    throw std::runtime_error("No connection to redis server, please start one on host'"
                             + std::string(REDIS_HOST) + "' and port "
                             + std::to_string(REDIS_PORT));
  }
  freeReplyObject(reply);
}

redisReply *LoggerDB::pull_data(const uint128_t &key) {
  redisReply* reply;
  std::string skey=uint128_to_string(key);
  reply = (redisReply *)redisCommand(redis, "SMEMBERS %s", skey.c_str());
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
  return reply;
}

std::vector<double> LoggerDB::pull(const uint128_t &key) {
  redisReply* reply=pull_data(key);
  unsigned long *sz = (unsigned long *)(reply->element)[0]->str;
  double *raw=(double *)(sz+2);
  std::vector<double> packedContainer;
  packedContainer.assign(raw, raw + *sz);
  freeReplyObject(reply);
  return packedContainer;
}

std::vector<double> LoggerDB::pull_key(const uint128_t &key) {
  redisReply* reply=pull_data(key);
  unsigned long *sz = (unsigned long *)(reply->element)[0]->str;
  unsigned long *key_l = sz+1;
  double *raw=(double *)(sz+2);
  std::vector<double> packedContainer;
  packedContainer.assign(raw, raw + *key_l);
  freeReplyObject(reply);
  return packedContainer;
}
//  Open a connection to redis (and retain the connection).
//  If the databse is already populated, zero it out.  Existing entries will screw up
//  CoEVP.
//
LoggerDB::LoggerDB() {
  std::cout << "Connecting to redis in LoggerDB..." << std::endl;
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
LoggerDB::~LoggerDB() {
  std::cout << "Closing redis in LoggerDB..." << std::endl;
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
