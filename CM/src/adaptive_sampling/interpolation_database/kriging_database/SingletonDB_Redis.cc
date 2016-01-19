#ifdef REDIS

#include "SingletonDB_Redis.h"

#ifndef REDIS_PORT
#define REDIS_PORT 6379
#endif
#ifndef NUTCRACKER_PORT
#define NUTCRACKER_PORT 6380
#endif

#include <iostream>
#include <string>
#include <assert.h>
#include <stdexcept>
#include <unistd.h>

static std::string uint128_to_string(const uint128_t &in){
   uint64_t *in64 = (uint64_t *)&in; 
   return std::to_string(*in64)+std::to_string(*(in64+1));
}
      
//  Will eventually be something like add_points
void  SingletonDB_Redis::push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length) {
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
                             + std::string(hostBuffer) + "' and port "
                             + std::to_string(REDIS_PORT));
  }
  freeReplyObject(reply);
}

void  SingletonDB_Redis::erase(const uint128_t &key){
  redisReply* reply;
  std::string skey=uint128_to_string(key);
  reply = (redisReply *)redisCommand(redis, "DEL %s", skey.c_str());

  if (!reply) {
    throw std::runtime_error("No connection to redis server, please start one on host'"
                             + std::string(hostBuffer) + "' and port "
                             + std::to_string(REDIS_PORT));
  }
  freeReplyObject(reply);
}

redisReply *SingletonDB_Redis::pull_data(const uint128_t &key) {
  redisReply* reply;
  std::string skey=uint128_to_string(key);
  reply = (redisReply *)redisCommand(redis, "SMEMBERS %s", skey.c_str());
  if (!reply) {
    throw std::runtime_error("No connection to redis server, please start one on host'"
                             + std::string(hostBuffer) + "' and port "
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

std::vector<double> SingletonDB_Redis::pull(const uint128_t &key) {
  redisReply* reply=pull_data(key);
  //TODO what should we do with the second elemet?
  unsigned long *sz = (unsigned long *)(reply->element)[0]->str;
  double *raw=(double *)(sz+2);
  std::vector<double> packedContainer;
  packedContainer.assign(raw, raw + *sz);
  freeReplyObject(reply);
  return packedContainer;
}

std::vector<double> SingletonDB_Redis::pull_key(const uint128_t &key) {
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
SingletonDB_Redis::SingletonDB_Redis(bool distributedRedis) {
  int port;
  gethostname(hostBuffer, 256);
  if(distributedRedis)
  {
    port = NUTCRACKER_PORT;
  }
  else
  {
    port = REDIS_PORT;
  }
  std::cout << "Connecting to redis in SingletonDB_Redis..." << std::endl;
  this->redisServerHandle = nullptr;
  this->nutcrackerServerHandle = nullptr;
  redis = redisConnect(hostBuffer, port);
  
  if (redis != NULL && redis->err) {
    //If needed, spawn nutcracker
    if(distributedRedis)
	{
		std::cout << "Enabling Twemproxy..." << std::endl;
		char cmdBuffer[256];
		sprintf(cmdBuffer, "%s -c ./%s/nutcracker.yml 2>&1", NUTCRACKER_SERVER, hostBuffer);
		this->nutcrackerServerHandle = popen(cmdBuffer, "r");
		bool serverNotReady = true;
		char buffer[256];
		while(serverNotReady)
		{
			char * crashedIfNull = fgets(buffer, 256, this->nutcrackerServerHandle);
			std::string lineStr(buffer);
			size_t funLine = lineStr.find("it's time to dig another one");
			if(crashedIfNull == nullptr || funLine != std::string::npos)
			{
				serverNotReady = false;
			}
		}
	}
    //Attempt to spawn redis-server
    std::cout << "Attempting to spawn redis in SingletonDB_Redis..." << std::endl;
    char cmdBuffer[256];
    sprintf(cmdBuffer, "%s --port %d --dbfilename %s.rdb", REDIS_SERVER, REDIS_PORT, hostBuffer);
    this->redisServerHandle = popen(cmdBuffer, "r");
    bool serverNotReady = true;
    char buffer[256];
    char portLine[256];
    sprintf(portLine, " port %d", REDIS_PORT);
    while(serverNotReady)
    {
      char * crashedIfNull = fgets(buffer, 256, this->redisServerHandle);
      std::string lineStr(buffer);
      size_t portIndex = lineStr.find(portLine);
      if(crashedIfNull == nullptr || portIndex != std::string::npos)
      {
        serverNotReady = false;
      }
    }
    //Try again
    redis = redisConnect(hostBuffer, port);
    if (redis != NULL && redis->err) {
		std::string hostString(hostBuffer);
      throw std::runtime_error("Error connecting to redis, please start one on host'"
                             + hostString  + "' and port "
                             + std::to_string(port));
    }
  }
  if(!distributedRedis)
  {
    redisReply *reply = (redisReply *) redisCommand(redis, "DBSIZE");
    if (!reply) {
      throw std::runtime_error("No connection to redis server, please start one on host'"
                             + std::string(hostBuffer) + "' and port "
                             + std::to_string(port));
    }
    if (reply->type != REDIS_REPLY_INTEGER){
      throw std::runtime_error("Wrong redis return type");
    }
    if(reply->integer > 0) {
      std::cout << "...existing database found (" << reply->integer << " keys), erasing it" << std::endl;
      redisCommand(redis, "FLUSHDB");    // Never supposed to fail :-)
    }
  }

}

//  Shutdown redis databse and print some simple info about the accumulated database.
//
SingletonDB_Redis::~SingletonDB_Redis() {
  std::cout << "Closing redis in SingletonDB_Redis..." << std::endl;
  redisReply *reply = (redisReply *) redisCommand(redis, "INFO");
  if (!reply) {
    throw std::runtime_error("No connection to redis server, please start one on host'"
                             + std::string(hostBuffer) + "' and port "
                             + std::to_string(REDIS_PORT));
  }
  if (reply->type != REDIS_REPLY_STRING){
    throw std::runtime_error("Wrong redis return type");
  }
  std::cout << reply->str << std::endl;
  if(this->redisServerHandle != nullptr)
  {
    std::cout << "Closing redis server in SingletonDB_Redis..." << std::endl;
    redisReply * reply = (redisReply *) redisCommand(redis, "SHUTDOWN");
    pclose(this->redisServerHandle);
  }
  if(this->nutcrackerServerHandle != nullptr)
  {
    std::cout << "Closing nutcracker server in SingletonDB_Redis..." << std::endl;
    pclose(this->nutcrackerServerHandle);
  }
  redisFree(redis);
}

#endif
