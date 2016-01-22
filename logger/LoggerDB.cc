//
//  Logging to a REDIS Databases
//
//  Provides functionality to:
//    - connect to a logging-specific REDIS database
//    - log CoEVP-specific events to the database
//    
//  Doesn't provide any capability to analyze the collected logs.
//  Python utility code will be use for that task.
//

#include "LoggerDB.h"
#include <unistd.h>
#include <stdexcept>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>


LoggerDB::LoggerDB(std::string db_node)
  : isDistributed(false), id(0) {
  std::cout << "Attempting to connect to REDIS logging database on " << db_node << std::endl;
  char buffer[256];
  gethostname(buffer, 256);
  hostname = std::string(buffer);
  std::cout << "from node(" << hostname << ")/id(" << id << ")" << std::endl;

  connectDB(db_node);
}



LoggerDB::LoggerDB(std::string db_node, std::string my_node, int my_rank)
  : isDistributed(true), hostname(my_node), id(my_rank) {
  std::cout << "Attempting to connect to REDIS logging database on " << db_node << std::endl;
  std::cout << "from node("<< my_node << ")/id(" << my_rank << ")" << std::endl;

  connectDB(db_node);
}


void  LoggerDB::connectDB(std::string db_node) {
  //  Separate the hostname from the port (input must be of form "host:nnn")
  std::vector<std::string> elems;
  std::stringstream ss(db_node);
  std::string item;
  while (std::getline(ss, item, ':')) {
    elems.push_back(item);
  }
  if (elems.size() != 2) {
    throw std::runtime_error("Invalid logging host:port (" + db_node + ") specified");
  }
  //  Try connecting to the logging database
  redis = redisConnect(elems[0].c_str(), std::stoi(elems[1]));
  if (redis != NULL && redis->err) {
    throw std::runtime_error("Error connecting to redis for logging, please start one on " + db_node);
  }
  redisReply *reply = (redisReply *) redisCommand(redis, "DBSIZE");
  if (!reply) {
    throw std::runtime_error("No connection to redis server for logging, please start one on " + db_node);
  }
  if (reply->type != REDIS_REPLY_INTEGER){
    throw std::runtime_error("Wrong redis return type (when opening for logging)");
  }
  if(reply->integer > 0) {
    std::cout << "...existing database found (" << reply->integer << " keys), erasing it" << std::endl;
    redisCommand(redis, "FLUSHDB");    // Never supposed to fail :-)
  }
  std::cout << "...connected for logging" << std::endl;
}


LoggerDB::~LoggerDB() {
  std::cout << "Closing redis for logging..." << std::endl;
  redisReply *reply = (redisReply *) redisCommand(redis, "INFO");
  if (!reply) {
    throw std::runtime_error("No connection to redis server, something went wrong along the way!");
  }
  if (reply->type != REDIS_REPLY_STRING){
    throw std::runtime_error("Wrong redis return type when closing for logging");
  }
  reply = (redisReply *) redisCommand(redis, "SHUTDOWN");
  redisFree(redis);
}


void  LoggerDB::logInfo(std::string txt) {
  std::cout << log_keywords[LOG_INFO] << "  " << hostname << "/" << id
            << "  " << txt << std::endl;
  redisReply *reply = (redisReply *)redisCommand(redis, "SET %s %s",
                                                 log_keywords[LOG_INFO].c_str(),
                                                 "See if this crusty REDIS thing works");
  if (!reply) {
    std::cerr << "No connection to redis for logging...continuing" << std::endl;
  }
  freeReplyObject(reply);
}


//  "Fast" timer (no name lookups)
//  Starts a default timer.
void  LoggerDB::startTimer(void) {
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_beg);
}


//  "Fast" timer (no name lookups)
//  Stops the default timer and logs an TIMER event.
//  Assumes that you've started the default timer previously.
void  LoggerDB::logStopTimer(std::string txt) {
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_end);
  float  et = (ts_end.tv_sec - ts_beg.tv_sec) + (ts_end.tv_nsec - ts_beg.tv_nsec) / 1e9;
  std::string key = makeKey(LOG_TIMER, txt);
  std::string val = makeVal(et);
  redisReply *reply = (redisReply *)redisCommand(redis, "SADD %s %s", key.c_str(), val.c_str());
  if (!reply) {
    std::cerr << "No connection to redis for logging...continuing" << std::endl;
  }
  freeReplyObject(reply);
  std::cout << key << " ----- " << val << std::endl;
}


std::string  LoggerDB::makeKey(enum LogKeyword keyword, std::string txt) {
  std::string  key = log_keywords[keyword];
  key += ':' + hostname + ':' + std::to_string(id) + ':' + txt;
  return key;
}


std::string  LoggerDB::makeVal(float et) {
  std::string  val = std::to_string(et);
  std::time_t result = std::time(nullptr);
  val += std::string(" sec   ") + std::asctime(std::localtime(&result));
  //  Remove the expected newline provided by asctime
  if (!val.empty() && val[val.length()-1] == '\n') {
    val.erase(val.length()-1);
  }
  return val;
}


