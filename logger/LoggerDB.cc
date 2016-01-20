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


LoggerDB::LoggerDB(std::string db_node, int port)
  : isDistributed(false), id(0) {
  std::cout << "Attempting to connect to REDIS logging database on:" << std::endl;
  std::cout << "  node: " << db_node << std::endl;
  std::cout << "  port: " << port    << std::endl;
  char buffer[256];
  gethostname(buffer, 256);
  hostname = std::string(buffer);
  std::cout << "from node(" << hostname << ")/id(" << id << ")" << std::endl;

  connectDB(db_node, port);
}



LoggerDB::LoggerDB(std::string db_node, int port, std::string my_node, int my_rank)
  : isDistributed(true), hostname(my_node), id(my_rank) {
  std::cout << "Attempting to connect to REDIS logging database on:" << std::endl;
  std::cout << "  node: " << db_node << std::endl;
  std::cout << "  port: " << port    << std::endl;
  std::cout << "from node("<< my_node << ")/id(" << my_rank << ")" << std::endl;

  connectDB(db_node, port);
}


void  LoggerDB::connectDB(std::string db_node, int port) {
  redis = redisConnect(db_node.c_str(), port);
  if (redis != NULL && redis->err) {
    throw std::runtime_error("Error connecting to redis for logging, please start one on host'"
                             + db_node + "' and port " + std::to_string(port));
  }
  redisReply *reply = (redisReply *) redisCommand(redis, "DBSIZE");
  if (!reply) {
    throw std::runtime_error("No connection to redis server for logging, please start one on host'"
                             + db_node + "' and port " + std::to_string(port));
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
  std::cout << log_keywords[LOG_TIMER] << "  " << hostname << "/" << id
              << "  " << txt << "  " << et << " sec" << std::endl;
}


