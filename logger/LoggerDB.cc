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
#include <iostream>
#include <stdexcept>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>


LoggerDB::LoggerDB(std::string db_node)
  : isDistributed(false), id(0), isLogging(true) {
  std::cout << "Attempting to connect to REDIS logging database on " << db_node << std::endl;
  char buffer[256];
  gethostname(buffer, 256);
  hostname = std::string(buffer);
  std::cout << "from node(" << hostname << ")/id(" << id << ")" << std::endl;

  connectDB(db_node);
}



LoggerDB::LoggerDB(std::string db_node, std::string my_node, int my_rank)
  : isDistributed(true), hostname(my_node), id(my_rank), isLogging(true) {
  std::cout << "Attempting to connect to REDIS logging database on " << db_node << std::endl;
  std::cout << "from node("<< my_node << ")/id(" << my_rank << ")" << std::endl;

  connectDB(db_node);
}


void  LoggerDB::connectDB(std::string db_node) {
  //  Get time we start logging (for stats on logging itself);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(loggingTimer.ts_beg));
  
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
  freeReplyObject(reply);
  std::cout << "...connected for logging" << std::endl;
}


LoggerDB::~LoggerDB() {
  if (isDistributed && (id==0)) {             // if MPI and rank 0 (mild kludge)...
    //  Write all the timers to REDIS
    for(auto &it : timers) {
      std::string key = it.first;
      std::string val = makeVal(it.second->et_secs, it.second->timesIncremented);
      redisReply *reply =
        (redisReply *)redisCommand(redis, "SET %s %s", key.c_str(), val.c_str());
      if (!reply) {
        std::cerr << "No connection to redis for logging...continuing" << std::endl;
      }
      freeReplyObject(reply);
    }
      
    //  Write all the counters to REDIS
    for(auto &it : counters) {
      std::string key = it.first;
      std::string val = makeVal(it.second);
      redisReply *reply =
        (redisReply *)redisCommand(redis, "SET %s %s", key.c_str(), val.c_str());
      if (!reply) {
        std::cerr << "No connection to redis for logging...continuing" << std::endl;
      }
      freeReplyObject(reply);
    }
      
#if 0
    //  Print stats for logging itself
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(loggingTimer.ts_end));
    timespec  ts_beg = loggingTimer.ts_beg;
    timespec  ts_end = loggingTimer.ts_end;
    float  et = (ts_end.tv_sec - ts_beg.tv_sec) + (ts_end.tv_nsec - ts_beg.tv_nsec) / 1e9;
    std::cout << timerCount   << " timers to REDIS in "   << et << " secs" << std::endl;
    std::cout << counterCount << " counters to REDIS in " << et << " secs" << std::endl;
    std::cout << "Timers per second:   " << (float)timerCount/et    << std::endl;
    std::cout << "Counters per second: " << (float)counterCount/et << std::endl;
#endif
    
    std::cout << "Clearing " << timers.size()   << " timers" << std::endl;
    timers.clear();
    std::cout << "Clearing " << counters.size() << " counters" << std::endl;
    counters.clear();
    
    std::cout << "Should be shutting down REDIS server..." << std::endl;
    std::cout << "Currently isn't compiled in to ease testing..." << std::endl;
    std::cout << "Make sure to put back in for 'production'..." << std::endl;
  }
#if 0
  std::cout << "Closing redis for logging..." << std::endl;
  redisReply *reply = (redisReply *) redisCommand(redis, "INFO");
  if (!reply) {
    throw std::runtime_error("No connection to redis server, something went wrong along the way!");
  }
  if (reply->type != REDIS_REPLY_STRING){
    throw std::runtime_error("Wrong redis return type when closing for logging");
  }
  reply = (redisReply *) redisCommand(redis, "SHUTDOWN");
  freeReplyObject(reply);
  redisFree(redis);
#endif
}


void  LoggerDB::logInfo(std::string txt) {
  if (!isLogging) return;
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


//  Starts a timer tied to a keyword (txt). If the timer already exists,
//  the function restarts an existing timer. It is the user's responsibility
//  to match keywords on subsequnet calls to logStopTimer.
void  LoggerDB::logStartTimer(std::string txt) {
  if (!isLogging) return;
  TimeVals  *tv;
  
  std::string key = makeKey(LOG_TIMER, txt);    // Add node, id, timestep
  std::map<std::string, TimeVals *>::iterator it = timers.find(key);
  
  if (it != timers.end()) {     // already exists
    tv = it->second;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(tv->ts_beg));
    return;
  } else {                      // create it
    tv = new TimeVals;
    tv->timesIncremented = 0;
    tv->et_secs = 0.0;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(tv->ts_beg));
    timers.insert(std::pair<std::string, TimeVals *>(key, tv));
  }
}


//  Increments a previously started timer (that has been tied to a key)
//  It is the user's responsibility to ensure keys match between starts
//  and stops.  We don't anticipate starting many timers (~100 per run),
//  so we wait until the LoggerDB destructor is called to delete them all.
//
//  TODO: this is very ugly and inefficient--fix it.
void  LoggerDB::logIncrTimer(std::string txt) {
  if (!isLogging) return;

  std::string key = makeKey(LOG_TIMER, txt);    // Add node, id, timestep
  std::map<std::string, struct TimeVals *>::iterator it = timers.find(key);

  if (it == timers.end()) {
    std::cerr << "Incrementing timer (" << txt << ") that hasn't been started?" << std::endl;
    return;
  } else {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(it->second->ts_end));
    timespec  beg = it->second->ts_beg;
    timespec  end = it->second->ts_end;
    float     et  = (end.tv_sec - beg.tv_sec) + (end.tv_nsec - beg.tv_nsec) / 1e9;
    it->second->et_secs += et;
    it->second->timesIncremented++;
  }
}


//  Increments a counter tied to a keyword (txt). If the counter already exists,
//  the function just increments it. If not, one is created (initialized to i).
//  In all case, on every call, the counter is written to the database. The
//  databse will overwrite a counter with an existing key. It is the user's
//  responsibility to match keywords on subsequnet calls.
//
//  We don't anticipate starting many counters (~100 per run), so we wait
//  until the LoggerDB object is deleted to free them all.
//
//  [TODO]  This is kind of kludgy. If it gets to be too much for the database,
//  counters can be cached and written to databse at the end of run. Possible
//  future optimization.
void  LoggerDB::logIncrCount(std::string txt, int i) {
  if (!isLogging) return;

  std::string key = makeKey(LOG_COUNT, txt);    // Add node, id, timestep
  std::map<std::string, int>::iterator it = counters.find(key);
  
  if (it != counters.end()) {     // already exists
    it->second += i;
  } else {                        // create it
    counters.insert(std::pair<std::string, int>(key, i));
  }
  return;
}


void  LoggerDB::incrTimeStep(void) {
  //  [TODO] Write the counters and timers to the database
step++;
}


std::string  LoggerDB::makeKey(enum LogKeyword keyword, std::string txt) {
  std::string  key = log_keywords[keyword];
  key += ',' + hostname + ',' + std::to_string(id) + ','
    + std::to_string(step) + ',' + txt;
  return key;
}


std::string  LoggerDB::makeVal(float et) {
  std::string  val = std::to_string(et);
  val += ',';
  //  No timestamps for now
#if 0
  std::time_t result = std::time(nullptr);
  val += std::string(" sec") + ':' + std::asctime(std::localtime(&result));
  //  Remove the expected newline provided by asctime
  if (!val.empty() && val[val.length()-1] == '\n') {
    val.erase(val.length()-1);
  }
#endif
  return val;
}


std::string  LoggerDB::makeVal(float et, int c) {
  std::string  val = std::to_string(et);
  val += ',' + std::to_string(c) + ',';
  //  No timestamps for now
#if 0
  std::time_t result = std::time(nullptr);
  val += std::asctime(std::localtime(&result));
  //  Remove the expected newline provided by asctime
  if (!val.empty() && val[val.length()-1] == '\n') {
    val.erase(val.length()-1);
  }
#endif
  return val;
}


std::string  LoggerDB::makeVal(int i) {
  std::string  val = std::to_string(i);
  val += ',';
  //  No timestamps for now
#if 0
  std::time_t result = std::time(nullptr);
  val += ':' + std::asctime(std::localtime(&result));
  //  Remove the expected newline provided by asctime
  if (!val.empty() && val[val.length()-1] == '\n') {
    val.erase(val.length()-1);
  }
#endif
  return val;
}

