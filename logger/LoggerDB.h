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

#ifndef LOGGERDB_H
#define LOGGERDB_H

#include <string>
#include <map>
#include <time.h>
#include "hiredis.h"
#include "Logger.h"


class LoggerDB : public Logger {
 public:
  LoggerDB(std::string db_node);
  LoggerDB(std::string db_node, std::string my_node, int my_rank);

  ~LoggerDB();
  
  virtual void  logInfo(std::string txt);
  //  Timers
  virtual void  logStartTimer(std::string);
  virtual void  logStopTimer(std::string txt);
  //  Counters
  virtual void  logCountIncr(std::string, int i=1);
  virtual void  logCount(std::string);

  virtual void  incrTimeStep(void);
  
 protected:
  std::string   hostname;
  int           id;
  int           step = 1;
  bool          isDistributed;
  redisContext *redis;

  //  Timers
  struct  TimeVals {
    timespec  ts_beg;
    timespec  ts_end;
  };
  std::map<std::string, TimeVals *> timers;
  //  Counters
  int  cnt;
  std::map<std::string, int> counters;
  
  void  connectDB(std::string db_name);
  std::string  makeKey(enum LogKeyword keyword, std::string txt);
  std::string  makeVal(float et);
  std::string  makeVal(int i);
};

#endif  // LOGGERDB_H
