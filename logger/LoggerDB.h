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
#include <time.h>
#include "hiredis.h"
#include "Logger.h"


class LoggerDB : public Logger {
 public:
  LoggerDB(std::string db_node);
  LoggerDB(std::string db_node, std::string my_node, int my_rank);

  ~LoggerDB();
  
  virtual void  logInfo(std::string txt);
  virtual void  startTimer(void);
  virtual void  logStopTimer(std::string txt);
  virtual void  incrTimeStep(void);
  
 protected:
  std::string   hostname;
  int           id;
  int           step = 1;
  bool          isDistributed;
  timespec      ts_beg, ts_end;  // these need to be indexed by keyword
  redisContext *redis;

  void  connectDB(std::string db_name);
  std::string  makeKey(enum LogKeyword keyword, std::string txt);
  std::string  makeVal(float et);
};

#endif  // LOGGERDB_H
