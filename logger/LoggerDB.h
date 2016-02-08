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
  virtual void  logIncrTimer(std::string txt);
  //  Counters
  virtual void  logIncrCount(std::string, int i=1);

  virtual void  incrTimeStep(void);
  
 protected:
  std::string   hostname;            //  Physical host (perhaps multiple IDs on it)
  int           id;                  //  Rank if MPI (Charm??)
  int           step = 1;            //  Track current time step
  bool          isDistributed;       //  If MPI (or Charm?)
  bool          isLogging = false;   //  This mirrors comman'd line switch
  redisContext *redis;

  //  Timers
  struct  TimeVals {
    timespec  ts_beg;
    timespec  ts_end;
    float     et_secs;
    int       timesIncremented;
  };
  std::map<std::string, TimeVals *> timers;

  //  Counters
  std::map<std::string, int> counters;
  
  //  Stats for the logger itself.
  TimeVals  loggingTimer;

  void         connectDB(std::string db_name);
  std::string  makeKey(enum LogKeyword keyword, std::string txt);
  std::string  makeVal(float et, int c);
  std::string  makeVal(float et);
  std::string  makeVal(int i);
};

#endif  // LOGGERDB_H
