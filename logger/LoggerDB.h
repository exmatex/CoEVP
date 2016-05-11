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

  virtual void  incrTimeStep(bool writeAtTimestepUpdate=true);
  
 protected:
  std::string   hostname;            //  physical host (perhaps multiple IDs on it)
  int           id;                  //  e.g. rank if MPI (Charm??)
  int           step = 1;            //  track current time step
  bool          isDistributed;       //  if MPI (or Charm?)
  bool          isLogging = false;   //  this mirrors command line switch
  redisContext *redis;               //  database connection

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
  
  //  Tracks how much time we've spent in logging.
  TimeVals  loggingTimer;

  void         connectDB(std::string db_name);
  void         writeAllTimers(void);
  void         writeAllCounters(void);
  std::string  makeKey(enum LogKeyword keyword, std::string txt);
  std::string  makeVal(float et, int c);
  std::string  makeVal(float et);
  std::string  makeVal(int i);
};

#endif  // LOGGERDB_H
