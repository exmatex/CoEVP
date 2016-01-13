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
//  Uses singleton design pattern (as in SingletonDB). Refactor all?
// http://stackoverflow.com/questions/1008019/c-singleton-design-pattern

#ifndef LOGGERDB_H
#define LOGGERDB_H

#include "Logger.h"

class LoggerDB : public Logger {
 public:
  virtual void  logInfo(std::string txt);
};

#endif  // LOGGERDB_H
