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

#include "LoggerDB.h"

void  LoggerDB::logInfo(std::string txt) {
  std::cout << "INFO: " << txt << std::endl;
}


