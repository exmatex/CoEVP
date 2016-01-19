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


LoggerDB::LoggerDB(std::string db_node, std::string port)
  : isDistributed(false) {
  std::cout << "Attempting to connect to REDIS logging database on:" << std::endl;
  std::cout << "  node: " << db_node << std::endl;
  std::cout << "  port: " << port    << std::endl;
}


LoggerDB::LoggerDB(std::string db_node, std::string port, std::string my_node, int my_rank)
  : isDistributed(true), my_node(my_node), my_rank(my_rank) {
  std::cout << "Attempting to connect to REDIS logging database on:" << std::endl;
  std::cout << "  node: " << db_node << std::endl;
  std::cout << "  port: " << port    << std::endl;
  std::cout << "from node("<< my_node << ")/rank(" << my_rank << ")" << std::endl;
}


void  LoggerDB::logInfo(std::string txt) {
  if (isDistributed) {
    std::cout << log_keywords[LOG_INFO] << "  " << my_node << " " << my_rank
              << "  " << txt << std::endl;
  } else {
    std::cout << log_keywords[LOG_INFO] << "  " << txt << std::endl;
  }
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
  if (isDistributed) {
    std::cout << log_keywords[LOG_TIMER] << "  " << my_node << " " << my_rank
              << "  " << txt << "  " << et << " sec" << std::endl;
  } else {
    std::cout << log_keywords[LOG_TIMER]
              << "  " << txt << "  " << et << " sec" << std::endl;
  }
}


