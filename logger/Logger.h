#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <iostream>

class Logger {
 public:
  virtual void  logInfo(std::string txt) = 0;
  virtual void  startTimer(void) = 0;
  virtual void  logStopTimer(std::string txt) = 0;
  virtual void  incrTimeStep(void) = 0;
  virtual      ~Logger() {}

 protected:
  enum LogKeyword {
    LOG_INFO = 0,
    LOG_TIMER,
    LOG_COUNT,
    LOG_EVENT,
    NUM_LOG_KEYWORDS
  };
  std::string log_keywords[NUM_LOG_KEYWORDS] = {
    "INFO",
    "TIMER",
    "COUNT",
    "EVENT"
  };
};

class NullLogger : public Logger
{
 public:
  virtual void  logInfo(std::string txt)                           { /* NO OP */ }
  virtual void  startTimer(void)                                   { /* NO OP */ }
  virtual void  logStopTimer(std::string txt)                      { /* NO OP */ }
  virtual void  incrTimeStep(void)                                 { /* NO OP */ }
};

#endif  // LOGGER_H
