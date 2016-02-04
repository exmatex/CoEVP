#ifndef LOGGER_H
#define LOGGER_H

#include <string>


class Logger {
 public:
  virtual void  logInfo(std::string) = 0;
  //  Timers
  virtual void  logStartTimer(std::string) = 0;
  virtual void  logStopTimer(std::string) = 0;
  //  Counters
  virtual void  logCountIncr(std::string, int i=1) = 0;
  virtual void  logCount(std::string) = 0;

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
  virtual void  logStartTimer(std::string txt)                     { /* NO OP */ }
  virtual void  logStopTimer(std::string txt)                      { /* NO OP */ }
  virtual void  logCountIncr(std::string, int i=1)                 { /* NO OP */ }
  virtual void  logCount(std::string)                              { /* NO OP */ }
  virtual void  incrTimeStep(void)                                 { /* NO OP */ }
};

#endif  // LOGGER_H
