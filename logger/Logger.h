#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <iostream>

class Logger {
 public:
  virtual void  logInfo(std::string txt) = 0;
  virtual      ~Logger() {}
};

class NullLogger : public Logger
{
 public:
  virtual void  logInfo(std::string txt) { /* NOOP */ }
};

#endif  // LOGGER_H
