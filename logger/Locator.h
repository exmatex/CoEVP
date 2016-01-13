#ifndef LOCATOR_H
#define LOCATOR_H

class Logger;
class NullLogger;

class Locator {
 public:
  static void     initialize() { service_ = &nullService_; }
  static Logger&  getLogger() { return *service_; }
  static void     provide(Logger *service) {
    if (service == NULL) {
      // Revert to null service.
      service_ = &nullService_;
    }
    else {
      service_ = service;
    }
  }
  
 private:
  static Logger*     service_;
  static NullLogger  nullService_;
};

#endif   // LOCATOR_H
