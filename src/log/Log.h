// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOGGING_LOG_H_
#define LOGGING_LOG_H_

#include <time.h>
#include <sys/timeb.h>
#include <iomanip>
#include <iostream>

#define DEBUG 3
#define INFO 2
#define WARN 1
#define ERROR 0

#define LOG(x) Log::log(x)

class Log {
 public:
  static std::ostream& log(uint8_t lvl) {
    if (lvl == ERROR) {
      return logPrefix(std::cerr, lvl) << " : ";
    } else {
      return logPrefix(std::cout, lvl) << " : ";
    }
  }

  static std::ostream& logPrefix(std::ostream& os, uint8_t lvl) {
    struct timeb tb;
    char tl[26];

    ftime(&tb);
    ctime_r(&tb.time, tl);
    tl[19] = '.';
    tl[20] = 0;

    os << "[" << tl << std::setfill('0') << std::setw(3)
       << tb.millitm << "] ";
    switch (lvl) {
      case DEBUG:
        return os << "DEBUG";
      case INFO:
        return os << "INFO";
      case WARN:
        return os << "WARN";
      case ERROR:
        return os << "ERROR";
    }
  }

  static void setLogLvl(uint8_t l) {
    _logLvl = l;
  }

 private:
  static uint8_t _logLvl;
};

uint8_t Log::_logLvl = 5;

#endif  // LOGGING_LOG_H_
