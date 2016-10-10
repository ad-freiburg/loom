// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOGGING_LOG_H_
#define LOGGING_LOG_H_

#include <time.h>
#include <sys/timeb.h>
#include <iomanip>
#include <iostream>

#define VDEBUG 4
#define DEBUG 3
#define INFO 2
#define WARN 1
#define ERROR 0

#ifndef LOGLEVEL
#define LOGLEVEL 2
#endif

/**
 * code mostly taken from the ad_utility logger
 */

// compiler will optimize statement away if x <= LOGLEVEL
#define LOG(x) if (x > LOGLEVEL) ; else Log::log<x>()

using std::setfill;
using std::setw;

class Log {
 public:
  template<char LEVEL>
  static std::ostream& log() {
    if (LEVEL == ERROR) {
      return getLogTime(std::cerr) << getLogPrefix<LEVEL>() << " : ";
    } else {
      return getLogTime(std::cout) << getLogPrefix<LEVEL>() << " : ";
    }
  }

 private:
  static std::ostream& getLogTime(std::ostream& os) {
    struct timeb tb;
    char tl[26];

    ftime(&tb);
    ctime_r(&tb.time, tl);
    tl[19] = '.';
    tl[20] = 0;

    return os << "[" << tl << setfill('0') << setw(3) << tb.millitm << "] ";
  }

  template<char LEVEL>
  static const char* getLogPrefix() {
    switch (LEVEL) {
      case DEBUG:
        return "DEBUG";
      case INFO:
        return "INFO ";
      case WARN:
        return "WARN ";
      case ERROR:
        return "ERROR";
      default:
        return "???  ";
    }
  }
};

#endif  // LOGGING_LOG_H_
