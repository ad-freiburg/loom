// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PBUTIL_LOG_LOG_H_
#define PBUTIL_LOG_LOG_H_

#include <sys/timeb.h>
#include <time.h>
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

// compiler will optimize statement away if x <= LOGLEVEL
#define LOG(x) if (x > LOGLEVEL) ; else pbutil::Log::log<x>()

using std::setfill;
using std::setw;

namespace pbutil {

class Log {
 public:
  template <char L>
  static std::ostream& log() {
    return getLogHead(L == ERROR ? std::cerr : std::cout) << getLogName<L>();
  }

 private:
  static std::ostream& getLogHead(std::ostream& os) {
    struct timeb tb;
    char tl[26];
    ftime(&tb);
    ctime_r(&tb.time, tl);
    tl[19] = '.';
    tl[20] = 0;
    return os << "[" << tl << setfill('0') << setw(3) << tb.millitm << "] ";
  }

  template <char LEVEL>
  static const char* getLogName() {
    switch (LEVEL) {
      case DEBUG:
        return "DEBUG: ";
      case INFO:
        return "INFO : ";
      case WARN:
        return "WARN : ";
      case ERROR:
        return "ERROR: ";
      default:
        return "(\?\?) : ";
    }
  }
};
}

#endif  // PBUTIL_LOG_LOG_H_
