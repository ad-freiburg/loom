// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_LOG_LOG_H_
#define UTIL_LOG_LOG_H_

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

#define VDEBUG 4
#define DEBUG 3
#define INFO 2
#define WARN 1
#define ERROR 0

#ifndef LOGLEVEL
#define LOGLEVEL 2
#endif

// compiler will optimize statement away if x > LOGLEVEL
#define LOG(x) if (x <= LOGLEVEL) util::Log<x>().log()
#define LOGTO(x, os) if (x <= LOGLEVEL) util::Log<x>(&os).log()

using std::setfill;
using std::setw;
using namespace std::chrono;

namespace util {

static const char* LOGS[] = {"ERROR", "WARN ", "INFO ", "DEBUG", "DEBUG"};

template <char LVL>
class Log {
 public:
  Log() { if (LVL < INFO) os = &std::cerr; else os = &std::cout; }
  Log(std::ostream* s) { os = s; }
  ~Log() { buf << std::endl; (*os) << buf.str(); }
  std::ostream& log() { return ts() << LOGS[(size_t)LVL] << ": "; }

 private:
  std::ostream* os;
  std::ostringstream buf;
  std::ostream& ts() {
    char tl[20];
    auto n = system_clock::now();
    time_t tt = system_clock::to_time_t(n);
    int m = duration_cast<milliseconds>(n-time_point_cast<seconds>(n)).count();
    struct tm t = *localtime(&tt);
    strftime(tl, 20, "%Y-%m-%d %H:%M:%S", &t);
    return buf << "[" << tl << "." << setfill('0') << setw(3) << m << "] ";
  }
};
}

#endif  // UTIL_LOG_LOG_H_
