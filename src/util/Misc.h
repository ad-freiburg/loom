// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_MISC_H_
#define UTIL_MISC_H_

#include <chrono>

#define UNUSED(expr) do { (void)(expr); } while (0)
#define TIME() std::chrono::high_resolution_clock::now()
#define TOOK(t1, t2) (std::chrono::duration_cast<microseconds>(t2 - t1).count() / 1000.0)

namespace util {

// cached first 10 powers of 10
static int pow10[10] = {
    1,           10,          100,           1000,          10000,
    100000,      1000000,     10000000,      100000000,     1000000000};

// _____________________________________________________________________________
inline uint64_t factorial(uint64_t n) {
  if (n == 1) return n;
  return n * factorial(n - 1);
}

// _____________________________________________________________________________
inline uint64_t atoul(const char* p) {
  uint64_t ret = 0;

  while (*p) {
    ret = ret * 10 + (*p++ - '0');
  }

  return ret;
}

// _____________________________________________________________________________
inline float atof(const char* p, uint8_t mn) {
  // this atof implementation works only on "normal" float strings like
  // 56.445 or -345.00, but should be faster than std::atof
  float ret = 0.0;
  bool neg = false;
  if (*p == '-') {
    neg = true;
    p++;
  }

  while (*p >= '0' && *p <= '9') {
    ret = ret * 10.0 + (*p - '0');
    p++;
  }

  if (*p == '.') {
    p++;
    float f = 0;
    uint8_t n = 0;

    for (; n < mn && *p >= '0' && *p <= '9'; n++, p++) {
      f = f * 10.0 + (*p - '0');
    }

    if (n < 11)
      ret += f / pow10[n];
    else
      ret += f / std::pow(10, n);
  }

  if (neg) return -ret;
  return ret;
}

// _____________________________________________________________________________
inline double atof(const char* p) { return atof(p, 38); }

}  // namespace util

#endif  // UTIL_MISC_H_
