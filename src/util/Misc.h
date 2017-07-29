// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_MISC_H_
#define UTIL_MISC_H_

namespace util {

// _____________________________________________________________________________
inline uint64_t factorial(uint64_t n) {
  if (n == 1) return n;
  return n * factorial(n - 1);
}

}  // namespace util


#endif  // UTIL_MISC_H_
