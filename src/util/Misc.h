// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_MISC_H_
#define UTIL_MISC_H_

#include <cmath>
#include <cstring>
#include <chrono>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

#define UNUSED(expr) do { (void)(expr); } while (0)
#define TIME() std::chrono::high_resolution_clock::now()
#define TOOK(t1, t2) (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.0)
#define T_START(n)  auto _tstart_##n = std::chrono::high_resolution_clock::now()
#define T_STOP(n) (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - _tstart_##n).count() / 1000.0)

#define _TEST3(s, o, e) if (!(s o e)) {  std::cerr << "\n" << __FILE__ << ":" << __LINE__ << ": Test failed!\n  Expected " << #s << " " << #o " " << (e) << ", got " << (s) << std::endl;  exit(1);}
#define _TEST2(s, e) _TEST3(s, ==, o)
#define _TEST1(s) _TEST3(static_cast<bool>(s), ==, true)

#define _GET_TEST_MACRO(_1,_2,_3,NAME,...) NAME
#define TEST(...) _GET_TEST_MACRO(__VA_ARGS__, _TEST3, _TEST2, _TEST1)(__VA_ARGS__)

#define TODO(msg) std::cerr << "\n" __FILE__ << ":" << __LINE__ << ": TODO: " << #msg << std::endl;

namespace util {

// cached first 10 powers of 10
static int pow10[10] = {
    1,           10,          100,           1000,          10000,
    100000,      1000000,     10000000,      100000000,     1000000000};

// _____________________________________________________________________________
inline uint64_t factorial(uint64_t n) {
  if (n < 2) return 1;
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
inline bool isFloatingPoint(const std::string& str) {
  std::stringstream ss(str);
  double f;
  ss >> std::noskipws >> f;
  return ss.eof() && ! ss.fail();
}

// _____________________________________________________________________________
inline double atof(const char* p, uint8_t mn) {
  // this atof implementation works only on "normal" float strings like
  // 56.445 or -345.00, but should be faster than std::atof
  double ret = 0.0;
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
    double f = 0;
    uint8_t n = 0;

    for (; n < mn && *p >= '0' && *p <= '9'; n++, p++) {
      f = f * 10.0 + (*p - '0');
    }

    if (n < 10)
      ret += f / pow10[n];
    else
      ret += f / std::pow(10, n);
  }

  if (neg) return -ret;
  return ret;
}

// _____________________________________________________________________________
inline double atof(const char* p) { return atof(p, 38); }

// _____________________________________________________________________________
inline std::string getHomeDir() {
  // parse implicit paths
  const char* homedir = 0;
  char* buf = 0;

  if ((homedir = getenv("HOME")) == 0) {
    homedir = "";
    struct passwd pwd;
    struct passwd* result;
    size_t bufsize;
    bufsize = sysconf(_SC_GETPW_R_SIZE_MAX);
    if (bufsize == static_cast<size_t>(-1)) bufsize = 0x4000;
    buf = static_cast<char*>(malloc(bufsize));
    if (buf != 0) {
      getpwuid_r(getuid(), &pwd, buf, bufsize, &result);
      if (result != NULL) homedir = result->pw_dir;
    }
  }

  std::string ret(homedir);
  if (buf) free(buf);

  return ret;
}

// _____________________________________________________________________________
inline std::string getTmpDir() {
  // first, check if an env variable is set
  const char* tmpdir = getenv("TMPDIR");
  if (tmpdir && std::strlen(tmpdir)) return std::string(tmpdir);

  // second, check if /tmp is writable
  if (access("/tmp/", W_OK) == 0) return "/tmp";

  // third, check if the cwd is writable
  if (access(".", W_OK) == 0)  return ".";

  // lastly, return the users home directory as a fallback
  return getHomeDir();
}

// _____________________________________________________________________________
class approx {
 public:
  explicit approx(double magnitude)
      : _epsilon{std::numeric_limits<float>::epsilon() * 100},
        _magnitude{magnitude} {}

  friend bool operator==(double lhs, approx const& rhs) {
    return std::abs(lhs - rhs._magnitude) < rhs._epsilon;
  }

  friend bool operator==(approx const& lhs, double rhs) {
    return operator==(rhs, lhs);
  }

  friend bool operator!=(double lhs, approx const& rhs) {
    return !operator==(lhs, rhs);
  }

  friend bool operator!=(approx const& lhs, double rhs) {
    return !operator==(rhs, lhs);
  }

  friend bool operator<=(double lhs, approx const& rhs) {
    return lhs < rhs._magnitude || lhs == rhs;
  }

  friend bool operator<=(approx const& lhs, double rhs) {
    return lhs._magnitude < rhs || lhs == rhs;
  }

  friend bool operator>=(double lhs, approx const& rhs) {
    return lhs > rhs._magnitude || lhs == rhs;
  }

  friend bool operator>=(approx const& lhs, double rhs) {
    return lhs._magnitude > rhs || lhs == rhs;
  }

  friend std::ostream& operator<< (std::ostream &out, const approx &a) {
    out << "~" << a._magnitude;
    return out;
  }

 private:
  double _epsilon;
  double _magnitude;
};

}  // namespace util

#endif  // UTIL_MISC_H_
