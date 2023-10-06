// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_MISC_H_
#define UTIL_MISC_H_

#include <cmath>
#include <cstring>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <map>
#include <thread>
#include "3rdparty/dtoa_milo.h"
#include "util/String.h"

#define UNUSED(expr) do { (void)(expr); } while (0)
#define TIME() std::chrono::high_resolution_clock::now()
#define TOOK(t1, t2) (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.0)
#define T_START(n)  auto _tstart_##n = std::chrono::high_resolution_clock::now()
#define T_STOP(n) (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - _tstart_##n).count() / 1000.0)

#define _TEST3(s, o, e) if (!(s o e)) {  std::cerr << "\n" << __FILE__ << ":" << __LINE__ << ": Test failed!\n  Expected " << #s << " " << #o " " << (e) << ", got " << (s) << std::endl;  exit(1);}
#define _TEST2(s, e) _TEST3(s, ==, e)
#define _TEST1(s) _TEST3(static_cast<bool>(s), ==, true)

#define _GET_TEST_MACRO(_1,_2,_3,NAME,...) NAME
#define TEST(...) _GET_TEST_MACRO(__VA_ARGS__, _TEST3, _TEST2, _TEST1, UNUSED)(__VA_ARGS__)

#define TODO(msg) std::cerr << "\n" __FILE__ << ":" << __LINE__ << ": TODO: " << #msg << std::endl;

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>
#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>
#endif
#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

namespace util {

const static std::map<std::string, std::string> HTML_COLOR_NAMES = {
	{"aliceblue","F0F8FF"},
	{"antiquewhite","FAEBD7"},
	{"aqua","00FFFF"},
	{"aquamarine","7FFFD4"},
	{"azure","F0FFFF"},
	{"beige","F5F5DC"},
	{"bisque","FFE4C4"},
	{"black","000000"},
	{"blanchedalmond","FFEBCD"},
	{"blue","0000FF"},
	{"blueviolet","8A2BE2"},
	{"brown","A52A2A"},
	{"burlywood","DEB887"},
	{"cadetblue","5F9EA0"},
	{"chartreuse","7FFF00"},
	{"chocolate","D2691E"},
	{"coral","FF7F50"},
	{"cornflowerblue","6495ED"},
	{"cornsilk","FFF8DC"},
	{"crimson","DC143C"},
	{"cyan","00FFFF"},
	{"darkblue","00008B"},
	{"darkcyan","008B8B"},
	{"darkgoldenrod","B8860B"},
	{"darkgray","A9A9A9"},
	{"darkgreen","006400"},
	{"darkgrey","A9A9A9"},
	{"darkkhaki","BDB76B"},
	{"darkmagenta","8B008B"},
	{"darkolivegreen","556B2F"},
	{"darkorange","FF8C00"},
	{"darkorchid","9932CC"},
	{"darkred","8B0000"},
	{"darksalmon","E9967A"},
	{"darkseagreen","8FBC8F"},
	{"darkslateblue","483D8B"},
	{"darkslategray","2F4F4F"},
	{"darkslategrey","2F4F4F"},
	{"darkturquoise","00CED1"},
	{"darkviolet","9400D3"},
	{"deeppink","FF1493"},
	{"deepskyblue","00BFFF"},
	{"dimgray","696969"},
	{"dimgrey","696969"},
	{"dodgerblue","1E90FF"},
	{"firebrick","B22222"},
	{"floralwhite","FFFAF0"},
	{"forestgreen","228B22"},
	{"fuchsia","FF00FF"},
	{"gainsboro","DCDCDC"},
	{"ghostwhite","F8F8FF"},
	{"gold","FFD700"},
	{"goldenrod","DAA520"},
	{"gray","808080"},
	{"green","008000"},
	{"greenyellow","ADFF2F"},
	{"grey","808080"},
	{"honeydew","F0FFF0"},
	{"hotpink","FF69B4"},
	{"indianred","CD5C5C"},
	{"indigo","4B0082"},
	{"ivory","FFFFF0"},
	{"khaki","F0E68C"},
	{"lavender","E6E6FA"},
	{"lavenderblush","FFF0F5"},
	{"lawngreen","7CFC00"},
	{"lemonchiffon","FFFACD"},
	{"lightblue","ADD8E6"},
	{"lightcoral","F08080"},
	{"lightcyan","E0FFFF"},
	{"lightgoldenrodyellow","FAFAD2"},
	{"lightgray","D3D3D3"},
	{"lightgreen","90EE90"},
	{"lightgrey","D3D3D3"},
	{"lightpink","FFB6C1"},
	{"lightsalmon","FFA07A"},
	{"lightseagreen","20B2AA"},
	{"lightskyblue","87CEFA"},
	{"lightslategray","778899"},
	{"lightslategrey","778899"},
	{"lightsteelblue","B0C4DE"},
	{"lightyellow","FFFFE0"},
	{"lime","00FF00"},
	{"limegreen","32CD32"},
	{"linen","FAF0E6"},
	{"magenta","FF00FF"},
	{"maroon","800000"},
	{"mediumaquamarine","66CDAA"},
	{"mediumblue","0000CD"},
	{"mediumorchid","BA55D3"},
	{"mediumpurple","9370DB"},
	{"mediumseagreen","3CB371"},
	{"mediumslateblue","7B68EE"},
	{"mediumspringgreen","00FA9A"},
	{"mediumturquoise","48D1CC"},
	{"mediumvioletred","C71585"},
	{"midnightblue","191970"},
	{"mintcream","F5FFFA"},
	{"mistyrose","FFE4E1"},
	{"moccasin","FFE4B5"},
	{"navajowhite","FFDEAD"},
	{"navy","000080"},
	{"oldlace","FDF5E6"},
	{"olive","808000"},
	{"olivedrab","6B8E23"},
	{"orange","FFA500"},
	{"orangered","FF4500"},
	{"orchid","DA70D6"},
	{"palegoldenrod","EEE8AA"},
	{"palegreen","98FB98"},
	{"paleturquoise","AFEEEE"},
	{"palevioletred","DB7093"},
	{"papayawhip","FFEFD5"},
	{"peachpuff","FFDAB9"},
	{"peru","CD853F"},
	{"pink","FFC0CB"},
	{"plum","DDA0DD"},
	{"powderblue","B0E0E6"},
	{"purple","800080"},
	{"red","FF0000"},
	{"rosybrown","BC8F8F"},
	{"royalblue","4169E1"},
	{"saddlebrown","8B4513"},
	{"salmon","FA8072"},
	{"sandybrown","F4A460"},
	{"seagreen","2E8B57"},
	{"seashell","FFF5EE"},
	{"sienna","A0522D"},
	{"silver","C0C0C0"},
	{"skyblue","87CEEB"},
	{"slateblue","6A5ACD"},
	{"slategray","708090"},
	{"slategrey","708090"},
	{"snow","FFFAFA"},
	{"springgreen","00FF7F"},
	{"steelblue","4682B4"},
	{"tan","D2B48C"},
	{"teal","008080"},
	{"thistle","D8BFD8"},
	{"tomato","FF6347"},
	{"turquoise","40E0D0"},
	{"violet","EE82EE"},
	{"wheat","F5DEB3"},
	{"white","FFFFFF"},
	{"whitesmoke","F5F5F5"},
	{"yellow","FFFF00"},
	{"yellowgreen","9ACD32"}
};

struct hashPair {
  template <class T1, class T2>
  size_t operator()(const std::pair<T1, T2>& p) const {
    auto h1 = std::hash<T1>{}(p.first);
    auto h2 = std::hash<T2>{}(p.second);
    return h1 ^  h2;
  }
};

template<typename Key, typename Val, Val Def>
class SparseMatrix {
 public:
  Val get(const Key& x, const Key& y) const {
    auto a = _m.find(std::pair<Key, Key>(x, y));
    if (a == _m.end()) return Def;
    return a->second;
  }

  void set(Key x, Key y, Val v) {
    _m[std::pair<Key, Key>(x, y)] = v;
  }

	const std::map<std::pair<Key, Key>, Val>& vals() const {
    return _m;
  }

 private:
	std::map<std::pair<Key, Key>, Val> _m;
};

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
inline std::string formatFloat(double f, size_t digits) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(digits) << f;
	std::string ret = ss.str();

  if (ret.find('.') != std::string::npos) {
		auto p = ret.find_last_not_of('0');
    if (ret[p] == '.') return ret.substr(0, p);
    return ret.substr(0, p + 1);
  }

  return ret;
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
template <typename V>
int merge(V* lst, V* tmpLst, size_t l, size_t m, size_t r) {
  size_t ret = 0;

  size_t lp = l;
  size_t rp = m;
  size_t outp = l;

  while (lp < m && rp < r + 1) {
    if (lst[lp] <= lst[rp]) {
      // if left element is smaller or equal, add it to return list,
      // increase left pointer
      tmpLst[outp] = lst[lp];
      lp++;
    } else {
      // if left element is bigger, add the right element, add it to ret,
      // increase right pointer
      tmpLst[outp] = lst[rp];
      rp++;

      // if the left element was bigger, everything to the right in the
      // left list is also bigger, and all these m - i elements were
      // initially in the wrong order! Count these inversions.
      ret += m - lp;
    }

    outp++;
  }

  // fill in remaining values
  if (lp < m) std::memcpy(tmpLst + outp, lst + lp, (m - lp) * sizeof(V));
  if (rp <= r) std::memcpy(tmpLst + outp, lst + rp, ((r + 1) - rp) * sizeof(V));

  // copy to output
  std::memcpy(lst + l, tmpLst + l, ((r + 1) - l) * sizeof(V));

  return ret;
}

// _____________________________________________________________________________
template <typename V>
size_t mergeInvCount(V* lst, V* tmpLst, size_t l, size_t r) {
  size_t ret = 0;
  if (l < r) {
    size_t m = (r + l) / 2;

    ret += mergeInvCount(lst, tmpLst, l, m);
    ret += mergeInvCount(lst, tmpLst, m + 1, r);

    ret += merge(lst, tmpLst, l, m + 1, r);
  }
  return ret;
}

// _____________________________________________________________________________
template <typename V>
size_t inversions(const std::vector<V>& v) {
  if (v.size() < 2) return 0;  // no inversions possible

  // unroll some simple cases
  if (v.size() == 2) return v[1] < v[0];
  if (v.size() == 3) return (v[0] > v[1]) + (v[0] > v[2]) + (v[1] > v[2]);

  auto tmpLst = new V[v.size()];
  auto lst = new V[v.size()];

  for (size_t i = 0; i < v.size(); i++) lst[i] = v[i];

  size_t ret = mergeInvCount<V>(lst, tmpLst, 0, v.size() - 1);
  delete[] tmpLst;
  delete[] lst;
  return ret;
}


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
inline std::string getTmpFName(std::string dir, std::string name,
    std::string postf) {
  if (postf.size()) postf = "-" + postf;
  if (dir == "<tmp>") dir = util::getTmpDir();
  if (dir.size() && dir.back() != '/') dir = dir + "/";

  std::string f = dir + name + postf;

  size_t c = 0;

  while (access(f.c_str(), F_OK) != -1) {
    c++;
    if (c > 10000) {
      // giving up...
      std::cerr << "Could not find temporary file name!" << std::endl;
      exit(1);
    }
    std::stringstream ss;
    ss << dir << name << postf << "-" << std::rand();
    f = ss.str().c_str();
  }

  return f;
}

// ___________________________________________________________________________
inline std::string rgbToHex(int r, int g, int b) {
  char hexcol[16];
  snprintf(hexcol, sizeof hexcol, "%02x%02x%02x", r, g, b);
  return hexcol;
}

// ___________________________________________________________________________
inline void hsvToRgb(float* r, float* g, float* b, float h, float s, float v) {
  int i;
  float f, p, q, t;

  if (s == 0) {
    *r = *g = *b = v;
    return;
  }

  h /= 60;
  i = floor(h);
  f = h - i;
  p = v * (1 - s);
  q = v * (1 - s * f);
  t = v * (1 - s * (1 - f));

  switch (i) {
    case 0:
      *r = v;
      *g = t;
      *b = p;
      break;
    case 1:
      *r = q;
      *g = v;
      *b = p;
      break;
    case 2:
      *r = p;
      *g = v;
      *b = t;
      break;
    case 3:
      *r = p;
      *g = q;
      *b = v;
      break;
    case 4:
      *r = t;
      *g = p;
      *b = v;
      break;
    default:
      *r = v;
      *g = p;
      *b = q;
      break;
  }
}

// ___________________________________________________________________________
inline std::string normHtmlColor(const std::string& col) {
	auto i = HTML_COLOR_NAMES.find(toLower(col));
  if (i != HTML_COLOR_NAMES.end()) return i->second;
  return col;
}

// ___________________________________________________________________________
inline std::string randomHtmlColor() {
  double goldenRatio = 0.618033988749895;
  double h = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  h += goldenRatio;
  h = fmod(h, 1.0);
  float r, g, b;
  hsvToRgb(&r, &g, &b, h * 360, 0.95, 0.95);
  return rgbToHex(r * 256, g * 256, b * 256);
}

// _____________________________________________________________________________
inline char* readableSize(double size, size_t n, char* buf) {
  int i = 0;
  const char* units[] = {"B", "kB", "MB", "GB", "TB", "PB"};
  while (size > 1024 && i < 5) {
    size /= 1024;
    i++;
  }
  snprintf(buf, n, "%.*f %s", i, size, units[i]);
  return buf;
}

// _____________________________________________________________________________
inline std::string readableSize(double size) {
  char buffer[30];
  return readableSize(size, 30, buffer);
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

/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */

// _____________________________________________________________________________
inline size_t getPeakRSS() {
#if defined(_WIN32)
  /* Windows -------------------------------------------------- */
  PROCESS_MEMORY_COUNTERS info;
  GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
  return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || \
    (defined(__sun__) || defined(__sun) ||        \
     defined(sun) && (defined(__SVR4) || defined(__svr4__)))
  /* AIX and Solaris ------------------------------------------ */
  struct psinfo psinfo;
  int fd = -1;
  if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
    return (size_t)0L; /* Can't open? */
  if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo)) {
    close(fd);
    return (size_t)0L; /* Can't read? */
  }
  close(fd);
  return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || \
    (defined(__APPLE__) && defined(__MACH__))
  /* BSD, Linux, and OSX -------------------------------------- */
  struct rusage rusage;
  getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
  return (size_t)rusage.ru_maxrss;
#else
  return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
  /* Unknown OS ----------------------------------------------- */
  return (size_t)0L; /* Unsupported. */
#endif
}

/*
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */

// _____________________________________________________________________________
inline size_t getCurrentRSS() {
#if defined(_WIN32)
  /* Windows -------------------------------------------------- */
  PROCESS_MEMORY_COUNTERS info;
  GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
  return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
  /* OSX ------------------------------------------------------ */
  struct mach_task_basic_info info;
  mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info,
                &infoCount) != KERN_SUCCESS)
    return (size_t)0L; /* Can't access? */
  return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || \
    defined(__gnu_linux__)
  /* Linux ---------------------------------------------------- */
  long rss = 0L;
  FILE* fp = NULL;
  if ((fp = fopen("/proc/self/statm", "r")) == NULL)
    return (size_t)0L; /* Can't open? */
  if (fscanf(fp, "%*s%ld", &rss) != 1) {
    fclose(fp);
    return (size_t)0L; /* Can't read? */
  }
  fclose(fp);
  return (size_t)rss * (size_t)sysconf(_SC_PAGESIZE);

#else
  /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
  return (size_t)0L; /* Unsupported. */
#endif
}

}  // namespace util

#endif  // UTIL_MISC_H_
