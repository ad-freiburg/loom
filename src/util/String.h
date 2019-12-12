// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_STRING_H_
#define UTIL_STRING_H_

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace util {

// _____________________________________________________________________________
inline std::string urlDecode(const std::string& encoded) {
  std::string decoded;
  for (size_t i = 0; i < encoded.size(); ++i) {
    char c = encoded[i];
    if (c == '%') {
      std::string ah = encoded.substr(i + 1, 2);
      char* nonProced = 0;
      char hexVal = strtol(ah.c_str(), &nonProced, 16);

      if (ah.find_first_of("+-") > 1 && ah.size() - strlen(nonProced) == 2) {
        c = hexVal;
        i += 2;
      }
    } else if (c == '+') {
      c = ' ';
    }
    decoded += c;
  }
  return decoded;
}

// _____________________________________________________________________________
inline std::string jsonStringEscape(const std::string& unesc) {
  // modified code from
  // http://stackoverflow.com/questions/7724448/simple-json-string-escape-for-c
  std::ostringstream o;
  for (auto c = unesc.cbegin(); c != unesc.cend(); c++) {
    switch (*c) {
      case '"':
        o << "\\\"";
        break;
      case '\\':
        o << "\\\\";
        break;
      case '\b':
        o << "\\b";
        break;
      case '\f':
        o << "\\f";
        break;
      case '\n':
        o << "\\n";
        break;
      case '\r':
        o << "\\r";
        break;
      case '\t':
        o << "\\t";
        break;
      default:
        if ('\x00' <= *c && *c <= '\x1f') {
          o << "\\u" << std::hex << std::setw(4) << std::setfill('0')
            << static_cast<int>(*c);
        } else {
          o << *c;
        }
    }
  }
  return o.str();
}

// _____________________________________________________________________________
inline bool replace(std::string& subj, const std::string& from,
                    const std::string& to) {
  if (from.empty()) return false;
  size_t start_pos = subj.find(from);
  if (start_pos != std::string::npos) {
    subj.replace(start_pos, from.length(), to);
    return true;
  }

  return false;
}

// _____________________________________________________________________________
inline bool replaceAll(std::string& subj, const std::string& from,
                       const std::string& to) {
  if (from.empty()) return false;
  bool found = false;
  size_t s = subj.find(from, 0);
  for (; s != std::string::npos; s = subj.find(from, s + to.length())) {
    found = true;
    subj.replace(s, from.length(), to);
  }

  return found;
}

// _____________________________________________________________________________
inline std::string unixBasename(const std::string& pathname) {
  return {std::find_if(pathname.rbegin(), pathname.rend(),
                       [](char c) { return c == '/'; })
              .base(),
          pathname.end()};
}

// _____________________________________________________________________________
template <typename T>
inline std::string toString(T obj) {
  std::stringstream ss;
  ss << obj;
  return ss.str();
}

// _____________________________________________________________________________
inline std::vector<std::string> split(std::string in, char sep) {
  std::stringstream ss(in);
  std::vector<std::string> ret(1);
  while (std::getline(ss, ret.back(), sep)) {
    ret.push_back("");
  }
  ret.pop_back();
  return ret;
}

// _____________________________________________________________________________
inline std::string ltrim(std::string str) {
  str.erase(0, str.find_first_not_of(" \t\n\v\f\r"));
  return str;
}

// _____________________________________________________________________________
inline std::string rtrim(std::string str) {
  str.erase(str.find_last_not_of(" \t\n\v\f\r") + 1);
  return str;
}

// _____________________________________________________________________________
inline std::string trim(std::string str) { return ltrim(rtrim(str)); }

// _____________________________________________________________________________
inline size_t editDist(const std::string& s1, const std::string& s2) {
  // https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
  size_t len1 = s1.size();
  size_t len2 = s2.size();
  std::vector<size_t> cur(len2 + 1);
  std::vector<size_t> prev(len2 + 1);

  for (size_t i = 0; i < prev.size(); i++) prev[i] = i;

  for (size_t i = 0; i < len1; i++) {
    cur[0] = i + 1;
    for (size_t j = 0; j < len2; j++) {
      cur[j + 1] =
          std::min(prev[1 + j] + 1,
                   std::min(cur[j] + 1, prev[j] + (s1[i] == s2[j] ? 0 : 1)));
    }
    std::swap(cur, prev);
  }

  return prev[len2];
}

// _____________________________________________________________________________
inline size_t prefixEditDist(const std::string& prefix, const std::string& s,
                             size_t deltaMax) {
  // https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
  size_t len1 = prefix.size();
  size_t len2 = std::min(s.size(), prefix.size() + deltaMax + 1);
  std::vector<size_t> d((len1 + 1) * (len2 + 1));

  d[0] = 0;
  for (size_t i = 1; i <= len1; ++i) d[i * (len2 + 1)] = i;
  for (size_t i = 1; i <= len2; ++i) d[i] = i;

  for (size_t i = 1; i <= len1; i++) {
    for (size_t j = 1; j <= len2; j++) {
      d[i * (len2 + 1) + j] = std::min(std::min(d[(i - 1) * (len2 + 1) + j] + 1,
                                                d[i * (len2 + 1) + j - 1] + 1),
                                       d[(i - 1) * (len2 + 1) + j - 1] +
                                           (prefix[i - 1] == s[j - 1] ? 0 : 1));
    }
  }

  // take min of last row
  size_t deltaMin = std::max(std::max(deltaMax + 1, prefix.size()), s.size());
  for (size_t i = 0; i <= len2; i++) {
    if (d[len1 * (len2 + 1) + i] < deltaMin)
      deltaMin = d[len1 * (len2 + 1) + i];
  }

  return deltaMin;
}

// _____________________________________________________________________________
inline size_t prefixEditDist(const std::string& prefix, const std::string& s) {
  return prefixEditDist(prefix, s, s.size());
}

// _____________________________________________________________________________
inline std::string toUpper(std::string str) {
  std::transform(str.begin(), str.end(), str.begin(), toupper);
  return str;
}

// _____________________________________________________________________________
inline std::string toLower(std::string str) {
  std::transform(str.begin(), str.end(), str.begin(), tolower);
  return str;
}

// _____________________________________________________________________________
template <class Iter>
inline std::string implode(Iter begin, const Iter& end, const char* del) {
  std::stringstream ss;
  size_t i = 0;
  while (begin != end) {
    if (i != 0) ss << del;
    ss << *begin;
    begin++;
    i++;
  }

  return ss.str();
}

// _____________________________________________________________________________
inline std::string normalizeWhiteSpace(const std::string& input) {
  std::string ret;
  bool ws = false;
  for (size_t i = 0; i < input.size(); i++) {
    if (std::isspace(input[i])) {
      if (!ws) {
        ret += " ";
        ws = true;
      }
      continue;
    } else {
      ws = false;
      ret += input[i];
    }
  }
  return ret;
}

// _____________________________________________________________________________
template <typename T>
inline std::string implode(const std::vector<T>& vec, const char* del) {
  return implode(vec.begin(), vec.end(), del);
}
}

#endif  // UTIL_STRING_H_
