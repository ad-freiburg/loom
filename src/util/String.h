// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_STRING_H_
#define UTIL_STRING_H_

#include <algorithm>
#include <cstring>
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
inline std::string jsonStringEscape(const std::string& unescaped) {
  std::string escaped;
  for (size_t i = 0; i < unescaped.size(); ++i) {
    if (unescaped[i] == '"' || unescaped[i] == '\\') {
      escaped += "\\";
    }
    if (iscntrl(unescaped[i])) {
      escaped += " ";
    }
    escaped += unescaped[i];
  }
  return escaped;
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
template <typename T>
inline std::string implode(const std::vector<T>& vec, const char* del) {
  return implode(vec.begin(), vec.end(), del);
}
}

#endif  // UTIL_STRING_H_
