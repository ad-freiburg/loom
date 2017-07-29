// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_STRING_H_
#define UTIL_STRING_H_

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
}

#endif  // UTIL_STRING_H_
