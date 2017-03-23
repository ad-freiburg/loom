// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_STRING_H_
#define UTIL_STRING_H_

#include <vector>
#include <string>
#include <cstring>

namespace pbutil {

// _____________________________________________________________________________
inline std::string urlDecode(const std::string& encoded) {
  std::string decoded;
  for (size_t i = 0; i < encoded.size(); ++i) {
    char c = encoded[i];
    if (c == '%') {
      std::string ah = encoded.substr(i+1, 2);
      char* nonProced = 0;
      char hexVal = strtol(ah.c_str(), &nonProced, 16);

      if (ah.find_first_of("+-") > 1 && ah.size() - strlen(nonProced) == 2) {
        c = hexVal;
        i += 2;
      }
    } else if (c == '+') { c = ' '; }
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
  bool found = false;
  if (from.empty()) return found;
  size_t s = subj.find(from, 0);
  for (; s != std::string::npos; s = subj.find(from, s + to.length())) {
    found = true;
    subj.replace(s, from.length(), to);
  }

  return found;
}

}


#endif  // UTIL_STRING_H_
