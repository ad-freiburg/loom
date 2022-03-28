// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_STRING_H_
#define UTIL_STRING_H_

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <codecvt>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <locale>
#include <set>
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
template <class String>
inline size_t prefixEditDist(const String& prefix, const String& s,
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
template <class String>
inline size_t prefixEditDist(const String& prefix, const String& s) {
  return prefixEditDist(prefix, s, s.size());
}

// _____________________________________________________________________________
inline size_t prefixEditDist(const char* prefix, const char* s) {
  return prefixEditDist<std::string>(std::string(prefix), std::string(s));
}

// _____________________________________________________________________________
inline size_t prefixEditDist(const char* prefix, const char* s,
                             size_t deltaMax) {
  return prefixEditDist<std::string>(std::string(prefix), std::string(s),
                                     deltaMax);
}

// _____________________________________________________________________________
inline size_t prefixEditDist(const char* prefix, const std::string& s) {
  return prefixEditDist<std::string>(std::string(prefix), s);
}

// _____________________________________________________________________________
inline size_t prefixEditDist(const char* prefix, const std::string& s,
                             size_t deltaMax) {
  return prefixEditDist<std::string>(std::string(prefix), s, deltaMax);
}

// _____________________________________________________________________________
inline size_t prefixEditDist(const std::string& prefix, const char* s) {
  return prefixEditDist<std::string>(prefix, std::string(s));
}

// _____________________________________________________________________________
inline size_t prefixEditDist(const std::string& prefix, const char* s,
                             size_t deltaMax) {
  return prefixEditDist<std::string>(prefix, std::string(s), deltaMax);
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
template <class T>
inline std::string implode(const std::vector<T>& vec, const char* del) {
  return implode(vec.begin(), vec.end(), del);
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
inline std::wstring toWStr(const std::string& str) {
  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
  return converter.from_bytes(str);
}

// _____________________________________________________________________________
inline std::string toNStr(const std::wstring& wstr) {
  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
  return converter.to_bytes(wstr);
}

// _____________________________________________________________________________
inline std::vector<std::string> tokenize(const std::string& str) {
  std::vector<std::string> ret;
  std::wstring wStr = toWStr(str);

  std::wstring cur;
  for (size_t i = 0; i < wStr.size(); ++i) {
    if (!std::iswalnum(wStr[i])) {
      if (cur.size()) ret.push_back(toNStr(cur));
      cur = L"";
      continue;
    }
    cur += wStr[i];
  }
  if (cur.size()) ret.push_back(toNStr(cur));

  return ret;
}

// _____________________________________________________________________________
inline double jaccardSimi(const std::string& a, const std::string& b) {
  if (a == b) return 1;

  std::set<std::string> sa, sb;

  auto toksA = tokenize(a);
  auto toksB = tokenize(b);

  // 0 if both are empty
  if (toksA.size() == 0 && toksB.size() == 0) return 0;

  sa.insert(toksA.begin(), toksA.end());
  sb.insert(toksB.begin(), toksB.end());

  std::set<std::string> isect;

  std::set_intersection(sa.begin(), sa.end(), sb.begin(), sb.end(),
                        std::inserter(isect, isect.begin()));

  double sInter = isect.size();
  double s1 = sa.size();
  double s2 = sb.size();
  return sInter / (s1 + s2 - sInter);
}

// _____________________________________________________________________________
inline double btsSimiInner(const std::vector<std::string>& toks,
                           const std::string& b, double best) {
  std::set<std::string> toksSet;
  toksSet.insert(toks.begin(), toks.end());
  std::vector<std::string> toksUniqSorted;
  toksUniqSorted.insert(toksUniqSorted.begin(), toksSet.begin(), toksSet.end());

  assert(toksUniqSorted.size() <= 8);

  for (uint8_t v = 1; v <= pow(2, toksUniqSorted.size()); v++) {
    std::bitset<8> bs(v);
    std::vector<std::string> cur(bs.count());

    size_t i = 0;
    for (size_t j = 0; j < toksUniqSorted.size(); j++) {
      if (bs[j]) {
        cur[i] = toksUniqSorted[j];
        i++;
      }
    }

    double tmp = util::implode(cur, " ").size();

    // ed between the two string will always be at least their length
    // difference - if this is already too big, skip it right now
    double dt = 1 - (fabs(tmp - b.size()) * 1.0) / (fmax(tmp, b.size()) * 1.0);

    if (dt <= best) continue;

    // cur is guaranteed to be sorted now
    do {
      const auto& comb = util::implode(cur, " ");

      double d =
          1 - ((editDist(comb, b) * 1.0) / (fmax(comb.size(), b.size()) * 1.0));

      if (fabs(d - 1) < 0.0001) return 1;

      if (d > best) best = d;
    } while (std::next_permutation(cur.begin(), cur.end()));
  }

  return best;
}

// _____________________________________________________________________________
inline double btsSimi(std::string a, std::string b) {
  // this follows the implementation for the station similarity paper in
  // https://github.com/ad-freiburg/statsimi/
  if (a == b) return 1;

  std::set<std::string> sa, sb;

  auto toksA = tokenize(a);
  auto toksB = tokenize(b);

  // fallback to jaccard if the token set is too large
  if (toksA.size() > 6 || toksB.size() > 6) {
    return jaccardSimi(a, b);
  }

  if (toksA.size() > toksB.size()) {
    std::swap(a, b);
    std::swap(toksA, toksB);
  }

  // this is already our best known value - simply the edit
  // distance similarity between the strings
  double best = 1 - (editDist(a, b) * 1.0) / std::fmax(a.size(), b.size());

  if (fabs(best) < 0.0001) return 0;

  best = btsSimiInner(toksA, b, best);

  if (fabs(best - 1) < 0.0001) return 1;

  return btsSimiInner(toksB, a, best);
}
}  // namespace util

#endif  // UTIL_STRING_H_
