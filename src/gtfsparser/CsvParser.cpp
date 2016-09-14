// Copyright 2013, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Hannah Bast <bast@informatik.uni-freiburg.de>,
//          Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <functional>
#include <string>
#include <cstring>
#include <algorithm>
#include <iostream>
#include "CsvParser.h"

using gtfsparser::CsvParser;
using std::string;
using std::remove;

// _____________________________________________________________________________
CsvParser::CsvParser() {}

CsvParser::CsvParser(std::istream* stream) : _curLine(0), _stream(stream) {
  readNextLine();
  parseHeader();
}

// _____________________________________________________________________________
bool CsvParser::readNextLine() {
  if (!_stream->good()) return false;
  getline(*_stream, _currentLine);
  _curLine++;
  // Remove new line characters
  _currentLine.erase(remove(_currentLine.begin(), _currentLine.end(), '\r'),
                     _currentLine.end());
  _currentLine.erase(remove(_currentLine.begin(), _currentLine.end(), '\n'),
                     _currentLine.end());
  if (_currentLine.empty()) return readNextLine();  // skip empty lines
  size_t pos = 0;
  size_t lastPos = pos;
  bool firstChar = false;
  bool esc = false;
  int esc_quotes_found = 0;
  _currentItems.clear();
  _currentModItems.clear();

  if (!_stream->eof() || !_currentLine.empty()) {
    if (static_cast<int>(_currentLine[0]) == -17 &&
        static_cast<int>(_currentLine[1]) == -69 &&
        static_cast<int>(_currentLine[2]) == -65) {
      pos = 3;
      lastPos = pos;
    }
    while (pos < _currentLine.size()) {
      if (!firstChar && std::isspace(_currentLine[pos]))  {
          pos++;
          lastPos = pos;
          continue;
      }

      firstChar = true;

      if (_currentLine[pos] == '"') {
        if (!esc) {
          esc = true;
          pos++;
          lastPos = pos;
          continue;
        } else {
          if (pos < _currentLine.size()-1 && _currentLine[pos+1] == '"') {
            pos++;
            esc_quotes_found++;
          } else {
            // we end this field here, because of the closing quotes
            // see CSV spec at http://tools.ietf.org/html/rfc4180#page-2
            _currentLine[pos] = 0;
            esc = false;
          }
        }
      }

      if ((esc || _currentLine[pos] != ',') && pos < _currentLine.size()-1) {
        pos++;
        continue;
      }

      if (_currentLine[pos] == ',') {
        _currentLine[pos] = 0;
        firstChar = false;
      }

      if (!esc_quotes_found) {
        _currentItems.push_back(
          inlineRightTrim(_currentLine.c_str() + lastPos));
      } else {
        size_t p = -1;

        // we have to modify to string (that is, we have to replace
        // characters. create a copy of this item
        // on the line-wise modified vector
        _currentModItems.push_back(_currentLine.c_str() + lastPos);

        while (esc_quotes_found) {
          p = _currentModItems.back().find("\"\"", p+1);
          if (p != string::npos) {
            _currentModItems.back().replace(p, 2, "\"");
          }

          esc_quotes_found--;
        }

        // pointers will point to strings on the modified string vector
        _currentItems.push_back(_currentModItems.back().c_str());
      }

      lastPos = ++pos;
    }
    return true;
  }
  return false;
}

// _____________________________________________________________________________
const char* CsvParser::getTString(const size_t i) const {
  if (i >= _currentItems.size()) return "";
  return _currentItems[i];
}

// _____________________________________________________________________________
double CsvParser::getDouble(const size_t i) const {
  if (i >= _currentItems.size() || !isDouble(_currentItems[i]))
    throw CsvParserException("expected float number", i, getFieldName(i),
                              _curLine);
  return atof(_currentItems[i]);
}

// _____________________________________________________________________________
bool CsvParser::lineIsEmpty(string* line) const {
  strtrim(line);
  return line->empty();
}

// _____________________________________________________________________________
bool CsvParser::lineIsEmpty(const char* line) const {
  size_t i = 0;
  while (line[i]) {
    if (!isspace(line[i])) return false;
    i++;
  }

  return true;
}

// _____________________________________________________________________________
int32_t CsvParser::getLong(const size_t i) const {
  if (i >= _currentItems.size() || !isLong(_currentItems[i]))
    throw CsvParserException("expected integer number", i, getFieldName(i),
                              _curLine);
  return atol(_currentItems[i]);
}

// _____________________________________________________________________________
bool CsvParser::hasItem(const string& fieldName) const {
  return _headerMap.find(fieldName) != _headerMap.end();
}

// _____________________________________________________________________________
bool CsvParser::fieldIsEmpty(const string& fieldName) const {
  return strlen(getTString(fieldName)) == 0;
}

// _____________________________________________________________________________
const char* CsvParser::getTString(const std::string& fieldName) const {
  return getTString(getFieldIndex(fieldName));
}

// _____________________________________________________________________________
double CsvParser::getDouble(const std::string& fieldName) const {
  return getDouble(getFieldIndex(fieldName));
}

// _____________________________________________________________________________
int32_t CsvParser::getLong(const std::string& fieldName) const {
  return getDouble(getFieldIndex(fieldName));
}

// _____________________________________________________________________________
size_t CsvParser::getNumColumns() const {
  return _currentItems.size();
}

// _____________________________________________________________________________
size_t CsvParser::getFieldIndex(const string& fieldName) const {
  if (_headerMap.find(fieldName) == _headerMap.end())
    throw CsvParserException("field " + fieldName + " does not exist.", -1,
                              fieldName, _curLine);
  return _headerMap.find(fieldName)->second;
}

// _____________________________________________________________________________
int32_t CsvParser::getCurLine() const {
  return _curLine;
}

// _____________________________________________________________________________
const string CsvParser::getFieldName(size_t i) const {
  if (i < _headerVec.size()) return _headerVec[i].c_str();
  return "(no field name given)";
}

// _____________________________________________________________________________
void CsvParser::parseHeader() {
  _headerMap.clear();
  for (size_t i = 0; i < getNumColumns(); ++i) {
    string s = getTString(i);
    s.erase(remove_if(s.begin(), s.end(), isspace), s.end());
    _headerMap[s] = i;
    _headerVec.push_back(s);
  }
}

// _____________________________________________________________________________
const char* CsvParser::inlineRightTrim(const char* t) const {
  char* s = const_cast<char*>(t);
  char *end = s + std::strlen(s)-1;
  while (std::isspace(*end)) *end-- = 0;
  return s;
}

// ___________________________________________________________________________
inline bool CsvParser::isDouble(string line, bool notEmpty) const {
  strtrim(&line);
  if (line.size() == 0 && notEmpty) return false;
  return isDouble(line);
}

// ___________________________________________________________________________
inline void CsvParser::ltrim(string* s) const {
  s->erase(s->begin(), std::find_if(s->begin(), s->end(),
  std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// ___________________________________________________________________________
inline void CsvParser::rtrim(string* s) const {
  s->erase(std::find_if(s->rbegin(), s->rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s->end());
}

// ___________________________________________________________________________
inline void CsvParser::strtrim(string* s) const {
  ltrim(s);
  rtrim(s);
}

// ___________________________________________________________________________
inline bool CsvParser::isLong(string line) const {
  strtrim(&line);
  char* p;
  strtol(line.c_str(), &p, 10);
  return *p == 0;
}

// ___________________________________________________________________________
inline bool CsvParser::isLong(string line, bool notEmpty) const {
  strtrim(&line);
  if (line.size() == 0 && notEmpty) return false;
  return isLong(line);
}

// ___________________________________________________________________________
inline bool CsvParser::isDouble(string line) const {
  strtrim(&line);
  char* p;
  strtod(line.c_str(), &p);
  return *p == 0;
}
