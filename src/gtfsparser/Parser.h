// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_PARSER_H_
#define GTFSPARSER_PARSER_H_

#include <stdint.h>
#include <boost/filesystem.hpp>
#include <unordered_map>
#include <istream>
#include <string>
#include <vector>
#include <exception>
#include <sstream>
#include <iostream>
#include "gtfs/Feed.h"
#include "CsvParser.h"

using std::string;

// A GTFS parser

namespace gtfsparser {

class ParserException : public std::exception {
 public:
  ParserException(std::string msg,
    std::string field_name, int64_t line, std::string file_name)
    : _msg(msg), _field_name(field_name),
      _line(line), _file_name(file_name) {}
  ParserException(std::string msg,
    std::string field_name, int64_t line)
    : _msg(msg), _field_name(field_name),
      _line(line), _file_name("?") {}
  ~ParserException() throw() {}

  virtual const char* what() const throw() {
    std::stringstream ss;
    ss << _file_name << ":";
    if (_line > -1) ss << _line << ":";
    if (!_field_name.empty()) ss << " in field '" << _field_name << "',";
    ss << " " << _msg;
    _what_msg = ss.str();
    return _what_msg.c_str();
  }

  virtual uint64_t getLine() const throw() {
    return _line;
  }

  void setFileName(const std::string& fn) {
    _file_name = fn;
  }

 private:
  mutable std::string _what_msg;
  std::string _msg;
  uint16_t _column;
  std::string _field_name;
  int64_t _line;
  std::string _file_name;
};

class Parser {
 public:
  // Default initialization.
  Parser() {}

  // parse a zip/folder into a GtfsFeed
  bool parse(gtfs::Feed* targetFeed, std::string path) const;

 private:
  void parseAgency(gtfs::Feed* targetFeed, std::istream*) const;
  void parseStops(gtfs::Feed* targetFeed, std::istream*) const;
  void parseRoutes(gtfs::Feed* targetFeed, std::istream*) const;
  void parseTrips(gtfs::Feed* targetFeed, std::istream*) const;
  void parseStopTimes(gtfs::Feed* targetFeed, std::istream*) const;
  void parseCalendar(gtfs::Feed* targetFeed, std::istream*) const;
  void parseCalendarDates(gtfs::Feed* targetFeed, std::istream*) const;
  void parseFareAttributes(gtfs::Feed* targetFeed, std::istream*) const;
  void parseFareRules(gtfs::Feed* targetFeed, std::istream*) const;
  void parseShapes(gtfs::Feed* targetFeed, std::istream*) const;
  void parseFrequencies(gtfs::Feed* targetFeed, std::istream*) const;
  void parseTransfers(gtfs::Feed* targetFeed, std::istream*) const;
  void parseFeedInfo(gtfs::Feed* targetFeed, std::istream*) const;

  std::string getString(const CsvParser& csv,
    const std::string& field) const;
  std::string getString(const CsvParser& csv,
    const std::string& fld, const std::string& def) const;

  double getDouble(const CsvParser& csv,
    const std::string& field) const;
  double getDouble(const CsvParser& csv,
    const std::string& fld, double def) const;

  int64_t getRangeInteger(const CsvParser& csv,
    const std::string& field, int64_t minv, int64_t maxv) const;
  int64_t getRangeInteger(const CsvParser& csv,
    const std::string& field, int64_t minv, int64_t maxv, int64_t def) const;

  uint32_t getColorFromHexString(const CsvParser& csv,
    const std::string& field, const std::string& def) const;

  gtfs::ServiceDate getServiceDate(const CsvParser& csv,
    const std::string& field) const;

  gtfs::Time getTime(const CsvParser& csv,
    const std::string& field) const;

  gtfs::Route::TYPE getRouteType(const CsvParser& csv,
    const std::string& field, int64_t t) const;
  void fileNotFound(boost::filesystem::path file) const;
};
}  // namespace gtfsparser

#endif  // GTFSPARSER_PARSER_H_
