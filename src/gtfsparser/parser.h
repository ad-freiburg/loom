// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_PARSER_H_
#define GTFSPARSER_PARSER_H_

#include <unordered_map>
#include <stdint.h>
#include <istream>
#include <string>
#include <vector>
#include <exception>
#include <sstream>
#include <iostream>
#include <boost/filesystem.hpp>
#include "gtfs/feed.h"
#include "csvparser.h"

using std::string;
using namespace boost;

// A GTFS parser

namespace gtfsparser {

class ParserException : public std::exception {
 public:
  ParserException(std::string msg,
    std::string field_name, uint64_t line, std::string file_name)
    : _msg(msg), _field_name(field_name),
      _line(line), _file_name(file_name) {}
  ~ParserException() throw() {}

  virtual const char* what() const throw() {
    std::stringstream ss;
    ss << _msg
      << " for field '" << _field_name << "'"
      << " in "
      << _file_name << ":" << getLine();
    _what_msg = ss.str();
    return _what_msg.c_str();
  };

  virtual uint64_t getLine() const throw() {
    return _line;
  };

 private:
  mutable std::string _what_msg;
  std::string _msg;
  uint16_t _column;
  std::string _field_name;
  uint64_t _line;
  std::string _file_name;
};

class Parser {
 public:
  // Default initialization.
  Parser() {};

  // parse a zip/folder into a GtfsFeed
  bool parse(gtfs::Feed* targetFeed, std::string path) const;

 private:
  void parseAgency(gtfs::Feed* targetFeed, std::string path) const;
  void parseStops(gtfs::Feed* targetFeed, std::string path) const;
  void parseRoutes(gtfs::Feed* targetFeed, std::string path) const;
  void parseTrips(gtfs::Feed* targetFeed, std::string path) const;
  void parseStopTimes(gtfs::Feed* targetFeed, std::string path) const;
  void parseCalendar(gtfs::Feed* targetFeed, std::string path) const;
  void parseCalendarDates(gtfs::Feed* targetFeed, std::string path) const;
  void parseFareAttributes(gtfs::Feed* targetFeed, std::string path) const;
  void parseFareRules(gtfs::Feed* targetFeed, std::string path) const;
  void parseShapes(gtfs::Feed* targetFeed, std::string path) const;
  void parseFrequencies(gtfs::Feed* targetFeed, std::string path) const;
  void parseTransfers(gtfs::Feed* targetFeed, std::string path) const;
  void parseFeedInfo(gtfs::Feed* targetFeed, std::string path) const;


  std::string getString(const filesystem::path& fn, const CsvParser& csv,
    const std::string& field) const;
  std::string getString(const filesystem::path& fn, const CsvParser& csv,
    const std::string& fld, const std::string& def) const;

  double getDouble(const filesystem::path& fn, const CsvParser& csv,
    const std::string& field) const;
  double getDouble(const filesystem::path& fn, const CsvParser& csv,
    const std::string& fld, double def) const;

  int64_t getRangeInteger(const filesystem::path& fn, const CsvParser& csv,
    const std::string& field, int64_t minv, int64_t maxv) const;
  int64_t getRangeInteger(const filesystem::path& fn, const CsvParser& csv,
    const std::string& field, int64_t minv, int64_t maxv, int64_t def) const;

  uint32_t getColorFromHexString(const filesystem::path& fn, const CsvParser& csv,
    const std::string& field, const std::string& def) const;


  void fileNotFound(boost::filesystem::path file) const;
};
}
#endif  // GTFSPARSER_PARSER_H_
