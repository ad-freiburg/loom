// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_GTFS_SERVICE_H_
#define GTFSPARSER_GTFS_SERVICE_H_

#include <map>
#include <ctime>
#include <vector>
#include <string>

using std::exception;
using std::string;


namespace gtfsparser {
namespace gtfs {

class ServiceDate {
 public:
  ServiceDate(uint8_t day, uint8_t month, uint16_t year)
  : _yyyymmdd(year * 10000 + month * 100 + day) {}

  explicit ServiceDate(uint32_t yyyymmdd) : _yyyymmdd(yyyymmdd) {}

  uint32_t getYYYYMMDD() const { return _yyyymmdd;}

  uint16_t getYear() const { return _yyyymmdd / 10000; }
  uint8_t getMonth() const { return (_yyyymmdd - (getYear() * 10000)) / 100; }
  uint8_t getDay() const {
    return _yyyymmdd - (getYear() * 10000) - (getMonth() * 100);
  }

  void setDay(uint8_t day) { _yyyymmdd = _yyyymmdd - getDay() + day; }

  void setMonth(uint8_t month) {
    _yyyymmdd = _yyyymmdd - getMonth() * 100 + month * 100;
  }

  void setYear(uint16_t year) {
    _yyyymmdd = _yyyymmdd - getYear() * 10000 + year * 10000;
  }

  // returns a time struct of this date at 12:00
  tm getTimeStrc() const {
    tm ret = {};
    ret.tm_year = getYear() - 1900;
    ret.tm_mon = getMonth() -1;
    ret.tm_mday = getDay();
    ret.tm_hour = 12;
    mktime(&ret);
    return ret;
  }

 private:
  uint32_t _yyyymmdd;
};

bool operator>(const ServiceDate& lh, const ServiceDate& rh);
bool operator<(const ServiceDate& lh, const ServiceDate& rh);
bool operator==(const ServiceDate& lh, const ServiceDate& rh);
bool operator!=(const ServiceDate& lh, const ServiceDate& rh);
bool operator>=(const ServiceDate& lh, const ServiceDate& rh);
bool operator<=(const ServiceDate& lh, const ServiceDate& rh);
ServiceDate operator+(const ServiceDate& lh, int i);
ServiceDate operator-(const ServiceDate& lh, int i);
ServiceDate operator++(ServiceDate& lh);
ServiceDate operator--(ServiceDate& lh);

class Service {
 public:
  enum SERVICE_DAY : uint8_t {
    NEVER = 0,        // 0000000
    MONDAYS  = 1,     // 0000001
    TUESDAYS = 2,     // 0000010
    WEDNESDAYS = 4,   // 0000100
    THURSDAYS = 8,    // 0001000
    FRIDAYS = 16,     // 0010000
    SATURDAYS = 32,   // 0100000
    SUNDAYS = 64,     // 1000000
    WEEKDAYS = 31,    // 0011111 (shorthand)
    WEEKENDS = 96,    // 1100000 (shorthand)
    ALL_WEEK = 127    // 1111111 (shorthand)
  };

  enum EXCEPTION_TYPE : uint8_t {
    NOT_SET = 0,
    SERVICE_ADDED = 1,
    SERVICE_REMOVED = 2
  };

  explicit Service(const string& id);
  Service(const string& id, uint8_t serviceDays, ServiceDate start,
    ServiceDate end);

  std::string getId() const;
  const std::map<ServiceDate, Service::EXCEPTION_TYPE>& getExceptions() const;

  void addException(const ServiceDate& d, EXCEPTION_TYPE t);

  bool isActiveOn(const ServiceDate& d) const;

  static SERVICE_DAY getServiceDay(const ServiceDate& d);

  EXCEPTION_TYPE getExceptionOn(const ServiceDate& d) const;

 private:
  string _id;
  uint8_t _serviceDays;
  std::map<ServiceDate, Service::EXCEPTION_TYPE> _exceptions;
  ServiceDate _exceptionsBegin, _exceptionsEnd;
};

}  // namespace gtfs
}  // namespace gtfsparser

#endif  // GTFSPARSER_GTFS_SERVICE_H_
