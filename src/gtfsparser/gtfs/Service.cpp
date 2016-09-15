// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <map>
#include <string>
#include <iostream>
#include "Service.h"

using gtfsparser::gtfs::Service;
using gtfsparser::gtfs::ServiceDate;

// _____________________________________________________________________________
Service::Service(const std::string& id)
: _id(id),
  _serviceDays(Service::SERVICE_DAY::NEVER),
  _exceptionsBegin(0),
  _exceptionsEnd(0)
{}

// _____________________________________________________________________________
Service::Service(const std::string& id, uint8_t serviceDays, ServiceDate start,
  ServiceDate end)
: _id(id), _serviceDays(serviceDays),
  _exceptionsBegin(start), _exceptionsEnd(end) {}

// _____________________________________________________________________________
std::string Service::getId() const {
  return _id;
}

// _____________________________________________________________________________
const std::map<ServiceDate, Service::EXCEPTION_TYPE>& Service::getExceptions()
const {
  return _exceptions;
}

// _____________________________________________________________________________
void Service::addException(const ServiceDate& d, Service::EXCEPTION_TYPE t) {
  _exceptions[d] = t;
}

// _____________________________________________________________________________
bool Service::isActiveOn(const ServiceDate& d) const {
  return ((d >= _exceptionsBegin && d <= _exceptionsEnd)
    && (_serviceDays & getServiceDay(d))
    && getExceptionOn(d) != EXCEPTION_TYPE::SERVICE_REMOVED)
    || getExceptionOn(d) == EXCEPTION_TYPE::SERVICE_ADDED;
}

// _____________________________________________________________________________
Service::SERVICE_DAY Service::getServiceDay(const ServiceDate& d) {
  return static_cast<SERVICE_DAY>(
    1 << (((d.getTimeStrc().tm_wday + 6) % 7)));
}

// _____________________________________________________________________________
Service::EXCEPTION_TYPE Service::getExceptionOn(const ServiceDate& d) const {
  auto ex = _exceptions.find(d);
  if (ex != _exceptions.end()) return ex->second;
  return EXCEPTION_TYPE::NOT_SET;
}

// _____________________________________________________________________________
bool gtfsparser::gtfs::operator>(const ServiceDate& lh, const ServiceDate& rh) {
  return lh.getYYYYMMDD() > rh.getYYYYMMDD();
}

// _____________________________________________________________________________
bool gtfsparser::gtfs::operator<(const ServiceDate& lh, const ServiceDate& rh) {
  return rh > lh;
}

// _____________________________________________________________________________
bool gtfsparser::gtfs::operator==(const ServiceDate& lh,
    const ServiceDate& rh) {
  return !(lh > rh) && !(lh < rh);
}

// _____________________________________________________________________________
bool gtfsparser::gtfs::operator!=(const ServiceDate& lh,
    const ServiceDate& rh) {
  return !(lh == rh);
}

// _____________________________________________________________________________
bool gtfsparser::gtfs::operator>=(const ServiceDate& lh,
    const ServiceDate& rh) {
  return lh > rh || lh == rh;
}

// _____________________________________________________________________________
bool gtfsparser::gtfs::operator<=(const ServiceDate& lh,
    const ServiceDate& rh) {
  return rh > lh || lh == rh;
}

// _____________________________________________________________________________
ServiceDate gtfsparser::gtfs::operator+(const ServiceDate& lh, int i) {
  tm tStrc = lh.getTimeStrc();
  tStrc.tm_mday += i;
  mktime(&tStrc);

  return ServiceDate(tStrc.tm_mday, tStrc.tm_mon + 1, tStrc.tm_year + 1900);
}

// _____________________________________________________________________________
ServiceDate gtfsparser::gtfs::operator-(const ServiceDate& lh, int i) {
  return lh + (-i);
}

// _____________________________________________________________________________
ServiceDate gtfsparser::gtfs::operator--(ServiceDate& lh) {
  ServiceDate ret = lh - 1;
  lh.setDay(ret.getDay());
  lh.setMonth(ret.getMonth());
  lh.setYear(ret.getYear());
  return ret;
}

// _____________________________________________________________________________
ServiceDate gtfsparser::gtfs::operator++(ServiceDate& lh) {
  ServiceDate ret = lh + 1;
  lh.setDay(ret.getDay());
  lh.setMonth(ret.getMonth());
  lh.setYear(ret.getYear());
  return ret;
}
