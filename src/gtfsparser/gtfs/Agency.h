// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_GTFS_AGENCY_H_
#define GTFSPARSER_GTFS_AGENCY_H_

#include <string>

using std::exception;
using std::string;

namespace gtfsparser {
namespace gtfs {

class Agency {
 public:
  Agency() {}

  Agency(string id, string name, string url, string timezone, string lang,
      string phone, string fare_url, string agency_email)
  : _id(id), _name(name), _url(url), _timezone(timezone), _lang(lang),
      _phone(phone), _fare_url(fare_url), _agency_email(agency_email) {}

  Agency(const char* id, const char* name, const char* url,
      const char* timezone, const char* lang, const char* phone,
      const char* fare_url, const char* agency_email)
  : _id(id), _name(name), _url(url), _timezone(timezone), _lang(lang),
      _phone(phone), _fare_url(fare_url), _agency_email(agency_email) {}

  std::string getId() const {
    return _id;
  }

  std::string getName() const {
    return _name;
  }

  std::string getUrl() const {
    return _url;
  }

  std::string getTimezone() const {
    return _timezone;
  }

  std::string getLang() const {
    return _lang;
  }

  std::string getPhone() const {
    return _phone;
  }

  std::string getFareUrl() const {
    return _fare_url;
  }

  std::string getAgencyEmail() const {
    return _agency_email;
  }

  // TODO(patrick): implement setters

 private:
  std::string _id, _name, _url, _timezone, _lang, _phone, _fare_url,
      _agency_email;
};

}  // namespace gtfs
}  // namespace gtfsparser

#endif  // GTFSPARSER_GTFS_AGENCY_H_
