// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_GTFS_STOP_H_
#define GTFSPARSER_GTFS_STOP_H_

#include <stdint.h>
#include <string>

using std::exception;
using std::string;


namespace gtfsparser {
namespace gtfs {

class Stop {

 public:
  enum LOCATION_TYPE : uint8_t {
    STOP_INSIDE_STATION = 0,
    STOP_OUTSIDE_STATION = 1,
    STATION = 2
  };

  enum WHEELCHAIR_BOARDING : uint8_t {
    NO_INFORMATION = 0,
    BOARDING_POSSIBLE = 1,
    BOARDING_NOT_POSSIBLE = 2
  };

  Stop() {};

  Stop(const string& id, const string& code, const string& name,
    const string& desc, double lat,
    double lng, string zone_id, const string& stop_url,
    Stop::LOCATION_TYPE location_type,
    Stop* parent_station, const string& stop_timezone,
    Stop::WHEELCHAIR_BOARDING wheelchair_boarding)
  : _id(id), _code(code), _name(name), _desc(desc), _zone_id(zone_id),
    _stop_url(stop_url),
    _stop_timezone(stop_timezone),
    _parent_station(parent_station),
    _lat(lat), _lng(lng), _wheelchair_boarding(wheelchair_boarding),
    _location_type(location_type) {}

  std::string getId() const {
    return _id;
  }

  std::string getCode() const {
    return _code;
  }

  std::string getName() const {
    return _name;
  }

  std::string getDesc() const {
    return _desc;
  }

  double getLat() const {
    return _lat;
  }

  double getLng() const {
    return _lng;
  }

  std::string getZoneId() const {
    return _zone_id;
  }

  std::string getStopUrl() const {
    return _stop_url;
  }

  Stop::LOCATION_TYPE getLocationType() const {
    return _location_type;
  }

  Stop* getParentStation() const {
    return _parent_station;
  }

  void setParentStation(Stop* p) {
    _parent_station = p;
  }

  std::string getStopTimezone() const {
    return _stop_timezone;
  }

  Stop::WHEELCHAIR_BOARDING getWheelchairBoarding() const {
    return _wheelchair_boarding;
  }

  // TODO: implement setters


 private:
  string _id, _code, _name, _desc, _zone_id, _stop_url,
    _stop_timezone;
  Stop* _parent_station;
  double _lat, _lng;
  Stop::WHEELCHAIR_BOARDING _wheelchair_boarding;
  Stop::LOCATION_TYPE _location_type;
};

}}

#endif  // GTFSPARSER_GTFS_STOP_H_
