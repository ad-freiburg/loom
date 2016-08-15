// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_GTFS_ROUTE_H_
#define GTFSPARSER_GTFS_ROUTE_H_

#include <stdint.h>
#include <string>
#include <sstream>
#include "Agency.h"

using std::exception;
using std::string;


namespace gtfsparser {
namespace gtfs {

class Route {

 public:
  enum TYPE: uint8_t {
    TRAM = 0,
    SUBWAY = 1,
    RAIL = 2,
    BUS = 3,
    FERRY = 4,
    CABLE_CAR = 5,
    GONDOLA = 6,
    FUNICULAR = 7
  };

  Route() {};

  Route(const string& id, Agency* agency, const string& short_name,
    const string& long_name, const string& desc, Route::TYPE type,
    const string& url, uint32_t color, uint32_t text_color)
  : _id(id), _agency(agency), _short_name(short_name), _long_name(long_name),
    _desc(desc),  _type(type), _url(url), _color(color),
    _text_color(text_color) {}

  const std::string& getId() const {
    return _id;
  }

  const Agency* getAgency() const {
    return _agency;
  }

  Agency* getAgency() {
    return _agency;
  }

  const std::string& getShortName() const {
    return _short_name;
  }

  const std::string& getLongName() const {
    return _long_name;
  }

  const std::string& getDesc() const {
    return _desc;
  }

  Route::TYPE getType() const {
    return _type;
  }

  const std::string& getUrl() const {
    return _url;
  }

  uint32_t getColor() const {
    return _color;
  }

  std::string getColorString() const {
    return getHexColorString(_color);
  }

  uint32_t getTextColor() const {
    return _text_color;
  }

  std::string getTextColorString() const {
    return getHexColorString(_text_color);
  }



  // TODO: implement setters


 private:
  string _id;
  Agency* _agency;
  string _short_name;
  string _long_name;
  string _desc;
  Route::TYPE _type;
  string _url;
  uint32_t _color;
  uint32_t _text_color;

  std::string getHexColorString(uint32_t color) const {
    // using stringstream here, because it doesnt add "0x" to the front
    std::stringstream ss;
    ss << std::hex << color;
    return ss.str();
  }

};

}}

#endif  // GTFSPARSER_GTFS_ROUTE_H_
