
// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_GTFS_STOPTIME_H_
#define GTFSPARSER_GTFS_STOPTIME_H_

#include <iomanip>
#include <vector>
#include <string>
#include "Stop.h"

using std::exception;
using std::string;


namespace gtfsparser {
namespace gtfs {

struct Time {
  Time(uint8_t h, uint8_t m, uint8_t s) : m(m), s(s), h(h) {}
  std::string toString() const {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) <<
      static_cast<int>(h) << ":" <<
      std::setfill('0') << std::setw(2) <<
      static_cast<int>(m) << ":" <<
      std::setfill('0') << std::setw(2) <<
      static_cast<int>(s);
    return ss.str();
  }
  uint8_t m : 6;
  uint8_t s : 6;
  uint8_t h : 8;
};

inline bool operator>(const Time& lh, const Time& rh) {
  return  lh.h * 3600 + lh.m * 60 + lh.s > rh.h * 3600 + rh.m * 60 + rh.s;
}
inline bool operator<(const Time& lh, const Time& rh) { return rh > lh; }
inline bool operator==(const Time& lh, const Time& rh) { return !(rh > lh) && !(rh < lh); }
inline bool operator!=(const Time& lh, const Time& rh) { return !(rh == lh); }
inline bool operator>=(const Time& lh, const Time& rh) { return lh > rh || lh == rh; }
inline bool operator<=(const Time& lh, const Time& rh) { return lh < rh || lh == rh; }


class StopTime {
 public:
  enum PU_DO_TYPE : uint8_t {
    REGULAR = 0,
    NEVER = 1,
    MUST_PHONE_AGENCY = 2,
    MUST_COORDINATE_W_DRIVER = 3
  };

  StopTime(const Time& at, const Time& dt, Stop* s, uint32_t seq,
    const std::string& hs,
    PU_DO_TYPE put, PU_DO_TYPE dot, float distTrav, bool isTp)
  : _at(at), _dt(dt), _s(s), _sequence(seq),
    _headsign(hs), _pickupType(put), _dropOffType(dot),
    _shapeDistTravelled(distTrav), _isTimepoint(isTp) {}

  const Time& getArrivalTime() const { return _at; }
  const Time& getDepartureTime() const { return _dt; }

  Stop* getStop() const { return _s; }
  const std::string& getHeadsign() const { return _headsign; }
  PU_DO_TYPE getPickupType() const { return _pickupType; }
  PU_DO_TYPE getDropOffType() const { return _dropOffType; }
  float getShapeDistanceTravelled() const { return _shapeDistTravelled; }
  bool isTimepoint() const { return _isTimepoint; }
  uint16_t getSeq() const { return _sequence; }
 private:
  Time _at;
  Time _dt;

  Stop* _s;
  uint32_t _sequence;
  std::string _headsign;
  PU_DO_TYPE _pickupType;
  PU_DO_TYPE _dropOffType;
  float _shapeDistTravelled;
  bool _isTimepoint;
};

struct StopTimeCompare {
  bool operator() (const StopTime& lh, const StopTime& rh) const {
    return lh.getSeq() < rh.getSeq();
  }
};

}}

#endif  // GTFSPARSER_GTFS_STOPTIME_H_
