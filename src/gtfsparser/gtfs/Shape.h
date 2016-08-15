// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_GTFS_SHAPE_H_
#define GTFSPARSER_GTFS_SHAPE_H_

#include <set>
#include <vector>
#include <string>

using std::exception;
using std::string;


namespace gtfsparser {
namespace gtfs {

struct ShapePoint {
  ShapePoint(double lat, double ln, double dist, uint32_t seq)
  : lat(lat), lng(ln), travelDist(dist), seq(seq) {}
  ShapePoint()
  : lat(0), lng(0), travelDist(-1), seq(0) {}
  double lat, lng, travelDist;
  uint32_t seq;
};

struct ShapePointCompare {
  bool operator() (const ShapePoint& lh, const ShapePoint& rh) const {
    return lh.seq < rh.seq;
  }
};

typedef std::set<ShapePoint, ShapePointCompare> ShapePoints;

class Shape {

 public:
  Shape() {}

  Shape(const string& id) : _id(id) {}

  std::string getId() const {
    return _id;
  }

  const ShapePoints& getPoints() const {
    return _shapePoints;
  }

  bool addPoint(ShapePoint p) {
    return _shapePoints.insert(p).second;
  }

 private:
  string _id;
  ShapePoints _shapePoints;
};

}}

#endif  // GTFSPARSER_GTFS_SHAPE_H_
