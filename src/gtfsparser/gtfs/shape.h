// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_GTFS_SHAPE_H_
#define GTFSPARSER_GTFS_SHAPE_H_

#include <vector>
#include <string>

using std::exception;
using std::string;


namespace gtfsparser {
namespace gtfs {

struct ShapePoint {
  double lat, lng, travelDist;
};

typedef std::vector<ShapePoint> ShapePoints;

class Shape {

 public:
  Shape() {}

  Shape(const string& id) {}

  std::string getId() const {
    return _id;
  }

  const ShapePoints& getPoints() const {
    return _shapePoints;
  }

  void addPoint(ShapePoint p) {
    _shapePoints.push_back(p);
  }

 private:
  string _id;
  ShapePoints _shapePoints;
};

}}

#endif  // GTFSPARSER_GTFS_SHAPE_H_
