// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_GRID_H_
#define UTIL_GEO_GRID_H_

#include <set>
#include <vector>
#include "util/geo/Geo.h"

namespace util {
namespace geo {

template <typename V, typename G>
class Grid {
 public:
  // initialization of a point grid with cell width w and cell height h
  // that covers the area of bounding box bbox
  Grid(double w, double h, const Box& bbox);

  // add object t to this grid
  void add(G geom, V val);

  void get(const Box& btbox, std::set<V>* s) const;
  void get(const G& geom, double d, std::set<V>* s) const;
  void remove(V val);

  void getNeighbors(const V& val, double d, std::set<V>* s) const;

  std::set<std::pair<size_t, size_t> > getCells(const V& val) const;

 private:
  double _width;
  double _height;

  double _cellWidth;
  double _cellHeight;

  Box _bb;

  size_t _counter;

  size_t _xWidth;
  size_t _yHeight;

  std::vector<std::vector<std::set<V> > > _grid;
  std::map<V, std::set<std::pair<size_t, size_t> > > _index;

  Box getBox(size_t x, size_t y) const;


  size_t getCellXFromX(double lon) const;
  size_t getCellYFromY(double lat) const;
};

#include "util/geo/Grid.tpp"

}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_GRID_H_
