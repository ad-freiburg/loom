// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_GRID_H_
#define UTIL_GEO_GRID_H_

#include <set>
#include <vector>
#include <map>
#include "util/geo/Geo.h"

namespace util {
namespace geo {

class GridException : public std::runtime_error {
 public:
  GridException(std::string const& msg) : std::runtime_error(msg) {}
};

template <typename V, template <typename> class G, typename T>
class Grid {
 public:
  // initialization of a point grid with cell width w and cell height h
  // that covers the area of bounding box bbox
  Grid(double w, double h, const Box<T>& bbox);

  // initialization of a point grid with cell width w and cell height h
  // that covers the area of bounding box bbox
  // optional parameters specifies whether a value->cell index
  // should be kept (true by default!)
  Grid(double w, double h, const Box<T>& bbox, bool buildValIdx);

  // the empty grid
  Grid();
  // the empty grid
  Grid(bool buildValIdx);

  // add object t to this grid
  void add(G<T> geom, V val);
  void add(size_t x, size_t y, V val);

  void get(const Box<T>& btbox, std::set<V>* s) const;
  void get(const G<T>& geom, double d, std::set<V>* s) const;
  void get(size_t x, size_t y, std::set<V>* s) const;
  void remove(V val);

  void getNeighbors(const V& val, double d, std::set<V>* s) const;
  void getCellNeighbors(const V& val, size_t d, std::set<V>* s) const;
  void getCellNeighbors(size_t x, size_t y, size_t xPerm, size_t yPerm,
                        std::set<V>* s) const;

  std::set<std::pair<size_t, size_t> > getCells(const V& val) const;

  size_t getXWidth() const;
  size_t getYHeight() const;

 private:
  double _width;
  double _height;

  double _cellWidth;
  double _cellHeight;

  Box<T> _bb;

  size_t _counter;

  size_t _xWidth;
  size_t _yHeight;

  bool _hasValIdx;

  std::vector<std::vector<std::set<V> > > _grid;
  std::map<V, std::set<std::pair<size_t, size_t> > > _index;
  std::set<V> _removed;

  Box<T> getBox(size_t x, size_t y) const;

  size_t getCellXFromX(double lon) const;
  size_t getCellYFromY(double lat) const;
};

#include "util/geo/Grid.tpp"

}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_GRID_H_
