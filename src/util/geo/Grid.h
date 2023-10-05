// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_GRID_H_
#define UTIL_GEO_GRID_H_

#include <map>
#include <set>
#include <vector>
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
  Grid(const Grid<V, G, T>&) = delete;
  Grid(Grid<V, G, T>&& o)
      : _width(o._width),
        _height(o._height),
        _cellWidth(o._cellWidth),
        _cellHeight(o._cellHeight),
        _bb(o._bb),
        _xWidth(o._xWidth),
        _yHeight(o._yHeight),
        _hasValIdx(o._hasValIdx),
        _grid(o._grid),
        _index(o._index),
        _removed(o._removed) {
    o._grid = 0;
  }

  Grid<V, G, T>& operator=(Grid<V, G, T>&& o) {
    _width = o._width;
    _height = o._height;
    _cellWidth = o._cellWidth;
    _cellHeight = o._cellHeight;
    _bb = o._bb;
    _xWidth = o._xWidth;
    _yHeight = o._yHeight;
    _hasValIdx = o._hasValIdx;
    _grid = o._grid;
    _index = std::move(o._index);
    _removed = std::move(o._removed);
    o._grid = 0;

    return *this;
  };

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

  ~Grid() {
    if (!_grid) return;
    for (size_t i = 0; i < _xWidth; i++) {
      delete[] _grid[i];
    }
    delete[] _grid;
  }

  // add object t to this grid
  void add(G<T> geom, V val);
  void add(size_t x, size_t y, V val);

  void get(const Box<T>& btbox, std::set<V>* s) const;

  template <template <typename> class GG>
  void get(const GG<T>& geom, double d, std::set<V>* s) const;

  template <template <typename> class GG>
  void get(const std::vector<GG<T>>& geom, double d, std::set<V>* s) const;

  void get(size_t x, size_t y, std::set<V>* s) const;
  void remove(V val);

  const std::set<V>& getCell(size_t x, size_t y) const;

  void getNeighbors(const V& val, double d, std::set<V>* s) const;
  void getCellNeighbors(const V& val, size_t d, std::set<V>* s) const;
  void getCellNeighbors(size_t x, size_t y, size_t xPerm, size_t yPerm,
                        std::set<V>* s) const;

  std::set<std::pair<size_t, size_t> > getCells(const V& val) const;

  size_t getXWidth() const;
  size_t getYHeight() const;

  size_t getCellXFromX(double lon) const;
  size_t getCellYFromY(double lat) const;

  Box<T> getBox(size_t x, size_t y) const;
  Box<T> getBBox() const { return _bb; };

 private:
  double _width;
  double _height;

  double _cellWidth;
  double _cellHeight;

  Box<T> _bb;

  size_t _xWidth;
  size_t _yHeight;

  bool _hasValIdx;

  // raw 2d array, less memory overhead
  std::set<V>** _grid;
  std::map<V, std::set<std::pair<size_t, size_t> > > _index;
  std::set<V> _removed;
};

#include "util/geo/Grid.tpp"

}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_GRID_H_
