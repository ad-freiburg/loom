// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosip@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
Grid<V, G, T>::Grid(bool bldIdx)
    : _width(0),
      _height(0),
      _cellWidth(0),
      _cellHeight(0),
      _xWidth(0),
      _yHeight(0),
      _hasValIdx(bldIdx) {}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
Grid<V, G, T>::Grid() : Grid<V, G, T>(true) {}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
Grid<V, G, T>::Grid(double w, double h, const Box<T>& bbox)
    : Grid<V, G, T>(w, h, bbox, true) {}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
Grid<V, G, T>::Grid(double w, double h, const Box<T>& bbox, bool bValIdx)
    : _cellWidth(fabs(w)),
      _cellHeight(fabs(h)),
      _bb(bbox),
      _hasValIdx(bValIdx) {
  _width =
      bbox.max_corner().template get<0>() - bbox.min_corner().template get<0>();
  _height =
      bbox.max_corner().template get<1>() - bbox.min_corner().template get<1>();

  if (_width < 0 || _height < 0) {
    _width = 0;
    _height = 0;
    _xWidth = 0;
    _yHeight = 0;
    return;
  }

  _xWidth = ceil(_width / _cellWidth);
  _yHeight = ceil(_height / _cellHeight);

  // resize rows
  _grid.resize(_xWidth);

  // resize columns
  for (size_t i = 0; i < _xWidth; i++) {
    _grid[i].resize(_yHeight);
  }
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
void Grid<V, G, T>::add(G<T> geom, V val) {
  Box<T> box = getBoundingBox(geom);
  size_t swX = getCellXFromX(box.min_corner().template get<0>());
  size_t swY = getCellYFromY(box.min_corner().template get<1>());

  size_t neX = getCellXFromX(box.max_corner().template get<0>());
  size_t neY = getCellYFromY(box.max_corner().template get<1>());

  for (size_t x = swX; x <= neX && x < _grid.size(); x++) {
    for (size_t y = swY; y <= neY && y < _grid[x].size(); y++) {
      if (intersects(geom, getBox(x, y))) {
        add(x, y, val);
      }
    }
  }
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
void Grid<V, G, T>::add(size_t x, size_t y, V val) {
  _grid[x][y].insert(val);
  if (_hasValIdx) _index[val].insert(std::pair<size_t, size_t>(x, y));
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
void Grid<V, G, T>::get(const Box<T>& box, std::set<V>* s) const {
  size_t swX = getCellXFromX(box.min_corner().template get<0>());
  size_t swY = getCellYFromY(box.min_corner().template get<1>());

  size_t neX = getCellXFromX(box.max_corner().template get<0>());
  size_t neY = getCellYFromY(box.max_corner().template get<1>());

  for (size_t x = swX; x <= neX && x >= 0 && x < _xWidth; x++)
    for (size_t y = swY; y <= neY && y >= 0 && y < _yHeight; y++) get(x, y, s);
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
void Grid<V, G, T>::get(const G<T>& geom, double d, std::set<V>* s) const {
  Box<T> a = getBoundingBox(geom);
  Box<T> b(Point<T>(a.min_corner().template get<0>() - d,
                    a.min_corner().template get<1>() - d),
           Point<T>(a.max_corner().template get<0>() + d,
                    a.max_corner().template get<1>() + d));
  return get(b, s);
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
void Grid<V, G, T>::get(size_t x, size_t y, std::set<V>* s) const {
  if (_hasValIdx) {
    s->insert(_grid[x][y].begin(), _grid[x][y].end());
  } else {
    // if we dont have a value index, we have a set of deleted nodes.
    // in this case, only insert if not deleted
    std::copy_if(_grid[x][y].begin(), _grid[x][y].end(),
                 std::inserter(*s, s->end()),
                 [&](const V& v) { return Grid<V, G, T>::_removed.count(v) == 0; });
  }
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
void Grid<V, G, T>::remove(V val) {
  if (_hasValIdx) {
    auto i = _index.find(val);
    if (i == _index.end()) return;

    for (auto pair : i->second) {
      _grid[pair.first][pair.second].erase(
          _grid[pair.first][pair.second].find(val));
    }

    _index.erase(i);
  } else {
    _removed.insert(val);
  }
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
void Grid<V, G, T>::getNeighbors(const V& val, double d, std::set<V>* s) const {
  if (!_hasValIdx) throw GridException("No value index build!");
  auto it = _index.find(val);
  if (it == _index.end()) return;

  size_t xPerm = ceil(d / _cellWidth);
  size_t yPerm = ceil(d / _cellHeight);

  for (auto pair : it->second) {
    getCellNeighbors(pair.first, pair.second, xPerm, yPerm, s);
  }
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
void Grid<V, G, T>::getCellNeighbors(const V& val, size_t d,
                                     std::set<V>* s) const {
  if (!_hasValIdx) throw GridException("No value index build!");
  auto it = _index.find(val);
  if (it == _index.end()) return;

  for (auto pair : it->second) {
    getCellNeighbors(pair.first, pair.second, d, d, s);
  }
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
void Grid<V, G, T>::getCellNeighbors(size_t cx, size_t cy, size_t xPerm,
                                     size_t yPerm, std::set<V>* s) const {
  size_t swX = xPerm > cx ? 0 : cx - xPerm;
  size_t swY = yPerm > cy ? 0 : cy - yPerm;

  size_t neX = xPerm + cx + 1 > _xWidth ? _xWidth : cx + xPerm + 1;
  size_t neY = yPerm + cy + 1 > _yHeight ? _yHeight : cy + yPerm + 1;

  for (size_t x = swX; x < neX; x++) {
    for (size_t y = swY; y < neY; y++) {
      s->insert(_grid[x][y].begin(), _grid[x][y].end());
    }
  }
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
std::set<std::pair<size_t, size_t> > Grid<V, G, T>::getCells(
    const V& val) const {
  if (!_hasValIdx) throw GridException("No value index build!");
  return _index.find(val)->second;
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
Box<T> Grid<V, G, T>::getBox(size_t x, size_t y) const {
  Point<T> sw(_bb.min_corner().template get<0>() + x * _cellWidth,
              _bb.min_corner().template get<1>() + y * _cellHeight);
  Point<T> ne(_bb.min_corner().template get<0>() + (x + 1) * _cellWidth,
              _bb.min_corner().template get<1>() + (y + 1) * _cellHeight);
  return Box<T>(sw, ne);
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
size_t Grid<V, G, T>::getCellXFromX(double x) const {
  float dist = x - _bb.min_corner().template get<0>();
  if (dist < 0) dist = 0;
  return floor(dist / _cellWidth);
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
size_t Grid<V, G, T>::getCellYFromY(double y) const {
  float dist = y - _bb.min_corner().template get<1>();
  if (dist < 0) dist = 0;
  return floor(dist / _cellHeight);
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
size_t Grid<V, G, T>::getXWidth() const {
  return _xWidth;
}

// _____________________________________________________________________________
template <typename V, template <typename> typename G, typename T>
size_t Grid<V, G, T>::getYHeight() const {
  return _yHeight;
}
