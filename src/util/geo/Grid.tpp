// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosip@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename V, typename G>
Grid<V, G>::Grid(double w, double h, const Box& bbox)
    : _cellWidth(fabs(w)),
      _cellHeight(fabs(h)),
      _bb(bbox) {
  _width = bbox.max_corner().get<0>() - bbox.min_corner().get<0>();
  _height = bbox.max_corner().get<1>() - bbox.min_corner().get<1>();

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
template <typename V, typename G>
void Grid<V, G>::add(G geom, V val) {
  Box box = getBoundingBox(geom);
  size_t swX = getCellXFromX(box.min_corner().get<0>());
  size_t swY = getCellYFromY(box.min_corner().get<1>());

  size_t neX = getCellXFromX(box.max_corner().get<0>());
  size_t neY = getCellYFromY(box.max_corner().get<1>());

  for (size_t x = swX; x <= neX && x < _grid.size(); x++) {
    for (size_t y = swY; y <= neY && y < _grid[x].size(); y++) {
      if (intersects(geom, getBox(x, y))) {
		add(x, y, val);
      }
    }
  }
}

// _____________________________________________________________________________
template <typename V, typename G>
void Grid<V, G>::add(size_t x, size_t y, V val) {
  _grid[x][y].insert(val);
  _index[val].insert(std::pair<size_t, size_t>(x, y));
}

// _____________________________________________________________________________
template <typename V, typename G>
void Grid<V, G>::get(const Box& box, std::set<V>* s) const {
  size_t swX = getCellXFromX(box.min_corner().get<0>());
  size_t swY = getCellYFromY(box.min_corner().get<1>());

  size_t neX = getCellXFromX(box.max_corner().get<0>());
  size_t neY = getCellYFromY(box.max_corner().get<1>());

  for (size_t x = swX; x <= neX && x >= 0 && x < _xWidth; x++)
    for (size_t y = swY; y <= neY && y >= 0 && y < _yHeight; y++)
	  get(x, y, s);
}

// _____________________________________________________________________________
template <typename V, typename G>
void Grid<V, G>::get(const G& geom, double d, std::set<V>* s) const {
  Box a = getBoundingBox(geom);
  Box b(Point(a.min_corner().get<0>() - d, a.min_corner().get<1>() - d),
    Point(a.max_corner().get<0>() + d, a.max_corner().get<1>() + d));
  return get(b, s);
}

// _____________________________________________________________________________
template <typename V, typename G>
void Grid<V, G>::get(size_t x, size_t y, std::set<V>* s) const {
  s->insert(_grid[x][y].begin(), _grid[x][y].end());
}

// _____________________________________________________________________________
template <typename V, typename G>
void Grid<V, G>::remove(V val) {
  auto i = _index.find(val);
  if (i == _index.end()) return;

  for (auto pair : i->second) {
    _grid[pair.first][pair.second].erase(_grid[pair.first][pair.second].find(val));
  }

  _index.erase(i);
}

// _____________________________________________________________________________
template <typename V, typename G>
void Grid<V, G>::getNeighbors(const V& val, double d, std::set<V>* s) const {
  auto it = _index.find(val);
  if (it == _index.end()) return;

  size_t xPerm = ceil(d / _cellWidth);
  size_t yPerm = ceil(d / _cellHeight);

  for (auto pair : it->second) {
	getCellNeighbors(pair.first, pair.second, xPerm, yPerm, s);
  }
}

// _____________________________________________________________________________
template <typename V, typename G>
void Grid<V, G>::getCellNeighbors(const V& val, size_t d, std::set<V>* s) const {
  auto it = _index.find(val);
  if (it == _index.end()) return;

  for (auto pair : it->second) {
	getCellNeighbors(pair.first, pair.second, d, d, s);
  }
}

// _____________________________________________________________________________
template <typename V, typename G>
void Grid<V, G>::getCellNeighbors(size_t cx, size_t cy, size_t xPerm, size_t yPerm, std::set<V>* s) const {
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
template <typename V, typename G>
std::set<std::pair<size_t, size_t> > Grid<V, G>::getCells(const V& val) const {
  return _index.find(val)->second;
}

// _____________________________________________________________________________
template <typename V, typename G>
Box Grid<V, G>::getBox(size_t x, size_t y) const {
  Point sw(_bb.min_corner().get<0>() + x * _cellWidth,
           _bb.min_corner().get<1>() + y * _cellHeight);
  Point ne(_bb.min_corner().get<0>() + (x + 1) * _cellWidth,
           _bb.min_corner().get<1>() + (y + 1) * _cellHeight);
  return Box(sw, ne);
}

// _____________________________________________________________________________
template <typename V, typename G>
size_t Grid<V, G>::getCellXFromX(double x) const {
  float dist = x - _bb.min_corner().get<0>();
  if (dist < 0) dist = 0;
  return floor(dist / _cellWidth);
}

// _____________________________________________________________________________
template <typename V, typename G>
size_t Grid<V, G>::getCellYFromY(double y) const {
  float dist = y - _bb.min_corner().get<1>();
  if (dist < 0) dist = 0;
  return floor(dist / _cellHeight);
}

// _____________________________________________________________________________
template <typename V, typename G>
size_t Grid<V, G>::getXWidth() const {
  return _xWidth;
}

// _____________________________________________________________________________
template <typename V, typename G>
size_t Grid<V, G>::getYHeight() const {
  return _yHeight;
}
