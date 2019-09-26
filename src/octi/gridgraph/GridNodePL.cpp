// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/gridgraph/GridNodePL.h"

using util::geo::Point;
using namespace octi::gridgraph;

// _____________________________________________________________________________
GridNodePL::GridNodePL(Point<double> pos)
    : _pos(pos),
      _parent(0),
      _closed(false),
      _sink(false),
      _station(false),
      _settled(false) {}

// _____________________________________________________________________________
const Point<double>* GridNodePL::getGeom() const { return &_pos; }

// _____________________________________________________________________________
util::json::Dict GridNodePL::getAttrs() const {
  util::json::Dict obj;
  std::vector<std::string> routes;

  return obj;
}

// _____________________________________________________________________________
GridNode* GridNodePL::getParent() const { return _parent; }

// _____________________________________________________________________________
void GridNodePL::setParent(GridNode* n) { _parent = n; }

// _____________________________________________________________________________
GridNode* GridNodePL::getPort(size_t i) const { return _ports[i]; }

// _____________________________________________________________________________
void GridNodePL::setPort(size_t p, GridNode* n) { _ports[p] = n; }

// _____________________________________________________________________________
void GridNodePL::setXY(size_t x, size_t y) {
  _x = x;
  _y = y;
}

// _____________________________________________________________________________
size_t GridNodePL::getX() const { return _parent->pl()._x; }

// _____________________________________________________________________________
size_t GridNodePL::getY() const { return _parent->pl()._y; }

// _____________________________________________________________________________
void GridNodePL::setClosed(bool c) { _closed = c; }

// _____________________________________________________________________________
bool GridNodePL::isClosed() const { return _closed; }

// _____________________________________________________________________________
void GridNodePL::setSettled(bool c) { _settled = c; }

// _____________________________________________________________________________
bool GridNodePL::isSettled() const { return _settled; }

// _____________________________________________________________________________
void GridNodePL::setSink() { _sink = true; }

// _____________________________________________________________________________
bool GridNodePL::isSink() const { return _sink; }

// _____________________________________________________________________________
void GridNodePL::setStation() { _station = true; }
