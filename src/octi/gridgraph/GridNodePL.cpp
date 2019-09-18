// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/gridgraph/GridNodePL.h"

using util::geo::Point;
using namespace octi::gridgraph;

// _____________________________________________________________________________
GridNodePL::GridNodePL(Point<double> pos)
    : _pos(pos), _parent(0), _closed(false) {}

// _____________________________________________________________________________
const Point<double>* GridNodePL::getGeom() const { return &_pos; }

// _____________________________________________________________________________
util::json::Dict GridNodePL::getAttrs() const { return util::json::Dict(); }

// _____________________________________________________________________________
GridNode* GridNodePL::getParent() const { return _parent; }

// _____________________________________________________________________________
void GridNodePL::setParent(GridNode* n) { _parent = n; }

// _____________________________________________________________________________
GridNode* GridNodePL::getStepMother() const { return _parent; }

// _____________________________________________________________________________
void GridNodePL::setStepMother(GridNode* n) { _parent = n; }

// _____________________________________________________________________________
GridNode* GridNodePL::getPort(size_t i) const {
  return _ports[i];
}

// _____________________________________________________________________________
void GridNodePL::setPort(size_t p, GridNode* n) { _ports[p] = n; }

// _____________________________________________________________________________
const std::vector<GridNode*>& GridNodePL::getMetaPorts(size_t i) const {
  return _metaPorts[i];
}

// _____________________________________________________________________________
void GridNodePL::addMetaPort(size_t i, GridNode* nd) {
  _metaPorts[i].push_back(nd);
}

// _____________________________________________________________________________
void GridNodePL::clearMetaPorts() {
  _metaPorts->clear();
}

// _____________________________________________________________________________
const std::vector<GridNode*>& GridNodePL::getStepChilds() const {
  return _stepChilds;
}

// _____________________________________________________________________________
void GridNodePL::addStepChild(GridNode* nd) {
  _stepChilds.push_back(nd);
}


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
