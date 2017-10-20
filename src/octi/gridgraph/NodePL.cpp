// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/gridgraph/NodePL.h"

using util::geo::Point;
using namespace octi::gridgraph;

// _____________________________________________________________________________
NodePL::NodePL(Point pos) : _pos(pos), _parent(0), _closed(false) {}

// _____________________________________________________________________________
const Point* NodePL::getGeom() const { return &_pos; }

// _____________________________________________________________________________
void NodePL::getAttrs(json::object_t& obj) const {}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* NodePL::getParent() const { return _parent; }

// _____________________________________________________________________________
void NodePL::setParent(Node<NodePL, EdgePL>* n) { _parent = n; }

// _____________________________________________________________________________
Node<NodePL, EdgePL>* NodePL::getPort(size_t i) const { return _ports[i]; }

// _____________________________________________________________________________
void NodePL::setPort(size_t p, Node<NodePL, EdgePL>* n) { _ports[p] = n; }

// _____________________________________________________________________________
void NodePL::setXY(size_t x, size_t y) {
  _x = x;
  _y = y;
}

// _____________________________________________________________________________
size_t NodePL::getX() const { return _parent->pl()._x; }

// _____________________________________________________________________________
size_t NodePL::getY() const { return _parent->pl()._y; }

// _____________________________________________________________________________
void NodePL::setClosed(bool c) { _closed = c; }

// _____________________________________________________________________________
bool NodePL::isClosed() const { return _closed; }

