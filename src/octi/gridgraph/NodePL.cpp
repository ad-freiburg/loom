// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/gridgraph/NodePL.h"

using util::geo::Point;
using namespace octi::gridgraph;

// _____________________________________________________________________________
NodePL::NodePL(Point pos) : _pos(pos) {
}

// _____________________________________________________________________________
const Point* NodePL::getGeom() const {
  return &_pos;
}

// _____________________________________________________________________________
void NodePL::getAttrs(json::object_t& obj) const {

}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* NodePL::getNorthPort() const {
  return _ports[0];
}

// _____________________________________________________________________________
void NodePL::setNorthPort(Node<NodePL, EdgePL>* n) {
  _ports[0] = n;
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* NodePL::getSouthPort() const {
  return _ports[4];
}

// _____________________________________________________________________________
void NodePL::setSouthPort(Node<NodePL, EdgePL>* n) {
  _ports[4] = n;
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* NodePL::getWestPort() const {
  return _ports[6];
}

// _____________________________________________________________________________
void NodePL::setWestPort(Node<NodePL, EdgePL>* n) {
  _ports[6] = n;
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* NodePL::getEastPort() const {
  return _ports[2];
}

// _____________________________________________________________________________
void NodePL::setEastPort(Node<NodePL, EdgePL>* n) {
  _ports[2] = n;
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* NodePL::getSouthWestPort() const {
  return _ports[5];
}

// _____________________________________________________________________________
void NodePL::setSouthWestPort(Node<NodePL, EdgePL>* n) {
  _ports[5] = n;
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* NodePL::getSouthEastPort() const {
  return _ports[3];
}

// _____________________________________________________________________________
void NodePL::setSouthEastPort(Node<NodePL, EdgePL>* n) {
  _ports[3] = n;
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* NodePL::getNorthWestPort() const {
  return _ports[7];
}

// _____________________________________________________________________________
void NodePL::setNorthWestPort(Node<NodePL, EdgePL>* n) {
  _ports[7] = n;
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* NodePL::getNorthEastPort() const {
  return _ports[1];
}

// _____________________________________________________________________________
void NodePL::setNorthEastPort(Node<NodePL, EdgePL>* n) {
  _ports[1] = n;
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* NodePL::getPort(size_t i) const {
  return _ports[i];
}