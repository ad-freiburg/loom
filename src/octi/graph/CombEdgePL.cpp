// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "util/geo/PolyLine.h"
#include "octi/graph/CombEdgePL.h"

using util::geo::PolyLine;
using namespace octi::graph;

// _____________________________________________________________________________
CombEdgePL::CombEdgePL(Edge* child) {
  _childs.push_back(child);
  _geom = PolyLine(*child->getFrom()->pl().getGeom(), *child->getTo()->pl().getGeom());
}

// _____________________________________________________________________________
const util::geo::Line* CombEdgePL::getGeom() const {
  return &_geom.getLine();
}

// _____________________________________________________________________________
const PolyLine& CombEdgePL::getPolyLine() const {
  return _geom;
}

// _____________________________________________________________________________
void CombEdgePL::setPolyLine(const PolyLine& p) {
  _geom = p;
}

// _____________________________________________________________________________
void CombEdgePL::getAttrs(json::object_t& obj) const {
  obj["included_edges"] = json::array();
  obj["generation"] = _generation;

  for (auto e : _childs) {
    json::object_t ret;
    e->pl().getAttrs(ret);

    json edge = json::object();
    edge["from"] = e->getFrom();
    edge["to"] = e->getTo();
    edge["lines"] = ret;

    obj["included_edges"].push_back(edge);
  }
}

// _____________________________________________________________________________
std::vector<Edge*>& CombEdgePL::getChilds() {
  return _childs;
}

// _____________________________________________________________________________
void CombEdgePL::setGeneration(size_t g) {
  _generation = g;
}

// _____________________________________________________________________________
size_t CombEdgePL::getGeneration() const {
  return _generation;
}

