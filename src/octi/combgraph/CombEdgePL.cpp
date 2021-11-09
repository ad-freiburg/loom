// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/combgraph/CombEdgePL.h"
#include "util/geo/PolyLine.h"

using util::geo::PolyLine;
using octi::combgraph::CombEdgePL;

// _____________________________________________________________________________
CombEdgePL::CombEdgePL(shared::linegraph::LineEdge* child) : _maxLineNum(0) {
  _childs.push_back(child);
  _geom = PolyLine<double>(*child->getFrom()->pl().getGeom(),
                           *child->getTo()->pl().getGeom());
}

// _____________________________________________________________________________
const util::geo::Line<double>* CombEdgePL::getGeom() const {
  return &_geom.getLine();
}

// _____________________________________________________________________________
const PolyLine<double>& CombEdgePL::getPolyLine() const { return _geom; }

// _____________________________________________________________________________
void CombEdgePL::setPolyLine(const PolyLine<double>& p) { _geom = p; }

// _____________________________________________________________________________
util::json::Dict CombEdgePL::getAttrs() const {
  util::json::Dict obj;
  util::json::Array inclEArr;

  for (auto e : _childs) {
    auto ret = e->pl().getAttrs();

    auto edge = util::json::Dict();
    edge["from"] = e->getFrom();
    edge["to"] = e->getTo();
    edge["lines"] = ret;

    inclEArr.push_back(edge);
  }
  obj["included_edges"] = inclEArr;

  return obj;
}

// _____________________________________________________________________________
const std::vector<shared::linegraph::LineEdge*>& CombEdgePL::getChilds() const {
  return _childs;
}

// _____________________________________________________________________________
std::vector<shared::linegraph::LineEdge*>& CombEdgePL::getChilds() {
  return _childs;
}
