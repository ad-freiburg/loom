// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "shared/linegraph/LineGraph.h"
#include "shared/linegraph/LineNodePL.h"
#include "shared/linegraph/NodeFront.h"

using shared::linegraph::LineNodePL;
using shared::linegraph::NodeFront;
using shared::linegraph::Station;
using util::geo::DPoint;
using util::geo::Point;

// _____________________________________________________________________________
LineNodePL::LineNodePL(Point<double> pos) : _pos(pos) {}

// _____________________________________________________________________________
const Point<double>* LineNodePL::getGeom() const { return &_pos; }

// _____________________________________________________________________________
void LineNodePL::clearConnExc() { _connEx.clear(); }

// _____________________________________________________________________________
void LineNodePL::setGeom(const Point<double>& p) { _pos = p; }

// _____________________________________________________________________________
std::string LineNodePL::toString() const {
  std::stringstream ret;

  ret << "<nd " << this;

  ret << " @ " << _pos.getX() << "," << _pos.getY();

  char sep = ' ';
  for (auto st : _is) {
    ret << sep;
    ret << "\"" << st.name << "\" (" << st.id << ")";
    sep = ',';
  }

  ret << ">";

  return ret.str();
}

// _____________________________________________________________________________
util::json::Dict LineNodePL::getAttrs() const {
  util::json::Dict obj;
  if (_is.size() > 0) {
    obj["station_id"] = _is.begin()->id;
    obj["station_label"] = _is.begin()->name;
  }

  auto arr = util::json::Array();

  for (const auto& ro : _connEx) {
    for (const auto& exFr : ro.second) {
      for (const auto* exTo : exFr.second) {
        util::json::Dict ex;
        ex["route"] = util::toString(ro.first->id());
        if (exFr.first == exTo) continue;
        auto shrd = LineGraph::sharedNode(exFr.first, exTo);
        if (!shrd) continue;
        auto nd1 = exFr.first->getOtherNd(shrd);
        auto nd2 = exTo->getOtherNd(shrd);
        ex["edge1_node"] = util::toString(nd1);
        ex["edge2_node"] = util::toString(nd2);
        arr.push_back(ex);
      }
    }
  }

  if (arr.size()) obj["excluded_line_conns"] = arr;

  auto nonServedArr = util::json::Array();

  for (const auto& no : _notServed) {
    nonServedArr.push_back(util::toString(no->id()));
  }

  if (nonServedArr.size()) obj["not_serving"] = nonServedArr;

  return obj;
}

// _____________________________________________________________________________
void LineNodePL::addStop(const Station& i) { _is.push_back(i); }

// _____________________________________________________________________________
const std::vector<Station>& LineNodePL::stops() const { return _is; }

// _____________________________________________________________________________
void LineNodePL::clearStops() { _is.clear(); }

// _____________________________________________________________________________
void LineNodePL::addConnExc(const Line* r, const LineEdge* edgeA,
                            const LineEdge* edgeB) {
  _connEx[r][edgeA].insert(edgeB);
  // index the other direction also, will lead to faster lookups later on
  _connEx[r][edgeB].insert(edgeA);
}

// _____________________________________________________________________________
void LineNodePL::delConnExc(const Line* r, const LineEdge* edgeA,
                            const LineEdge* edgeB) {
  _connEx[r][edgeA].erase(edgeB);
  // index the other direction also, will lead to faster lookups later on
  _connEx[r][edgeB].erase(edgeA);
}

// _____________________________________________________________________________
bool LineNodePL::connOccurs(const Line* r, const LineEdge* edgeA,
                            const LineEdge* edgeB) const {
  const auto& i = _connEx.find(r);
  if (_connEx.find(r) == _connEx.end()) return true;

  const auto& ii = i->second.find(edgeA);
  if (ii == i->second.end()) return true;

  return ii->second.find(edgeB) == ii->second.end();
}

// _____________________________________________________________________________
const NodeFront* LineNodePL::frontFor(const LineEdge* e) const {
  auto i = _edgToNf.find(e);
  if (i == _edgToNf.end()) return 0;
  return &_nodeFronts[i->second];
}

// _____________________________________________________________________________
const std::vector<NodeFront>& LineNodePL::fronts() const { return _nodeFronts; }

// _____________________________________________________________________________
void LineNodePL::delFrontFor(const LineEdge* e) {
  auto i = _edgToNf.find(e);
  if (i == _edgToNf.end()) return;
  size_t idx = i->second;
  _edgToNf.erase(e);

  _edgToNf[_nodeFronts.back().edge] = idx;
  _nodeFronts[idx] = _nodeFronts.back();
  _nodeFronts.pop_back();
}

// _____________________________________________________________________________
std::vector<NodeFront>& LineNodePL::fronts() { return _nodeFronts; }

// _____________________________________________________________________________
void LineNodePL::addFront(const NodeFront& f) {
  assert(!frontFor(f.edge));
  _edgToNf[f.edge] = _nodeFronts.size();
  _nodeFronts.push_back(f);
}

// _____________________________________________________________________________
void LineNodePL::addLineNotServed(const Line* r) { _notServed.insert(r); }

// _____________________________________________________________________________
bool LineNodePL::lineServed(const Line* r) const {
  return !_notServed.count(r);
}
