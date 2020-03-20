// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "shared/linegraph/Line.h"
#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineGraph.h"
#include "shared/linegraph/LineNodePL.h"
#include "shared/style/LineStyle.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"

using util::geo::PolyLine;
using shared::linegraph::LineNode;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineOcc;

// _____________________________________________________________________________
LineEdgePL::LineEdgePL() : _dontContract(false) {}

// _____________________________________________________________________________
LineEdgePL::LineEdgePL(const PolyLine<double>& p) : _p(p), _dontContract(false) {}

// _____________________________________________________________________________
const util::geo::Line<double>* LineEdgePL::getGeom() const {
  return &_p.getLine();
}

// _____________________________________________________________________________
void LineEdgePL::setGeom(const util::geo::Line<double>& l) { _p = l; }

// _____________________________________________________________________________
const PolyLine<double>& LineEdgePL::getPolyline() const { return _p; }

// _____________________________________________________________________________
void LineEdgePL::setPolyline(const PolyLine<double>& p) { _p = p; }

// _____________________________________________________________________________
void LineEdgePL::addLine(const Line* r, const LineNode* dir,
                          util::Nullable<shared::style::LineStyle> ls) {
  LineOcc occ(r, dir, ls);
  auto f = _lines.find(occ);
  if (f != _lines.end()) {
    const auto& prev = *f;
    // the route is already present in both directions, ignore newly inserted
    if (prev.direction == 0) return;

    // the route is already present in the same direction, ignore newly inserted
    if (prev.direction == dir) return;

    // the route is already present in the other direction, make two-way
    if (prev.direction != dir) {
      occ.direction = 0;
      _lines.erase(f);
    }
  }
  _lines.insert(occ);
}

// _____________________________________________________________________________
void LineEdgePL::addLine(const Line* r, const LineNode* dir) {
  addLine(r, dir, util::Nullable<shared::style::LineStyle>());
}

// _____________________________________________________________________________
void LineEdgePL::delLine(const Line* r) {
  LineOcc occ(r, 0);
  _lines.erase(occ);
}

// _____________________________________________________________________________
const std::set<LineOcc>& LineEdgePL::getLines() const { return _lines; }

// _____________________________________________________________________________
std::set<LineOcc>& LineEdgePL::getLines() { return _lines; }

// _____________________________________________________________________________
util::json::Dict LineEdgePL::getAttrs() const {
  util::json::Dict obj;
  auto arr = util::json::Array();

  std::string dbg_lines = "";
  bool first = true;

  for (auto r : getLines()) {
    auto route = util::json::Dict();
    route["id"] = r.line->id();
    route["label"] = r.line->label();
    route["color"] = r.line->color();

    if (r.direction != 0) {
      route["direction"] = util::toString(r.direction);
      dbg_lines += (first ? "" : "$") + r.line->label() + ">";
    } else {
      dbg_lines += (first ? "" : "$") + r.line->label();
    }

    arr.push_back(route);
    first = false;
  }

  obj["lines"] = arr;
  obj["dbg_lines"] = dbg_lines;

  return obj;
}

// _____________________________________________________________________________
bool LineEdgePL::hasLine(const Line* r) const {
  return _lines.count(LineOcc(r, 0)) > 0;
}

// _____________________________________________________________________________
const LineOcc& LineEdgePL::lineOcc(const Line* r) const {
  return *_lines.find(LineOcc(r, 0));
}

// _____________________________________________________________________________
const LineOcc& LineEdgePL::lineOccAtPos(size_t i) const {
  auto it = _lines.begin();
  for (size_t j = 0; j < i; j++) it++;
  return *it;
}

// _____________________________________________________________________________
size_t LineEdgePL::linePosUnder(const Line* r,
                                const std::vector<size_t> ordering) const {
  size_t i = 0;
  for (const LineOcc& lo : _lines) {
    if (lo.line == r) {
      assert(std::find(ordering.begin(), ordering.end(), i) != ordering.end());
      size_t pos =
          std::find(ordering.begin(), ordering.end(), i) - ordering.begin();
      return pos;
    }
    i++;
  }
  assert(false);
  return -1;
}

// _____________________________________________________________________________
size_t LineEdgePL::linePos(const Line* r) const {
  size_t i = 0;
  for (const LineOcc& ro : _lines) {
    if (ro.line == r) return i;
    i++;
  }
  return -1;
}
