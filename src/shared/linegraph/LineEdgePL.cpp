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

using shared::linegraph::LineEdgePL;
using shared::linegraph::LineNode;
using shared::linegraph::LineOcc;
using util::geo::PolyLine;

// _____________________________________________________________________________
LineEdgePL::LineEdgePL() : _dontContract(false) {}

// _____________________________________________________________________________
LineEdgePL::LineEdgePL(const PolyLine<double>& p)
    : _dontContract(false), _p(p) {}

// _____________________________________________________________________________
LineEdgePL::LineEdgePL(const LineEdgePL& other)
    : _lineToIdx(other._lineToIdx),
      _lines(other._lines),
      _dontContract(other._dontContract),
      _comp(other._comp),
      _p(other._p) {}

// _____________________________________________________________________________
LineEdgePL::LineEdgePL(LineEdgePL&& other)
    : _lineToIdx(std::move(other._lineToIdx)),
      _lines(std::move(other._lines)),
      _dontContract(other._dontContract),
      _comp(other._comp),
      _p(std::move(other._p)) {}

// _____________________________________________________________________________
const util::geo::Line<double>* LineEdgePL::getGeom() const {
  return &_p.getLine();
}

// _____________________________________________________________________________
void LineEdgePL::setGeom(const util::geo::Line<double>& l) {
  _p = util::geo::PolyLine<double>(l);
}

// _____________________________________________________________________________
const PolyLine<double>& LineEdgePL::getPolyline() const { return _p; }

// _____________________________________________________________________________
PolyLine<double>& LineEdgePL::getPolyline() { return _p; }

// _____________________________________________________________________________
void LineEdgePL::setPolyline(const PolyLine<double>& p) { _p = p; }

// _____________________________________________________________________________
void LineEdgePL::addLine(const Line* r, const LineNode* dir,
                         util::Nullable<shared::style::LineStyle> ls) {
  auto f = _lineToIdx.find(r);
  if (f != _lineToIdx.end()) {
    size_t prevIdx = f->second;
    const auto& prev = _lines[prevIdx];
    // the route is already present in both directions, ignore newly inserted
    if (prev.direction == 0) return;

    // the route is already present in the same direction, ignore newly inserted
    if (prev.direction == dir) return;

    // the route is already present in the other direction, make two-way
    if (prev.direction != dir) {
      _lines[prevIdx].direction = 0;
      return;
    }
  }
  _lineToIdx[r] = _lines.size();
  LineOcc occ(r, dir, ls);
  _lines.push_back(occ);
}

// _____________________________________________________________________________
void LineEdgePL::addLine(const Line* r, const LineNode* dir) {
  addLine(r, dir, util::Nullable<shared::style::LineStyle>());
}

// _____________________________________________________________________________
void LineEdgePL::delLine(const Line* r) {
  _lineToIdx[_lines.back().line] = _lineToIdx.find(r)->second;
  _lines[_lineToIdx.find(r)->second] = _lines.back();
  _lines.resize(_lines.size() - 1);
  _lineToIdx.erase(r);
}

// _____________________________________________________________________________
const std::vector<LineOcc>& LineEdgePL::getLines() const { return _lines; }

// _____________________________________________________________________________
util::json::Dict LineEdgePL::getAttrs() const {
  util::json::Dict obj;
  auto arr = util::json::Array();
  std::string dbg_lines = "";

  for (auto r : getLines()) {
    auto line = util::json::Dict();
    line["id"] = r.line->id();
    line["label"] = r.line->label();
    line["color"] = r.line->color();
    if (!r.style.isNull()) {
      if (r.style.get().getCss().size()) line["style"] = r.style.get().getCss();
      if (r.style.get().getOutlineCss().size())
        line["outline-style"] = r.style.get().getOutlineCss();
    }

    if (r.direction != 0) {
      line["direction"] = util::toString(r.direction);
      dbg_lines += (!arr.size() ? "" : ",") + r.line->label();
    } else {
      dbg_lines += (!arr.size() ? "" : ",") + r.line->label();
    }

    arr.push_back(line);
  }

  obj["lines"] = arr;
  obj["dbg_lines"] = dbg_lines;

  if (_comp != std::numeric_limits<uint32_t>::max()) obj["component"] = _comp;

  return obj;
}

// _____________________________________________________________________________
bool LineEdgePL::hasLine(const Line* l) const { return _lineToIdx.count(l); }

// _____________________________________________________________________________
const LineOcc& LineEdgePL::lineOcc(const Line* l) const {
  return _lines[_lineToIdx.find(l)->second];
}

// _____________________________________________________________________________
const LineOcc& LineEdgePL::lineOccAtPos(size_t i) const { return _lines[i]; }

// _____________________________________________________________________________
void LineEdgePL::updateLineOcc(const LineOcc& occ) {
  _lines[_lineToIdx.find(occ.line)->second] = occ;
}

// _____________________________________________________________________________
void LineEdgePL::writePermutation(const std::vector<size_t> order) {
  std::vector<LineOcc> linesNew(_lines.size());
  for (size_t i = 0; i < order.size(); i++) {
    linesNew[i] = _lines[order[i]];
    _lineToIdx[_lines[order[i]].line] = i;
  }
  _lines = linesNew;
}

// _____________________________________________________________________________
size_t LineEdgePL::linePos(const Line* r) const {
  auto it = _lineToIdx.find(r);
  if (it == _lineToIdx.end()) return -1;
  return it->second;
}
