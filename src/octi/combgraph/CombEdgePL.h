// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_COMBGRAPH_COMBEDGEPL_H_
#define OCTI_COMBGRAPH_COMBEDGEPL_H_

#include <set>
#include "shared/linegraph/LineGraph.h"
#include "util/geo/GeoGraph.h"

using util::geo::PolyLine;

namespace octi {
namespace combgraph {

class CombEdgePL : util::geograph::GeoEdgePL<double> {
 public:
  CombEdgePL(shared::linegraph::LineEdge* child);

  const util::geo::Line<double>* getGeom() const;
  util::json::Dict getAttrs() const;

  const std::vector<shared::linegraph::LineEdge*>& getChilds() const;
  std::vector<shared::linegraph::LineEdge*>& getChilds();

  const PolyLine<double>& getPolyLine() const;
  void setPolyLine(const PolyLine<double>& p);

  size_t getNumLines() const { return _maxLineNum; }
  void setNumLines(size_t numLines) { _maxLineNum = numLines; }

 private:
  std::vector<shared::linegraph::LineEdge*> _childs;

  size_t _maxLineNum;

  PolyLine<double> _geom;
};
}
}

#endif  // OCTI_COMBGRAPH_COMBEDGEPL_H_
