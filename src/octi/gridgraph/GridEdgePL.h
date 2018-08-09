// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRIDGRAPH_GRIDEDGEPL_H_
#define OCTI_GRIDGRAPH_GRIDEDGEPL_H_

#include <set>
#include "octi/combgraph/CombEdgePL.h"
#include "octi/combgraph/CombNodePL.h"
#include "util/geo/GeoGraph.h"
#include "util/geo/PolyLine.h"
#include "util/graph/Node.h"

using util::geo::PolyLine;

namespace octi {
namespace gridgraph {

class GridEdgePL : util::geograph::GeoEdgePL<double> {
 public:
  GridEdgePL(const PolyLine<double>& p, double c, bool secondar);
  GridEdgePL(const PolyLine<double>& p, double c, bool secondar, bool closed);

  void addRoute(std::string);
  const std::set<util::graph::Edge<octi::combgraph::CombNodePL,
                                   octi::combgraph::CombEdgePL>*>&
  getResEdges() const;

  const util::geo::Line<double>* getGeom() const;
  util::json::Dict getAttrs() const;

  void setCost(double c);
  double cost() const;
  double rawCost() const;
  void addResidentEdge(util::graph::Edge<octi::combgraph::CombNodePL,
                                         octi::combgraph::CombEdgePL>* e);
  bool isSecondary() const;

  void close();
  void open();
  bool closed() const;
  void setVisited(int i);

 private:
  PolyLine<double> _pl;
  double _c;

  bool _isSecondary;
  bool _closed;

  std::set<util::graph::Edge<octi::combgraph::CombNodePL,
                             octi::combgraph::CombEdgePL>*>
      _resEdges;
  int _visited;
};
}
}

#endif  // OCTI_GRIDGRAPH_GRIDEDGEPL_H_
