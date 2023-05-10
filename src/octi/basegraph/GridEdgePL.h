// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.  // Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_GRIDEDGEPL_H_
#define OCTI_BASEGRAPH_GRIDEDGEPL_H_

#include <set>
#include "octi/combgraph/CombEdgePL.h"
#include "octi/combgraph/CombNodePL.h"
#include "util/geo/GeoGraph.h"
#include "util/geo/PolyLine.h"
#include "util/graph/Node.h"

using util::geo::PolyLine;

namespace octi {
namespace basegraph {

class GridEdgePL : util::geograph::GeoEdgePL<double> {
 public:
  GridEdgePL(double c, bool secondar, bool sink);

  const util::geo::Line<double>* getGeom() const;
  util::json::Dict getAttrs() const;

  void setCost(double c);
  double cost() const;
  double rawCost() const;
  bool isSecondary() const;

  void close();
  void softClose();
  void open();
  bool closed() const;

  size_t resEdgs() const;

  void block();
  void unblock();

  void reset();

  void delResEdg();
  void addResEdge();

  void setId(size_t id);
  size_t getId() const;

 private:
  float _c;

  bool _isSecondary : 1;
  bool _isSink : 1;

  bool _closed : 1;
  bool _softClosed : 1;

  // edges are blocked if they would cross a settled edge
  bool _blocked : 1;

  uint8_t _resEdgs : 8;

  uint32_t _id;
};
}
}

#endif  // OCTI_BASEGRAPH_GRIDEDGEPL_H_
