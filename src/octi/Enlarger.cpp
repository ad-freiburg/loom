// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/Enlarger.h"
#include "octi/combgraph/CombGraph.h"
#include "util/geo/QuadTree.h"
#include "util/log/Log.h"

using octi::Enlarger;
using octi::combgraph::CombEdge;
using octi::combgraph::CombGraph;
using octi::combgraph::CombNode;
using util::geo::DBox;
using util::geo::QuadNode;
using util::geo::QuadTree;
using util::geo::QuadValue;

using util::DEBUG;

// _____________________________________________________________________________
void Enlarger::enlarge(CombGraph& cg, double theta) const {
  LOGTO(DEBUG, std::cerr) << "Enlarging input graph with theta=" << theta
                          << std::endl;

  auto bbox = util::geo::getBoundingRect(cg.getBBox());

  struct SplitFunc : util::geo::SplitFunc<const CombNode*, double> {
    SplitFunc(double theta) : _theta(theta) {}
    virtual bool operator()(
        const QuadNode<double>& nd,
        const QuadValue<const CombNode*, double>& newVal) const {
      UNUSED(newVal);
      double l = nd.bbox.getUpperRight().getX() - nd.bbox.getLowerLeft().getX();

      return nd.numEls > 0 && l > _theta;
    }
    double _theta;
  } sFunc(theta);

  QuadTree<const CombNode*, double> qt(1024, sFunc, bbox);

  // write nodes to quadtree
  for (auto cNd : cg.getNds()) {
    qt.insert(cNd, *cNd->pl().getGeom());
  }

  qt.print(std::cout);

  exit(1);
}
