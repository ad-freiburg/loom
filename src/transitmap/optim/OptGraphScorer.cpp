// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/OptGraph.h"
#include "transitmap/graph/Penalties.h"
#include "transitmap/optim/OptGraphScorer.h"

using namespace transitmapper;
using namespace optim;
using transitmapper::graph::TransitGraph;
using transitmapper::graph::Node;
using transitmapper::graph::Edge;
using transitmapper::graph::Route;
using transitmapper::graph::InnerGeometry;
using transitmapper::graph::IDENTITY_PENALTIES;

// _____________________________________________________________________________
double OptGraphScorer::getScore(const std::set<OptNode*>& g,
                                const OrderingConfig& c) const {
  double ret = 0;

  for (auto n : g) {
    ret += getScore(n, c);
  }

  return ret;
}
