// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_OPTGRAPHSCORER_H_
#define TRANSITMAP_OPTIM_OPTGRAPHSCORER_H_

#include <string>
#include <vector>
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/optim/OptGraph.h"
#include "transitmap/graph/Penalties.h"

using transitmapper::graph::Penalties;
using transitmapper::graph::OrderingConfig;

namespace transitmapper {
namespace optim {

class OptGraphScorer {
 public:
  OptGraphScorer(const Penalties& pens) : _pens(pens) {}

  double getScore(const std::set<OptNode*>& g, const OrderingConfig& c) const;

 private:
  Penalties _pens;
};
}
}

#endif  // TRANSITMAP_OPTIM_OPTGRAPHSCORER_H_
