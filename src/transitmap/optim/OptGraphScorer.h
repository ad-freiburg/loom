// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_OPTGRAPHSCORER_H_
#define TRANSITMAP_OPTIM_OPTGRAPHSCORER_H_

#include <string>
#include <vector>
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/Penalties.h"
#include "transitmap/optim/OptGraph.h"

using transitmapper::graph::Penalties;
using transitmapper::graph::OrderingConfig;

namespace transitmapper {
namespace optim {

class OptGraphScorer {
 public:
  OptGraphScorer(const Scorer* scorer) : _scorer(scorer) {}

  double getCrossingScore(OptGraph* og, const std::set<OptNode*>& g,
                  const OptOrderingConfig& c) const;
  double getCrossingScore(OptGraph* og, OptNode* n, const OptOrderingConfig& c) const;

  double getCrossingScore(OptGraph* og, OptEdge* e,
                  const OptOrderingConfig& c) const;
  double getSplittingScore(OptGraph* og, OptEdge* e, const OptOrderingConfig& c) const;

  double getSplittingScore(OptGraph* og, const std::set<OptNode*>& g,
                  const OptOrderingConfig& c) const;
  double getSplittingScore(OptGraph* og, OptNode* n, const OptOrderingConfig& c) const;

 private:
  const Scorer* _scorer;

  std::pair<size_t, size_t> getNumCrossings(OptGraph* og, OptNode* n,
                                            const OptOrderingConfig& c) const;
  size_t getNumSeparations(OptGraph* og, OptNode* n,
                                            const OptOrderingConfig& c) const;
};
}
}

#endif  // TRANSITMAP_OPTIM_OPTGRAPHSCORER_H_
