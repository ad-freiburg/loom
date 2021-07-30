// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_OPTGRAPHSCORER_H_
#define LOOM_OPTIM_OPTGRAPHSCORER_H_

#include <string>
#include <vector>
#include "loom/optim/OptGraph.h"

namespace loom {
namespace optim {

class OptGraphScorer {
 public:
  OptGraphScorer(const Scorer* scorer) : _scorer(scorer) {}

  double getCrossingScore(const std::set<OptNode*>& g,
                          const OptOrderCfg& c) const;
  double getCrossingScore(OptNode* n, const OptOrderCfg& c) const;

  double getCrossingScore(OptEdge* e, const OptOrderCfg& c) const;
  double getSplittingScore(OptEdge* e, const OptOrderCfg& c) const;

  double getSplittingScore(const std::set<OptNode*>& g,
                           const OptOrderCfg& c) const;
  double getSplittingScore(OptNode* n, const OptOrderCfg& c) const;

 private:
  const Scorer* _scorer;

  std::pair<size_t, size_t> getNumCrossings(OptNode* n,
                                            const OptOrderCfg& c) const;
  size_t getNumSeparations(OptNode* n, const OptOrderCfg& c) const;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_OPTGRAPHSCORER_H_
