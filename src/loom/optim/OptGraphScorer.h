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

  double getCrossingScore(const OptGraph* g, const OptOrderCfg& c) const;
  double getCrossingScore(const std::set<OptNode*>& g,
                          const OptOrderCfg& c) const;
  double getCrossingScore(OptNode* n, const OptOrderCfg& c) const;

  double getCrossingScore(OptEdge* e, const OptOrderCfg& c) const;

  double getSplittingScore(const OptGraph* g, const OptOrderCfg& c) const;

  double getSplittingScore(OptEdge* e, const OptOrderCfg& c) const;

  double getSplittingScore(const std::set<OptNode*>& g,
                           const OptOrderCfg& c) const;
  double getSplittingScore(OptNode* n, const OptOrderCfg& c) const;

  std::pair<size_t, size_t> getNumCrossings(OptNode* n,
                                            const OptOrderCfg& c) const;

  size_t getNumSeparations(OptNode* n, const OptOrderCfg& c) const;

  std::pair<size_t, size_t> getNumCrossings(const OptGraph* g,
                                            const OptOrderCfg& c) const;

  size_t getNumSeparations(const OptGraph* g, const OptOrderCfg& c) const;

 private:
  const Scorer* _scorer;

};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_OPTGRAPHSCORER_H_
