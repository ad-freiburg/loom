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
  OptGraphScorer(const shared::rendergraph::Penalties& pens) : _pens(pens) {}

  double getTotalScore(const OptGraph* g, const OptOrderCfg& c) const;
  double getTotalScore(const std::set<OptNode*>& g, const OptOrderCfg& c) const;
  double getTotalScore(OptNode* n, const OptOrderCfg& c) const;
  double getTotalScore(OptEdge* n, const OptOrderCfg& c) const;

  double getCrossingScore(const OptGraph* g, const OptOrderCfg& c) const;
  double getCrossingScore(const std::set<OptNode*>& g,
                          const OptOrderCfg& c) const;
  double getCrossingScore(OptNode* n, const OptOrderCfg& c) const;
  double getCrossingScore(OptEdge* e, const OptOrderCfg& c) const;

  double getSeparationScore(const OptGraph* g, const OptOrderCfg& c) const;

  double getSeparationScore(OptEdge* e, const OptOrderCfg& c) const;

  double getSeparationScore(const std::set<OptNode*>& g,
                           const OptOrderCfg& c) const;
  double getSeparationScore(OptNode* n, const OptOrderCfg& c) const;

  std::pair<std::pair<size_t, size_t>, size_t> getNumCrossSeps(
      OptNode* n, const OptOrderCfg& c) const;

  std::pair<size_t, size_t> getNumCrossings(OptNode* n,
                                            const OptOrderCfg& c) const;

  size_t getNumSeparations(OptNode* n, const OptOrderCfg& c) const;

  std::pair<size_t, size_t> getNumCrossings(const OptGraph* g,
                                            const OptOrderCfg& c) const;

  std::pair<size_t, size_t> getNumCrossings(const std::set<OptNode*>& g,
                                            const OptOrderCfg& c) const;

  size_t getNumSeparations(const OptGraph* g, const OptOrderCfg& c) const;
  size_t getNumSeparations(const std::set<OptNode*>& g,
                           const OptOrderCfg& c) const;

  bool optimizeSep() const;

  double getSeparationPen(const OptNode* n) const;
  double getCrossingPenSameSeg(const OptNode* n) const;
  double getCrossingPenDiffSeg(const OptNode* n) const;

  const shared::rendergraph::Penalties& getPens() const { return _pens;}

 private:
  shared::rendergraph::Penalties _pens;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_OPTGRAPHSCORER_H_
