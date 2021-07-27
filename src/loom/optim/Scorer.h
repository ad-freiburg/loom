// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_SCORER_H_
#define LOOM_OPTIM_SCORER_H_

#include <string>
#include <vector>
#include "shared/rendergraph/OrderCfg.h"
#include "shared/rendergraph/Penalties.h"
#include "shared/rendergraph/RenderGraph.h"

namespace loom {
namespace optim {

class Scorer {
 public:
  Scorer(const shared::rendergraph::RenderGraph* g,
         const shared::rendergraph::Penalties& pens)
      : _g(g), _pens(pens) {}

  const shared::rendergraph::RenderGraph* getGraph() const { return _g; }

  double getScore(const shared::rendergraph::OrderCfg& c) const;

  double getCrossScore(const shared::rendergraph::OrderCfg& c) const;

  double getSeparationScore(const shared::rendergraph::OrderCfg& c) const;

  size_t getNumCrossings(const shared::rendergraph::OrderCfg& c) const;

  size_t getNumSeparations(const shared::rendergraph::OrderCfg& c) const;

  double getNumPossSolutions() const;
  size_t getMaxCrossPenalty() const;
  size_t getMaxSplitPenalty() const;

  double getScore(const shared::linegraph::LineNode* n,
                  const shared::rendergraph::OrderCfg& cfg) const;
  size_t getNumCrossings(const shared::linegraph::LineNode* n,
                         const shared::rendergraph::OrderCfg& c) const;
  size_t getNumSeparations(const shared::linegraph::LineNode* n,
                           const shared::rendergraph::OrderCfg& c) const;
  double getSeparationScore(const shared::linegraph::LineNode* n,
                            const shared::rendergraph::OrderCfg& c,
                            const shared::rendergraph::Penalties& pens) const;
  double getCrossingScore(const shared::linegraph::LineNode* n,
                          const shared::rendergraph::OrderCfg& c,
                          const shared::rendergraph::Penalties& pens) const;

  double getSeparationScore(const shared::linegraph::LineNode* n,
                            const shared::rendergraph::OrderCfg& c) const;
  double getCrossingScore(const shared::linegraph::LineNode* n,
                          const shared::rendergraph::OrderCfg& c) const;

  int getCrossingPenaltySameSeg(
      const shared::linegraph::LineNode* n,
      const shared::rendergraph::Penalties& pens) const;
  int getCrossingPenaltyDiffSeg(
      const shared::linegraph::LineNode* n,
      const shared::rendergraph::Penalties& pens) const;
  int getSplittingPenalty(const shared::linegraph::LineNode* n,
                          const shared::rendergraph::Penalties& pens) const;

  int getCrossingPenaltySameSeg(const shared::linegraph::LineNode* n) const;
  int getCrossingPenaltyDiffSeg(const shared::linegraph::LineNode* n) const;
  int getSplittingPenalty(const shared::linegraph::LineNode* n) const;

 private:
  const shared::rendergraph::RenderGraph* _g;
  shared::rendergraph::Penalties _pens;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_GRAPH_ROUTE_H_
