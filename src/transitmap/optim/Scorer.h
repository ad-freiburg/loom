// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_SCORER_H_
#define TRANSITMAP_OPTIM_SCORER_H_

#include <string>
#include <vector>
#include "transitmap/graph/OrderCfg.h"
#include "transitmap/graph/Penalties.h"
#include "transitmap/graph/RenderGraph.h"

namespace transitmapper {
namespace optim {

class Scorer {
 public:
  Scorer(const graph::RenderGraph* g, const graph::Penalties& pens)
      : _g(g), _pens(pens) {}

  const graph::RenderGraph* getGraph() const { return _g; }

  double getScore() const;
  double getScore(const graph::OrderCfg& c) const;

  double getCrossScore() const;
  double getCrossScore(const graph::OrderCfg& c) const;

  double getSeparationScore() const;
  double getSeparationScore(const graph::OrderCfg& c) const;

  size_t getNumCrossings() const;
  size_t getNumCrossings(const graph::OrderCfg& c) const;

  size_t getNumSeparations() const;
  size_t getNumSeparations(const graph::OrderCfg& c) const;

  double getNumPossSolutions() const;
  size_t getMaxCrossPenalty() const;
  size_t getMaxSplitPenalty() const;

  double getScore(const shared::linegraph::LineNode* n,
                  const graph::OrderCfg& cfg) const;
  size_t getNumCrossings(const shared::linegraph::LineNode* n,
                         const graph::OrderCfg& c) const;
  size_t getNumSeparations(const shared::linegraph::LineNode* n,
                           const graph::OrderCfg& c) const;
  double getSeparationScore(const shared::linegraph::LineNode* n,
                            const graph::OrderCfg& c,
                            const graph::Penalties& pens) const;
  double getCrossingScore(const shared::linegraph::LineNode* n,
                          const graph::OrderCfg& c,
                          const graph::Penalties& pens) const;

  double getSeparationScore(const shared::linegraph::LineNode* n,
                            const graph::OrderCfg& c) const;
  double getCrossingScore(const shared::linegraph::LineNode* n,
                          const graph::OrderCfg& c) const;

  int getCrossingPenaltySameSeg(const shared::linegraph::LineNode* n,
                                const graph::Penalties& pens) const;
  int getCrossingPenaltyDiffSeg(const shared::linegraph::LineNode* n,
                                const graph::Penalties& pens) const;
  int getSplittingPenalty(const shared::linegraph::LineNode* n,
                          const graph::Penalties& pens) const;

  int getCrossingPenaltySameSeg(const shared::linegraph::LineNode* n) const;
  int getCrossingPenaltyDiffSeg(const shared::linegraph::LineNode* n) const;
  int getSplittingPenalty(const shared::linegraph::LineNode* n) const;

 private:
  const graph::RenderGraph* _g;
  graph::Penalties _pens;
};
}
}

#endif  // TRANSITMAP_GRAPH_ROUTE_H_
