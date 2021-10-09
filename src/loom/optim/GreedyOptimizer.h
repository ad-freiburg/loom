// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_GREEDYOPTIMIZER_H_
#define LOOM_OPTIM_GREEDYOPTIMIZER_H_

#include <map>
#include "loom/config/LoomConfig.h"
#include "loom/optim/ExhaustiveOptimizer.h"
#include "loom/optim/ILPEdgeOrderOptimizer.h"
#include "loom/optim/NullOptimizer.h"
#include "loom/optim/OptGraph.h"
#include "loom/optim/OptGraphScorer.h"
#include "shared/rendergraph/OrderCfg.h"
#include "shared/rendergraph/RenderGraph.h"

namespace loom {
namespace optim {

typedef std::set<const OptEdge*> SettledEdgs;
typedef std::map<
    std::pair<const shared::linegraph::Line*, const shared::linegraph::Line*>,
    std::pair<bool, double>>
    Cmp;

struct LineCmp {
  LineCmp(const Cmp& cmp, bool rev) : _map(cmp), _rev(rev){};

  bool operator()(const shared::linegraph::Line* a,
                  const shared::linegraph::Line* b) {
    return _map.find({a, b})->second.first ^ _rev;
  }

  Cmp _map;
  bool _rev;
};

class GreedyOptimizer : public ExhaustiveOptimizer {
 public:
  GreedyOptimizer(const config::Config* cfg,
                  const shared::rendergraph::Penalties& pens, bool lookAhead)
      : ExhaustiveOptimizer(cfg, pens), _lookAhead(lookAhead){};

  virtual int optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                           shared::rendergraph::HierarOrderCfg* c,
                           size_t depth) const;

  void getFlatConfig(const std::set<OptNode*>& g,
                     OptOrderCfg* cfg) const;

 private:
  bool _lookAhead;

  const OptEdge* getNextEdge(const std::set<OptNode*>& g,
                             SettledEdgs* settled) const;
  const OptEdge* getInitialEdge(const std::set<OptNode*>& g) const;

  std::pair<bool, double> guess(const shared::linegraph::Line* a,
                                const shared::linegraph::Line* b,
                                const OptEdge* start, const OptNode* refNd,
                                const OptOrderCfg& cfg) const;
  std::pair<int, double> smallerThanAt(const shared::linegraph::Line* a,
                                       const shared::linegraph::Line* b,
                                       const OptEdge* e, const OptNode* nd,
                                       const OptEdge* ignore,
                                       const OptOrderCfg& cfg) const;

  const OptEdge* eligibleNextEdge(const OptEdge* start, const OptNode* nd,
                                  const shared::linegraph::Line* a,
                                  const shared::linegraph::Line* b) const;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_GREEDYOPTIMIZER_H_
