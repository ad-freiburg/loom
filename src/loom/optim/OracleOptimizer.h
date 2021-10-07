// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_ORACLEOPTIMIZER_H_
#define LOOM_OPTIM_ORACLEOPTIMIZER_H_

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

struct LineCmp {
  LineCmp(const OptEdge* e, const OptOrderCfg& cfg,
          const OptGraphScorer& scorer)
      : _e(e), _cfg(cfg), _scorer(scorer){};

  bool operator()(const shared::linegraph::Line* a,
                  const shared::linegraph::Line* b) {
    auto leftPair = smallerThanAt(a, b, _e, _e->getFrom(), 0);
    auto rightPair = smallerThanAt(a, b, _e, _e->getTo(), 0);

    int left = leftPair.first;
    int right = rightPair.first;

    double leftCost = leftPair.second;
    double rightCost = rightPair.second;

    bool rev = !_e->pl().lnEdgParts.front().dir;

    if (left != 0 && right != 0) {
      if (left == -right) {
        // both sides induce the same, and are not undecided
        if (left == 1) return false ^ rev;
        return true ^ rev;
      } else {
        // both sides induce different orderings, take the order of the side
        // where a crossing is more expensive
        if (leftCost > rightCost) {
          if (left == 1) return false ^ rev;
          return true ^ rev;
        } else {
          if (right == 1) return true ^ rev;
          return false ^ rev;
        }
      }
    } else if (left != 0 && right == 0) {
      // left side induces
      if (left == 1) return false ^ rev;
      return true ^ rev;
    } else if (left == 0 && right != 0) {
      // right side induces
      if (right == 1) return true ^ rev;
      return false ^ rev;
    } else
      return oracle(a, b);
  }

  std::pair<int, double> smallerThanAt(const shared::linegraph::Line* a,
                                       const shared::linegraph::Line* b,
                                       const OptEdge* e,
                                       const OptNode* nd,
                                       const OptEdge* ignore) const;

  bool oracle(const shared::linegraph::Line* a,
              const shared::linegraph::Line* b) const;

  const OptEdge* eligibleNextEdge(const OptEdge* start, const OptNode* nd,
                                  const shared::linegraph::Line* a,
                                  const shared::linegraph::Line* b) const;

  const OptEdge* _e;
  const OptOrderCfg& _cfg;
  const OptGraphScorer& _scorer;
};

class OracleOptimizer : public ExhaustiveOptimizer {
 public:
  OracleOptimizer(const config::Config* cfg,
                  const shared::rendergraph::Penalties& pens)
      : ExhaustiveOptimizer(cfg, pens){};

  virtual int optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                           shared::rendergraph::HierarOrderCfg* c,
                           size_t depth) const;

 private:
  const OptEdge* getNextEdge(const std::set<OptNode*>& g,
                             SettledEdgs* settled) const;
  const OptEdge* getInitialEdge(const std::set<OptNode*>& g) const;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_ORACLEOPTIMIZER_H_
