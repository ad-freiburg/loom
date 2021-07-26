// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "loom/config/TransitMapConfig.h"
#include "loom/graph/OrderCfg.h"
#include "loom/graph/RenderGraph.h"
#include "loom/optim/OptGraph.h"

#ifndef LOOM_OPTIM_OPTIMIZER_H_
#define LOOM_OPTIM_OPTIMIZER_H_

namespace loom {
namespace optim {

typedef std::pair<const OptLO, const OptLO> LinePair;
typedef std::pair<size_t, size_t> PosCom;
typedef std::pair<PosCom, PosCom> PosComPair;
typedef std::pair<OptEdge*, OptEdge*> EdgePair;

class Optimizer {
 public:
  Optimizer(const config::Config* cfg, const Scorer* scorer)
      : _cfg(cfg), _scorer(scorer){};

  virtual int optimize(graph::RenderGraph* tg) const;
  int optimizeComp(OptGraph* g, const std::set<OptNode*>& cmp,
                   graph::HierarOrderCfg* c) const;
  virtual int optimizeComp(OptGraph* g, const std::set<OptNode*>& cmp,
                           graph::HierarOrderCfg* c, size_t depth) const = 0;

  static std::vector<LinePair> getLinePairs(OptEdge* segment);
  static std::vector<LinePair> getLinePairs(OptEdge* segment, bool unique);

  static bool crosses(OptGraph* g, OptNode* node, OptEdge* segmentA,
                      OptEdge* segmentB, PosComPair postcomb);

  static bool crosses(OptGraph* g, OptNode* node, OptEdge* segmentA,
                      EdgePair segments, PosCom postcomb);
  static util::geo::DPoint getPos(OptGraph* g, OptNode* n, OptEdge* segment,
                                  size_t p);

  static std::vector<EdgePair> getEdgePartnerPairs(OptNode* node,
                                                   OptEdge* segmentA,
                                                   const LinePair& linepair);

  static std::vector<OptEdge*> getEdgePartners(OptNode* node, OptEdge* segmentA,
                                               const LinePair& linepair);
  static size_t maxCard(const std::set<OptNode*>& g);
  static double solutionSpaceSize(const std::set<OptNode*>& g);
  static double numEdges(const std::set<OptNode*>& g);

 protected:
  const config::Config* _cfg;
  const Scorer* _scorer;

  static std::string prefix(size_t depth);
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_OPTIMIZER_H_
