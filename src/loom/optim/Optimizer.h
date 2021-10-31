// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "loom/config/LoomConfig.h"
#include "loom/optim/OptGraph.h"
#include "loom/optim/OptGraphScorer.h"
#include "shared/rendergraph/OrderCfg.h"
#include "shared/rendergraph/RenderGraph.h"

#ifndef LOOM_OPTIM_OPTIMIZER_H_
#define LOOM_OPTIM_OPTIMIZER_H_

namespace loom {
namespace optim {

typedef std::pair<const OptLO, const OptLO> LinePair;
typedef std::pair<size_t, size_t> PosCom;
typedef std::pair<PosCom, PosCom> PosComPair;
typedef std::pair<OptEdge*, OptEdge*> EdgePair;

struct OptResStats {
  size_t numNodesOrig, numStationsOrig, numEdgesOrig, maxLineCardOrig, numLinesOrig, maxDegOrig;
  size_t numStations, numNodes, numEdges, maxLineCard, nonTrivialComponents, numCompsSolSpaceOne, maxNumNodesPerComp, maxNumEdgesPerComp, maxCardPerComp, numCompsOrig, maxNumRowsPerComp, maxNumColsPerComp;
  size_t runs;
  double avgSolveTime, avgIterations, avgScore, avgCross, avgSameSegCross, avgDiffSegCross, avgSeps, solutionSpaceSize, solutionSpaceSizeOrig, maxCompSolSpace, simplificationTime;

  // best score for multiple runs
  size_t sameSegCrossings;
  size_t diffSegCrossings;
  size_t separations;
  double score;
};

class Optimizer {
 public:
  Optimizer(const config::Config* cfg,
            const shared::rendergraph::Penalties& pens)
      : _cfg(cfg), _scorer(pens){};

  virtual OptResStats optimize(shared::rendergraph::RenderGraph* rg) const;
  double optimizeComp(OptGraph* g, const std::set<OptNode*>& cmp,
                   shared::rendergraph::HierarOrderCfg* c,
                   OptResStats& stats) const;
  virtual double optimizeComp(OptGraph* g, const std::set<OptNode*>& cmp,
                           shared::rendergraph::HierarOrderCfg* c,
                           size_t depth, OptResStats& stats) const = 0;

  static std::vector<LinePair> getLinePairs(OptEdge* segment);
  static std::vector<LinePair> getLinePairs(OptEdge* segment, bool unique);

  static bool crosses(OptNode* node, OptEdge* segmentA, OptEdge* segmentB,
                      PosComPair postcomb);

  static bool crosses(OptNode* node, OptEdge* segmentA, EdgePair segments,
                      PosCom postcomb);

  static bool separates(PosComPair postcomb);

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
  const OptGraphScorer _scorer;

  static std::string prefix(size_t depth);

 private:
  static OptOrderCfg getOptOrderCfg(
      const shared::rendergraph::OrderCfg&,
      const std::map<const shared::linegraph::LineNode*, OptNode*>& ndMap,
      const OptGraph* g);
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_OPTIMIZER_H_
