// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_ILPEDGEORDEROPTIMIZER_H_
#define TRANSITMAP_OPTIM_ILPEDGEORDEROPTIMIZER_H_

#include <glpk.h>
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/Route.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/optim/ILPOptimizer.h"
#include "transitmap/optim/OptGraph.h"
#include "transitmap/optim/Optimizer.h"
#include "transitmap/optim/Scorer.h"

using std::exception;
using std::string;

namespace transitmapper {
namespace optim {

using namespace graph;

typedef std::pair<const Route*, const Route*> LinePair;
typedef std::pair<size_t, size_t> PosCom;
typedef std::pair<PosCom, PosCom> PosComPair;
typedef std::pair<OptEdge*, OptEdge*> EdgePair;

class ILPEdgeOrderOptimizer : public ILPOptimizer {
 public:
  ILPEdgeOrderOptimizer(TransitGraph* g, const config::Config* cfg,
                        const Scorer* scorer)
      : ILPOptimizer(g, cfg, scorer){};

 private:
  virtual glp_prob* createProblem(const std::set<OptNode*>& g) const;

  virtual void getConfigurationFromSolution(glp_prob* lp, OrderingConfig* c,
                                            const std::set<OptNode*>& g) const;

  void writeCrossingOracle(const std::set<OptNode*>& g, VariableMatrix* vm,
                           glp_prob* lp) const;

  void writeDiffSegConstraintsImpr(const std::set<OptNode*>& g, VariableMatrix* vm,
                                   glp_prob* lp) const;
};
}  // namespace optim
}  // namespace transitmapper

#endif  // TRANSITMAP_OPTIM_ILPEDGEORDEROPTIMIZER_H_
