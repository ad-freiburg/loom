// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_ILPEDGEORDEROPTIMIZER_H_
#define TRANSITMAP_OPTIM_ILPEDGEORDEROPTIMIZER_H_

#include <glpk.h>
#include "./../config/TransitMapConfig.h"
#include "./../graph/OrderingConfiguration.h"
#include "./../graph/Route.h"
#include "./../graph/TransitGraph.h"
#include "./ILPOptimizer.h"
#include "./OptGraph.h"
#include "./Optimizer.h"

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
  ILPEdgeOrderOptimizer(TransitGraph* g, const config::Config* cfg)
      : ILPOptimizer(g, cfg){};

 private:
  virtual glp_prob* createProblem(const OptGraph& g) const;

  virtual void getConfigurationFromSolution(glp_prob* lp, Configuration* c,
                                            const OptGraph& g) const;

  void writeCrossingOracle(const OptGraph& g, int* ia, int* ja, double* res,
                           size_t* c, glp_prob* lp) const;

  void writeDiffSegConstraintsImpr(const OptGraph& g, int* ia, int* ja,
                                   double* res, size_t* c, glp_prob* lp) const;
};
}  // namespace optim
}  // namespace transitmapper

#endif  // TRANSITMAP_OPTIM_ILPEDGEORDEROPTIMIZER_H_
