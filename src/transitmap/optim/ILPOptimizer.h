// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_ILPOPTIMIZER_H_
#define TRANSITMAP_OPTIM_ILPOPTIMIZER_H_

#include <glpk.h>
#include "./../config/TransitMapConfig.h"
#include "./../graph/OrderingConfiguration.h"
#include "./../graph/Route.h"
#include "./../graph/TransitGraph.h"
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

class ILPOptimizer : public Optimizer {
 public:
  ILPOptimizer(TransitGraph* g, const config::Config* cfg)
      : _g(g), _cfg(cfg){};

  int optimize() const;

 protected:
  TransitGraph* _g;
  const config::Config* _cfg;

  virtual glp_prob* createProblem(const OptGraph& g) const;

  void solveProblem(glp_prob* lp) const;

  virtual void getConfigurationFromSolution(glp_prob* lp, Configuration* c,
                                    const OptGraph& g) const;

  std::string getILPVarName(OptEdge* e, const Route* r, size_t p) const;

  void writeSameSegConstraints(const OptGraph& g, int* ia, int* ja, double* res,
                               size_t* c, glp_prob* lp) const;

  void writeDiffSegConstraints(const OptGraph& g, int* ia, int* ja, double* res,
                               size_t* c, glp_prob* lp) const;

  bool printHumanReadable(glp_prob* lp, const std::string& path) const;
  double getConstraintCoeff(glp_prob* lp, size_t constraint, size_t col) const;

  std::vector<LinePair> getLinePairs(OptEdge* segment) const;
  std::vector<OptEdge*> getEdgePartners(OptNode* node, OptEdge* segmentA,
                                        const LinePair& linepair) const;
  std::vector<EdgePair> getEdgePartnerPairs(OptNode* node, OptEdge* segmentA,
                                            const LinePair& linepair) const;
  std::vector<PosComPair> getPositionCombinations(OptEdge* a, OptEdge* b) const;
  std::vector<PosCom> getPositionCombinations(OptEdge* a) const;

  bool crosses(OptNode* node, OptEdge* segmentA, OptEdge* segmentB,
               PosComPair postcomb) const;

  bool crosses(OptNode* node, OptEdge* segmentA, EdgePair segments,
               PosCom postcomb) const;

  Point getPos(OptNode* n, OptEdge* segment, size_t p) const;
};
}
}

#endif  // TRANSITMAP_OPTIM_ILPOPTIMIZER_H_
