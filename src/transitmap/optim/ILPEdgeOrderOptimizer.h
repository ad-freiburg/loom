// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_ILPEDGEORDEROPTIMIZER_H_
#define TRANSITMAP_OPTIM_ILPEDGEORDEROPTIMIZER_H_

#include "./../graph/OrderingConfiguration.h"
#include "./../graph/TransitGraph.h"
#include <glpk.h>

using std::exception;
using std::string;

namespace transitmapper {
namespace optim {

using namespace graph;

class ILPEdgeOrderOptimizer {
 public:
  ILPEdgeOrderOptimizer(TransitGraph* g) : _g(g) {};

  void optimize();
 private:
  TransitGraph* _g;

  glp_prob* createProblem() const;
  void solveProblem(glp_prob* lp) const;

  void getConfigurationFromSoluation(glp_prob* lp,
      Configuration* c) const;
  std::string getILPVarName(const EdgeTripGeom* seg,
    const gtfs::Route* r, size_t p) const;
};

}}

#endif  // TRANSITMAP_OPTIM_ILPEDGEORDEROPTIMIZER_H_

