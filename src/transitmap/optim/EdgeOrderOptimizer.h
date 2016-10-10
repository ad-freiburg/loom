// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_EDGEORDEROPTIMIZER_H_
#define TRANSITMAP_OPTIM_EDGEORDEROPTIMIZER_H_

#include "./../graph/OrderingConfiguration.h"
#include "./../graph/TransitGraph.h"

using std::exception;
using std::string;

namespace transitmapper {
namespace optim {

using namespace graph;

class EdgeOrderOptimizer {
 public:
  EdgeOrderOptimizer(TransitGraph* g) : _g(g) {};

  void optimize(size_t numRuns);
 private:
  TransitGraph* _g;

  bool doOptimStep(Configuration* c);

  void applyConfig(const Configuration& c);

  void getPermutations(size_t n, std::vector<std::vector<size_t> >* ret) const;

  void generateRandConfig(Configuration* c) const;
};

}}

#endif  // TRANSITMAP_OPTIM_EDGEORDEROPTIMIZER_H_

