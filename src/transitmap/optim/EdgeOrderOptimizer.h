// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

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

  void optimize();
 private:
  TransitGraph* _g;

  bool doOptimStep(Configuration* c);

  void applyConfig(const Configuration& c);

  std::vector<std::vector<size_t> > getPermutations(std::vector<size_t> order) const;

  void generateRandConfig(Configuration* c) const;
  void getConfig(Configuration* c) const;
};

}}

#endif  // TRANSITMAP_OPTIM_EDGEORDEROPTIMIZER_H_

