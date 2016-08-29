// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_EDGEORDEROPTIMIZER_H_
#define TRANSITMAP_OPTIM_EDGEORDEROPTIMIZER_H_

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

  void doOptimStep();
};

}}

#endif  // TRANSITMAP_OPTIM_EDGEORDEROPTIMIZER_H_

