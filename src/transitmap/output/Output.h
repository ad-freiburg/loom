// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OUTPUT_OUTPUT_H_
#define TRANSITMAP_OUTPUT_OUTPUT_H_

#include "../graph/TransitGraph.h"

namespace transitmapper {
namespace output {

class Output {
 public:
  virtual ~Output() {};

  // print the outputGraph
  virtual void print(const graph::TransitGraph& outputGraph) = 0;
};

}}

#endif  // TRANSITMAP_OUTPUT_OUTPUT_H_
