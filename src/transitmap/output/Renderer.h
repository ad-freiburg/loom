// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OUTPUT_RENDERER_H_
#define TRANSITMAP_OUTPUT_RENDERER_H_

#include "transitmap/graph/RenderGraph.h"

namespace transitmapper {
namespace output {

class Renderer {
 public:
  virtual ~Renderer() {};

  // print the outputGraph
  virtual void print(const graph::RenderGraph& outG) = 0;
};

}}

#endif  // TRANSITMAP_OUTPUT_RENDERER_H_
