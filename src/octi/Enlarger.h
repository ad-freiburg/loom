// Copyright 2020, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_ENLARGER_H_
#define OCTI_ENLARGER_H_

#include "shared/linegraph/LineGraph.h"
#include "octi/combgraph/CombGraph.h"

namespace octi {

class Enlarger {
 public:
  Enlarger() {}
  void enlarge(octi::combgraph::CombGraph& lg, double theta) const;
};

}  // namespace octi

#endif  // OCTI_ENLARGER_H_
