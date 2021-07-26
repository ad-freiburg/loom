// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_RENDERGRAPH_PENALTIES_H_
#define SHARED_RENDERGRAPH_PENALTIES_H_

namespace transitmapper {
namespace graph {

struct Penalties {
  double inStatCrossPenDegTwo;
  double inStatSplitPenDegTwo;
  double sameSegCrossPen;
  double diffSegCrossPen;
  double splitPen;
  double inStatCrossPenSameSeg;
  double inStatCrossPenDiffSeg;
  double inStatSplitPen;
  bool crossAdjPen;
  bool splitAdjPen;
};

const Penalties IDENTITY_PENALTIES = Penalties {1, 1, 1, 1, 1, 1, 1, 1, 0, 0};

}  // namespace rendergraph
}  // namespace shared

#endif  // SHARED_RENDERGRAPH_PENALTIES_H_
