// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/gridgraph/GridGraph.h"

using namespace octi::gridgraph;

// _____________________________________________________________________________
GridGraph::GridGraph(util::geo::Box bbox, double cellSize) : _bbox(bbox), _cellSize(cellSize) {

}
