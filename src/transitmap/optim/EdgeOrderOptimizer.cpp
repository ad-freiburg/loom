// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include "./EdgeOrderOptimizer.h"


using namespace transitmapper;
using namespace optim;
using namespace graph;

// _____________________________________________________________________________
void EdgeOrderOptimizer::optimize() {
  

}

// _____________________________________________________________________________
void EdgeOrderOptimizer::doOptimeStep() {
  double currentScore = _g->getScore();

  for (graph::Node* n : _g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (const graph::EdgeTripGeom& g : *e->getEdgeTripGeoms()) {




      }
    }
  }
}
