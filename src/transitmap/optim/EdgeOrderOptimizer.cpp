// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include "./EdgeOrderOptimizer.h"


using namespace transitmapper;
using namespace optim;
using namespace graph;

// _____________________________________________________________________________
void EdgeOrderOptimizer::optimize() {
  srand(time(0));

  Configuration bestConfig;
  double bestScore = DBL_MAX;
  for (size_t i = 0; i < 100; i++) {
    std::cout << "Round " << i << std::endl;
    Configuration c;
    generateRandConfig(&c);
    applyConfig(c);

    while (doOptimStep(&c)) {};

    double newScore = _g->getScore();
    if (newScore < bestScore) {
      bestScore = newScore;
      bestConfig = c;
    }
  }

  applyConfig(bestConfig);
}

// _____________________________________________________________________________
bool EdgeOrderOptimizer::doOptimStep(Configuration* c) {
  double currentScore = _g->getScore();
  std::pair<EdgeTripGeom*, std::vector<size_t> > bestCand;

  for (graph::Node* n : *_g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (graph::EdgeTripGeom& g : *e->getEdgeTripGeoms()) {
        if (g.getTripOrdering().size() == 1) continue;
        Ordering origOrdering = g.getTripOrdering();
        std::vector<Ordering > permutations = getPermutations(origOrdering);

        /**
        //take random permutation
        size_t r = rand() % permutations.size();
        auto& perm = permutations[r];
        **/

        for (auto& perm : permutations) {
          g.setTripOrdering(perm);
          double newScore = _g->getScore();
          if (newScore < currentScore) {
            bestCand = std::pair<EdgeTripGeom*, Ordering>(&g, perm);
            currentScore = newScore;
          }
        }

        // set the ordering back
        g.setTripOrdering(origOrdering);
      }
    }
  }

  if (currentScore < _g->getScore()) {
    bestCand.first->setTripOrdering(bestCand.second);
    (*c)[bestCand.first] = bestCand.second;
    return true;
  } else {
    return false;
  }
}

// _____________________________________________________________________________
void EdgeOrderOptimizer::generateRandConfig(Configuration* c) const {
  for (graph::Node* n : *_g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (graph::EdgeTripGeom& g : *e->getEdgeTripGeoms()) {
        if (g.getTripOrdering().size() == 1) continue;
        Ordering origOrdering = g.getTripOrdering();
        std::vector<Ordering > permutations = getPermutations(origOrdering);
        size_t i = rand() % permutations.size();
        (*c)[&g] = permutations[i];
      }
    }
  }
}

// _____________________________________________________________________________
std::vector<std::vector<size_t> >
EdgeOrderOptimizer::getPermutations(std::vector<size_t> order)
const {
  std::sort(order.begin(), order.end());
  std::vector<std::vector<size_t> > ret;

  do {
    ret.push_back(order);
  } while (std::next_permutation(order.begin(), order.end()));

  return ret;
}

// _____________________________________________________________________________
void EdgeOrderOptimizer::applyConfig(const Configuration& c) {
  for (auto a : c) {
    a.first->setTripOrdering(a.second);
  }
}
