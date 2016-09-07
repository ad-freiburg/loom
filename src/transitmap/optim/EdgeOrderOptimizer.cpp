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
  double bestScore = _g->getScore();

  for (size_t i = 0; i < 100; i++) {
    std::cout << "Round " << i << std::endl;
    Configuration c;
    getConfig(&c);

    size_t s = 0;
    while (doOptimStep(&c)) {
      s++;
      std::cout << "Step: " << s << std::endl;
    };

    applyConfig(c);
    double newScore = _g->getScore();
    if (newScore < bestScore) {
      std::cout << "New best score: " << newScore << std::endl;
      bestScore = newScore;
      bestConfig = c;
    }

    generateRandConfig(&c);
    applyConfig(c);
  }

  applyConfig(bestConfig);
}

// _____________________________________________________________________________
bool EdgeOrderOptimizer::doOptimStep(Configuration* c) {
  std::pair<EdgeTripGeom*, std::vector<size_t> > bestCand;
  bestCand.first = 0;

  for (graph::Node* n : *_g->getNodes()) {
    // the area score is the score of the node plus the score of all its
    // neighbors
    for (graph::Edge* e : n->getAdjListOut()) {
      for (graph::EdgeTripGeom& g : *e->getEdgeTripGeoms()) {
        if (g.getTripOrdering().size() == 1) continue;

        double oldAreaScore = n->getAreaScore();
        Ordering origOrdering = g.getTripOrdering();
        std::vector<Ordering > permutations = getPermutations(origOrdering);

        #pragma omp parallel
        for (size_t i = 0; i < permutations.size(); ++i) {
          double newAreaScore = n->getAreaScore(g, permutations[i]);
          if (newAreaScore < oldAreaScore) {
            std::cout << "TEST" << std::endl;
            #pragma omp critical
            {
              bestCand = std::pair<EdgeTripGeom*, Ordering>(&g, permutations[i]);
              oldAreaScore = newAreaScore;
            }
          }
        }
      }
    }
  }

  if (bestCand.first) {
    bestCand.first->setTripOrdering(bestCand.second);
    (*c)[bestCand.first] = bestCand.second;
    return true;
  } else {
    // no improvement could be made
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
void EdgeOrderOptimizer::getConfig(Configuration* c) const {
  for (graph::Node* n : *_g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (graph::EdgeTripGeom& g : *e->getEdgeTripGeoms()) {
        if (g.getTripOrdering().size() == 1) continue;
        Ordering origOrdering = g.getTripOrdering();
        (*c)[&g] = origOrdering;
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
