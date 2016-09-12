// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include "easylogging/easylogging.h"
#include "./EdgeOrderOptimizer.h"


using namespace transitmapper;
using namespace optim;
using namespace graph;

// _____________________________________________________________________________
void EdgeOrderOptimizer::optimize() {
  srand(time(0));

  Configuration startConfig;
  getConfig(&startConfig);

  size_t numRuns = 4;

  std::vector<std::pair<double, Configuration> > results;
  results.resize(numRuns);

  for (size_t i = 0; i < numRuns; i++) {
    LOG(INFO) << "Round " << i;

    Configuration c;

    if (i == 0) {
      c = startConfig;
    } else {
      generateRandConfig(&c);
    }

    size_t step = 0;

    while (doOptimStep(&c)) {
      std::cout << step++ << std::endl;
    };

    double newScore = _g->getScore(c);
    LOG(INFO) << newScore;
    results[i] = std::pair<double, Configuration>(newScore, c);
  }

  // take best one
  size_t best;
  double bestScore = DBL_MAX;
  for (size_t i = 0; i < numRuns; i++) {
    if (results[i].first < bestScore) {
      best = i;
      bestScore = results[i].first;
    }
  }

  std::cout << "result: " << best << " score: " << bestScore << std::endl;
  applyConfig(results[best].second);
}

// _____________________________________________________________________________
bool EdgeOrderOptimizer::doOptimStep(Configuration* c) {
  std::pair<EdgeTripGeom*, std::vector<size_t> > bestCand;
  bestCand.first = 0;

  bool simple = false; // TRUE for simple, false for steepest hill climbing
  bool abort = false;

  double bestAreaScoreImprov = 0;

  for (graph::Node* n : *_g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (graph::EdgeTripGeom& g : *e->getEdgeTripGeoms()) {
        if (abort) goto exitloop;
        if (g.getTripOrdering().size() == 1) continue;

        double oldAreaScore = n->getAreaScore(*c);
        std::vector<Ordering > permutations;
        getPermutations(g.getTripOrdering(), &permutations);

        #pragma omp parallel for
        for (size_t i = 0; i < permutations.size(); ++i) {
          if (!abort) {
            double newAreaScore = n->getAreaScore(*c, &g, &permutations[i]);

            double diff = oldAreaScore - newAreaScore;

            if (diff - bestAreaScoreImprov > 0.000001) {
              #pragma omp critical
              {
                bestCand = std::pair<EdgeTripGeom*, Ordering>(&g, permutations[i]);
                bestAreaScoreImprov = diff;
              }

              if (simple) {
                abort = true;
              }
            }
          }
        }
      }
    }
  }

exitloop:

  if (bestCand.first) {
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
        Ordering origOrdering = g.getTripOrdering();
        std::vector<Ordering > permutations;
        getPermutations(origOrdering, &permutations);
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
        Ordering origOrdering = g.getTripOrdering();
        (*c)[&g] = origOrdering;
      }
    }
  }
}

// _____________________________________________________________________________
void
EdgeOrderOptimizer::getPermutations(std::vector<size_t> order, std::vector<std::vector<size_t> >* ret)
const {
  //std::sort(order.begin(), order.end());

  do {
    ret->push_back(order);
  } while (std::next_permutation(order.begin(), order.end()));
}

// _____________________________________________________________________________
void EdgeOrderOptimizer::applyConfig(const Configuration& c) {
  for (auto a : c) {
    a.first->setTripOrdering(a.second);
  }
}
