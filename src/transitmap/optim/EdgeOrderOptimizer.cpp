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

  size_t numRuns = 100;

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

    while (doOptimStep(&c)) {
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
  double s = _g->getScore(results[best].second);
  applyConfig(results[best].second);
  double ss = _g->getScore(results[best].second);
  std::cout << s << " vs1 " << ss << std::endl;
  assert(s == ss);
  std::cout << s << " vs " << _g->getScore() << std::endl;
}

// _____________________________________________________________________________
bool EdgeOrderOptimizer::doOptimStep(Configuration* c) {
  std::pair<EdgeTripGeom*, std::vector<size_t> > bestCand;
  bestCand.first = 0;

  double oldScore = _g->getScore(*c);

  for (graph::Node* n : *_g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (graph::EdgeTripGeom& g : *e->getEdgeTripGeoms()) {
        if (g.getTripOrdering().size() == 1) continue;

        double oldAreaScore = n->getAreaScore(*c);
        Ordering origOrdering = g.getTripOrdering();
        std::vector<Ordering > permutations = getPermutations(origOrdering);

        for (size_t i = 0; i < permutations.size(); ++i) {
          double newAreaScore = n->getAreaScore(*c, g, permutations[i]);

          // for testin
          double diff = oldAreaScore - newAreaScore;
          Ordering old = (*c)[&g];
          (*c)[&g] = permutations[i];
          std::cout << newAreaScore << " " <<  n->getAreaScore(*c) << std::endl;
          assert(newAreaScore == n->getAreaScore(*c));
          double diff2 = oldScore - _g->getScore(*c);

          //std::cout << "Scores: " << diff << " vs " << diff2 << std::endl;
          //assert(fabs(diff - diff2) < 0.001);

          (*c)[&g] = old;

          if (diff2 > 0) { //newAreaScore < oldAreaScore) {
            bestCand = std::pair<EdgeTripGeom*, Ordering>(&g, permutations[i]);
            oldAreaScore = newAreaScore;
          }
        }
      }
    }
  }

  if (bestCand.first) {
    bestCand.first->setTripOrdering(bestCand.second);
    (*c)[bestCand.first] = bestCand.second;
    assert(_g->getScore(*c) < oldScore);
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
