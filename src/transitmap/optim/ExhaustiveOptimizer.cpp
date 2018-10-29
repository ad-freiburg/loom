// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/ExhaustiveOptimizer.h"

using namespace transitmapper;
using namespace optim;
using namespace transitmapper::graph;
using transitmapper::optim::ExhaustiveOptimizer;

// _____________________________________________________________________________
int ExhaustiveOptimizer::optimize(TransitGraph* tg) const {
  // TODO

  return 0;
}

// _____________________________________________________________________________
int ExhaustiveOptimizer::optimize(const std::set<OptNode*>& g,
                                  HierarchOrderingConfig* hc) const {
  OptOrderingConfig best;
  initialConfig(g, &best);

  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;

      // double score = _scorer.getScore(
    }
  }

  writeHierarch(&best, hc);
  return 0;
}

// _____________________________________________________________________________
void ExhaustiveOptimizer::initialConfig(const std::set<OptNode*>& g,
                                        OptOrderingConfig* cfg) const {
  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      (*cfg)[e] = std::vector<size_t>(e->pl().getCardinality());
      for (size_t i = 0; i < (*cfg)[e].size(); i++) (*cfg)[e][i] = i;
    }
  }
}

// _____________________________________________________________________________
void ExhaustiveOptimizer::writeHierarch(OptOrderingConfig* cfg, HierarchOrderingConfig* hc) const {
  for (auto ep : *cfg) {
    auto e = ep.first;

    for (auto etgp : e->pl().etgs) {
      if (etgp.wasCut) continue;
      for (auto tp : ep.second) {
        auto r = ep.first->pl().getRoutes()[tp];
        for (size_t p = 0; p < etgp.etg->getCardinality(); p++) {
          auto ro = (*etgp.etg->getRoutes())[p];
          if (!(r == ro)) continue;

          if (std::find(e->pl().getRoutes().begin(),
                        e->pl().getRoutes().end(),
                        ro) == e->pl().getRoutes().end())
            continue;

          if (ro.route->relativeTo()) continue;

          if (!(etgp.dir ^ e->pl().etgs.front().dir)) {
            (*hc)[etgp.etg][etgp.order].insert(
                (*hc)[etgp.etg][etgp.order].begin(), p);
          } else {
            (*hc)[etgp.etg][etgp.order].push_back(p);
          }
        }
      }

      std::cout << std::endl;
      std::cout << "From ";

      if (etgp.etg->getFrom()->getStops().size())
      std::cout << etgp.etg->getFrom()->getStops().front().name;
      std::cout << " to ";
      if (etgp.etg->getTo()->getStops().size())
      std::cout << etgp.etg->getTo()->getStops().front().name;
      std::cout << " " << etgp.order << ": [";
      for (auto p : (*hc)[etgp.etg][etgp.order]) std::cout << p << ", ";
      std::cout << "] ";
      std::cout << "(" << ep.second.size() << ")";
      std::cout << "(" << ep.first->pl().getCardinality() << ")";
    }
  }

  std::cout << std::endl;
}
