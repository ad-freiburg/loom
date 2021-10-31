// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "loom/optim/NullOptimizer.h"
#include "util/Misc.h"
#include "util/log/Log.h"

using loom::optim::NullOptimizer;
using shared::rendergraph::HierarOrderCfg;

// _____________________________________________________________________________
double NullOptimizer::optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                                HierarOrderCfg* hc, size_t depth,
                                OptResStats& stats) const {
  LOGTO(DEBUG, std::cerr) << prefix(depth)
                          << "(NullOptimizer) Optimizing component with "
                          << g.size() << " nodes.";
  T_START(1);

  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      for (auto lnEdgPart : e->pl().lnEdgParts) {
        if (lnEdgPart.wasCut) continue;
        for (auto ro : e->pl().getLines()) {
          for (auto rel : ro.relatives) {
            // retrieve the original line pos
            size_t p = lnEdgPart.lnEdg->pl().linePos(rel);
            if (!(lnEdgPart.dir ^ e->pl().lnEdgParts.front().dir)) {
              (*hc)[lnEdgPart.lnEdg][lnEdgPart.order].insert(
                  (*hc)[lnEdgPart.lnEdg][lnEdgPart.order].begin(), p);
            } else {
              (*hc)[lnEdgPart.lnEdg][lnEdgPart.order].push_back(p);
            }
          }
        }
      }
    }
  }
  return T_STOP(1);
}
