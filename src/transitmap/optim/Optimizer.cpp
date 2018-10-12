// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/Optimizer.h"

using transitmapper::optim::Optimizer;

// _____________________________________________________________________________
void Optimizer::expandRelatives(TransitGraph* g, OrderingConfig* c) {
  std::set<std::pair<graph::Edge*, const Route*>> visited;
  for (graph::Node* n : *g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (const auto& ra : *(e->getRoutes())) {
        if (ra.route->relativeTo()) {
          const Route* ref = ra.route->relativeTo();
          expandRelativesFor(c, ref, e, e->getRoutesRelTo(ref), visited);
        }
      }
    }
  }
}

// _____________________________________________________________________________
void Optimizer::expandRelativesFor(
    OrderingConfig* c, const Route* ref, graph::Edge* start,
    const std::set<const Route*>& rs,
    std::set<std::pair<graph::Edge*, const Route*>>& visited) {
  std::stack<std::pair<graph::Edge*, graph::Edge*>> todo;

  todo.push(std::pair<graph::Edge*, graph::Edge*>(0, start));
  while (!todo.empty()) {
    auto cur = todo.top();
    todo.pop();

    if (!visited.insert(std::pair<graph::Edge*, const Route*>(cur.second, ref))
             .second)
      continue;

    for (auto r : rs) {
      size_t index = cur.second->getRouteWithPos(ref).second;
      auto it =
          std::find((*c)[cur.second].begin(), (*c)[cur.second].end(), index);

      size_t p = cur.second->getRouteWithPos(r).second;

      assert(it != (*c)[cur.second].end());

      if (cur.first != 0 &&
          ((cur.first->getTo() == cur.second->getTo() ||
            cur.first->getFrom() == cur.second->getFrom()) ^
           (cur.first->getRouteWithPosUnder(r, (*c)[cur.first]).second >
            cur.first->getRouteWithPosUnder(ref, (*c)[cur.first]).second))) {
        (*c)[cur.second].insert(it + 1, p);
      } else {
        (*c)[cur.second].insert(it, p);
      }
    }

    for (const Node* n : {cur.second->getFrom(), cur.second->getTo()}) {
      if (cur.first != 0 &&
          (cur.first->getTo() == n || cur.first->getFrom() == n)) {
        continue;
      }

      for (auto e : n->getAdjListIn()) {
        if (e->containsRoute(ref) &&
            visited.find(std::pair<graph::Edge*, const Route*>(e, ref)) ==
                visited.end()) {
          todo.push(std::pair<graph::Edge*, graph::Edge*>(cur.second, e));
        }
      }

      for (auto e : n->getAdjListOut()) {
        if (e->containsRoute(ref) &&
            visited.find(std::pair<graph::Edge*, const Route*>(e, ref)) ==
                visited.end()) {
          todo.push(std::pair<graph::Edge*, graph::Edge*>(cur.second, e));
        }
      }
    }
  }
}
