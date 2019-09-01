// Copyright 2016
// Author: Patrick Brosi

#ifndef TOPO_TEST_TOPOTESTUTIL_H_
#define TOPO_TEST_TOPOTESTUTIL_H_

#include "shared/transitgraph/TransitGraph.h"

using shared::transitgraph::TransitGraph;
using shared::transitgraph::TransitNode;
using shared::transitgraph::TransitEdge;

inline bool hasEdge(const TransitNode* n, const TransitEdge* f) {
  for (auto e : n->getAdjList())
    if (e == f) return true;
  std::cerr << "Could not find edge " << f << " adjacent to node " << n
            << std::endl;
  return false;
}

inline bool validExceptions(const TransitNode* n) {
  for (const auto& ro : n->pl().getConnExc()) {
    for (const auto& exFr : ro.second) {
      for (const auto& exTo : exFr.second) {
        if (!hasEdge(n, exFr.first)) return false;
        if (!hasEdge(n, exTo)) return false;
      }
    }
  }
  return true;
}

inline bool validExceptions(const TransitGraph* g) {
  for (auto nd : g->getNds()) {
    if (!validExceptions(nd)) return false;
  }
  return true;
}

#endif
