// Copyright 2016
// Author: Patrick Brosi

#ifndef TOPO_TEST_TOPOTESTUTIL_H_
#define TOPO_TEST_TOPOTESTUTIL_H_

#include "shared/linegraph/LineGraph.h"

using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using shared::linegraph::LineEdge;

inline bool hasEdge(const LineNode* n, const LineEdge* f) {
  for (auto e : n->getAdjList())
    if (e == f) return true;
  std::cerr << "Could not find edge " << f << " adjacent to node " << n
            << std::endl;
  return false;
}

inline bool validExceptions(const LineNode* n) {
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

inline bool validExceptions(const LineGraph* g) {
  for (auto nd : g->getNds()) {
    if (!validExceptions(nd)) return false;
  }
  return true;
}

#endif
