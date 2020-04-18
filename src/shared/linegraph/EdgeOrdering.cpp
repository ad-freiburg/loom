// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "shared/linegraph/EdgeOrdering.h"

using util::geo::Point;
using shared::linegraph::EdgeOrdering;
using shared::linegraph::LineEdge;

// _____________________________________________________________________________
void EdgeOrdering::add(LineEdge* e, double deg) {
  _edgeOrder.push_back(std::pair<LineEdge*, double>(e, deg));
  std::sort(_edgeOrder.begin(), _edgeOrder.end(), PairCmp());
}

// _____________________________________________________________________________
int64_t EdgeOrdering::dist(LineEdge* a, LineEdge* b) const {
  auto aIt = _edgeOrder.begin();
  auto bIt = _edgeOrder.begin();

  for (; aIt != _edgeOrder.end(); aIt++) {
    if (aIt->first == a) break;
  }

  for (; bIt != _edgeOrder.end(); bIt++) {
    if (bIt->first == b) break;
  }

  assert(aIt != _edgeOrder.end());
  assert(bIt != _edgeOrder.end());

  int32_t ap = std::distance(_edgeOrder.begin(), aIt);
  int32_t bp = std::distance(_edgeOrder.begin(), bIt);

  int32_t ret = ((ap > bp ? -1 : 1) * (abs(bp - ap)));

  if (ret < 0) ret = (_edgeOrder.size() + ret);

  return ret;
}

// _____________________________________________________________________________
const std::vector<std::pair<LineEdge*, double>>& EdgeOrdering::getOrderedSet()
    const {
  return _edgeOrder;
}

// _____________________________________________________________________________
bool EdgeOrdering::has(LineEdge* e) const {
  auto aIt = _edgeOrder.begin();

  for (; aIt != _edgeOrder.end(); aIt++) {
    if (aIt->first == e) return true;
  }

  return false;
}

// _____________________________________________________________________________
bool EdgeOrdering::equals(const EdgeOrdering& e) const {
  auto a = getOrderedSet().begin();
  auto b = e.getOrderedSet().begin();
  for (; b != e.getOrderedSet().end(); b++) {
    if (a->first == b->first) break;
  }

  for (; a != getOrderedSet().end(); a++) {
    if (a->first != b->first) return false;
    b++;
    if (b == e.getOrderedSet().end()) b = e.getOrderedSet().begin();
  }

  return true;
}
