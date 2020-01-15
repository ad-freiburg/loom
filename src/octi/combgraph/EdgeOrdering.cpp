// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/combgraph/CombNodePL.h"
#include "octi/combgraph/EdgeOrdering.h"

using util::geo::Point;
using octi::combgraph::EdgeOrdering;
using octi::combgraph::CombEdge;

// _____________________________________________________________________________
void EdgeOrdering::add(CombEdge* e, double deg) {
  _edgeOrder.push_back(std::pair<CombEdge*, double>(e, deg));
  std::sort(_edgeOrder.begin(), _edgeOrder.end(), PairCmp());
}

// _____________________________________________________________________________
int64_t EdgeOrdering::dist(CombEdge* a, CombEdge* b) const {
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
const std::vector<std::pair<CombEdge*, double>>& EdgeOrdering::getOrderedSet()
    const {
  return _edgeOrder;
}

// _____________________________________________________________________________
bool EdgeOrdering::has(CombEdge* e) const {
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

// _____________________________________________________________________________
std::string EdgeOrdering::toString(CombNode* origNode) const {
  std::stringstream ss;
  for (auto e : getOrderedSet()) {
    if (e.first->getOtherNd(origNode)->pl().getParent()->pl().stops().size()) {
      ss << e.first << "(to "
         << e.first->getOtherNd(origNode)
                ->pl()
                .getParent()
                ->pl()
                .stops()
                .front()
                .name
         << "), ";
    } else {
      ss << e.first << ", ";
    }
  }
  return ss.str();
}
