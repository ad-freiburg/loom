// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_map>
#include "loom/optim/OracleOptimizer.h"
#include "shared/linegraph/Line.h"
#include "util/log/Log.h"

using namespace loom;
using namespace optim;
using loom::optim::OracleOptimizer;
using loom::optim::SettledEdgs;
using shared::linegraph::Line;
using shared::rendergraph::HierarOrderCfg;

// _____________________________________________________________________________
int OracleOptimizer::optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                                  HierarOrderCfg* hc, size_t depth) const {
  UNUSED(depth);
  OptOrderCfg cur, null;

  const OptEdge* e = 0;
  SettledEdgs settled;

  while ((e = getNextEdge(g, &settled))) {
    std::cerr << "Optimizing " << e << " with " << e->pl().getCardinality()
              << " lines" << std::endl;
    for (auto lo : e->pl().getLines()) {
      std::cerr << lo.line->label() << ", ";
    }
    std::cerr << std::endl;

    // fill lines into empty config
    for (const auto& lo : e->pl().getLines()) {
      cur[e].push_back(lo.line);
    }

    auto cmp = LineCmp(e, cur, _optScorer);

    auto order = cur[e];
    std::sort(order.begin(), order.end(), cmp);
    cur[e] = order;

    settled.insert(e);
  }

  writeHierarch(&cur, hc);
  return 0;
}

// _____________________________________________________________________________
const OptEdge* OracleOptimizer::getNextEdge(const std::set<OptNode*>& g,
                                            SettledEdgs* settled) const {
  if (settled->size() == 0) return getInitialEdge(g);

  // else, use some unsettled edge adjacent to the settled set
  for (auto eid : *settled) {
    for (auto adj : eid->getFrom()->getAdjList()) {
      // if (adj->pl().getCardinality() < 2) continue;
      if (!settled->count(adj)) return adj;
    }
    for (auto adj : eid->getTo()->getAdjList()) {
      // if (adj->pl().getCardinality() < 2) continue;
      if (!settled->count(adj)) return adj;
    }
  }

  return 0;
}

// _____________________________________________________________________________
const OptEdge* OracleOptimizer::getInitialEdge(
    const std::set<OptNode*>& g) const {
  const OptEdge* ret = 0;

  for (auto n : g) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n || e->pl().getCardinality() < 2) continue;
      if (!ret || e->pl().getCardinality() > ret->pl().getCardinality() ||
          (e->pl().getCardinality() == ret->pl().getCardinality() &&
           e->getTo()->getDeg() + e->getFrom()->getDeg() >
               ret->getTo()->getDeg() + ret->getFrom()->getDeg())) {
        ret = e;
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
std::pair<int, double> LineCmp::smallerThanAt(const shared::linegraph::Line* a,
                                              const shared::linegraph::Line* b,
                                              const OptEdge* e,
                                              const OptNode* nd,
                                              const OptEdge* ign) const {
  // return -1 for false, 0 for undecided, 1 for true

  std::vector<size_t> positionsA;
  std::vector<size_t> positionsB;

  double cost = 0;

  size_t offset = 0;

  for (auto e : OptGraph::clockwEdges(e, nd)) {
    if (e == ign) continue;
    auto loA = e->pl().getLineOcc(a);
    auto loB = e->pl().getLineOcc(b);

    if (loA && loB) {
      auto eCfg = _cfg.find(e);
      if (eCfg != _cfg.end()) {
        bool rev = (e->getFrom() != nd) ^ e->pl().lnEdgParts.front().dir;
        size_t peaA = std::find(eCfg->second.begin(), eCfg->second.end(), a) -
                      eCfg->second.begin();
        size_t peaB = std::find(eCfg->second.begin(), eCfg->second.end(), b) -
                      eCfg->second.begin();
        if (rev) {
          positionsA.push_back(offset + peaA);
          positionsB.push_back(offset + peaB);
        } else {
          positionsA.push_back(offset + (e->pl().getCardinality() - 1 - peaA));
          positionsB.push_back(offset + (e->pl().getCardinality() - 1 - peaB));
        }
        cost += _scorer.getCrossingPenSameSeg(nd);
      }
    } else if (loA) {
      positionsA.push_back(offset);
      cost += _scorer.getCrossingPenDiffSeg(nd);
    } else if (loB) {
      positionsB.push_back(offset);
      cost += _scorer.getCrossingPenDiffSeg(nd);
    }

    offset += e->pl().getCardinality();
  }

  if (positionsA.size() == 0 || positionsB.size() == 0) return {0, 0};

  if (positionsA.back() < positionsB.front()) return {1, cost};
  if (positionsA.front() > positionsB.back()) return {-1, cost};
  return {0, 0};
}

// _____________________________________________________________________________
bool LineCmp::oracle(const shared::linegraph::Line* a,
                     const shared::linegraph::Line* b) const {
  // return a < b;
  int left = 0;
  int right = 0;

  double leftCost = 0;
  double rightCost = 0;

  std::cerr << "oracleing... " << a->label() << " vs " << b->label() << std::endl;

  // go to the left
  auto e = _e;
  auto curNd = e->getFrom();
  while (curNd) {
    std::cerr << "A" << std::endl;
    auto newECand = eligibleNextEdge(e, curNd, a, b);
    std::cerr << newECand << std::endl;
    if (!newECand) break;
    curNd = newECand->getOtherNd(curNd);
    auto i = smallerThanAt(a, b, newECand, curNd, e);
    if (i.first != 0) {
      left = i.first;
      leftCost = i.second;
      break;
    }
    e = newECand;
  }

  // go to the right
  e = _e;
  curNd = e->getTo();
  while (curNd) {
    std::cerr << "B" << std::endl;
    auto newECand = eligibleNextEdge(e, curNd, a, b);
    if (!newECand) break;
    curNd = newECand->getOtherNd(curNd);
    auto i = smallerThanAt(a, b, newECand, curNd, e);
    if (i.first != 0) {
      right = i.first;
      rightCost = i.second;
      break;
    }
    e = newECand;
  }

  // build return value

  bool rev = !_e->pl().lnEdgParts.front().dir;

  if (left != 0 && right != 0) {
    if (left == -right) {
      // both sides induce the same, and are not undecided
      if (left == 1) return false ^ rev;
      return true ^ rev;
    } else {
      // both sides induce different orderings, take the order of the side
      // where a crossing is more expensive
      if (leftCost > rightCost) {
        if (left == 1) return false ^ rev;
        return true ^ rev;
      } else {
        if (right == 1) return true ^ rev;
        return false ^ rev;
      }
    }
  } else if (left != 0 && right == 0) {
    // left side induces
    if (left == 1) return false ^ rev;
    return true ^ rev;
  } else if (left == 0 && right != 0) {
    // right side induces
    if (right == 1) return true ^ rev;
    return false ^ rev;
  }

  std::cerr << "(failure)" << std::endl;

  return a < b;
}

// _____________________________________________________________________________
const OptEdge* LineCmp::eligibleNextEdge(
    const OptEdge* start, const OptNode* nd, const shared::linegraph::Line* a,
    const shared::linegraph::Line* b) const {
  const OptEdge* newECand = 0;
  for (auto newE : nd->getAdjList()) {
    if (newE == start) continue;
    if (newE->pl().getLineOcc(a) && newE->pl().getLineOcc(b)) {
      if (newECand != 0) return 0;  // both branch, abort
      newECand = newE;
    }
  }
  return newECand;
}
