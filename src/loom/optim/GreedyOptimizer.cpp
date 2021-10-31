// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_map>
#include "loom/optim/GreedyOptimizer.h"
#include "shared/linegraph/Line.h"
#include "util/log/Log.h"

using namespace loom;
using namespace optim;
using loom::optim::GreedyOptimizer;
using loom::optim::SettledEdgs;
using shared::linegraph::Line;
using shared::rendergraph::HierarOrderCfg;

// _____________________________________________________________________________
double GreedyOptimizer::optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                                  HierarOrderCfg* hc, size_t depth,
                                  OptResStats& stats) const {
  LOGTO(DEBUG, std::cerr) << prefix(depth)
                          << "(GreedyOptimizer) Optimizing component with "
                          << g.size() << " nodes.";
  T_START(1);
  UNUSED(depth);
  OptOrderCfg cfg;

  getFlatConfig(g, &cfg);

  writeHierarch(&cfg, hc);
  return T_STOP(1);
}

// _____________________________________________________________________________
void GreedyOptimizer::getFlatConfig(const std::set<OptNode*>& g,
                                   OptOrderCfg* cfg) const {
  const OptEdge* e = 0;
  SettledEdgs settled;

  while ((e = getNextEdge(g, &settled))) {
    Cmp left, right;

    // build cmp function for left
    for (const auto& lo1 : e->pl().getLines()) {
      for (const auto& lo2 : e->pl().getLines()) {
        if (lo1.line == lo2.line) continue;
        left[{lo1.line, lo2.line}] =
            guess(lo1.line, lo2.line, e, e->getFrom(), *cfg);
        right[{lo1.line, lo2.line}] =
            guess(lo1.line, lo2.line, e, e->getTo(), *cfg);
      }
    }

    double costLeft = 0;
    double costRight = 0;

    // which one is cheaper?
    for (const auto& lo1 : e->pl().getLines()) {
      for (const auto& lo2 : e->pl().getLines()) {
        if (lo1.line == lo2.line) continue;
        if (left[{lo1.line, lo2.line}].first ==
            right[{lo1.line, lo2.line}].first) {
          costLeft += right[{lo1.line, lo2.line}].second;
          costRight += left[{lo1.line, lo2.line}].second;
        }
      }
    }

    // build cmp function for right

    LineCmp cmp({}, false);
    if (costLeft < costRight) {
      cmp = LineCmp(left, false);
    } else {
      cmp = LineCmp(right, true);
    }

    // fill lines into empty config
    for (const auto& lo : e->pl().getLines()) {
      (*cfg)[e].push_back(lo.line);
    }

    std::sort((*cfg)[e].begin(), (*cfg)[e].end(), cmp);

    settled.insert(e);
  }
}

// _____________________________________________________________________________
const OptEdge* GreedyOptimizer::getNextEdge(const std::set<OptNode*>& g,
                                            SettledEdgs* settled) const {
  if (settled->size() == 0) return getInitialEdge(g);

  // else, use some unsettled edge adjacent to the settled set
  for (auto eid : *settled) {
    for (auto adj : eid->getFrom()->getAdjList()) {
      if (!settled->count(adj)) return adj;
    }
    for (auto adj : eid->getTo()->getAdjList()) {
      if (!settled->count(adj)) return adj;
    }
  }

  return 0;
}

// _____________________________________________________________________________
const OptEdge* GreedyOptimizer::getInitialEdge(
    const std::set<OptNode*>& g) const {
  const OptEdge* ret = 0;

  for (auto n : g) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
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
std::pair<int, double> GreedyOptimizer::smallerThanAt(
    const shared::linegraph::Line* a, const shared::linegraph::Line* b,
    const OptEdge* start, const OptNode* nd, const OptEdge* ign,
    const OptOrderCfg& cfg) const {
  // return -1 for false, 0 for undecided, 1 for true
  std::vector<size_t> positionsA;
  std::vector<size_t> positionsB;

  double cost = 0;

  size_t offset = 0;

  for (auto e : OptGraph::clockwEdges(start, nd)) {
    if (e == ign) continue;
    auto loA = e->pl().getLineOcc(a);
    auto loB = e->pl().getLineOcc(b);

    if (loA && loB) {
      auto eCfg = cfg.find(e);
      if (eCfg != cfg.end()) {
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
std::pair<bool, double> GreedyOptimizer::guess(const shared::linegraph::Line* a,
                                               const shared::linegraph::Line* b,
                                               const OptEdge* start,
                                               const OptNode* refNd,
                                               const OptOrderCfg& cfg) const {
  int dec = 0;
  bool notRef = false;

  double cost = 0;

  // on ref node
  auto e = start;
  auto curNd = refNd;
  while (true) {
    auto i = smallerThanAt(a, b, e, curNd, e, cfg);
    if (i.first != 0) {
      dec = i.first;
      cost = i.second;
      break;
    }
    e = eligibleNextEdge(e, curNd, a, b);
    if (!e || e == start) break;
    curNd = e->getOtherNd(curNd);
  }

  if (!dec) {
    e = start;
    curNd = start->getOtherNd(refNd);
    while (true) {
      auto i = smallerThanAt(a, b, e, curNd, e, cfg);
      if (i.first != 0) {
        dec = i.first;
        cost = i.second;
        notRef = true;
        break;
      }
      e = eligibleNextEdge(e, curNd, a, b);
      if (!e || e == start) break;
      curNd = e->getOtherNd(curNd);
    }
  }

  // build return value
  bool rev = !start->pl().lnEdgParts.front().dir;

  if (dec == 1) return {notRef ^ rev, cost};
  if (dec == -1) return {!notRef ^ rev, cost};

  // undecided
  return {a < b, 0};
}

// _____________________________________________________________________________
const OptEdge* GreedyOptimizer::eligibleNextEdge(
    const OptEdge* start, const OptNode* nd, const shared::linegraph::Line* a,
    const shared::linegraph::Line* b) const {
  if (!_lookAhead) return 0;
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
