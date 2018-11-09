// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/Optimizer.h"

using transitmapper::optim::Optimizer;
using transitmapper::optim::LinePair;
using transitmapper::optim::PosComPair;
using transitmapper::optim::EdgePair;
using transitmapper::graph::NodeFront;

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

// _____________________________________________________________________________
std::vector<LinePair> Optimizer::getLinePairs(OptEdge* segment) {
  return getLinePairs(segment, false);
}

// _____________________________________________________________________________
std::vector<LinePair> Optimizer::getLinePairs(OptEdge* segment,
                                                 bool unique) {
  std::set<const Route*> processed;
  std::vector<LinePair> ret;
  for (auto& toA : segment->pl().getRoutes()) {
    if (toA.route->relativeTo()) continue;
    processed.insert(toA.route);
    for (auto& toB : segment->pl().getRoutes()) {
      if (unique && processed.count(toB.route)) continue;
      if (toB.route->relativeTo()) continue;
      if (toA.route == toB.route) continue;

      // this is to make sure that we always get the same line pairs in unique
      // mode -> take the smaller pointer first
      if (!unique || toA.route < toB.route)
        ret.push_back(LinePair(toA.route, toB.route));
      else
        ret.push_back(LinePair(toB.route, toA.route));
    }
  }
  return ret;
}

// _____________________________________________________________________________
bool Optimizer::crosses(OptNode* node, OptEdge* segmentA, OptEdge* segmentB,
                           PosComPair poscomb) {
  bool otherWayA =
      (segmentA->getFrom() != node) ^ segmentA->pl().etgs.front().dir;
  bool otherWayB =
      (segmentB->getFrom() != node) ^ segmentB->pl().etgs.front().dir;

  size_t cardA = segmentA->pl().etgs.front().etg->getCardinality(true);
  size_t cardB = segmentB->pl().etgs.front().etg->getCardinality(true);

  size_t posAinA =
      otherWayA ? cardA - 1 - poscomb.first.first : poscomb.first.first;
  size_t posAinB =
      otherWayB ? cardB - 1 - poscomb.first.second : poscomb.first.second;
  size_t posBinA =
      otherWayA ? cardA - 1 - poscomb.second.first : poscomb.second.first;
  size_t posBinB =
      otherWayB ? cardB - 1 - poscomb.second.second : poscomb.second.second;

  bool pCrossing = false;

  DPoint aInA = getPos(node, segmentA, posAinA);
  DPoint bInA = getPos(node, segmentA, posBinA);
  DPoint aInB = getPos(node, segmentB, posAinB);
  DPoint bInB = getPos(node, segmentB, posBinB);

  DLine a;
  a.push_back(aInA);
  a.push_back(aInB);

  DLine b;
  b.push_back(bInA);
  b.push_back(bInB);

  if (util::geo::intersects(aInA, aInB, bInA, bInB) ||
      util::geo::dist(a, b) < 1)
    pCrossing = true;

  // pCrossing = (posAinA < posBinA && posAinB < posBinB);

  return pCrossing;
}

// _____________________________________________________________________________
bool Optimizer::crosses(OptNode* node, OptEdge* segmentA, EdgePair segments,
                           PosCom postcomb) {
  bool otherWayA =
      (segmentA->getFrom() != node) ^ segmentA->pl().etgs.front().dir;
  bool otherWayB = (segments.first->getFrom() != node) ^
                   segments.first->pl().etgs.front().dir;
  bool otherWayC = (segments.second->getFrom() != node) ^
                   segments.second->pl().etgs.front().dir;

  size_t cardA = segmentA->pl().etgs.front().etg->getCardinality(true);
  size_t cardB = segments.first->pl().etgs.front().etg->getCardinality(true);
  size_t cardC = segments.second->pl().etgs.front().etg->getCardinality(true);

  size_t posAinA = otherWayA ? cardA - 1 - postcomb.first : postcomb.first;
  size_t posBinA = otherWayA ? cardA - 1 - postcomb.second : postcomb.second;

  DPoint aInA = getPos(node, segmentA, posAinA);
  DPoint bInA = getPos(node, segmentA, posBinA);

  for (size_t i = 0;
       i < segments.first->pl().etgs.front().etg->getCardinality(true); ++i) {
    for (size_t j = 0;
         j < segments.second->pl().etgs.front().etg->getCardinality(true);
         ++j) {
      size_t posAinB = otherWayB ? cardB - 1 - i : i;
      size_t posBinC = otherWayC ? cardC - 1 - j : j;

      DPoint aInB = getPos(node, segments.first, posAinB);
      DPoint bInC = getPos(node, segments.second, posBinC);

      DLine a;
      a.push_back(aInA);
      a.push_back(aInB);

      DLine b;
      b.push_back(bInA);
      b.push_back(bInC);

      if (util::geo::intersects(aInA, aInB, bInA, bInC) ||
          util::geo::dist(a, b) < 1)
        return true;
    }
  }
  return false;
}

// _____________________________________________________________________________
DPoint Optimizer::getPos(OptNode* n, OptEdge* segment, size_t p) {
  // look for correct nodefront
  const NodeFront* nf = 0;
  for (auto etg : segment->pl().etgs) {
    const NodeFront* test = n->pl().node->getNodeFrontFor(etg.etg);
    if (test) {
      nf = test;
      break;
    }
  }

  return nf->getTripPos(segment->pl().etgs.front().etg, p, false);
}

// _____________________________________________________________________________
std::vector<OptEdge*> Optimizer::getEdgePartners(
    OptNode* node, OptEdge* segmentA, const LinePair& linepair) {
  std::vector<OptEdge*> ret;

  graph::Edge* fromEtg = OptGraph::getAdjEdg(segmentA, node);
  const Node* dirA = fromEtg->getRoute(linepair.first)->direction;
  const Node* dirB = fromEtg->getRoute(linepair.second)->direction;

  for (OptEdge* segmentB : node->getAdjList()) {
    if (segmentB == segmentA) continue;

    if (OptGraph::getCtdRoutesIn(linepair.first, dirA, segmentA, segmentB)
            .size() &&
        OptGraph::getCtdRoutesIn(linepair.second, dirB, segmentA, segmentB)
            .size()) {
      ret.push_back(segmentB);
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::vector<EdgePair> Optimizer::getEdgePartnerPairs(
    OptNode* node, OptEdge* segmentA, const LinePair& linepair) {
  std::vector<EdgePair> ret;

  graph::Edge* fromEtg = OptGraph::getAdjEdg(segmentA, node);
  const Node* dirA = fromEtg->getRoute(linepair.first)->direction;
  const Node* dirB = fromEtg->getRoute(linepair.second)->direction;

  for (OptEdge* segmentB : node->getAdjList()) {
    if (segmentB == segmentA) continue;

    if (OptGraph::getCtdRoutesIn(linepair.first, dirA, segmentA, segmentB)
            .size()) {
      EdgePair curPair;
      curPair.first = segmentB;
      for (OptEdge* segmentC : node->getAdjList()) {
        if (segmentC == segmentA || segmentC == segmentB) continue;

        if (OptGraph::getCtdRoutesIn(linepair.second, dirB, segmentA, segmentC)
                .size()) {
          curPair.second = segmentC;
          ret.push_back(curPair);
        }
      }
    }
  }
  return ret;
}
