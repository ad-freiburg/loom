// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fstream>
#include "transitmap/optim/Optimizer.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Algorithm.h"
#include "util/log/Log.h"

using transitmapper::optim::Optimizer;
using transitmapper::optim::LinePair;
using transitmapper::optim::PosComPair;
using transitmapper::optim::EdgePair;
using transitmapper::graph::NodeFront;

// _____________________________________________________________________________
int Optimizer::optimize(TransitGraph* tg) const {
  // create optim graph
  OptGraph g(tg, _scorer);

  size_t maxC = maxCard(*g.getNds());
  double solSp = solutionSpaceSize(*g.getNds());
  LOG(DEBUG) << "Optimizing line graph of size " << tg->getNodes()->size()
             << " with max cardinality = " << maxC
             << " and solution space size = " << solSp;

  if (_cfg->untangleGraph) {
    // do full untangling
    LOG(DEBUG) << "Untangling graph...";
    T_START(1);
    for (size_t i = 0; i < maxC; i++) {
      g.untangle();
      g.simplify();
      g.split();
    }
    LOG(DEBUG) << "Done (" << T_STOP(1) << " ms)";
  } else if (_cfg->createCoreOptimGraph) {
    // only apply core graph rules
    T_START(1);
    LOG(DEBUG) << "Creating core optimization graph...";
    g.simplify();
    LOG(DEBUG) << "Done (" << T_STOP(1) << " ms)";
  }

  if (_cfg->outOptGraph) {
    util::geo::output::GeoGraphJsonOutput out;
    std::ofstream fstr(_cfg->dbgPath + "/optgraph.json");
    out.print(g, fstr);
  }

  if (_cfg->outputStats) {
    LOG(INFO) << "(stats) Stats for <optim> graph of '" << tg->getName() << "'";
    LOG(INFO) << "(stats)   Total node count: " << g.getNumNodes();
    LOG(INFO) << "(stats)   Total edge count: " << g.getNumEdges();
    // LOG(INFO) << "(stats)   Total unique route count: " << g.getNumRoutes();
    LOG(INFO) << "(stats)   Max edge route cardinality: "
              << maxCard(*g.getNds());
    LOG(INFO) << "(stats)   Solution space: "
              << solutionSpaceSize(*g.getNds());
  }

  // iterate over components and optimize all of them separately
  const auto& comps = util::graph::Algorithm::connectedComponents(g);

  size_t runs = _cfg->optimRuns;
  double tSum = 0;
  double iterSum = 0;
  double scoreSum = 0;
  double crossSum = 0;
  double sepSum = 0;

  for (size_t run = 0; run < runs; run++) {
    OrderingConfig c;
    HierarchOrderingConfig hc;

    T_START(1);
    size_t iters = 0;
    double maxCompSolSpace = 0;
    size_t maxCompC = 0;
    size_t maxNumNodes = 0;
    size_t maxNumEdges = 0;
    size_t numM1Comps = 0;

    for (const auto nds : comps) {
      if (_cfg->outputStats) {
        size_t maxC = maxCard(nds);
        double solSp = solutionSpaceSize(nds);
        if (maxC > maxCompC) maxCompC = maxC;
        if (solSp > maxCompSolSpace) maxCompSolSpace = solSp;
        if (solSp == 1) numM1Comps++;
        if (nds.size() > maxNumNodes) maxNumNodes = nds.size();
        if (numEdges(nds) > maxNumEdges) maxNumEdges = numEdges(nds);

        LOG(INFO) << " (stats) Optimizing subgraph of size " << nds.size()
                   << " with max cardinality = " << maxC
                   << " and solution space size = " << solSp;
      }
      iters += optimizeComp(nds, &hc);
    }

    double t = T_STOP(1);

    LOG(INFO) << " -- Optimization took " << t << " ms -- ";

    if (_cfg->outputStats) {
      LOG(INFO) << " (stats) Number of components: " << comps.size();
      LOG(INFO) << " (stats) Number of components with M=1: " << numM1Comps;
      LOG(INFO) << " (stats) Max number of nodes of all components: " << maxNumNodes;
      LOG(INFO) << " (stats) Max number of edges of all components: " << maxNumEdges;
      LOG(INFO) << " (stats) Max cardinality of all components: " << maxCompC;
      LOG(INFO) << " (stats) Max solution space size of all components: "
                 << maxCompSolSpace;
    }

    tSum += t;
    iterSum += iters;

    hc.writeFlatCfg(&c);

    Optimizer::expandRelatives(tg, &c);
    tg->setConfig(c);

    if (runs > 1) {
      if (_cfg->splittingOpt) scoreSum += _scorer->getScore();
      else scoreSum += _scorer->getCrossScore();
      crossSum += _scorer->getNumCrossings();
      sepSum += _scorer->getNumSeparations();
    }
  }

  if (runs > 1 && _cfg->outputStats) {
    LOG(INFO) << "";
    LOG(INFO) << "(multiple opt runs stats) avg time: " << tSum / (1.0 * runs) << " ms";
    LOG(INFO) << "(multiple opt runs stats) avg iters: " << iterSum / (1.0 * runs);
    LOG(INFO) << "(multiple opt runs stats) avg score: -- " << scoreSum / (1.0 * runs) << " --";
    LOG(INFO) << "(multiple opt runs stats) avg num crossings: -- " << crossSum / (1.0 * runs) << " --";
    LOG(INFO) << "(multiple opt runs stats) avg num separations: -- " << sepSum / (1.0 * runs) << " --";
    LOG(INFO) << "";
  }

  return 0;
}

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
std::vector<LinePair> Optimizer::getLinePairs(OptEdge* segment, bool unique) {
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

  // size_t cardA = segmentA->pl().etgs.front().etg->getCardinality(true);
  // size_t cardB = segmentB->pl().etgs.front().etg->getCardinality(true);

  size_t cardA = segmentA->pl().getCardinality();
  size_t cardB = segmentB->pl().getCardinality();

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

  // size_t cardA = segmentA->pl().etgs.front().etg->getCardinality(true);
  // size_t cardB = segments.first->pl().etgs.front().etg->getCardinality(true);
  // size_t cardC =
  // segments.second->pl().etgs.front().etg->getCardinality(true);

  size_t cardA = segmentA->pl().getCardinality();
  size_t cardB = segments.first->pl().getCardinality();
  size_t cardC = segments.second->pl().getCardinality();

  size_t posAinA = otherWayA ? cardA - 1 - postcomb.first : postcomb.first;
  size_t posBinA = otherWayA ? cardA - 1 - postcomb.second : postcomb.second;

  size_t posEdgeA =
      std::distance(node->pl().orderedEdges.begin(),
                    std::find(node->pl().orderedEdges.begin(),
                              node->pl().orderedEdges.end(), segments.first));
  size_t posEdgeB =
      std::distance(node->pl().orderedEdges.begin(),
                    std::find(node->pl().orderedEdges.begin(),
                              node->pl().orderedEdges.end(), segments.second));

  bool retB = false;

  if (posAinA > posBinA && posEdgeA < posEdgeB) retB = true;
  retB = false;

  DPoint aInA = getPos(node, segmentA, posAinA);
  DPoint bInA = getPos(node, segmentA, posBinA);

  bool ret = false;

  for (size_t i = 0; i < segments.first->pl().getCardinality(); ++i) {
    for (size_t j = 0; j < segments.second->pl().getCardinality(); ++j) {
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
          util::geo::dist(a, b) < 1) {
        ret = true;
        break;
      }
    }
  }

  // std::cout << "---" << std::endl;
  // std::cout << segmentA->pl().toStr() << std::endl;
  // std::cout << segments.first->pl().toStr() << std::endl;
  // std::cout << segments.second->pl().toStr() << std::endl;
  // std::cout << postcomb.first << ":" << postcomb.second << std::endl;
  // assert(ret == retB);
  return ret;
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
std::vector<OptEdge*> Optimizer::getEdgePartners(OptNode* node,
                                                 OptEdge* segmentA,
                                                 const LinePair& linepair) {
  std::vector<OptEdge*> ret;

  graph::Edge* fromEtg = OptGraph::getAdjEdg(segmentA, node);
  const Node* dirA = fromEtg->getRoute(linepair.first)->direction;
  const Node* dirB = fromEtg->getRoute(linepair.second)->direction;

  for (OptEdge* segmentB : node->getAdjList()) {
    if (segmentB == segmentA) continue;

    if (OptGraph::hasCtdRoutesIn(linepair.first, dirA, segmentA, segmentB) &&
        OptGraph::hasCtdRoutesIn(linepair.second, dirB, segmentA, segmentB)) {
      ret.push_back(segmentB);
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::vector<EdgePair> Optimizer::getEdgePartnerPairs(OptNode* node,
                                                     OptEdge* segmentA,
                                                     const LinePair& linepair) {
  std::vector<EdgePair> ret;

  graph::Edge* fromEtg = OptGraph::getAdjEdg(segmentA, node);
  const Node* dirA = fromEtg->getRoute(linepair.first)->direction;
  const Node* dirB = fromEtg->getRoute(linepair.second)->direction;

  for (OptEdge* segmentB : node->getAdjList()) {
    if (segmentB == segmentA) continue;

    if (OptGraph::hasCtdRoutesIn(linepair.first, dirA, segmentA, segmentB)) {
      EdgePair curPair;
      curPair.first = segmentB;
      for (OptEdge* segmentC : node->getAdjList()) {
        if (segmentC == segmentA || segmentC == segmentB) continue;

        if (OptGraph::hasCtdRoutesIn(linepair.second, dirB, segmentA,
                                     segmentC)) {
          curPair.second = segmentC;
          ret.push_back(curPair);
        }
      }
    }
  }
  return ret;
}

// _____________________________________________________________________________
size_t Optimizer::maxCard(const std::set<OptNode*>& g) {
  size_t ret = 0;
  for (const auto* n : g) {
    for (const auto* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (e->pl().getCardinality() > ret) ret = e->pl().getCardinality();
    }
  }

  return ret;
}

// _____________________________________________________________________________
double Optimizer::numEdges(const std::set<OptNode*>& g) {
  double ret = 0;
  for (const auto* n : g) {
    for (const auto* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      ret++;
    }
  }
  return ret;
}

// _____________________________________________________________________________
double Optimizer::solutionSpaceSize(const std::set<OptNode*>& g) {
  double ret = 1;
  for (const auto* n : g) {
    for (const auto* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      ret *= factorial(e->pl().getCardinality());
    }
  }
  return ret;
}
