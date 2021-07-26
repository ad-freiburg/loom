// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fstream>
#include "loom/optim/OptGraph.h"
#include "loom/optim/Optimizer.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Algorithm.h"
#include "util/log/Log.h"

using loom::optim::Optimizer;
using loom::optim::LinePair;
using loom::optim::PosComPair;
using loom::optim::EdgePair;
using loom::optim::OptGraph;
using loom::optim::OptNode;
using loom::optim::OptEdge;
using loom::graph::RenderGraph;
using loom::graph::HierarOrderCfg;
using shared::linegraph::NodeFront;
using shared::linegraph::LineNode;
using shared::linegraph::LineEdge;
using shared::linegraph::Line;
using util::geo::DPoint;
using util::geo::DLine;
using util::factorial;

// _____________________________________________________________________________
int Optimizer::optimize(RenderGraph* tg) const {
  // create optim graph
  OptGraph g(tg, _scorer);

  size_t maxC = maxCard(*g.getNds());
  double solSp = solutionSpaceSize(*g.getNds());
  LOG(DEBUG) << "Optimizing line graph of size " << tg->getNds()->size()
             << " with max cardinality = " << maxC
             << " and solution space size = " << solSp;

  g.partnerLines();

  if (_cfg->untangleGraph) {
    // do full untangling
    LOG(DEBUG) << "Untangling graph...";
    T_START(1);
    for (size_t i = 0; i < 2 * maxC; i++) {
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
    LOG(INFO) << "Outputting optimization graph to " << _cfg->dbgPath
              << "/optgraph.json";
    util::geo::output::GeoGraphJsonOutput out;
    std::ofstream fstr(_cfg->dbgPath + "/optgraph.json");
    out.print(g, fstr);
  }

  if (_cfg->outputStats) {
    LOG(INFO) << "(stats) Stats for <optim> graph of '" << tg << "'";
    LOG(INFO) << "(stats)   Total node count: " << g.getNumNodes();
    LOG(INFO) << "(stats)   Total edge count: " << g.getNumEdges();
    // LOG(INFO) << "(stats)   Total unique line count: " << g.getNumLines();
    LOG(INFO) << "(stats)   Max edge line cardinality: "
              << maxCard(*g.getNds());
    LOG(INFO) << "(stats)   Solution space: " << solutionSpaceSize(*g.getNds());
  }

  // iterate over components and optimize all of them separately
  const auto& comps = util::graph::Algorithm::connectedComponents(g);

  size_t runs = _cfg->optimRuns;
  double tSum = 0;
  double iterSum = 0;
  double scoreSum = 0;
  double crossSum = 0;
  double sepSum = 0;

  LOG(INFO) << "Optimization graph has " << comps.size() << " components.";

  for (size_t run = 0; run < runs; run++) {
    graph::OrderCfg c;
    graph::HierarOrderCfg hc;

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
      iters += optimizeComp(&g, nds, &hc);
    }

    double t = T_STOP(1);

    LOG(INFO) << " -- Optimization took " << t << " ms -- ";

    if (_cfg->outputStats) {
      LOG(INFO) << " (stats) Number of components: " << comps.size();
      LOG(INFO) << " (stats) Number of components with M=1: " << numM1Comps;
      LOG(INFO) << " (stats) Max number of nodes of all components: "
                << maxNumNodes;
      LOG(INFO) << " (stats) Max number of edges of all components: "
                << maxNumEdges;
      LOG(INFO) << " (stats) Max cardinality of all components: " << maxCompC;
      LOG(INFO) << " (stats) Max solution space size of all components: "
                << maxCompSolSpace;
    }

    tSum += t;
    iterSum += iters;

    hc.writeFlatCfg(&c);

    tg->setConfig(c);

    if (runs > 1) {
      if (_cfg->splittingOpt)
        scoreSum += _scorer->getScore();
      else
        scoreSum += _scorer->getCrossScore();
      crossSum += _scorer->getNumCrossings();
      sepSum += _scorer->getNumSeparations();
    }
  }

  if (runs > 1 && _cfg->outputStats) {
    LOG(INFO) << "";
    LOG(INFO) << "(multiple opt runs stats) avg time: " << tSum / (1.0 * runs)
              << " ms";
    LOG(INFO) << "(multiple opt runs stats) avg iters: "
              << iterSum / (1.0 * runs);
    LOG(INFO) << "(multiple opt runs stats) avg score: -- "
              << scoreSum / (1.0 * runs) << " --";
    LOG(INFO) << "(multiple opt runs stats) avg num crossings: -- "
              << crossSum / (1.0 * runs) << " --";
    LOG(INFO) << "(multiple opt runs stats) avg num separations: -- "
              << sepSum / (1.0 * runs) << " --";
    LOG(INFO) << "";
  }

  return 0;
}

// _____________________________________________________________________________
std::vector<LinePair> Optimizer::getLinePairs(OptEdge* segment) {
  return getLinePairs(segment, false);
}

// _____________________________________________________________________________
std::vector<LinePair> Optimizer::getLinePairs(OptEdge* segment, bool unique) {
  std::set<const Line*> processed;
  std::vector<LinePair> ret;
  for (auto& toA : segment->pl().getLines()) {
    processed.insert(toA.line);
    for (auto& toB : segment->pl().getLines()) {
      if (unique && processed.count(toB.line)) continue;
      if (toA.line == toB.line) continue;

      // this is to make sure that we always get the same line pairs in unique
      // mode -> take the smaller pointer first
      if (!unique || toA.line < toB.line)
        ret.push_back(LinePair(toA, toB));
      else
        ret.push_back(LinePair(toB, toA));
    }
  }
  return ret;
}

// _____________________________________________________________________________
bool Optimizer::crosses(OptGraph* g, OptNode* node, OptEdge* segmentA,
                        OptEdge* segmentB, PosComPair poscomb) {
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

  // NOTE: as soon as we use the non-geometric version of the crossing calc
  // here,
  // delete the getPos() method, as it is only used for this case. As soon as
  // this
  // is done, it is no longer necessary to give the OptGraph* g argument to this
  // (and other) functions, as g is only used to derive the original
  // RenderGraph
  // for the position calculation
  DPoint aInA = getPos(g, node, segmentA, posAinA);
  DPoint bInA = getPos(g, node, segmentA, posBinA);
  DPoint aInB = getPos(g, node, segmentB, posAinB);
  DPoint bInB = getPos(g, node, segmentB, posBinB);

  DLine a;
  a.push_back(aInA);
  a.push_back(aInB);

  DLine b;
  b.push_back(bInA);
  b.push_back(bInB);

  if (util::geo::intersects(aInA, aInB, bInA, bInB) ||
      util::geo::dist(a, b) < 1)
    pCrossing = true;

  // TODO(patrick): use this instead
  // pCrossing = (posAinA < posBinA && posAinB < posBinB);

  return pCrossing;
}

// _____________________________________________________________________________
bool Optimizer::crosses(OptGraph* g, OptNode* node, OptEdge* segmentA,
                        EdgePair segments, PosCom postcomb) {
  bool otherWayA =
      (segmentA->getFrom() != node) ^ segmentA->pl().etgs.front().dir;
  bool otherWayB = (segments.first->getFrom() != node) ^
                   segments.first->pl().etgs.front().dir;
  bool otherWayC = (segments.second->getFrom() != node) ^
                   segments.second->pl().etgs.front().dir;

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

  // bool retB = false;

  // if (posAinA > posBinA && posEdgeA < posEdgeB) retB = true;

  // return retB;

  DPoint aInA = getPos(g, node, segmentA, posAinA);
  DPoint bInA = getPos(g, node, segmentA, posBinA);

  bool ret = false;

  for (size_t i = 0; i < segments.first->pl().getCardinality(); ++i) {
    for (size_t j = 0; j < segments.second->pl().getCardinality(); ++j) {
      size_t posAinB = otherWayB ? cardB - 1 - i : i;
      size_t posBinC = otherWayC ? cardC - 1 - j : j;

      DPoint aInB = getPos(g, node, segments.first, posAinB);
      DPoint bInC = getPos(g, node, segments.second, posBinC);

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

  return ret;
}

// _____________________________________________________________________________
DPoint Optimizer::getPos(OptGraph* g, OptNode* n, OptEdge* segment, size_t p) {
  // look for correct nodefront
  const NodeFront* nf = 0;
  for (auto etg : segment->pl().etgs) {
    const NodeFront* test = n->pl().node->pl().frontFor(etg.etg);
    if (test) {
      nf = test;
      break;
    }
  }

  return g->getGraph()->linePosOn(*nf, segment->pl().etgs.front().etg, p, false,
                                  false);
}

// _____________________________________________________________________________
std::vector<OptEdge*> Optimizer::getEdgePartners(OptNode* node,
                                                 OptEdge* segmentA,
                                                 const LinePair& linepair) {
  std::vector<OptEdge*> ret;

  auto* fromEtg = OptGraph::getAdjEdg(segmentA, node);
  auto* dirA = fromEtg->pl().lineOcc(linepair.first.line).direction;
  auto* dirB = fromEtg->pl().lineOcc(linepair.second.line).direction;

  for (OptEdge* segmentB : node->getAdjList()) {
    if (segmentB == segmentA) continue;

    if (OptGraph::hasCtdLinesIn(linepair.first.line, dirA, segmentA,
                                segmentB) &&
        OptGraph::hasCtdLinesIn(linepair.second.line, dirB, segmentA,
                                segmentB)) {
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

  auto* fromEtg = OptGraph::getAdjEdg(segmentA, node);
  auto* dirA = fromEtg->pl().lineOcc(linepair.first.line).direction;
  auto* dirB = fromEtg->pl().lineOcc(linepair.second.line).direction;

  for (OptEdge* segmentB : node->getAdjList()) {
    if (segmentB == segmentA) continue;

    if (OptGraph::hasCtdLinesIn(linepair.first.line, dirA, segmentA,
                                segmentB)) {
      EdgePair curPair;
      curPair.first = segmentB;
      for (OptEdge* segmentC : node->getAdjList()) {
        if (segmentC == segmentA || segmentC == segmentB) continue;

        if (OptGraph::hasCtdLinesIn(linepair.second.line, dirB, segmentA,
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

// _____________________________________________________________________________
int Optimizer::optimizeComp(OptGraph* g, const std::set<OptNode*>& cmp,
                            HierarOrderCfg* c) const {
  return optimizeComp(g, cmp, c, 0);
}

// _____________________________________________________________________________
std::string Optimizer::prefix(size_t depth) {
  std::stringstream ret;
  for (size_t i = 0; i < depth * 2 + 1; i++) ret << " ";
  return ret.str();
}
