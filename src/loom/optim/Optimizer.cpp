// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fstream>
#include <numeric>
#include "loom/optim/NullOptimizer.h"
#include "loom/optim/OptGraph.h"
#include "loom/optim/OptGraphScorer.h"
#include "loom/optim/Optimizer.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Algorithm.h"
#include "util/log/Log.h"

using loom::optim::EdgePair;
using loom::optim::LinePair;
using loom::optim::NullOptimizer;
using loom::optim::OptEdge;
using loom::optim::OptGraph;
using loom::optim::OptGraphScorer;
using loom::optim::Optimizer;
using loom::optim::OptNode;
using loom::optim::OptOrderCfg;
using loom::optim::OptResStats;
using loom::optim::PosComPair;
using shared::linegraph::Line;
using shared::linegraph::LineNode;
using shared::rendergraph::HierarOrderCfg;
using shared::rendergraph::OrderCfg;
using shared::rendergraph::Ordering;
using shared::rendergraph::RenderGraph;
using util::factorial;
using util::geo::DLine;
using util::geo::DPoint;

// _____________________________________________________________________________
OptResStats Optimizer::optimize(RenderGraph* rg) const {
  // create optim graph
  OptGraph g(&_scorer);
  g.build(rg);

  OptResStats optResStats;

  size_t maxC = maxCard(*g.getNds());
  double solSp = solutionSpaceSize(*g.getNds());
  LOGTO(DEBUG, std::cerr) << "Optimizing line graph of size "
                          << rg->getNds()->size()
                          << " with max cardinality = " << maxC
                          << " and solution space size = " << solSp;

  g.partnerLines();

  if (_cfg->untangleGraph) {
    // do full untangling
    LOGTO(DEBUG, std::cerr) << "Untangling graph...";
    T_START(1);

    for (size_t i = 0; i <= maxC + 10; i++) {
      g.untangle();
      g.contractDeg2Nds();
      g.splitSingleLineEdgs();
      g.terminusDetach();
    }

    LOGTO(DEBUG, std::cerr) << "Done (" << T_STOP(1) << " ms)";
  } else if (_cfg->pruneGraph) {
    // only apply core graph rules
    T_START(1);
    LOGTO(DEBUG, std::cerr) << "Creating core optimization graph...";
    g.contractDeg2Nds();
    LOGTO(DEBUG, std::cerr) << "Done (" << T_STOP(1) << " ms)";
  }

  if (_cfg->outOptGraph) {
    LOGTO(INFO, std::cerr) << "Outputting optimization graph to "
                           << _cfg->dbgPath << "/optgraph.json";
    util::geo::output::GeoGraphJsonOutput out;
    std::ofstream fstr(_cfg->dbgPath + "/optgraph.json");
    out.print(g, fstr);
  }

  optResStats.numNodes = g.getNumNodes();
  optResStats.numEdges = g.getNumEdges();
  optResStats.maxLineCard = maxCard(*g.getNds());
  optResStats.solutionSpaceSize = solutionSpaceSize(*g.getNds());

  if (_cfg->outputStats) {
    LOGTO(INFO, std::cerr) << "(stats) Stats for <optim> graph of '" << rg
                           << "'";
    LOGTO(INFO, std::cerr) << "(stats)   Total node count: "
                           << optResStats.numNodes;
    LOGTO(INFO, std::cerr) << "(stats)   Total edge count: "
                           << optResStats.numEdges;
    LOGTO(INFO, std::cerr) << "(stats)   Max edge line cardinality: "
                           << optResStats.maxLineCard;
    LOGTO(INFO, std::cerr) << "(stats)   Solution space: "
                           << optResStats.solutionSpaceSize;
  }

  // iterate over components and optimize all of them separately
  const auto& comps = util::graph::Algorithm::connectedComponents(g);

  size_t nonTrivialComponents = 0;

  for (const auto& nds : comps) {
    // skip trivial components
    if (nds.size() < 3) continue;
    nonTrivialComponents++;
  }

  size_t runs = _cfg->optimRuns;
  double tSum = 0;
  double iterSum = 0;
  double scoreSum = 0;
  double crossSumSame = 0;
  double crossSumDiff = 0;
  double sepSum = 0;

  LOGTO(INFO, std::cerr) << "Optimization graph has " << nonTrivialComponents
                         << " nontrivial components.";

  // for trivial cases
  const NullOptimizer nullOpt(_cfg, _scorer.getPens());

  double bestScore = std::numeric_limits<double>::infinity();
  OrderCfg bestCfg;

  for (size_t run = 0; run < runs; run++) {
    OrderCfg c;
    HierarOrderCfg hc;

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

        // skip trivial components
        if (nds.size() > 2) {
          if (maxC > maxCompC) maxCompC = maxC;
          if (solSp > maxCompSolSpace) maxCompSolSpace = solSp;
          if (solSp == 1) numM1Comps++;
          if (nds.size() > maxNumNodes) maxNumNodes = nds.size();
          if (numEdges(nds) > maxNumEdges) maxNumEdges = numEdges(nds);

          LOGTO(INFO, std::cerr)
              << " (stats) Optimizing subgraph of size " << nds.size()
              << " with max cardinality = " << maxC
              << " and solution space size = " << solSp;
        }
      }

      // this is the implementation of the single edge pruning described in the
      // publication - simple skip such components
      if (nds.size() > 2) {
        iters += optimizeComp(&g, nds, &hc);
      } else {
        nullOpt.optimizeComp(&g, nds, &hc, 0);
      }
    }

    double t = T_STOP(1);

    optResStats.nonTrivialComponents = nonTrivialComponents;
    optResStats.numCompsSolSpaceOne = numM1Comps;
    optResStats.maxNumNodesPerComp = maxNumNodes;
    optResStats.maxNumEdgesPerComp = maxNumEdges;
    optResStats.maxCardPerComp = maxCompC;
    optResStats.maxCompSolSpace = maxCompSolSpace;

    if (_cfg->outputStats) {
      LOGTO(INFO, std::cerr) << "(stats) Number of nontrivial components: "
                             << optResStats.nonTrivialComponents;
      LOGTO(INFO, std::cerr)
          << "(stats) Number of nontrivial components with sol space size 1: "
          << optResStats.numCompsSolSpaceOne;
      LOGTO(INFO, std::cerr)
          << "(stats) Max number of nodes of all nontrivial components: "
          << optResStats.maxNumNodesPerComp;
      LOGTO(INFO, std::cerr)
          << "(stats) Max number of edges of all nontrivial components: "
          << optResStats.maxNumEdgesPerComp;
      LOGTO(INFO, std::cerr)
          << "(stats) Max cardinality of all nontrivial components: "
          << optResStats.maxCardPerComp;
      LOGTO(INFO, std::cerr)
          << "(stats) Max solution space size of all nontrivial components: "
          << optResStats.maxCompSolSpace;
    }

    hc.writeFlatCfg(&c);

    // fill in missing edges (which may have been pruned in the optim graph)
    // use the input ordering for these edges
    for (auto n : *rg->getNds()) {
      for (auto e : n->getAdjList()) {
        if (e->getFrom() != n) continue;
        if (c.find(e) == c.end()) {
          Ordering o(e->pl().getLines().size());
          std::iota(o.begin(), o.end(), 0);
          c[e] = o;
        }
      }
    }

    tSum += t;
    iterSum += iters;

    OptGraph gg(&_scorer);
    auto ndMap = gg.build(rg);

    auto optCfg = getOptOrderCfg(c, ndMap, &gg);

    double score = _scorer.getCrossingScore(&gg, optCfg);

    if (_scorer.optimizeSep())
      score += _scorer.getSeparationScore(&gg, optCfg);

    scoreSum += score;

    auto crossings = _scorer.getNumCrossings(&gg, optCfg);
    crossSumSame += crossings.first;
    crossSumDiff += crossings.second;

    size_t separations = _scorer.getNumSeparations(&gg, optCfg);

    sepSum += separations;

    if (score < bestScore) {
      bestCfg = c;
      bestScore = score;

      optResStats.score = score;
      optResStats.sameSegCrossings = crossings.first;
      optResStats.diffSegCrossings = crossings.second;
      optResStats.separations = separations;
    }
  }

  rg->writePermutation(bestCfg);

  optResStats.runs = runs;
  optResStats.avgSolveTime = tSum / (1.0 * runs);
  optResStats.avgIterations = iterSum / (1.0 * runs);
  optResStats.avgScore = scoreSum / (1.0 * runs);
  optResStats.avgSameSegCross = crossSumSame / (1.0 * runs);
  optResStats.avgDiffSegCross = crossSumDiff / (1.0 * runs);
  optResStats.avgSeps = sepSum / (1.0 * runs);

  if (_cfg->outputStats) {
    LOGTO(INFO, std::cerr) << "";
    LOGTO(INFO, std::cerr) << "(stats) After " << optResStats.runs
                           << " run(s):";
    LOGTO(INFO, std::cerr) << "(stats) avg solve time: "
                           << optResStats.avgSolveTime << " ms";
    LOGTO(INFO, std::cerr) << "(stats) avg iters: "
                           << optResStats.avgIterations;
    LOGTO(INFO, std::cerr) << "(stats) avg score: -- " << optResStats.avgScore
                           << " --";
    LOGTO(INFO, std::cerr) << "(stats) avg num crossings: -- "
                           << optResStats.avgSameSegCross +
                                  optResStats.avgDiffSegCross
                           << " -- (" << optResStats.avgSameSegCross
                           << " same, " << optResStats.avgDiffSegCross
                           << " diff)";
    LOGTO(INFO, std::cerr) << "(stats) avg num separations: -- "
                           << optResStats.avgSeps << " --";
    LOGTO(INFO, std::cerr) << "";
  }

  return optResStats;
}

// _____________________________________________________________________________
std::vector<LinePair> Optimizer::getLinePairs(OptEdge* segment) {
  return getLinePairs(segment, false);
}

// _____________________________________________________________________________
std::vector<LinePair> Optimizer::getLinePairs(OptEdge* segment, bool unique) {
  std::vector<LinePair> ret;
  for (size_t a = 0; a < segment->pl().getLines().size(); a++) {
    for (size_t b = a + 1; b < segment->pl().getLines().size(); b++) {
      const auto& loA = segment->pl().getLines()[a];
      const auto& loB = segment->pl().getLines()[b];

      if (!unique) {
        ret.push_back(LinePair(loA, loB));
        ret.push_back(LinePair(loB, loA));
        continue;
      }

      if (loA.line < loB.line) {
        ret.push_back(LinePair(loA, loB));
      } else {
        ret.push_back(LinePair(loB, loA));
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
bool Optimizer::crosses(OptNode* node, OptEdge* segA, OptEdge* segB,
                        PosComPair poscom) {
  bool revA = (segA->getFrom() != node) ^ segA->pl().lnEdgParts.front().dir;
  bool revB = (segB->getFrom() != node) ^ segB->pl().lnEdgParts.front().dir;

  size_t carA = segA->pl().getCardinality();
  size_t carB = segB->pl().getCardinality();

  size_t pAinA = revA ? carA - 1 - poscom.first.first : poscom.first.first;
  size_t pAinB = revB ? carB - 1 - poscom.first.second : poscom.first.second;
  size_t pBinA = revA ? carA - 1 - poscom.second.first : poscom.second.first;
  size_t pBinB = revB ? carB - 1 - poscom.second.second : poscom.second.second;

  return !((pAinA < pBinA) ^ (pAinB < pBinB));
}

// _____________________________________________________________________________
bool Optimizer::crosses(OptNode* node, OptEdge* segA, EdgePair segments,
                        PosCom poscomb) {
  bool revA = (segA->getFrom() != node) ^ segA->pl().lnEdgParts.front().dir;

  size_t carA = segA->pl().getCardinality();

  size_t pAinA = revA ? carA - 1 - poscomb.first : poscomb.first;
  size_t pBinA = revA ? carA - 1 - poscomb.second : poscomb.second;

  size_t pSegA = node->pl().circOrder(segA);
  size_t pEdgeA = node->pl().circOrder(segments.first);
  size_t pEdgeB = node->pl().circOrder(segments.second);

  // make position relative to segment A
  if (pEdgeA > pSegA)
    pEdgeA = pEdgeA - pSegA;
  else
    pEdgeA = (node->getDeg() - pSegA) + pEdgeA;

  if (pEdgeB > pSegA)
    pEdgeB = pEdgeB - pSegA;
  else
    pEdgeB = (node->getDeg() - pSegA) + pEdgeB;

  return !((pAinA > pBinA) ^ (pEdgeA < pEdgeB));
}

// _____________________________________________________________________________
bool Optimizer::separates(PosComPair poscom) {
  size_t pAinA = poscom.first.first;
  size_t pAinB = poscom.first.second;
  size_t pBinA = poscom.second.first;
  size_t pBinB = poscom.second.second;

  return ((pAinA > pBinA ? pAinA - pBinA : pBinA - pAinA) == 1 &&
          (pAinB > pBinB ? pAinB - pBinB : pBinB - pAinB) != 1) ||
         ((pAinA > pBinA ? pAinA - pBinA : pBinA - pAinA) != 1 &&
          (pAinB > pBinB ? pAinB - pBinB : pBinB - pAinB) == 1);
}

// _____________________________________________________________________________
std::vector<OptEdge*> Optimizer::getEdgePartners(OptNode* node, OptEdge* segA,
                                                 const LinePair& linepair) {
  std::vector<OptEdge*> ret;

  auto* fromLnEdg = OptGraph::getAdjEdg(segA, node);
  auto* dirA = fromLnEdg->pl().lineOcc(linepair.first.line).direction;
  auto* dirB = fromLnEdg->pl().lineOcc(linepair.second.line).direction;

  for (OptEdge* segB : node->getAdjList()) {
    if (segB == segA) continue;

    if (OptGraph::hasCtdLineIn(linepair.first.line, dirA, segA, segB) &&
        OptGraph::hasCtdLineIn(linepair.second.line, dirB, segA, segB)) {
      ret.push_back(segB);
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::vector<EdgePair> Optimizer::getEdgePartnerPairs(OptNode* node,
                                                     OptEdge* segA,
                                                     const LinePair& linepair) {
  std::vector<EdgePair> ret;
  if (node->getDeg() < 3) return ret;

  auto fromLnEdg = OptGraph::getAdjEdg(segA, node);
  auto dirA = fromLnEdg->pl().lineOcc(linepair.first.line).direction;
  auto dirB = fromLnEdg->pl().lineOcc(linepair.second.line).direction;

  for (OptEdge* segB : node->getAdjList()) {
    if (segB == segA) continue;

    if (OptGraph::hasCtdLineIn(linepair.first.line, dirA, segA, segB)) {
      EdgePair curPair;
      curPair.first = segB;
      for (OptEdge* segmentC : node->getAdjList()) {
        if (segmentC == segA || segmentC == segB) continue;

        if (OptGraph::hasCtdLineIn(linepair.second.line, dirB, segA,
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
OptOrderCfg Optimizer::getOptOrderCfg(
    const shared::rendergraph::OrderCfg& cfg,
    const std::map<const LineNode*, OptNode*>& ndMap, const OptGraph* g) {
  OptOrderCfg ret;
  for (auto i : cfg) {
    auto e = i.first;
    auto order = i.second;

    auto opNdFr = ndMap.find(e->getFrom())->second;
    auto opNdTo = ndMap.find(e->getTo())->second;
    auto opEdg = g->getEdg(opNdFr, opNdTo);

    for (auto pos = order.rbegin(); pos != order.rend(); pos++) {
      auto lo = e->pl().lineOccAtPos(*pos);
      ret[opEdg].push_back(lo.line);
    }
  }

  return ret;
}

// _____________________________________________________________________________
std::string Optimizer::prefix(size_t depth) {
  std::stringstream ret;
  for (size_t i = 0; i < depth * 2 + 1; i++) ret << " ";
  return ret.str();
}
