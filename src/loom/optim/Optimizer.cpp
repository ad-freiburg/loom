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

using loom::optim::EdgePair;
using loom::optim::LinePair;
using loom::optim::OptEdge;
using loom::optim::OptGraph;
using loom::optim::Optimizer;
using loom::optim::OptNode;
using loom::optim::PosComPair;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineNode;
using shared::linegraph::NodeFront;
using shared::rendergraph::HierarOrderCfg;
using shared::rendergraph::OrderCfg;
using shared::rendergraph::RenderGraph;
using util::factorial;
using util::geo::DLine;
using util::geo::DPoint;

// _____________________________________________________________________________
int Optimizer::optimize(RenderGraph* rg) const {
  // create optim graph
  OptGraph g(rg, _scorer);

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
    for (size_t i = 0; i < 2 * maxC; i++) {
      g.untangle();
      g.simplify();
      g.split();
    }
    LOGTO(DEBUG, std::cerr) << "Done (" << T_STOP(1) << " ms)";
  } else if (_cfg->createCoreOptimGraph) {
    // only apply core graph rules
    T_START(1);
    LOGTO(DEBUG, std::cerr) << "Creating core optimization graph...";
    g.simplify();
    LOGTO(DEBUG, std::cerr) << "Done (" << T_STOP(1) << " ms)";
  }

  if (_cfg->outOptGraph) {
    LOGTO(INFO, std::cerr) << "Outputting optimization graph to "
                           << _cfg->dbgPath << "/optgraph.json";
    util::geo::output::GeoGraphJsonOutput out;
    std::ofstream fstr(_cfg->dbgPath + "/optgraph.json");
    out.print(g, fstr);
  }

  if (_cfg->outputStats) {
    LOGTO(INFO, std::cerr) << "(stats) Stats for <optim> graph of '" << rg
                           << "'";
    LOGTO(INFO, std::cerr) << "(stats)   Total node count: " << g.getNumNodes();
    LOGTO(INFO, std::cerr) << "(stats)   Total edge count: " << g.getNumEdges();
    // LOGTO(INFO,std::cerr) << "(stats)   Total unique line count: " <<
    // g.getNumLines();
    LOGTO(INFO, std::cerr) << "(stats)   Max edge line cardinality: "
                           << maxCard(*g.getNds());
    LOGTO(INFO, std::cerr) << "(stats)   Solution space: "
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

  LOGTO(INFO, std::cerr) << "Optimization graph has " << comps.size()
                         << " components.";

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
      iters += optimizeComp(&g, nds, &hc);
    }

    double t = T_STOP(1);

    LOGTO(INFO, std::cerr) << " -- Optimization took " << t << " ms -- ";

    if (_cfg->outputStats) {
      LOGTO(INFO, std::cerr)
          << " (stats) Number of components: " << comps.size();
      LOGTO(INFO, std::cerr)
          << " (stats) Number of components with M=1: " << numM1Comps;
      LOGTO(INFO, std::cerr)
          << " (stats) Max number of nodes of all components: " << maxNumNodes;
      LOGTO(INFO, std::cerr)
          << " (stats) Max number of edges of all components: " << maxNumEdges;
      LOGTO(INFO, std::cerr)
          << " (stats) Max cardinality of all components: " << maxCompC;
      LOGTO(INFO, std::cerr)
          << " (stats) Max solution space size of all components: "
          << maxCompSolSpace;
    }

    tSum += t;
    iterSum += iters;

    hc.writeFlatCfg(&c);

    rg->writePermutation(c);

    if (runs > 1) {
      if (_cfg->splittingOpt)
        scoreSum += _scorer->getScore(c);
      else
        scoreSum += _scorer->getCrossScore(c);
      crossSum += _scorer->getNumCrossings(c);
      sepSum += _scorer->getNumSeparations(c);
    }
  }

  if (runs > 1 && _cfg->outputStats) {
    LOGTO(INFO, std::cerr) << "";
    LOGTO(INFO, std::cerr) << "(multiple opt runs stats) avg time: "
                           << tSum / (1.0 * runs) << " ms";
    LOGTO(INFO, std::cerr) << "(multiple opt runs stats) avg iters: "
                           << iterSum / (1.0 * runs);
    LOGTO(INFO, std::cerr) << "(multiple opt runs stats) avg score: -- "
                           << scoreSum / (1.0 * runs) << " --";
    LOGTO(INFO, std::cerr) << "(multiple opt runs stats) avg num crossings: -- "
                           << crossSum / (1.0 * runs) << " --";
    LOGTO(INFO, std::cerr)
        << "(multiple opt runs stats) avg num separations: -- "
        << sepSum / (1.0 * runs) << " --";
    LOGTO(INFO, std::cerr) << "";
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
bool Optimizer::crosses(OptNode* node, OptEdge* segmentA, OptEdge* segmentB,
                        PosComPair poscom) {
  bool revA = (segmentA->getFrom() != node) ^ segmentA->pl().etgs.front().dir;
  bool revB = (segmentB->getFrom() != node) ^ segmentB->pl().etgs.front().dir;

  size_t carA = segmentA->pl().getCardinality();
  size_t carB = segmentB->pl().getCardinality();

  size_t pAinA = revA ? carA - 1 - poscom.first.first : poscom.first.first;
  size_t pAinB = revB ? carB - 1 - poscom.first.second : poscom.first.second;
  size_t pBinA = revA ? carA - 1 - poscom.second.first : poscom.second.first;
  size_t pBinB = revB ? carB - 1 - poscom.second.second : poscom.second.second;

  return pAinA < pBinA && pAinB < pBinB;
}

// _____________________________________________________________________________
bool Optimizer::crosses(OptNode* node, OptEdge* segmentA, EdgePair segments,
                        PosCom poscomb) {
  bool revA = (segmentA->getFrom() != node) ^ segmentA->pl().etgs.front().dir;

  size_t carA = segmentA->pl().getCardinality();

  size_t pAinA = revA ? carA - 1 - poscomb.first : poscomb.first;
  size_t pBinA = revA ? carA - 1 - poscomb.second : poscomb.second;

  size_t pSegmentA = std::find(node->pl().orderedEdges.begin(),
                               node->pl().orderedEdges.end(), segmentA) -
                     node->pl().orderedEdges.begin();

  size_t pEdgeA = std::find(node->pl().orderedEdges.begin(),
                            node->pl().orderedEdges.end(), segments.first) -
                  node->pl().orderedEdges.begin();

  size_t pEdgeB = std::find(node->pl().orderedEdges.begin(),
                            node->pl().orderedEdges.end(), segments.second) -
                  node->pl().orderedEdges.begin();

  // make position relative to segment A
  if (pEdgeA > pSegmentA)
    pEdgeA = pEdgeA - pSegmentA;
  else
    pEdgeA = (node->getDeg() - pSegmentA) + pEdgeA;

  if (pEdgeB > pSegmentA)
    pEdgeB = pEdgeB - pSegmentA;
  else
    pEdgeB = (node->getDeg() - pSegmentA) + pEdgeB;

  return pAinA > pBinA && pEdgeA < pEdgeB;
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

  auto fromEtg = OptGraph::getAdjEdg(segmentA, node);
  auto dirA = fromEtg->pl().lineOcc(linepair.first.line).direction;
  auto dirB = fromEtg->pl().lineOcc(linepair.second.line).direction;

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
