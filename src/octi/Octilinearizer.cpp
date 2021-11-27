// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <thread>
#include "ilp/ILPGridOptimizer.h"
#include "octi/Octilinearizer.h"
#include "octi/basegraph/BaseGraph.h"
#include "octi/basegraph/ConvexHullOctiGridGraph.h"
#include "octi/basegraph/GridGraph.h"
#include "octi/basegraph/HexGridGraph.h"
#include "octi/basegraph/NodeCost.h"
#include "octi/basegraph/OctiGridGraph.h"
#include "octi/basegraph/OctiHananGraph.h"
#include "octi/basegraph/OctiQuadTree.h"
#include "octi/basegraph/OrthoRadialGraph.h"
#include "octi/basegraph/PseudoOrthoRadialGraph.h"
#include "octi/combgraph/Drawing.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/BiDijkstra.h"
#include "util/graph/Dijkstra.h"
#include "util/log/Log.h"

using namespace octi;
using namespace basegraph;
using namespace util;

using combgraph::EdgeOrdering;
using octi::basegraph::BaseGraph;
using octi::combgraph::Drawing;
using octi::config::OrderMethod;
using octi::ilp::ILPStats;
using util::geo::DBox;
using util::geo::dist;
using util::geo::DPoint;
using util::geo::DPolygon;
using util::geo::len;
using util::geo::MultiPoint;
using util::geo::Polygon;
using util::graph::BiDijkstra;
using util::graph::Dijkstra;

// _____________________________________________________________________________
Score Octilinearizer::drawILP(
    const CombGraph& cg, const util::geo::DBox& box, LineGraph* outTg,
    BaseGraph** retGg, Drawing* dOut, const Penalties& pens, double gridSize,
    double borderRad, double maxGrDist, OrderMethod orderMethod, bool noSolve,
    double enfGeoPen, size_t hananIters, int timeLim,
    const std::string& cacheDir, double cacheThreshold, int numThreads,
    octi::ilp::ILPStats* stats, const std::string& solverStr,
    const std::string& path) {
  BaseGraph* gg;
  Drawing drawing;

  // always set density penality to 0, cannot by used in ILP and prevents proper
  // presolve by our approximate approach
  Penalties pensCpy = pens;
  pensCpy.densityPen = 0;

  LOGTO(DEBUG, std::cerr) << "Presolving...";
  try {
    // presolve using heuristical approach to get a first feasible solution
    LineGraph tmpOutTg;
    // important: always use restrLocSearch here!
    auto score = draw(cg, box, &tmpOutTg, &gg, &drawing, pensCpy, gridSize,
                      borderRad, maxGrDist, orderMethod, true, enfGeoPen,
                      hananIters, {}, 100, std::numeric_limits<size_t>::max());
    if (score.violations) throw NoEmbeddingFoundExc();
    LOGTO(DEBUG, std::cerr) << "Presolving finished.";
  } catch (const NoEmbeddingFoundExc& exc) {
    LOGTO(DEBUG, std::cerr) << "Presolve was not successful.";
    gg = newBaseGraph(box, cg, gridSize, borderRad, hananIters, pensCpy);
    gg->init();
    drawing = Drawing(gg);
  }

  GeoPensMap enfGeoPens;
  const GeoPensMap* geoPens = 0;

  if (enfGeoPen) {
    LOGTO(DEBUG, std::cerr) << "Writing geopens... ";
    auto edges = getOrdering(cg, OrderMethod::NUM_LINES);
    T_START(geopens);
    for (auto cmbEdg : edges) {
      gg->writeGeoCoursePens(cmbEdg, &enfGeoPens, enfGeoPen);
    }
    LOGTO(DEBUG, std::cerr) << "Done. (" << T_STOP(geopens) << "ms)";
    geoPens = &enfGeoPens;
  }

  // TODO
  // if (obstacles.size()) {
  // std::cerr << "Writing obstacles... ";
  // T_START(obstacles);
  // for (auto gg : ggs) for (const auto& obst : obstacles)
  // gg->addObstacle(obst);
  // std::cerr << " done (" << T_STOP(obstacles) << "ms)" << std::endl;
  // }

  ilp::ILPGridOptimizer ilpoptim;

  *stats =
      ilpoptim.optimize(gg, cg, &drawing, maxGrDist, noSolve, geoPens, timeLim,
                        cacheDir, cacheThreshold, numThreads, solverStr, path);

  drawing.getLineGraph(outTg);
  *retGg = gg;
  *dOut = drawing;

  Score a;
  a.full = stats->score;

  return a;
}

// _____________________________________________________________________________
Score Octilinearizer::draw(const CombGraph& cg, const DBox& box,
                           LineGraph* outTg, BaseGraph** retGg, Drawing* dOut,
                           const Penalties& pens, double gridSize,
                           double borderRad, double maxGrDist,
                           OrderMethod orderMethod, bool restrLocSearch,
                           double enfGeoPen, size_t hananIters,
                           const std::vector<Polygon<double>>& obstacles,
                           size_t locSearchIters, size_t abortAfter) {
  size_t jobs = 4;
  std::vector<BaseGraph*> ggs(jobs);

  LOGTO(DEBUG, std::cerr) << "Creating grid graphs... ";
  T_START(ggraph);
#pragma omp parallel for
  for (size_t i = 0; i < jobs; i++) {
    ggs[i] = newBaseGraph(box, cg, gridSize, borderRad, hananIters, pens);
    ggs[i]->init();
  }

  LOGTO(DEBUG, std::cerr) << "Done. (" << T_STOP(ggraph) << "ms)";

  LOGTO(DEBUG, std::cerr) << "Grid graph has " << ggs[0]->getNds()->size()
                          << " nodes";

  size_t LOCAL_SEARCH_ITERS = locSearchIters;
  double CONVERGENCE_THRESHOLD = 0.05;

  GeoPensMap enfGeoPens;
  const GeoPensMap* geoPens = 0;

  // ordering is irrelevant, this is a just a shortcut to get all edges
  auto edges = getOrdering(cg, OrderMethod::NUM_LINES);

  if (enfGeoPen > 0) {
    LOGTO(DEBUG, std::cerr) << "Writing geopens... ";
    T_START(geopens);
    for (auto cmbEdg : edges) {
      ggs[0]->writeGeoCoursePens(cmbEdg, &enfGeoPens, enfGeoPen);
    }
    LOGTO(DEBUG, std::cerr) << "Done. (" << T_STOP(geopens) << "ms)";
    geoPens = &enfGeoPens;
  }

  if (obstacles.size()) {
    LOGTO(DEBUG, std::cerr) << "Writing obstacles... ";
    T_START(obstacles);
    for (auto gg : ggs)
      for (const auto& obst : obstacles) gg->addObstacle(obst);
    LOGTO(DEBUG, std::cerr) << "Done. (" << T_STOP(obstacles) << "ms)";
  }

  // this is the best drawing
  Drawing drawing(ggs[0]);

  // try our default edge ordering first, without any randomization

  std::vector<OrderMethod> methods = {
      OrderMethod::NUM_LINES,     OrderMethod::LENGTH,
      OrderMethod::ADJ_ND_DEGREE, OrderMethod::ADJ_ND_LDEGREE,
      OrderMethod::GROWTH_DEG,    OrderMethod::GROWTH_LDEG};

  if (orderMethod != OrderMethod::ALL) {
    methods = {orderMethod};
  }

  std::vector<std::vector<OrderMethod>> batches(jobs);
  for (size_t i = 0; i < methods.size(); i++) {
    batches[i % jobs].push_back(methods[i]);
  }

  LOGTO(DEBUG, std::cerr) << "Searching initial drawing... ";

#pragma omp parallel for
  for (size_t btch = 0; btch < jobs; btch++) {
    for (OrderMethod meth : batches[btch]) {
      T_START(draw);
      Drawing drawingCp(ggs[btch]);

      // get a randomized ordering
      std::vector<CombEdge*> iterOrder = getOrdering(cg, meth);

      double bestScoreSoFar = 0;

#pragma omp critical
      { bestScoreSoFar = drawing.score(); }

      auto status = draw(iterOrder, ggs[btch], &drawingCp, bestScoreSoFar,
                         maxGrDist, geoPens, abortAfter);

      drawingCp.eraseFromGrid(ggs[btch]);

      statLine(status, std::string("Try ") + std::to_string(meth), drawingCp,
               T_STOP(draw), "*");

#pragma omp critical
      {
        if (status == DRAWN && drawingCp.score() < drawing.score()) {
          drawing = drawingCp;
        } else {
          drawingCp.crumble();
        }
      }
    }
  }

  if (drawing.score() == INF) throw NoEmbeddingFoundExc();

  LOGTO(DEBUG, std::cerr) << "Done.";

  for (size_t i = 0; i < jobs; i++) drawing.applyToGrid(ggs[i]);

  size_t iters = 0;

  LOGTO(DEBUG, std::cerr) << "Initial score: " << drawing.score() << " ("
                          << drawing.violations() << " topology violations).";
  LOGTO(DEBUG, std::cerr) << "Starting local search...";

  // dont use local search if abortAfter is set
  if (abortAfter != std::numeric_limits<size_t>::max()) LOCAL_SEARCH_ITERS = 0;

  std::vector<std::vector<CombNode*>> batchesLoc(jobs);
  size_t c = 0;
  for (auto nd : cg.getNds()) {
    if (nd->getDeg() == 0) continue;
    batchesLoc[c % jobs].push_back(nd);
    c++;
  }

  for (; iters < LOCAL_SEARCH_ITERS; iters++) {
    T_START(iter);
    std::vector<Drawing> bestFrIters(jobs);

#pragma omp parallel for
    for (size_t btch = 0; btch < jobs; btch++) {
      for (auto a : batchesLoc[btch]) {
        Drawing drawingCp = drawing;

        // use the batches grid graph
        drawingCp.setBaseGraph(ggs[btch]);

        // reverting a
        std::vector<CombEdge*> test;
        for (auto ce : a->getAdjList()) {
          test.push_back(ce);

          drawingCp.eraseFromGrid(ce, ggs[btch]);
          drawingCp.erase(ce);
        }

        drawingCp.erase(a);
        ggs[btch]->unSettleNd(a);

        for (size_t pos = 0; pos < ggs[btch]->maxDeg() + 1; pos++) {
          SettledPos p;

          auto n = ggs[btch]->neigh(drawing.getGrNd(a), pos);
          if (!n) continue;

          p[a] = n;

          if (restrLocSearch) {
            // dont try positions outside the move radius for consistency with
            // ILP approach
            double gridD = dist(*a->pl().getGeom(), *n->pl().getGeom());
            double maxDis = ggs[btch]->getCellSize() * maxGrDist;
            if (gridD >= maxDis) continue;
          }

          Drawing run = drawingCp;

          // we can use bestFromIter.score() as the limit for the shortest
          // path computation, as we can already do at least as good.
          auto error =
              draw(test, p, ggs[btch], &run, bestFrIters[btch].score(),
                   maxGrDist, geoPens, std::numeric_limits<size_t>::max());

          if (!error && bestFrIters[btch].score() > run.score()) {
            bestFrIters[btch] = run;
          }

          // reset grid
          for (auto ce : a->getAdjList()) run.eraseFromGrid(ce, ggs[btch]);
          if (ggs[btch]->isSettled(a)) ggs[btch]->unSettleNd(a);
        }

        ggs[btch]->settleNd(const_cast<GridNode*>(ggs[btch]->getGrNdById(
                                drawing.getGrNd(a)->pl().getId())),
                            a);

        // re-settle edges
        for (auto ce : a->getAdjList()) drawing.applyToGrid(ce, ggs[btch]);
      }
    }

    size_t bestCore;
    double bestScore = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < jobs; i++) {
      if (bestFrIters[i].score() < bestScore) {
        bestScore = bestFrIters[i].score();
        bestCore = i;
      }
    }

    double imp = (drawing.score() - bestFrIters[bestCore].score());
    LOGTO(DEBUG, std::cerr)
        << " ++ Iter " << iters << ", prev " << drawing.score() << ", next "
        << bestFrIters[bestCore].score() << " (" << (imp >= 0 ? "+" : "") << imp
        << ", " << T_STOP(iter) << " ms)";

    for (size_t i = 0; i < jobs; i++) {
      drawing.eraseFromGrid(ggs[i]);
      bestFrIters[bestCore].applyToGrid(ggs[i]);
    }
    drawing = bestFrIters[bestCore];

    if (imp < CONVERGENCE_THRESHOLD) break;
  }

  drawing.getLineGraph(outTg);
  auto fullScore = drawing.fullScore();
  LOGTO(DEBUG, std::cerr) << "Topo violations: " << drawing.violations()
                          << ", hop costs: " << fullScore.hop
                          << ", bend costs: " << fullScore.bend
                          << ", mv costs: " << fullScore.move
                          << ", dense costs: " << fullScore.dense;

  *retGg = ggs[0];
  *dOut = drawing;

  // the drawing might still have another internal grid graph, make sure they
  // match (this is important for drawILP)
  dOut->setBaseGraph(ggs[0]);
  fullScore.iters = iters;
  return fullScore;
}

// _____________________________________________________________________________
void Octilinearizer::settleRes(GridNode* frGrNd, GridNode* toGrNd,
                               BaseGraph* gg, CombNode* from, CombNode* to,
                               const GrEdgList& res, CombEdge* e,
                               size_t rndrOrder) {
  gg->settleNd(toGrNd, to);
  gg->settleNd(frGrNd, from);

  // balance edges
  for (auto f : res) {
    if (f->pl().isSecondary()) continue;
    gg->settleEdg(f->getFrom()->pl().getParent(), f->getTo()->pl().getParent(),
                  e, rndrOrder);
  }
}

// _____________________________________________________________________________
void Octilinearizer::writeNdCosts(GridNode* n, CombNode* origNode, CombEdge* e,
                                  BaseGraph* g) {
  NodeCost c;
  c += g->topoBlockPen(n, origNode, e);
  c += g->spacingPen(n, origNode, e);
  c += g->nodeBendPen(n, origNode, e);

  g->addCostVec(n, c);
}

// _____________________________________________________________________________
Undrawable Octilinearizer::draw(const std::vector<CombEdge*>& order,
                                BaseGraph* gg, Drawing* drawing, double cutoff,
                                double maxGrDist, const GeoPensMap* geoPensMap,
                                size_t abortAfter) {
  SettledPos emptyPos;
  return draw(order, emptyPos, gg, drawing, cutoff, maxGrDist, geoPensMap,
              abortAfter);
}

// _____________________________________________________________________________
Undrawable Octilinearizer::draw(const std::vector<CombEdge*>& ord,
                                const SettledPos& settled, BaseGraph* gg,
                                Drawing* drawing, double globCutoff,
                                double maxGrDist, const GeoPensMap* geoPensMap,
                                size_t abortAfter) {
  SettledPos retPos;

  size_t i = 0;

  for (auto cmbEdg : ord) {
    double cutoff = globCutoff - drawing->score();
    i++;
    if (drawing->score() == std::numeric_limits<double>::infinity()) {
      cutoff = drawing->score();
    }
    bool rev = false;
    auto frCmbNd = cmbEdg->getFrom();
    auto toCmbNd = cmbEdg->getTo();

    std::set<GridNode*> frGrNds, toGrNds;

    std::tie(frGrNds, toGrNds) =
        getRtPair(frCmbNd, toCmbNd, settled, gg, maxGrDist);

    if (frGrNds.size() == 0 || toGrNds.size() == 0) return NO_CANDS;

    if (toGrNds.size() > frGrNds.size()) {
      auto tmp = frCmbNd;
      frCmbNd = toCmbNd;
      toCmbNd = tmp;
      auto tmp2 = frGrNds;
      frGrNds = toGrNds;
      toGrNds = tmp2;
      rev = true;
    }

    // if we open node sinks, we have to offset their cost by the highest
    // possible turn cost + 1 to not distort turn penalties
    double costOffsetFrom = 0;
    double costOffsetTo = 0;

    // open the source nodes
    for (auto n : frGrNds) {
      if (gg->isSettled(frCmbNd)) {
        // only count displacement penalty ONCE
        gg->openSinkFr(n, 0);
      } else {
        costOffsetFrom = (gg->getPens().p_45 - gg->getPens().p_135);
        gg->openSinkFr(n, costOffsetFrom + gg->ndMovePen(frCmbNd, n));
      }
    }

    // open the target nodes
    for (auto n : toGrNds) {
      if (gg->isSettled(toCmbNd)) {
        // only count displacement penalty ONCE
        gg->openSinkTo(n, 0);
      } else {
        costOffsetTo = (gg->getPens().p_45 - gg->getPens().p_135);
        gg->openSinkTo(n, costOffsetTo + gg->ndMovePen(toCmbNd, n));
      }
    }

    // IMPORTANT: node costs are only written to sinks if they are already
    // settled. There is no need to add node costs before, as they handle
    // relations between two or more adjacent edges. If the node has not
    // already been settled, such a relation does not exist.
    //
    // Even more importantly, is a node is settled, its turn edges have
    // already been closed.
    //
    // the size() == 1 check is important, because nd cost writing will
    // not work if the to node is not already settled!

    if (frGrNds.size() == 1 && gg->isSettled(frCmbNd)) {
      writeNdCosts(*frGrNds.begin(), frCmbNd, cmbEdg, gg);
    }

    if (toGrNds.size() == 1 && gg->isSettled(toCmbNd)) {
      writeNdCosts(*toGrNds.begin(), toCmbNd, cmbEdg, gg);
    }

    GrEdgList eL;
    GrNdList nL;
    GridNode* toGrNd = 0;
    GridNode* frGrNd = 0;

    auto heur = gg->getHeur(toGrNds);

    if (geoPensMap) {
      // init cost function with geo distance penalties
      auto cost = GridCostGeoPen(cutoff + costOffsetTo + costOffsetFrom,
                                 &geoPensMap->find(cmbEdg)->second);
      Dijkstra::shortestPath(frGrNds, toGrNds, cost, *heur, &eL, &nL);
    } else {
      auto cost = GridCost(cutoff + costOffsetTo + costOffsetFrom);

      Dijkstra::shortestPath(frGrNds, toGrNds, cost, *heur, &eL, &nL);
    }

    delete heur;

    if (!nL.size()) {
      // cleanup
      for (auto n : toGrNds) gg->closeSinkTo(n);
      for (auto n : frGrNds) gg->closeSinkFr(n);

      return NO_PATH;
    }

    toGrNd = nL.front();
    frGrNd = nL.back();

    // remove the cost offsets to not distort final costs
    eL.front()->pl().setCost(eL.front()->pl().cost() - costOffsetTo);
    eL.back()->pl().setCost(eL.back()->pl().cost() - costOffsetFrom);

    // draw
    drawing->draw(cmbEdg, eL, rev);

    if (i > abortAfter) {
      settleRes(frGrNd, toGrNd, gg, frCmbNd, toCmbNd, eL, cmbEdg, i);
      return DRAWN;
    }

    // close the source and target node
    for (auto n : toGrNds) gg->closeSinkTo(n);
    for (auto n : frGrNds) gg->closeSinkFr(n);

    settleRes(frGrNd, toGrNd, gg, frCmbNd, toCmbNd, eL, cmbEdg, i);
  }

  return DRAWN;
}

// _____________________________________________________________________________
std::vector<CombEdge*> Octilinearizer::getOrdering(
    const CombGraph& cg, config::OrderMethod method) const {
  if (method == OrderMethod::GROWTH_DEG)
    return getGrowthOrder<EdgeCmpDeg, NodeCmpDeg>(cg);
  else if (method == OrderMethod::GROWTH_LDEG)
    return getGrowthOrder<EdgeCmpLdeg, NodeCmpLdeg>(cg);

  std::vector<CombEdge*> retOrder;

  for (auto n : cg.getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      retOrder.push_back(e);
    }
  }

  if (method == OrderMethod::ADJ_ND_DEGREE) {
    EdgeCmpDeg cmp;
    std::sort(retOrder.begin(), retOrder.end(), cmp);
  } else if (method == OrderMethod::ADJ_ND_LDEGREE) {
    EdgeCmpLdeg cmp;
    std::sort(retOrder.begin(), retOrder.end(), cmp);
  } else if (method == OrderMethod::NUM_LINES) {
    EdgeCmpNumLines cmp;
    std::sort(retOrder.begin(), retOrder.end(), cmp);
  } else if (method == OrderMethod::LENGTH) {
    EdgeCmpLength cmp;
    std::sort(retOrder.begin(), retOrder.end(), cmp);
  }

  return retOrder;
}

// _____________________________________________________________________________
RtPair Octilinearizer::getRtPair(CombNode* frCmbNd, CombNode* toCmbNd,
                                 const SettledPos& preSettled, BaseGraph* gg,
                                 double maxGrDist) {
  // shortcut
  if (gg->getSettled(frCmbNd) && gg->getSettled(toCmbNd)) {
    return {getCands(frCmbNd, preSettled, gg, 0),
            getCands(toCmbNd, preSettled, gg, 0)};
  }

  std::set<GridNode*> frGrNds, toGrNds;

  size_t i = 0;

  while ((!frGrNds.size() || !toGrNds.size()) && i < 10) {
    auto frCands = getCands(frCmbNd, preSettled, gg, maxGrDist);
    auto toCands = getCands(toCmbNd, preSettled, gg, maxGrDist);

    std::set<GridNode*> isect;
    std::set_intersection(frCands.begin(), frCands.end(), toCands.begin(),
                          toCands.end(), std::inserter(isect, isect.begin()));

    std::set_difference(frCands.begin(), frCands.end(), isect.begin(),
                        isect.end(), std::inserter(frGrNds, frGrNds.begin()));
    std::set_difference(toCands.begin(), toCands.end(), isect.begin(),
                        isect.end(), std::inserter(toGrNds, toGrNds.begin()));

    // this effectively builds a Voronoi diagram
    for (auto iNd : isect) {
      if (util::geo::dist(*iNd->pl().getGeom(), *frCmbNd->pl().getGeom()) <
          util::geo::dist(*iNd->pl().getGeom(), *toCmbNd->pl().getGeom())) {
        frGrNds.insert(iNd);
      } else {
        toGrNds.insert(iNd);
      }
    }

    maxGrDist += i * 2;
    i++;
  }

  return {frGrNds, toGrNds};
}

// _____________________________________________________________________________
std::set<GridNode*> Octilinearizer::getCands(CombNode* cmbNd,
                                             const SettledPos& preSettled,
                                             BaseGraph* gg, size_t maxGrDist) {
  std::set<GridNode*> ret;

  const auto& settled = gg->getSettled(cmbNd);

  if (settled) {
    ret.insert(settled);
  } else if (preSettled.count(cmbNd)) {
    auto nd = preSettled.find(cmbNd)->second->pl().getParent();
    if (nd && !nd->pl().isClosed()) ret.insert(nd);
  } else {
    ret = gg->getGrNdCands(cmbNd, maxGrDist);
  }

  return ret;
}

// _____________________________________________________________________________
BaseGraph* Octilinearizer::newBaseGraph(const DBox& bbox, const CombGraph& cg,
                                        double cellSize, double spacer,
                                        size_t hananIters,
                                        const Penalties& pens) const {
  switch (_baseGraphType) {
    case OCTIGRID:
      return new OctiGridGraph(bbox, cellSize, spacer, pens);
    case CONVEXHULLOCTIGRID:
      return new ConvexHullOctiGridGraph(hull(cg), bbox, cellSize, spacer,
                                         pens);
    case GRID:
      return new GridGraph(bbox, cellSize, spacer, pens);
    case ORTHORADIAL:
      return new OrthoRadialGraph(bbox, cellSize, spacer, pens);
    case PSEUDOORTHORADIAL:
      return new PseudoOrthoRadialGraph(bbox, cellSize, spacer, pens);
    case OCTIHANANGRID:
      return new OctiHananGraph(bbox, cg, cellSize, spacer, hananIters, pens);
    case OCTIQUADTREE:
      return new OctiQuadTree(bbox, cg, cellSize, spacer, pens);
    case HEXGRID:
      return new HexGridGraph(bbox, cellSize, spacer, pens);
    default:
      return 0;
  }
}

// _____________________________________________________________________________
DPolygon Octilinearizer::hull(const CombGraph& cg) const {
  MultiPoint<double> points;
  for (auto nd : cg.getNds()) {
    points.push_back(*nd->pl().getGeom());
  }
  return util::geo::convexHull(points);
}

// _____________________________________________________________________________
void Octilinearizer::statLine(Undrawable status, const std::string& msg,
                              const Drawing& drawing, double ms,
                              const std::string& mark) const {
  switch (status) {
    case DRAWN:
      LOGTO(DEBUG, std::cerr) << " ++ " << msg << ", score " << drawing.score()
                              << ", (" << ms << " ms)" << mark;
      break;
    case NO_PATH:
      LOGTO(DEBUG, std::cerr) << " ++ " << msg << " score <inf>"
                              << " <no path>"
                              << " (" << ms << " ms)" << mark;
      break;
    case NO_CANDS:
      LOGTO(DEBUG, std::cerr) << " ++ " << msg << ", score <inf>"
                              << " <no cands>"
                              << " (" << ms << " ms)" << mark;
      break;
  }
}

// _____________________________________________________________________________
size_t Octilinearizer::maxNodeDeg() const {
  // TODO: this is currently at two locations, in the base graph class and here,
  // fix this
  switch (_baseGraphType) {
    case OCTIGRID:
      return 8;
    case CONVEXHULLOCTIGRID:
      return 8;
    case GRID:
      return 4;
    case ORTHORADIAL:
      return 4;
    case PSEUDOORTHORADIAL:
      return 4;
    case OCTIHANANGRID:
      return 8;
    case OCTIQUADTREE:
      return 8;
    case HEXGRID:
      return 6;
    default:
      return 8;
  }
}
