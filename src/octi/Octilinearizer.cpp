// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <thread>
#include "octi/Octilinearizer.h"
#include "octi/basegraph/BaseGraph.h"
#include "octi/basegraph/ConvexHullOctiGridGraph.h"
#include "octi/basegraph/GridGraph.h"
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

#include "ilp/ILPGridOptimizer.h"

using namespace octi;
using namespace basegraph;
using namespace util;

using combgraph::EdgeOrdering;
using octi::basegraph::BaseGraph;
using octi::combgraph::Drawing;
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
Score Octilinearizer::drawILP(const CombGraph& cg, const util::geo::DBox& box,
                              LineGraph* outTg, BaseGraph** retGg,
                              const Penalties& pens, double gridSize,
                              double borderRad, double maxGrDist, bool noSolve,
                              double enfGeoPen, int timeLim,
                              const std::string& solverStr,
                              const std::string& path) {
  BaseGraph* gg;

  LOGTO(DEBUG, std::cerr) << "Presolving...";
  try {
    // presolve using heuristical approach to get a first feasible solution
    LineGraph tmpOutTg;
    // important: always use restrLocSearch here!
    draw(cg, box, &tmpOutTg, &gg, pens, gridSize, borderRad, maxGrDist, true,
         enfGeoPen, {}, std::numeric_limits<size_t>::max());
    LOGTO(DEBUG, std::cerr) << "Presolving finished.";
  } catch (const NoEmbeddingFoundExc& exc) {
    LOGTO(DEBUG, std::cerr) << "Presolve was not sucessful.";
    gg = newBaseGraph(box, cg, gridSize, borderRad, pens);
    gg->init();
  }

  Drawing drawing(gg);
  GeoPensMap enfGeoPens;
  const GeoPensMap* geoPens = 0;

  if (enfGeoPen) {
    LOGTO(DEBUG, std::cerr) << "Writing geopens... ";
    auto initOrder = getOrdering(cg, false);
    T_START(geopens);
    for (auto cmbEdg : initOrder) {
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

  ilpoptim.optimize(gg, cg, &drawing, maxGrDist, noSolve, geoPens, timeLim,
                    solverStr, path);

  drawing.getLineGraph(outTg);

  *retGg = gg;

  return drawing.fullScore();
}

// _____________________________________________________________________________
Score Octilinearizer::draw(const CombGraph& cg, const DBox& box,
                           LineGraph* outTg, BaseGraph** retGg,
                           const Penalties& pens, double gridSize,
                           double borderRad, double maxGrDist,
                           bool restrLocSearch, double enfGeoPen,
                           const std::vector<Polygon<double>>& obstacles,
                           size_t abortAfter) {
  size_t jobs = 4;
  std::vector<BaseGraph*> ggs(jobs);

  auto gg = newBaseGraph(box, cg, gridSize, borderRad, pens);
  gg->init();

  // util::geo::output::GeoGraphJsonOutput out;
  // out.print(*gg, std::cout);
  // exit(0);

  LOGTO(DEBUG, std::cerr) << "Creating grid graphs... ";
  T_START(ggraph);
#pragma omp parallel for
  for (size_t i = 0; i < jobs; i++) {
    ggs[i] = newBaseGraph(box, cg, gridSize, borderRad, pens);
    ggs[i]->init();
  }
  LOGTO(DEBUG, std::cerr) << "Done. (" << T_STOP(ggraph) << "ms)";

  size_t INITIAL_TRIES = 100;
  size_t LOCAL_SEARCH_ITERS = 100;
  double CONVERGENCE_THRESHOLD = 0.05;

  GeoPensMap enfGeoPens;
  const GeoPensMap* geoPens = 0;

  // get a non-randomized initial ordering
  auto initOrder = getOrdering(cg, false);

  if (enfGeoPen > 0) {
    LOGTO(DEBUG, std::cerr) << "Writing geopens... ";
    T_START(geopens);
    for (auto cmbEdg : initOrder) {
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
  T_START(draw);
  auto status = draw(initOrder, ggs[0], &drawing, drawing.score(), maxGrDist,
                     geoPens, abortAfter);

  statLine(status, "Try <init>", drawing, T_STOP(draw), "");
  if (status != DRAWN) drawing.crumble();

  size_t a = 1;

  // TODO: better check topology violations
  if (drawing.score() == INF || drawing.violations()) {
    bool abort = false;
    // clear the initial grid
    drawing.eraseFromGrid(ggs[0]);
    // set a to 0 to re-add the best solution later to grid 0
    a = 0;

#pragma omp parallel for
    for (size_t btch = 0; btch < jobs; btch++) {
      for (size_t i = 0; i < INITIAL_TRIES / jobs; i++) {
        size_t tryNr = INITIAL_TRIES / jobs * btch + i;
        if (abort) continue;
        T_START(draw);
        Drawing drawingCp(ggs[btch]);

        // get a randomized ordering
        std::vector<CombEdge*> iterOrder = getOrdering(cg, true);

        auto status = draw(iterOrder, ggs[btch], &drawingCp, drawingCp.score(),
                           maxGrDist, geoPens, abortAfter);

        drawingCp.eraseFromGrid(ggs[btch]);

#pragma omp critical
        {
          if (status == DRAWN && drawingCp.violations() == 0) {
            // found solution without topology violations, immediately abort
            abort = true;
            statLine(status, std::string("Try ") + std::to_string(tryNr),
                     drawingCp, T_STOP(draw), "**");
          }
          if (status == DRAWN && drawingCp.score() < drawing.score()) {
            drawing = drawingCp;
            statLine(status, std::string("Try ") + std::to_string(tryNr),
                     drawingCp, T_STOP(draw), "*");
          } else {
            drawingCp.crumble();
            statLine(status, std::string("Try ") + std::to_string(tryNr),
                     drawingCp, T_STOP(draw), "");
          }
        }
      }
    }
  }

  if (drawing.score() == INF) throw NoEmbeddingFoundExc();

  LOGTO(DEBUG, std::cerr) << "Done.";

  // skipping grid 0 here if a remains set to 1
  for (size_t i = a; i < jobs; i++) drawing.applyToGrid(ggs[i]);

  size_t iters = 0;

  LOGTO(DEBUG, std::cerr) << "Initial score: " << drawing.score() << " ("
                          << drawing.violations() << " topology violations).";
  LOGTO(DEBUG, std::cerr) << "Starting local search...";

  // dont use local search if abortAfter is set
  if (abortAfter != std::numeric_limits<size_t>::max()) LOCAL_SEARCH_ITERS = 0;

  std::vector<std::vector<CombNode*>> batches(jobs);
  size_t c = 0;
  for (auto nd : cg.getNds()) {
    if (nd->getDeg() == 0) continue;
    batches[c % jobs].push_back(nd);
    c++;
  }

  for (; iters < LOCAL_SEARCH_ITERS; iters++) {
    T_START(iter);
    std::vector<Drawing> bestFrIters(jobs);

#pragma omp parallel for
    for (size_t btch = 0; btch < jobs; btch++) {
      for (auto a : batches[btch]) {
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

      // auto c = Dijkstra::shortestPath(frGrNds, toGrNds, cost, *heur, &eL,
      // &nL);
      Dijkstra::shortestPath(frGrNds, toGrNds, cost, *heur, &eL, &nL);

      // for (auto e : eL) std::cerr << "E " << e << " " <<
      // e->getFrom()->pl().getX() << "," << e->getFrom()->pl().getY() << " -> "
      // << e->getTo()->pl().getX() << "," << e->getTo()->pl().getY() << " " <<
      // e->pl().cost() << std::endl; eL.clear(); nL.clear();

      // eL.clear();
      // nL.clear();

      // auto c2 = Dijkstra::shortestPath(frGrNds, toGrNds, cost, &eL, &nL);
      // std::cerr << c << " vs " << c2 << " " << fabs(c - c2) << std::endl;
      // assert((std::isinf(c) && std::isinf(c2)) || (std::isnan(c) &&
      // std::isnan(c2)) || fabs(c - c2) < 0.001);

      // for (auto e : eL) std::cerr << "E " << e << " " <<
      // e->getFrom()->pl().getX() << "," << e->getFrom()->pl().getY() << " -> "
      // << e->getTo()->pl().getX() << "," << e->getTo()->pl().getY() << " " <<
      // e->pl().cost() << std::endl; assert((isinf(c) && isinf(c2)) || (fabs(c
      // - c2) < 0.0001)); exit(0);
    }

    delete heur;

    // if (i == 4) {
    // util::geo::output::GeoGraphJsonOutput out;
    // out.print(*gg, std::cout);
    // exit(0);
    // }

    if (!nL.size()) {
      // std::cerr << "FAILED TO FIND A ROUTE FROM " << std::endl;
      // for (auto fr : frGrNds)
      // std::cerr << fr->pl().getX() << "," << fr->pl().getY() << std::endl;
      // std::cerr << " TO" << std::endl;
      // for (auto to : toGrNds)
      // std::cerr << to->pl().getX() << "," << to->pl().getY() << std::endl;
      // std::cerr << "DEG FROM CMB NODE " << frCmbNd->getDeg() << std::endl;
      // std::cerr << "DEG TO CMB NODE " << toCmbNd->getDeg() << std::endl;
      // util::geo::output::GeoGraphJsonOutput out;
      // out.print(*gg, std::cout);
      // exit(0);

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
std::vector<CombEdge*> Octilinearizer::getOrdering(const CombGraph& cg,
                                                   bool randr) const {
  NodePQ globalPq, dangling;

  std::set<CombNode*> settled;
  std::vector<CombEdge*> order;

  for (auto n : cg.getNds()) globalPq.push(n);
  std::set<CombEdge*> done;

  while (!globalPq.empty()) {
    auto n = globalPq.top();
    globalPq.pop();
    dangling.push(n);

    while (!dangling.empty()) {
      auto n = dangling.top();
      dangling.pop();

      if (settled.find(n) != settled.end()) continue;

      auto odSet = n->pl().getEdgeOrdering().getOrderedSet();
      if (randr) std::random_shuffle(odSet.begin(), odSet.end());

      for (auto ee : odSet) {
        if (done.find(ee.first) != done.end()) continue;
        done.insert(ee.first);
        dangling.push(ee.first->getOtherNd(n));

        order.push_back(ee.first);
      }
      settled.insert(n);
    }
  }

  return order;
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
      return new OctiHananGraph(bbox, cg, cellSize, spacer, pens);
    case OCTIQUADTREE:
      return new OctiQuadTree(bbox, cg, cellSize, spacer, pens);
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
                              << ", (" << ms << " ms) " << mark;
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
