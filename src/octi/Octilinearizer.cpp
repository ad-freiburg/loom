// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <thread>
#include "octi/Octilinearizer.h"
#include "octi/combgraph/Drawing.h"
#include "octi/gridgraph/NodeCost.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Dijkstra.h"
#include "util/log/Log.h"

#include "ilp/ILPGridOptimizer.h"

using namespace octi;
using namespace gridgraph;

using combgraph::EdgeOrdering;
using octi::combgraph::Drawing;
using util::geo::dist;
using util::geo::DPoint;
using util::geo::len;
using util::graph::Dijkstra;

// _____________________________________________________________________________
void Octilinearizer::removeEdgesShorterThan(LineGraph* g, double d) {
start:
  for (auto n1 : *g->getNds()) {
    for (auto e1 : n1->getAdjList()) {
      if (!e1->pl().dontContract() && e1->pl().getPolyline().getLength() < d) {
        if (e1->getOtherNd(n1)->getAdjList().size() > 1 &&
            n1->getAdjList().size() > 1 &&
            (n1->pl().stops().size() == 0 ||
             e1->getOtherNd(n1)->pl().stops().size() == 0)) {
          auto otherP = e1->getOtherNd(n1)->pl().getGeom();
          auto newGeom =
              DPoint((n1->pl().getGeom()->getX() + otherP->getX()) / 2,
                     (n1->pl().getGeom()->getY() + otherP->getY()) / 2);
          LineNode* n = 0;

          if (e1->getTo()->pl().stops().size() > 0) {
            auto servedLines = g->servedLines(e1->getTo());
            n = g->mergeNds(e1->getFrom(), e1->getTo());

            for (auto l : g->servedLines(n)) {
              if (!servedLines.count(l)) {
                n->pl().addLineNotServed(l);
              }
            }

          } else if (e1->getFrom()->pl().stops().size() > 0) {
            auto servedLines = g->servedLines(e1->getFrom());
            n = g->mergeNds(e1->getTo(), e1->getFrom());

            for (auto l : g->servedLines(n)) {
              if (!servedLines.count(l)) {
                n->pl().addLineNotServed(l);
              }
            }

          } else {
            n = g->mergeNds(e1->getTo(), e1->getFrom());
          }

          n->pl().setGeom(newGeom);
          goto start;
        }
      }
    }
  }
}

// _____________________________________________________________________________
double Octilinearizer::drawILP(LineGraph* tg, LineGraph* outTg,
                               GridGraph** retGg, const Penalties& pens,
                               double gridSize, double borderRad, bool deg2heur,
                               double maxGrDist, bool noSolve, double enfGeoPen,
                               const std::string& path) {
  GridGraph* gg;
  removeEdgesShorterThan(tg, gridSize / 2);

  CombGraph cg(tg, deg2heur);
  auto box = tg->getBBox();
  box = util::geo::pad(box, gridSize + 1);

  std::cerr << "Presolving..." << std::endl;
  try {
    // presolve using heuristical approach to get a first feasible solution
    LineGraph tmpOutTg;
    // important: always use restrLocSearch here!
    draw(cg, box, &tmpOutTg, &gg, pens, gridSize, borderRad, deg2heur, maxGrDist,
         true, enfGeoPen, {});
    std::cerr << "Presolving finished." << std::endl;
  } catch (const NoEmbeddingFoundExc& exc) {
    std::cerr << "Presolve was not sucessful." << std::endl;
    gg = new GridGraph(box, gridSize, borderRad, pens);
  }
  Drawing drawing(gg);
  GeoPensMap enfGeoPens;
  const GeoPensMap* geoPens = 0;

  if (enfGeoPen) {
    std::cerr << "Writing geopens... ";
    auto initOrder = getOrdering(cg, false);
    T_START(geopens);
    for (auto cmbEdg : initOrder) {
      gg->writeGeoCoursePens(cmbEdg, &enfGeoPens, enfGeoPen);
    }
    std::cerr << " done (" << T_STOP(geopens) << "ms)" << std::endl;
    geoPens = &enfGeoPens;
  }

  // TODO
  // if (obstacles.size()) {
    // std::cerr << "Writing obstacles... ";
    // T_START(obstacles);
    // for (auto gg : ggs) for (const auto& obst : obstacles) gg->addObstacle(obst);
    // std::cerr << " done (" << T_STOP(obstacles) << "ms)" << std::endl;
  // }

  ilp::ILPGridOptimizer ilpoptim;

  double score = ilpoptim.optimize(gg, cg, &drawing, maxGrDist, noSolve, geoPens, path);

  drawing.getLineGraph(outTg);

  *retGg = gg;

  return score;
}

// _____________________________________________________________________________
double Octilinearizer::draw(
    LineGraph* tg, LineGraph* outTg, GridGraph** retGg, const Penalties& pens,
    double gridSize, double borderRad, bool deg2heur, double maxGrDist,
    bool restrLocSearch, double enfGeoPen,
    const std::vector<util::geo::Polygon<double>>& obstacles) {
  removeEdgesShorterThan(tg, gridSize / 2);

  // util::geo::output::GeoGraphJsonOutput out;
  // out.print(*tg, std::cout);
  // exit(1);

  CombGraph cg(tg, deg2heur);

  auto box = tg->getBBox();
  box = util::geo::pad(box, gridSize + 1);

  return draw(cg, box, outTg, retGg, pens, gridSize, borderRad, deg2heur,
              maxGrDist, restrLocSearch, enfGeoPen, obstacles);
}
// _____________________________________________________________________________
double Octilinearizer::draw(
    const CombGraph& cg, const util::geo::DBox& box, LineGraph* outTg,
    GridGraph** retGg, const Penalties& pens, double gridSize, double borderRad,
    bool deg2heur, double maxGrDist, bool restrLocSearch, double enfGeoPen,
    const std::vector<util::geo::Polygon<double>>& obstacles) {
  size_t jobs = 4;
  std::vector<GridGraph*> ggs(jobs);

  std::cerr << "Creating grid graphs... ";
  T_START(ggraph);
#pragma omp parallel for
  for (size_t i = 0; i < jobs; i++) {
    ggs[i] = new GridGraph(box, gridSize, borderRad, pens);
  }
  std::cerr << " done (" << T_STOP(ggraph) << "ms)" << std::endl;

  bool found = false;

  size_t tries = 100;
  size_t ITERS = 100;

  GeoPensMap enfGeoPens;
  const GeoPensMap* geoPens = 0;

  auto initOrder = getOrdering(cg, false);

  if (enfGeoPen > 0) {
    std::cerr << "Writing geopens... ";
    T_START(geopens);
    for (auto cmbEdg : initOrder) {
      ggs[0]->writeGeoCoursePens(cmbEdg, &enfGeoPens, enfGeoPen);
    }
    std::cerr << " done (" << T_STOP(geopens) << "ms)" << std::endl;
    geoPens = &enfGeoPens;
  }

  if (obstacles.size()) {
    std::cerr << "Writing obstacles... ";
    T_START(obstacles);
    for (auto gg : ggs) for (const auto& obst : obstacles) gg->addObstacle(obst);
    std::cerr << " done (" << T_STOP(obstacles) << "ms)" << std::endl;
  }

  Drawing drawing(ggs[0]);

  for (size_t i = 0; i < tries; i++) {
    T_START(draw);
    std::vector<CombEdge*> iterOrder;
    if (i != 0)
      iterOrder = getOrdering(cg, true);
    else
      iterOrder = initOrder;

    bool locFound =
        draw(iterOrder, ggs[0], &drawing, drawing.score(), maxGrDist, geoPens);

    if (locFound) {
      std::cerr << " ++ Try " << i << ", score " << drawing.score()
                << ", (took " << T_STOP(draw) << " ms)" << std::endl;
      found = true;
    } else {
      std::cerr << " ++ Try " << i << ", score <inf>"
                << ", next <not found>"
                << " (took " << T_STOP(draw) << " ms)" << std::endl;
    }

    drawing.eraseFromGrid(ggs[0]);
    if (found)
      break;
    else
      drawing.crumble();
  }

  if (!found) throw NoEmbeddingFoundExc();

  std::cerr << "Done.." << std::endl;

  for (size_t i = 0; i < jobs; i++) drawing.applyToGrid(ggs[i]);

  size_t iters = 0;

  std::cerr << "Iterating..." << std::endl;

  std::vector<std::vector<CombNode*>> batches(jobs);
  size_t c = 0;
  for (auto nd : cg.getNds()) {
    if (nd->getDeg() == 0) continue;
    batches[c % jobs].push_back(nd);
    c++;
  }

  for (; iters < ITERS; iters++) {
    T_START(iter);
    std::vector<Drawing> bestFrIters(jobs);

#pragma omp parallel for
    for (size_t btch = 0; btch < jobs; btch++) {
      for (auto a : batches[btch]) {
        Drawing drawingCp = drawing;

        // use the batches grid graph
        drawingCp.setGridGraph(ggs[btch]);

        size_t origX = drawing.getGrNd(a)->pl().getX();
        size_t origY = drawingCp.getGrNd(a)->pl().getY();

        // reverting a
        std::vector<CombEdge*> test;
        for (auto ce : a->getAdjList()) {
          test.push_back(ce);

          drawingCp.eraseFromGrid(ce, ggs[btch]);
          drawingCp.erase(ce);
        }

        drawingCp.erase(a);
        ggs[btch]->unSettleNd(a);

        for (size_t pos = 0; pos < 9; pos++) {
          SettledPos p;

          if (pos == 1) p[a] = {origX + 1, origY + 1};
          if (pos == 2) p[a] = {origX + 1, origY};
          if (pos == 3) p[a] = {origX + 1, origY - 1};

          if (pos == 7) p[a] = {origX - 1, origY + 1};
          if (pos == 6) p[a] = {origX - 1, origY};
          if (pos == 5) p[a] = {origX - 1, origY - 1};

          if (pos == 0) p[a] = {origX, origY + 1};
          if (pos == 8) p[a] = {origX, origY};
          if (pos == 4) p[a] = {origX, origY - 1};

          if (restrLocSearch && ggs[btch]->getNode(p[a].first, p[a].second)) {
            // dont try positions outside the move radius for consistency with
            // ILP approach
            double gridD = dist(
                *a->pl().getGeom(),
                *ggs[btch]->getNode(p[a].first, p[a].second)->pl().getGeom());
            double maxDis = ggs[btch]->getCellSize() * maxGrDist;
            if (gridD >= maxDis) continue;
          }

          Drawing run = drawingCp;

          // we can use bestFromIter.score() as the limit for the shortest
          // path computation, as we can already do at least as good.
          bool found = draw(test, p, ggs[btch], &run, bestFrIters[btch].score(),
                            maxGrDist, geoPens);

          if (found && bestFrIters[btch].score() > run.score()) {
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
    std::cerr << " ++ Iter " << iters << ", prev " << drawing.score()
              << ", next " << bestFrIters[bestCore].score() << " ("
              << (imp >= 0 ? "+" : "") << imp << ", took " << T_STOP(iter)
              << " ms)" << std::endl;

    for (size_t i = 0; i < jobs; i++) {
      drawing.eraseFromGrid(ggs[i]);
      bestFrIters[bestCore].applyToGrid(ggs[i]);
    }
    drawing = bestFrIters[bestCore];

    if (imp < 0.05) break;
  }

  drawing.getLineGraph(outTg);
  auto fullScore = drawing.fullScore();
  std::cerr << "Hop costs: " << fullScore.hop
            << ", bend costs: " << fullScore.bend
            << ", mv costs: " << fullScore.move
            << ", dense costs: " << fullScore.dense << std::endl;

  *retGg = ggs[0];
  return drawing.score();
}

// _____________________________________________________________________________
void Octilinearizer::settleRes(GridNode* frGrNd, GridNode* toGrNd,
                               GridGraph* gg, CombNode* from, CombNode* to,
                               const GrEdgList& res, CombEdge* e) {
  gg->settleNd(toGrNd, to);
  gg->settleNd(frGrNd, from);

  // balance edges
  for (auto f : res) {
    if (f->pl().isSecondary()) continue;
    gg->settleEdg(f->getFrom()->pl().getParent(), f->getTo()->pl().getParent(),
                  e);
  }
}

// _____________________________________________________________________________
void Octilinearizer::writeNdCosts(GridNode* n, CombNode* origNode, CombEdge* e,
                                  GridGraph* g) {
  NodeCost c;
  c += g->topoBlockPenalty(n, origNode, e);
  c += g->spacingPenalty(n, origNode, e);
  c += g->nodeBendPenalty(n, e);

  g->addCostVector(n, c);
}

// _____________________________________________________________________________
bool Octilinearizer::draw(const std::vector<CombEdge*>& order, GridGraph* gg,
                          Drawing* drawing, double cutoff, double maxGrDist,
                          const GeoPensMap* geoPensMap) {
  SettledPos emptyPos;
  return draw(order, emptyPos, gg, drawing, cutoff, maxGrDist, geoPensMap);
}

// _____________________________________________________________________________
bool Octilinearizer::draw(const std::vector<CombEdge*>& ord,
                          const SettledPos& settled, GridGraph* gg,
                          Drawing* drawing, double globCutoff, double maxGrDist,
                          const GeoPensMap* geoPensMap) {
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

    if (frGrNds.size() == 0 || toGrNds.size() == 0) {
      return false;
    }

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
        gg->openNodeSinkFr(n, 0);
      } else {
        costOffsetFrom = (gg->getPenalties().p_45 - gg->getPenalties().p_135);
        gg->openNodeSinkFr(n, costOffsetFrom + gg->ndMovePen(frCmbNd, n));
      }
    }

    // open the target nodes
    for (auto n : toGrNds) {
      if (gg->isSettled(toCmbNd)) {
        // only count displacement penalty ONCE
        gg->openNodeSinkTo(n, 0);
      } else {
        costOffsetTo = (gg->getPenalties().p_45 - gg->getPenalties().p_135);
        gg->openNodeSinkTo(n, costOffsetTo + gg->ndMovePen(toCmbNd, n));
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

    auto heur = GridHeur(gg, toGrNds);

    if (geoPensMap) {
      // init cost function with geo distance penalties
      auto cost = GridCostGeoPen(cutoff + costOffsetTo + costOffsetFrom,
                                 &geoPensMap->find(cmbEdg)->second);
      Dijkstra::shortestPath(frGrNds, toGrNds, cost, heur, &eL, &nL);
    } else {
      auto cost = GridCost(cutoff + costOffsetTo + costOffsetFrom);
      Dijkstra::shortestPath(frGrNds, toGrNds, cost, heur, &eL, &nL);
    }

    if (!nL.size()) {
      // cleanup
      for (auto n : toGrNds) gg->closeNodeSinkTo(n);
      for (auto n : frGrNds) gg->closeNodeSinkFr(n);

      return false;
    }

    toGrNd = nL.front();
    frGrNd = nL.back();

    // remove the cost offsets to not distort final costs
    eL.front()->pl().setCost(eL.front()->pl().cost() - costOffsetTo);
    eL.back()->pl().setCost(eL.back()->pl().cost() - costOffsetFrom);

    // draw
    drawing->draw(cmbEdg, eL, rev);

    // close the source and target node
    for (auto n : toGrNds) gg->closeNodeSinkTo(n);
    for (auto n : frGrNds) gg->closeNodeSinkFr(n);

    settleRes(frGrNd, toGrNd, gg, frCmbNd, toCmbNd, eL, cmbEdg);
  }

  return true;
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
                                 const SettledPos& preSettled, GridGraph* gg,
                                 double maxGrDist) {
  // shortcut
  if (gg->getSettled(frCmbNd) && gg->getSettled(toCmbNd)) {
    return {getCands(frCmbNd, preSettled, gg, 0),
            getCands(toCmbNd, preSettled, gg, 0)};
  }

  double maxDis = gg->getCellSize() * maxGrDist;

  std::set<GridNode*> frGrNds;
  std::set<GridNode*> toGrNds;

  size_t i = 0;

  while ((!frGrNds.size() || !toGrNds.size()) && i < 10) {
    std::set<GridNode*> frCands = getCands(frCmbNd, preSettled, gg, maxDis);
    std::set<GridNode*> toCands = getCands(toCmbNd, preSettled, gg, maxDis);

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

    maxDis += i * 2.0;
    i++;
  }

  return {frGrNds, toGrNds};
}

// _____________________________________________________________________________
std::set<GridNode*> Octilinearizer::getCands(CombNode* cmbNd,
                                             const SettledPos& preSettled,
                                             GridGraph* gg, double maxDis) {
  std::set<GridNode*> ret;

  if (gg->getSettled(cmbNd)) {
    ret.insert(gg->getSettled(cmbNd));
  } else if (preSettled.count(cmbNd)) {
    auto nd = gg->getNode(preSettled.find(cmbNd)->second.first,
                          preSettled.find(cmbNd)->second.second);
    if (nd && !nd->pl().isClosed()) ret.insert(nd);
  } else {
    ret = gg->getGrNdCands(cmbNd, maxDis);
  }

  return ret;
}
