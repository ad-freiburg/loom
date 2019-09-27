// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fstream>
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

using octi::combgraph::Drawing;
using combgraph::EdgeOrdering;
using util::graph::Dijkstra;
using util::graph::Dijkstra;
using util::geo::len;
using util::geo::dist;
using util::geo::DPoint;

// _____________________________________________________________________________
void Octilinearizer::removeEdgesShorterThan(TransitGraph* g, double d) {
start:
  for (auto n1 : *g->getNds()) {
    for (auto e1 : n1->getAdjList()) {
      if (e1->pl().getPolyline().getLength() < d) {
        if (e1->getOtherNd(n1)->getAdjList().size() > 1 &&
            n1->getAdjList().size() > 1 &&
            (n1->pl().getStops().size() == 0 ||
             e1->getOtherNd(n1)->pl().getStops().size() == 0)) {
          auto otherP = e1->getFrom()->pl().getGeom();

          TransitNode* n = 0;

          if (e1->getTo()->pl().getStops().size() > 0) {
            n = g->mergeNds(e1->getFrom(), e1->getTo());
          } else {
            n = g->mergeNds(e1->getTo(), e1->getFrom());
          }

          n->pl().setGeom(
              DPoint((n->pl().getGeom()->getX() + otherP->getX()) / 2,
                     (n->pl().getGeom()->getY() + otherP->getY()) / 2));
          goto start;
        }
      }
    }
  }
}

// _____________________________________________________________________________
double Octilinearizer::getMaxDis(CombNode* to, CombEdge* e, double gridSize) {
  return gridSize * 1.7;
}

// _____________________________________________________________________________
TransitGraph Octilinearizer::draw(TransitGraph* tg, GridGraph** retGg,
                                  const Penalties& pens, double gridSize,
                                  double borderRad) {
  std::cerr << "Removing short edges... ";
  T_START(remshortegs);
  removeEdgesShorterThan(tg, gridSize / 2);
  std::cerr << " done (" << T_STOP(remshortegs) << "ms)" << std::endl;

  std::cerr << "Building combination graph... ";
  T_START(combgraph);
  CombGraph cg(tg);
  std::cerr << " done (" << T_STOP(combgraph) << "ms)" << std::endl;

  auto box = tg->getBBox();

  auto gg = new GridGraph(box, gridSize, borderRad, pens);

  ///////////
  // ilp::ILPGridOptimizer ilpoptim;

  // ilpoptim.optimize(gg, cg);

  // util::geo::output::GeoGraphJsonOutput out;
  // std::ofstream of;
  // of.open("octi.json");
  // out.print(*gg, of);
  // of << std::flush;

  // exit(0);

  /////////////

  T_START(grid);
  std::cerr << "Build grid graph in " << T_STOP(grid) << " ms " << std::endl;

  NodePQ globalPq, dangling;

  std::set<CombNode*> settled;

  std::vector<CombEdge*> order;

  for (auto n : *cg.getNds()) globalPq.push(n);

  std::set<CombEdge*> done;
  int64_t gen = 0;

  while (!globalPq.empty()) {
    auto n = globalPq.top();
    globalPq.pop();
    dangling.push(n);

    while (!dangling.empty()) {
      auto n = dangling.top();
      dangling.pop();

      if (settled.find(n) != settled.end()) continue;

      for (auto ee : n->pl().getEdgeOrdering().getOrderedSet()) {
        if (done.find(ee.first) != done.end()) continue;
        done.insert(ee.first);
        dangling.push(ee.first->getOtherNd(n));

        order.push_back(ee.first);
      }
      settled.insert(n);
    }
  }

  Drawing drawing(gg);
  draw(order, gg, &drawing);

  std::cerr << " +++++++ Cost from graph: " << drawing.score() << std::endl;

  Drawing bestFromIter = drawing;

  for (auto a : *cg.getNds()) {
    Drawing drawingCp = drawing;
    size_t origX = drawing.getGrNd(a)->pl().getX();
    size_t origY = drawing.getGrNd(a)->pl().getX();

    std::cerr << drawingCp.score() << std::endl;

    // reverting a
    std::vector<CombEdge*> test;
    for (auto ce : a->getAdjList()) {
      auto es = drawingCp.getGrEdgs(ce);
      drawingCp.erase(ce);
      test.push_back(ce);

      for (auto e : es) {
        gg->unSettleEdg(e->getFrom()->pl().getParent(),
                        e->getTo()->pl().getParent());
      }
    }

    drawingCp.erase(a);
    gg->unSettleNd(a);

    std::map<GridEdge*, double> oldCosts;
    for (auto gn :  *gg->getNds()) {
      for (auto ge : gn->getAdjListOut()) {
        oldCosts[ge] = ge->pl().cost();
      }
    }

    for (size_t  pos = 0; pos < 8; pos++) {
      std::cerr << "Checking " << a << " at " << pos << std::endl;
      Drawing run = drawingCp;
      SettledPos p;

      if (pos == 0) p[a] = {origX, origY + 1};
      if (pos == 1) p[a] = {origX + 1, origY + 1};
      if (pos == 2) p[a] = {origX + 1,origY };
      if (pos == 3) p[a] = {origX + 1,  origY- 1};
      if (pos == 4) p[a] = {origX,  origY- 1};
      if (pos == 5) p[a] = {origX - 1, origY - 1};
      if (pos == 6) p[a] = {origX - 1, origY};
      if (pos == 7) p[a] = {origX - 1, origY + 1};

      bool found = draw(test, p, gg, &run);

      if (!found) {
        std::cerr << "No solution found..." << std::endl;
      } else {
        std::cerr << " ++++++ Cost from graph after: " << run.score() << std::endl;
        if (bestFromIter.score() > run.score()) {
          std::cerr << "  ----- NEWBEST -----" << std::endl;
          bestFromIter = run;
        }
      }

      // reset grid
      for (auto ce : a->getAdjList()) {
        if (!run.drawn(ce)) continue;
        auto es = run.getGrEdgs(ce);

        std::cerr << "Reverting " << ce << std::endl;

        for (auto e : es) {
          gg->unSettleEdg(e->getFrom()->pl().getParent(),
                          e->getTo()->pl().getParent());
        }
      }

      if (gg->isSettled(a)) gg->unSettleNd(a);
    }

    for (auto gn :  *gg->getNds()) {
      for (auto ge : gn->getAdjListOut()) {
        if (oldCosts[ge] == ge->pl().cost()) continue;
        if (!(fabs(oldCosts[ge] - ge->pl().cost()) < 0.0001)) {
          std::cerr << oldCosts[ge] << " vs " << ge->pl().cost() << std::endl;
          assert(fabs(oldCosts[ge] - ge->pl().cost()) < 0.0001);
        }
      }
    }

    gg->settleNd(const_cast<GridNode*>(drawing.getGrNd(a)), a);

    // RE-SETTLE EDGES!
    for (auto ce : a->getAdjList()) {
      auto es = drawing.getGrEdgs(ce);

      for (auto e : es) {
        if (e->pl().isSecondary()) continue;
        gg->settleEdg(e->getFrom()->pl().getParent(), e->getTo()->pl().getParent(),
                      ce);
      }
    }
  }

  // TODO:

  // drawing.remFromGrid(gg);
  // bestFromIter.applyToGrid(gg);

  std::cerr << " ++++++ Initial: " << drawing.score() << std::endl;
  std::cerr << " ++++++ Cost after one iteration: " << bestFromIter.score() << std::endl;

  util::geo::output::GeoGraphJsonOutput out;
  std::ofstream of;
  of.open("octi.json");
  out.print(*gg, of);
  of << std::flush;
  exit(1);
  auto best = order;

  for (size_t iter = 0; iter < 10; iter++) {
    Drawing bestBef = drawing;
    for (size_t i = 1; i < order.size(); i++) {
      std::random_shuffle(order.begin(), order.end());
      gg->reset();

      Drawing nextDrawing(gg);
      draw(order, gg, &nextDrawing);
      if (nextDrawing.score() < drawing.score()) {
        drawing = nextDrawing;
        best = order;
      }
    }

    if (fabs(bestBef.score() - drawing.score()) < 1) {
      break;
    }
  }

  // TODO: move nodes...

  gg->reset();

  // arbitrary, but fixed ordering of nodes
  // std::vector<CombNode*> orderedNds(cg.getNds()->begin(),
  // cg.getNds()->end());

  draw(best, gg, &drawing);

  // now we have the node positions in curPos

  // double bestNeighCost = curCost;
  // SettledPos bestNeighPos = curPos;
  // for (size_t i = 0; i < cg.getNds()->size() * 8; i++) {
  // break;
  // std::cerr << "Change node " << (i / 8) << " to " << ((i % 8)) << std::endl;

  // SettledPos newPos = getNeighbor(curPos, orderedNds, i);

  // std::cerr << "Redrawing.." << std::endl;
  // double neighCost = draw(best, gg);
  // std::cerr << "@ cost " << neighCost << std::endl;

  // if (neighCost < bestNeighCost) {
  // std::cerr << "MURR" << std::endl;
  // bestNeighPos = newPos;
  // neighCost = bestNeighCost;
  // }
  // }

  // std::cerr << "Best node position neighbor score was " << bestNeighCost
  // << std::endl;

  // curCost = draw(best, &curPos, gg);
  // std::cerr << curCost << std::endl;

  if (drawing.score() == std::numeric_limits<double>::infinity()) {
    LOG(ERROR) << "Could not find planar embedding for input graph.";
    exit(1);
  }

  std::cerr << "Cost: " << drawing.score() << std::endl;

  TransitGraph ret;
  drawing.getTransitGraph(&ret);

  *retGg = gg;

  return ret;
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
void Octilinearizer::writeNdCosts(GridNode* n, CombNode* origNode,
                                      CombEdge* e, GridGraph* g) {
  NodeCost c;
  c += g->topoBlockPenalty(n, origNode, e);
  c += g->nodeBendPenalty(n, e);

  g->addCostVector(n, c);
}

// _____________________________________________________________________________
bool Octilinearizer::draw(const std::vector<CombEdge*>& order, GridGraph* gg,
                          Drawing* drawing) {
  SettledPos emptyPos;
  return draw(order, emptyPos, gg, drawing);
}

// _____________________________________________________________________________
bool Octilinearizer::draw(const std::vector<CombEdge*>& ord,
                          const SettledPos& settled, GridGraph* gg,
                          Drawing* drawing) {
  size_t gen = 0;
  bool fail = false;

  SettledPos retPos;

  for (auto cmbEdg : ord) {
    // if (gen == 5) return true;
    auto frCmbNd = cmbEdg->getFrom();
    auto toCmbNd = cmbEdg->getTo();

    bool reversed = false;

    if (frCmbNd->getDeg() < toCmbNd->getDeg() ||
        (frCmbNd->getDeg() == toCmbNd->getDeg() &&
         frCmbNd->pl().getRouteNumber() < toCmbNd->pl().getRouteNumber()) ||
        (!gg->isSettled(frCmbNd) && gg->isSettled(toCmbNd))) {
      auto tmp = frCmbNd;
      frCmbNd = toCmbNd;
      toCmbNd = tmp;
      reversed = !reversed;
    }

    GridNode* frGrNd;
    std::set<GridNode*> toGrNds;

    // STEP 1
    // grid node selection
    if (settled.count(frCmbNd)) {
      frGrNd = gg->getNode(settled.find(frCmbNd)->second.first,
                           settled.find(frCmbNd)->second.second);
      if (!frGrNd || frGrNd->pl().isClosed()) {
        fail = true;
        break;
      }
    } else {
      frGrNd = gg->getGridNodeFrom(frCmbNd, gg->getCellSize() * 1.7, 0);
    }

    if (settled.count(toCmbNd)) {
      toGrNds.insert(gg->getNode(settled.find(toCmbNd)->second.first,
                                 settled.find(toCmbNd)->second.second));
      if (!(*toGrNds.begin()) || *toGrNds.begin() == frGrNd) {
        fail = true;
        break;
      }
      if ((*toGrNds.begin())->pl().isClosed()) {
        fail = true;
        break;
      }
    } else {
      // get surrounding displacement nodes
      double maxDis = getMaxDis(toCmbNd, cmbEdg, gg->getCellSize());
      toGrNds = gg->getGridNodesTo(toCmbNd, maxDis, frGrNd);

      // TODO: abort criteria
      while (!toGrNds.size()) {
        maxDis *= 2;
        toGrNds = gg->getGridNodesTo(toCmbNd, maxDis, frGrNd);
      }
    }

    if (!frGrNd || toGrNds.size() == 0) {
        fail = true;
        break;
    }

    for (auto to : toGrNds) assert(to != frGrNd);

    // END STEP 1

    // why not distance based? (TODO)
    double penPerGrid = 5;

    // open the target nodes
    for (auto n : toGrNds) {
      double gridD = dist(*n->pl().getGeom(), *toCmbNd->pl().getGeom());
      gridD = gridD / gg->getCellSize();

      if (gg->isSettled(toCmbNd)) {
        // only count displacement penalty ONCE
        gg->openNodeSink(n, 0);
      } else {
        gg->openNodeSink(n, gridD * penPerGrid);
      }
    }

    // open from source node
    if (gg->isSettled(frCmbNd)) {
      // only count displacement penalty ONCE
      gg->openNodeSink(frGrNd, 0);
    } else {
      double gridD = dist(*frGrNd->pl().getGeom(), *frCmbNd->pl().getGeom());
      gridD = gridD / gg->getCellSize();
      gg->openNodeSink(frGrNd, gridD * penPerGrid);
    }

    if (gg->isSettled(frCmbNd)) {
      writeNdCosts(frGrNd, frCmbNd, cmbEdg, gg);
    }

    if (toGrNds.size() == 1 && gg->isSettled(toCmbNd)) {
      // the size() == 1 check is important, because nd cost writing will
      // not work if the to node is not already settled!
      writeNdCosts(*toGrNds.begin(), toCmbNd, cmbEdg, gg);
    }

    GrEdgList res;
    GrNdList nList;
    GridNode* toGrNd = 0;
    Dijkstra::shortestPath(frGrNd, toGrNds, GridCost(),
                           GridHeur(gg, frGrNd, toGrNds), &res, &nList);

    if (nList.size()) toGrNd = nList.front();
    if (toGrNd == 0) {
      // cleanup
      for (auto n : toGrNds) gg->closeNodeSink(n);
      gg->closeNodeSink(frGrNd);

      fail = true;
      break;
    }

    // draw
    drawing->draw(cmbEdg, res, reversed);

    // close the target node
    for (auto n : toGrNds) gg->closeNodeSink(n);

    // close the start node
    gg->closeNodeSink(frGrNd);

    settleRes(frGrNd, toGrNd, gg, frCmbNd, toCmbNd, res, cmbEdg);

    gen++;
  }

  if (fail) {
    std::cerr << "FAIL!" << std::endl;
    // undraw...

    return false;
  }

  return true;
}
