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

  std::random_shuffle(order.begin(), order.end());
  Drawing curCost = draw(order, gg);
  auto best = order;

  for (size_t iter = 0; iter < 10; iter++) {
    // std::cerr << "+++ ITERATION " << iter << " +++ " << std::endl;
    Drawing bestBef = curCost;
    for (size_t i = 1; i < order.size(); i++) {
      std::random_shuffle(order.begin(), order.end());
      gg->reset();

      Drawing nextCost = draw(order, gg);
      // std::cerr << "Prev: " << curCost << " This: " << nextCost << std::endl;
      if (nextCost.score() < curCost.score()) {
        curCost = nextCost;
        best = order;
      }
    }

    if (fabs(bestBef.score() - curCost.score()) < 1) {
      break;
    }
  }

  // TODO: move nodes...

  gg->reset();
  SettledPos curPos;

  // arbitrary, but fixed ordering of nodes
  // std::vector<CombNode*> orderedNds(cg.getNds()->begin(),
  // cg.getNds()->end());

  curCost = draw(best, gg);

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

  if (curCost.score() == std::numeric_limits<double>::infinity()) {
    LOG(ERROR) << "Could not find planar embedding for input graph.";
    exit(1);
  }

  std::cerr << "Cost: " << curCost.score() << std::endl;

  TransitGraph ret;
  curCost.getTransitGraph(&ret);

  *retGg = gg;

  return ret;
}

// _____________________________________________________________________________
SettledPos Octilinearizer::getNeighbor(const SettledPos& pos,
                                       const std::vector<CombNode*>& ordered,
                                       size_t i) const {
  size_t n = i / 8;
  size_t p = (i % 8);
  auto nd = ordered[n];

  SettledPos ret = pos;

  // std::cerr << "Cur: " << ret[nd].first << ", " << ret[nd].second <<
  // std::endl;

  if (p == 0) {
    ret[nd].second += 1;
  }

  if (p == 1) {
    ret[nd].first += 1;
    ret[nd].second += 1;
  }

  if (p == 2) {
    ret[nd].first += 1;
  }

  if (p == 3) {
    ret[nd].first += 1;
    ret[nd].second -= 1;
  }

  if (p == 4) {
    ret[nd].second -= 1;
  }

  if (p == 5) {
    ret[nd].first -= 1;
    ret[nd].second -= 1;
  }

  if (p == 6) {
    ret[nd].first -= 1;
  }

  if (p == 7) {
    ret[nd].first -= 1;
    ret[nd].second += 1;
  }

  // std::cerr << "Neigh: " << ret[nd].first << ", " << ret[nd].second <<
  // std::endl;

  return ret;
}

// _____________________________________________________________________________
void Octilinearizer::settleRes(GridNode* frGrNd, GridNode* toGrNd,
                               GridGraph* gg, CombNode* from, CombNode* to,
                               const GrEdgList& res, CombEdge* e) {
  gg->settleGridNode(toGrNd, to);
  gg->settleGridNode(frGrNd, from);

  // add resident edges to grid graph
  for (auto ge : res) {
    assert(ge->pl().getResEdges().size() == 0);
    ge->pl().addResidentEdge(e);

    auto gf = gg->getEdg(ge->getTo(), ge->getFrom());
    assert(gf->pl().getResEdges().size() == 0);
    gf->pl().addResidentEdge(e);
  }

  // balance edges
  for (auto f : res) {
    if (f->pl().isSecondary()) continue;
    gg->balanceEdge(f->getFrom()->pl().getParent(),
                    f->getTo()->pl().getParent());
  }

  // write everything to the result graph
  // PolyLine<double> pl = buildPolylineFromRes(res);
  // if (e->getFrom() != from) pl.reverse();
  // e->pl().setPolyLine(pl);
}

// _____________________________________________________________________________
NodeCost Octilinearizer::writeNdCosts(GridNode* n, CombNode* origNode,
                                      CombEdge* e, GridGraph* g) {
  NodeCost c;
  c += g->topoBlockPenalty(n, origNode, e);
  c += g->nodeBendPenalty(n, e);

  return g->addCostVector(n, c);
}

// _____________________________________________________________________________
Drawing Octilinearizer::draw(const std::vector<CombEdge*>& ord, GridGraph* gg) {
  Drawing drawing;
  size_t gen = 0;

  SettledPos retPos;

  for (auto cmbEdg : ord) {
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

    frGrNd = gg->getGridNodeFrom(frCmbNd, gg->getCellSize() * 1.7, 0);

    // get surrounding displacement nodes
    double maxDis = getMaxDis(toCmbNd, cmbEdg, gg->getCellSize());
    toGrNds = gg->getGridNodesTo(toCmbNd, maxDis, frGrNd);

    // TODO: abort criteria
    while (!toGrNds.size()) {
      maxDis *= 2;
      toGrNds = gg->getGridNodesTo(toCmbNd, maxDis, frGrNd);
    }

    if (!frGrNd || toGrNds.size() == 0) {
      drawing = Drawing(std::numeric_limits<double>::infinity());
      break;
    }

    for (auto to : toGrNds) assert(to != frGrNd);

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

      gg->openNode(n);
    }

    // open from source node
    if (gg->isSettled(toCmbNd)) {
      // only count displacement penalty ONCE
      gg->openNodeSink(frGrNd, 0);
    } else {
      double gridD = dist(*frGrNd->pl().getGeom(), *frCmbNd->pl().getGeom());
      gridD = gridD / gg->getCellSize();
      gg->openNodeSink(frGrNd, gridD * penPerGrid);
    }

    gg->openNode(frGrNd);

    NodeCost addCFromInv = writeNdCosts(frGrNd, frCmbNd, cmbEdg, gg);
    NodeCost addCToInv;

    if (toGrNds.size() == 1) {
      addCToInv = writeNdCosts(*toGrNds.begin(), toCmbNd, cmbEdg, gg);
    }

    GrEdgList res;
    GrNdList nList;
    GridNode* toGrNd = 0;
    Dijkstra::shortestPath(frGrNd, toGrNds, GridCost(),
                           GridHeur(gg, frGrNd, toGrNds), &res, &nList);

    if (nList.size()) toGrNd = nList.front();
    if (toGrNd == 0) {
      drawing = Drawing(std::numeric_limits<double>::infinity());
      break;
    }

    // draw
    drawing.draw(cmbEdg, res);

    gg->removeCostVector(frGrNd, addCFromInv);
    gg->removeCostVector(*toGrNds.begin(), addCToInv);

    // close the target node
    for (auto n : toGrNds) gg->closeNodeSink(n);
    gg->closeNode(toGrNd);

    // close the start node
    gg->closeNodeSink(frGrNd);
    gg->closeNode(frGrNd);

    settleRes(frGrNd, toGrNd, gg, frCmbNd, toCmbNd, res, cmbEdg);

    gen++;
  }

  return drawing;
}
