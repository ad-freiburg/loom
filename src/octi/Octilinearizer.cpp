// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fstream>
#include "octi/Octilinearizer.h"
#include "octi/gridgraph/NodeCost.h"
#include "util/Misc.h"
#include "util/geo/BezierCurve.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Dijkstra.h"
#include "util/log/Log.h"

#include "ilp/ILPGridOptimizer.h"

using namespace octi;
using namespace gridgraph;

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

  std::random_shuffle(order.begin(), order.end() );
  double curCost = draw(order, gg);
  auto best = order;

  for (size_t iter = 0; iter < 10; iter++) {
    std::cerr << "+++ ITERATION " << iter << " +++ " << std::endl;
    double bestBef = curCost;
    for (size_t i = 1; i < order.size(); i++) {
      std::random_shuffle(order.begin(), order.end() );

      gg->reset();
      // util::geo::output::GeoGraphJsonOutput out;
      // std::ofstream of;
      // of.open("octi.json");
      // out.print(*gg, of);
      // of << std::flush;
      // exit(0);

      double nextCost = draw(order, gg);
      std::cerr << "Prev: " << curCost << " This: " << nextCost << std::endl;
      if (nextCost < curCost) {
        curCost = nextCost;
        best = order;
      }
    }

    if (fabs(bestBef - curCost) < 1) {
      std::cerr << "Aborting..." << std::endl;
      break;
    }
  }

  gg->reset();
  curCost = draw(best, gg);

  std::cerr << "Cost: " << curCost << std::endl;

  TransitGraph ret;
  cg.getTransitGraph(&ret);

  *retGg = gg;

  return ret;
}

// _____________________________________________________________________________
void Octilinearizer::settleRes(GridNode* frGrNd, GridNode* toGrNd,
                               GridGraph* gg, CombNode* from, CombNode* to,
                               const GrEdgList& res, CombEdge* e, size_t gen) {
  gg->settleGridNode(toGrNd, to);
  gg->settleGridNode(frGrNd, from);

  // write everything to the result graph
  PolyLine<double> pl = buildPolylineFromRes(res);
  if (e->getFrom() != from) pl.reverse();

  addResidentEdges(gg, e, res);

  for (auto f : res) {
    if (f->pl().isSecondary()) continue;
    gg->balanceEdge(f->getFrom()->pl().getParent(),
                    f->getTo()->pl().getParent());
  }

  e->pl().setPolyLine(pl);
  e->pl().setGeneration(gen);
}

// _____________________________________________________________________________
PolyLine<double> Octilinearizer::buildPolylineFromRes(
    const std::vector<GridEdge*>& res) {
  PolyLine<double> pl;
  for (auto revIt = res.rbegin(); revIt != res.rend(); revIt++) {
    auto f = *revIt;
    if (!f->pl().isSecondary()) {
      if (pl.getLine().size() > 0 &&
          dist(pl.getLine().back(), *f->getFrom()->pl().getGeom()) > 0) {
        BezierCurve<double> bc(pl.getLine().back(),
                               *f->getFrom()->pl().getParent()->pl().getGeom(),
                               *f->getFrom()->pl().getParent()->pl().getGeom(),
                               *f->getFrom()->pl().getGeom());

        for (auto p : bc.render(10).getLine()) pl << p;
      } else {
        pl << *f->getFrom()->pl().getParent()->pl().getGeom();
      }

      pl << *f->getFrom()->pl().getGeom();
      pl << *f->getTo()->pl().getGeom();
    }
  }

  if (res.size()) pl << *res.front()->getTo()->pl().getParent()->pl().getGeom();

  return pl;
}

// _____________________________________________________________________________
double Octilinearizer::getCostFromRes(const std::vector<GridEdge*>& res) {
  size_t c = 0;
  double cost = 0;
  for (auto f : res) {
    cost += c > 0 ? f->pl().cost() : 0;
    c++;
  }
  return cost;
}

// _____________________________________________________________________________
void Octilinearizer::addResidentEdges(GridGraph* g, CombEdge* ce,
                                      const std::vector<GridEdge*>& res) {
  for (auto e : res) {
    assert(e->pl().getResEdges().size() == 0);
    e->pl().addResidentEdge(ce);

    auto f = g->getEdg(e->getTo(), e->getFrom());
    assert(f->pl().getResEdges().size() == 0);
    f->pl().addResidentEdge(ce);
  }
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
double Octilinearizer::draw(const std::vector<CombEdge*>& order, GridGraph* gg) {
  double cost = 0;
  size_t gen = 0;

  for (auto cmbEdg : order) {
    // if (gen > 16) break;

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

    GridNode* frGrNd = gg->getGridNodeFrom(frCmbNd, gg->getCellSize() * 1.7, 0);
    // get surrounding displacement nodes
    double maxDis = getMaxDis(toCmbNd, cmbEdg,  gg->getCellSize());
    std::set<GridNode*> toGrNds = gg->getGridNodesTo(toCmbNd, maxDis, frGrNd);

    // TODO: abort criteria
    while (!toGrNds.size()) {
      maxDis *= 2;
      toGrNds = gg->getGridNodesTo(toCmbNd, maxDis, frGrNd);
    }

    if (!frGrNd) {
      LOG(ERROR) << "Could not sort in source node " << frCmbNd << " ("
                 << frCmbNd->pl().getParent()->pl().getStops().front().name
                 << ")";
      cost = std::numeric_limits<double>::infinity();
      break;
    }

    if (toGrNds.size() == 0) {
      LOG(ERROR) << "Could not sort in target node " << toCmbNd << " ("
                 << toCmbNd->pl().getParent()->pl().getStops().front().name
                 << ") with displacement distance " << maxDis;
      cost = std::numeric_limits<double>::infinity();
      break;
    }

    // why not distance based? (TODO)
    double movePenPerGrid = 5;

    // open the target nodes
    for (auto n : toGrNds) {
      double gridDist = floor(
          dist(*n->pl().getGeom(), *toCmbNd->pl().getGeom()) / gg->getCellSize());

      gg->openNodeSink(n, gridDist * movePenPerGrid);
      gg->openNode(n);
    }

    // open from source node
    gg->openNodeSink(frGrNd, 0);
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

    cost += getCostFromRes(res);

    if (nList.size()) toGrNd = nList.front();

    gg->removeCostVector(frGrNd, addCFromInv);
    gg->removeCostVector(*toGrNds.begin(), addCToInv);

    if (toGrNd == 0) {
      LOG(ERROR) << "Could not route to target node " << toCmbNd << " ("
                 << toCmbNd->pl().getParent()->pl().getStops().front().name
                 << ")";
      cost = std::numeric_limits<double>::infinity();
      break;
    }

    // close the target node
    for (auto n : toGrNds) gg->closeNodeSink(n);
    gg->closeNode(toGrNd);

    // close the start node
    gg->closeNodeSink(frGrNd);
    gg->closeNode(frGrNd);

    settleRes(frGrNd, toGrNd, gg, frCmbNd, toCmbNd, res, cmbEdg, gen);

    gen++;
  }

  return cost;
}
