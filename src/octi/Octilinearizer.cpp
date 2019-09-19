// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/Octilinearizer.h"
#include "octi/gridgraph/NodeCost.h"
#include "util/Misc.h"
#include "util/geo/BezierCurve.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Dijkstra.h"
#include "util/log/Log.h"

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
  // this is just a collection of displacement heuristics, should be
  // made configurable

  if (to->getAdjList().size() == 1) {
    return len(*e->pl().getGeom()) / 1.5;
  }

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

  T_START(grid);
  std::cerr << "Build grid graph in " << T_STOP(grid) << " ms " << std::endl;

  NodePQ globalPq, dangling;

  std::set<CombNode*> settled;
  std::map<CombNode*, DPoint> newPositions;

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
        auto cmbEdg = ee.first;
        if (done.find(cmbEdg) != done.end()) continue;
        done.insert(cmbEdg);

        dangling.push(cmbEdg->getOtherNd(n));

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

        GridNode* frGrNd = gg->getGridNodeFrom(frCmbNd, gridSize * 1.7);
        // get surrounding displacement nodes
        double maxDis = getMaxDis(toCmbNd, cmbEdg, gridSize);
        std::set<GridNode*> toGrNds = gg->getGridNodesTo(toCmbNd, frGrNd, maxDis);


        // TODO: abort criteria
        size_t i = 0;
        while (!toGrNds.size() && i < 5) {
          i++;
          maxDis *= 2;
          toGrNds = gg->getGridNodesTo(toCmbNd, frGrNd, maxDis);
        }

        if (!frGrNd) {
          LOG(ERROR) << "Could not sort in source node " << frCmbNd;
          break;
        }

        if (toGrNds.size() == 0) {
          LOG(ERROR) << "Could not sort in target node " << toCmbNd;
          break;
        }

        // why not distance based? (TODO)
        double movePenPerGrid = toGrNds.size() > 1 ? 10 : 0;

        bool frWasClosed = frGrNd->pl().isClosed();
        bool toWasClosed = false;

        // open the target nodes
        for (auto n : toGrNds) {
          assert(n != frGrNd);
          double gridDist = floor(
              dist(*n->pl().getGeom(), *toCmbNd->pl().getGeom()) / gridSize);

          double topoPen =
              cg.changesTopology(toCmbNd, *n->pl().getGeom(), newPositions) *
              50;

          gg->openNodeSink(n, gridDist * movePenPerGrid + topoPen);
        }

        // open from source node
        gg->openNodeSink(frGrNd, 0);
        gg->openNode(frGrNd);

        NodeCost addCFromInv = writeNdCosts(frGrNd, frCmbNd, cmbEdg, gg);
        NodeCost addCToInv;

        if (toGrNds.size() == 1) {
          assert(*toGrNds.begin() != frGrNd);
          toWasClosed = (*toGrNds.begin())->pl().isClosed();
          gg->openNode(*toGrNds.begin());
          addCToInv = writeNdCosts(*toGrNds.begin(), toCmbNd, cmbEdg, gg);
        }

        GrEdgList res;
        GrNdList nList;
        GridNode* toGrNd = 0;
        Dijkstra::shortestPath(frGrNd, toGrNds, GridCost(),
                               GridHeur(gg, frGrNd, toGrNds), &res, &nList);

        if (nList.size()) {
          if (toGrNds.size() == 1) assert(toGrNd = *toGrNds.begin());
          toGrNd = nList.front();
        }

        gg->removeCostVector(frGrNd, addCFromInv);
        gg->removeCostVector(*toGrNds.begin(), addCToInv);

        // close the target node
        for (auto n : toGrNds) gg->closeNodeSink(n);
        // close the start node
        gg->closeNodeSink(frGrNd);

        if (toGrNd == 0) {
          LOG(ERROR) << "Could not route to target node " << toCmbNd;

          // rollback
          if (frWasClosed) gg->closeNode(frGrNd);
          if (toGrNds.size() == 1 && toWasClosed) gg->closeNode(*toGrNds.begin());

          break;
        }

        newPositions[toCmbNd] = *toGrNd->pl().getGeom();
        newPositions[frCmbNd] = *frGrNd->pl().getGeom();

        gg->closeNode(toGrNd);
        gg->closeNode(frGrNd);

        settleRes(frGrNd, toGrNd, gg, frCmbNd, toCmbNd, res, cmbEdg, &cg);

        gen++;
      }

      settled.insert(n);
    }
  }

  TransitGraph ret;
  writeTrGraph(&cg, gg, &ret);

  *retGg = gg;

  return ret;
}

// _____________________________________________________________________________
void Octilinearizer::settleRes(GridNode* frGrNd, GridNode* toGrNd,
                               GridGraph* gg, CombNode* from, CombNode* to,
                               const GrEdgList& res, CombEdge* e,
                               CombGraph* cg) {
  gg->settleGridNode(toGrNd, to);
  gg->settleGridNode(frGrNd, from);

  gg->addResidentEdges(e, frGrNd, res, toGrNd);

  for (auto f : res) {
    if (f->pl().isSecondary()) continue;
    gg->balanceEdge(f->getFrom()->pl().getParent(),
                    f->getTo()->pl().getParent());
  }

  gg->splitAlong(res);
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
NodeCost Octilinearizer::writeNdCosts(GridNode* n, CombNode* origNode,
                                      CombEdge* e, GridGraph* g) {
  NodeCost c = g->spacingPenalty(n, origNode, e);
  // c += g->outDegDeviationPenalty(origNode, e);

  return g->addCostVector(n, c);
}

// _____________________________________________________________________________
void Octilinearizer::writeTrGraph(CombGraph* cg, GridGraph* gg,
                                  TransitGraph* tg) {
  std::map<std::pair<size_t, size_t>, TransitNode*> nodes;

  for (auto cn : *cg->getNds()) {
    for (auto ce : cn->getAdjList()) {
      if (ce->getFrom() != cn) continue;

      for (size_t i = 1; i < gg->getGridNodes(ce).size(); i++) {
        auto gn1 = gg->getGridNodes(ce)[i - 1];
        auto gn2 = gg->getGridNodes(ce)[i];

        if (nodes.find({gn1->pl().getX(), gn1->pl().getY()}) == nodes.end()) {
          nodes[{gn1->pl().getX(), gn1->pl().getY()}] =
              tg->addNd(*gn1->pl().getGeom());
        }

        if (nodes.find({gn2->pl().getX(), gn2->pl().getY()}) == nodes.end()) {
          nodes[{gn2->pl().getX(), gn2->pl().getY()}] =
              tg->addNd(*gn2->pl().getGeom());
        }

        auto e = tg->addEdg(nodes[{gn1->pl().getX(), gn1->pl().getY()}],
                   nodes[{gn2->pl().getX(), gn2->pl().getY()}]);

        for (auto ro : ce->pl().getChilds().front()->pl().getRoutes()) {
          e->pl().addRoute(ro.route, 0);
        }
      }
    }

  }
}
