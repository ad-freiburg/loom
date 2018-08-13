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

using transitgraph::EdgeOrdering;
using util::graph::Dijkstra;
using util::graph::Dijkstra;

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

          n->pl().setGeom(util::geo::DPoint(
              (n->pl().getGeom()->getX() + otherP->getX()) / 2,
              (n->pl().getGeom()->getY() + otherP->getY()) / 2));
          goto start;
        }
      }
    }
  }
}

// _____________________________________________________________________________
double Octilinearizer::getMaxDis(CombNode* to, CombEdge* e, double gridSize) {
  double tooMuch = gridSize * 4;

  if (to->getAdjList().size() == 1) {
    return util::geo::len(*e->pl().getGeom()) / 1.5;
  }

  if (to->getAdjList().size() > 1) {
    if (e->pl().getChilds().size() > 5 &&
        util::geo::len(*e->pl().getGeom()) / e->pl().getChilds().size() >
            tooMuch) {
      return ((util::geo::len(*e->pl().getGeom()) /
               e->pl().getChilds().size()) -
              tooMuch) *
             e->pl().getChilds().size();
    }
  }

  return gridSize * 1.7;
}

// _____________________________________________________________________________
TransitGraph Octilinearizer::draw(TransitGraph* tg, GridGraph** gg,
                                  const Penalties& pens) {
  double gridSize = 250;

  std::cerr << "Removing short edges... ";
  T_START(remshortegs);
  removeEdgesShorterThan(tg, gridSize / 2);
  std::cerr << " done (" << T_STOP(remshortegs) << "ms)" << std::endl;

  std::cerr << "Building combination graph... ";
  T_START(combgraph);
  CombGraph cg(tg);
  std::cerr << " done (" << T_STOP(combgraph) << "ms)" << std::endl;

  auto box = tg->getBBox();

  GridGraph* g;

  T_START(grid);
  g = new GridGraph(box, gridSize, pens);
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
        auto e = ee.first;
        if (done.find(e) != done.end()) continue;
        done.insert(e);

        dangling.push(e->getOtherNd(n));

        auto from = e->getFrom();
        auto to = e->getTo();

        bool reversed = false;

        if (from->getDeg() < to->getDeg() ||
            (from->getDeg() == to->getDeg() &&
             from->pl().getRouteNumber() < to->pl().getRouteNumber())) {
          auto tmp = from;
          from = to;
          to = tmp;
          reversed = !reversed;
        }

        // TODO: move into if clause above
        if (!g->isSettled(from) && g->isSettled(to)) {
          auto tmp = from;
          from = to;
          to = tmp;
          reversed = !reversed;
        }

        std::cerr << "\n\n ++ Generation " << gen << std::endl;
        std::cerr << "Edge from " << from->pl().toString() << " to "
                  << to->pl().toString() << std::endl;

        GridNode* startGridNd = g->getGridNodeFrom(from, gridSize * 1.7);

        if (!startGridNd) {
          LOG(ERROR) << "Could not sort in source node " << from << std::endl;
        }

        // get surrounding displacement nodes
        double maxDis = getMaxDis(to, e, gridSize);
        std::set<GridNode*> toGridNds = g->getGridNodesTo(to, maxDis);

        if (toGridNds.size() == 0) {
          LOG(ERROR) << "Could not sort in target node " << to << std::endl;
        }

        // why not distance based? (TODO)
        double movePenPerGrid = toGridNds.size() > 1 ? 10 : 0;

        // open the target nodes
        for (auto n : toGridNds) {
          double gridDist =
              floor(util::geo::dist(*n->pl().getGeom(), *to->pl().getGeom()) /
                    gridSize);

          double topoPen =
              cg.changesTopology(to, *n->pl().getGeom(), newPositions) * 50;

          g->openNodeSink(n, gridDist * movePenPerGrid + topoPen);
          g->openNode(n);
        }

        // open from source node
        g->openNodeSink(startGridNd, 0);
        g->openNode(startGridNd);

        NodeCost addCFromInv = writeNodeCosts(startGridNd, from, e, g);
        NodeCost addCToInv;

        if (toGridNds.size() == 1) {
          addCToInv = writeNodeCosts(*toGridNds.begin(), to, e, g);
        }

        Dijkstra::EList<GridNodePL, GridEdgePL> res;
        Dijkstra::NList<GridNodePL, GridEdgePL> nList;
        GridNode* toGridNd = 0;
        Dijkstra::shortestPath(startGridNd, toGridNds, GridCost(),
                               GridHeur(g, startGridNd, toGridNds), &res,
                               &nList);
        std::reverse(res.begin(), res.end());

        if (nList.size()) toGridNd = nList.front();

        g->removeCostVector(startGridNd, addCFromInv);
        g->removeCostVector(*toGridNds.begin(), addCToInv);

        if (toGridNd == 0) {
          LOG(ERROR) << "Could not sort in node " << to << std::endl;
          for (auto n : toGridNds) {
            g->closeNodeSink(n);
          }
          g->closeNodeSink(startGridNd);
          g->closeNode(startGridNd);
          continue;
        }

        newPositions[to] = *toGridNd->pl().getGeom();
        newPositions[from] = *startGridNd->pl().getGeom();
        g->settleGridNode(toGridNd, to);
        g->settleGridNode(startGridNd, from);

        // write everything to the result graph
        PolyLine<double> pl = buildPolylineFromRes(res);
        if (reversed) pl.reverse();

        addResidentEdges(g, e, res);

        // close the target node
        for (auto n : toGridNds) {
          g->closeNodeSink(n);
        }
        g->closeNode(toGridNd);
        g->closeNodeSink(startGridNd);
        g->closeNode(startGridNd);

        for (auto f : res) {
          if (f->pl().isSecondary()) continue;
          g->balanceEdge(f->getFrom()->pl().getParent(),
                         f->getTo()->pl().getParent());
        }

        e->pl().setPolyLine(pl);
        e->pl().setGeneration(gen);

        gen++;
      }

      settled.insert(n);
    }
  }

  TransitGraph ret;
  cg.getTransitGraph(&ret);

  *gg = g;

  return ret;
}

// _____________________________________________________________________________
PolyLine<double> Octilinearizer::buildPolylineFromRes(
    const std::vector<GridEdge*>& res) {
  PolyLine<double> pl;
  for (auto f : res) {
    if (!f->pl().isSecondary()) {
      if (pl.getLine().size() > 0) {
        BezierCurve<double> bc(pl.getLine().back(),
                               *f->getFrom()->pl().getParent()->pl().getGeom(),
                               *f->getFrom()->pl().getParent()->pl().getGeom(),
                               *f->getFrom()->pl().getGeom());

        for (auto p : bc.render(10).getLine()) {
          pl << p;
        }
      } else {
        pl << *f->getFrom()->pl().getParent()->pl().getGeom();
      }

      pl << *f->getFrom()->pl().getGeom();
      pl << *f->getTo()->pl().getGeom();
    }
  }

  if (res.size()) pl << *res.back()->getTo()->pl().getParent()->pl().getGeom();

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
void Octilinearizer::addResidentEdges(GridGraph* g, CombEdge* e,
                                      const std::vector<GridEdge*>& res) {
  for (auto f : res) {
    assert(f->pl().getResEdges().size() == 0);
    f->pl().addResidentEdge(e);
    assert(g->getEdg(f->getTo(), f->getFrom())->pl().getResEdges().size() == 0);
    g->getEdg(f->getTo(), f->getFrom())->pl().addResidentEdge(e);
  }
}

// _____________________________________________________________________________
NodeCost Octilinearizer::writeNodeCosts(GridNode* n, CombNode* origNode,
                                    CombEdge* e, GridGraph* g) {
  NodeCost c = g->spacingPenalty(n, origNode, e);
  c += g->topoBlockPenalty(n, origNode, e);
  NodeCost test = c;
  // c.normalize();
  c += g->outDegDeviationPenalty(origNode, e);

  return g->addCostVector(n, c);
}
