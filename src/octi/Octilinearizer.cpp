// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/Octilinearizer.h"
#include "util/Misc.h"
#include "util/geo/BezierCurve.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Dijkstra.h"
#include "util/log/Log.h"

using namespace octi;
using namespace gridgraph;

using transitgraph::EdgeOrdering;

// _____________________________________________________________________________
void Octilinearizer::writeEdgeOrdering(CombGraph* target) {
  for (auto n : *target->getNds()) {
    n->pl().setEdgeOrdering(getEdgeOrderingForNode(n));
  }
}

// _____________________________________________________________________________
EdgeOrdering Octilinearizer::getEdgeOrderingForNode(CombNode* n) const {
  return getEdgeOrderingForNode(n, true, std::map<CombNode*, DPoint>());
}

// _____________________________________________________________________________
EdgeOrdering Octilinearizer::getEdgeOrderingForNode(
    CombNode* n, bool useOrigNextNode,
    const std::map<CombNode*, DPoint>& newPos) const {
  EdgeOrdering order;
  for (auto e : n->getAdjList()) {
    auto r = e->pl().getChilds().front();
    DPoint a = *n->pl().getGeom();
    if (newPos.find(n) != newPos.end()) a = newPos.find(n)->second;

    DPoint b;
    if (useOrigNextNode) {
      b = *r->getOtherNd(n->pl().getParent())->pl().getGeom();
    } else {
      auto other = e->getOtherNd(n);
      if (e->pl().getGeom()->size() > 2) {
        if (e->getTo() == n) {
          b = e->pl().getGeom()->at(e->pl().getGeom()->size() - 2);
        } else {
          b = e->pl().getGeom()->at(1);
        }
      } else {
        b = *other->pl().getGeom();
        if (newPos.find(other) != newPos.end()) b = newPos.find(other)->second;
      }
    }

    // get the angles
    double deg = util::geo::angBetween(a, b);

    order.add(e, deg);
  }

  return order;
}

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
void Octilinearizer::normalizeCostVector(double* vec) const {
  double smallest = std::numeric_limits<double>::max();
  for (int i = 0; i < 8; i++) {
    if (vec[i] > -1 && vec[i] < smallest) {
      smallest = vec[i];
    }
  }

  for (int i = 0; i < 8; i++) {
    if (vec[i] > -1) {
      vec[i] -= smallest;
    }
  }
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

  writeEdgeOrdering(&cg);

  auto box = tg->getBBox();

  TransitGraph bestGraph;
  GridGraph* g;

  T_START(gridgraph);
  g = new GridGraph(box, gridSize, pens);
  std::cerr << "Build grid graph in " << T_STOP(gridgraph) << " ms "
            << std::endl;

  // comparator for nodes, based on degree
  struct NodeCompare {
    bool operator()(CombNode* a, CombNode* b) {
      // return a->getAdjList().size() < b->getAdjList().size();
      return a->getAdjList().size() < b->getAdjList().size() ||
             (a->getAdjList().size() == b->getAdjList().size() &&
              a->pl().getRouteNumber() < b->pl().getRouteNumber());
    }
  };

  std::priority_queue<CombNode*, std::vector<CombNode*>, NodeCompare> globalPq;
  std::priority_queue<CombNode*, std::vector<CombNode*>, NodeCompare> dangling;

  std::set<CombNode*> settled;
  std::map<CombNode*, DPoint> newPositions;

  for (auto n : *cg.getNds()) {
    globalPq.push(n);
  }

  std::set<CombEdge*> done;
  int64_t gen = 0;
  double totalCost = 0;

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

        if (from->getAdjList().size() < to->getAdjList().size() ||
            (from->getAdjList().size() == to->getAdjList().size() &&
             from->pl().getRouteNumber() < to->pl().getRouteNumber())) {
          auto tmp = from;
          from = to;
          to = tmp;
          reversed = !reversed;
        }

        // TODO: move into if clause above
        if (!g->isSettled(from) &&
            g->isSettled(
                to)) {  // m.find(from) == m.end() && m.find(to) != m.end()) {
          auto tmp = from;
          from = to;
          to = tmp;
          reversed = !reversed;
        }

        std::cerr << "\n\n ++ Generation " << gen << std::endl;
        std::cerr << "Edge from " << from->pl().toString() << " to "
                  << to->pl().toString() << std::endl;

        GridNode* fromGridNode = g->getGridNodeFrom(from, gridSize * 1.7);

        if (!fromGridNode) {
          LOG(ERROR) << "Could not sort in source node " << from << std::endl;
        }

        // get surrounding displacement nodes
        double maxDis = getMaxDis(to, e, gridSize);
        std::set<GridNode*> tos = g->getGridNodesTo(to, maxDis);

        if (tos.size() == 0) {
          LOG(ERROR) << "Could not sort in target node " << to << std::endl;
        }

        // why not distance based? (TODO)
        double movePenPerGrid = tos.size() > 1 ? 10 : 0;

        // open to target node
        for (auto n : tos) {
          double gridDist =
              floor(util::geo::dist(*n->pl().getGeom(), *to->pl().getGeom()) /
                    gridSize);

          double topoPen =
              changesTopology(to, *n->pl().getGeom(), newPositions) * 50;

          g->openNodeSink(n, gridDist * movePenPerGrid + topoPen);
          g->openNode(n);
        }

        // open from source node
        g->openNodeSink(fromGridNode, 0);
        g->openNode(fromGridNode);

        double addCFrom[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        double addCTo[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        double addCFromInv[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        double addCToInv[8] = {0, 0, 0, 0, 0, 0, 0, 0};

        g->spacingPenalty(fromGridNode, from, e, addCFrom);
        g->topoBlockPenalty(fromGridNode, from, e, addCFrom);
        normalizeCostVector(addCFrom);
        g->outDegDeviationPenalty(fromGridNode, from, e, addCFrom);

        g->addCostVector(fromGridNode, addCFrom, addCFromInv);

        if (tos.size() == 1) {
          g->spacingPenalty(*tos.begin(), to, e, addCTo);
          g->topoBlockPenalty(*tos.begin(), to, e, addCTo);
          normalizeCostVector(addCTo);
          g->outDegDeviationPenalty(*tos.begin(), to, e, addCTo);

          g->addCostVector(*tos.begin(), addCTo, addCToInv);
        }
        int STOP_GEN = -1;

        util::graph::Dijkstra::EList<GridNodePL, GridEdgePL> res;
        util::graph::Dijkstra::NList<GridNodePL, GridEdgePL> nList;
        GridNode* target = 0;
        util::graph::Dijkstra::shortestPath(fromGridNode, tos, GridCost(),
                                            GridHeur(g, fromGridNode, tos),
                                            &res, &nList);
        std::reverse(res.begin(), res.end());

        if (nList.size()) target = nList.front();

        g->removeCostVector(fromGridNode, addCFromInv);
        if (tos.size() == 1) {
          g->removeCostVector(*tos.begin(), addCToInv);
        }

        if (gen == STOP_GEN) {
          util::geo::output::GeoGraphJsonOutput out;
          out.print(*g, std::cout);
          exit(0);
        }

        if (target == 0) {
          LOG(ERROR) << "Could not sort in node " << to << std::endl;
          for (auto n : tos) {
            g->closeNodeSink(n);
          }
          g->closeNodeSink(fromGridNode);
          g->closeNode(fromGridNode);
          continue;
        }

        newPositions[to] = *target->pl().getGeom();
        newPositions[from] = *fromGridNode->pl().getGeom();
        g->settleGridNode(target, to);
        g->settleGridNode(fromGridNode, from);

        // write everything to the result graph
        PolyLine<double> pl;
        buildPolylineFromRes(res, pl);
        if (reversed) pl.reverse();

        double cost = getCostFromRes(res);
        addResidentEdges(g, e, res);

        double estimatedCost = g->heurCost(
            g->getNodeCoords(fromGridNode).first,
            g->getNodeCoords(fromGridNode).second,
            g->getNodeCoords(target).first, g->getNodeCoords(target).second);

        std::cerr << std::round(cost) << " vs " << std::round(estimatedCost)
                  << std::endl;
        assert(std::round(cost) >= std::round(estimatedCost));

        totalCost += cost;

        // close the target node
        for (auto n : tos) {
          g->closeNodeSink(n);
        }
        g->closeNode(target);
        g->closeNodeSink(fromGridNode);
        g->closeNode(fromGridNode);

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

  bestGraph = TransitGraph();
  cg.getTransitGraph(&bestGraph);

  *gg = g;

  return bestGraph;
}

// _____________________________________________________________________________
void Octilinearizer::buildPolylineFromRes(const std::vector<GridEdge*>& res,
                                          PolyLine<double>& pl) {
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
size_t Octilinearizer::changesTopology(
    CombNode* nOrig, DPoint p,
    const std::map<CombNode*, DPoint>& newPos) const {
  // collect the affected nodes
  size_t ret = 0;
  std::set<CombNode*> aff;
  auto newPosA = newPos;
  newPosA[nOrig] = p;
  for (auto e : nOrig->getAdjList()) {
    if (newPos.find(e->getFrom()) != newPos.end()) aff.insert(e->getFrom());
    if (newPos.find(e->getTo()) != newPos.end()) aff.insert(e->getTo());
  }

  for (auto n : aff) {
    if (!getEdgeOrderingForNode(n, false, newPosA)
             .equals(getEdgeOrderingForNode(n))) {
      ret += 1;
    }
  }

  return ret;
}
