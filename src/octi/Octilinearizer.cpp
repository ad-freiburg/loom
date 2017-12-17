// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <chrono>
#include "octi/Octilinearizer.h"
#include "util/geo/BezierCurve.h"
#include "util/geo/output/GeoJsonOutput.h"
#include "util/graph/Dijkstra.h"
#include "util/log/Log.h"

using namespace octi;
using namespace gridgraph;

// _____________________________________________________________________________
void Octilinearizer::rotate(CombGraph* g, Point center, double deg) {
  for (auto n : *g->getNodes()) {
    n->pl().getParent()->pl().setGeom(
        util::geo::rotate(*n->pl().getGeom(), deg, center));
    for (auto e : n->getAdjListOut()) {
      e->pl().setPolyLine(
          PolyLine(util::geo::rotate(*e->pl().getGeom(), deg, center)));
    }
  }
}

// _____________________________________________________________________________
void Octilinearizer::buildCombGraph(TransitGraph* source, CombGraph* target) {
  auto nodes = source->getNodes();

  std::map<TransitNode*, CombNode*> m;

  for (auto n : *nodes) {
    CombNode* cn = new CombNode(n);
    target->addNode(cn);
    m[n] = cn;
  }

  for (auto n : *nodes) {
    for (auto e : n->getAdjListOut()) {
      target->addEdge(m[e->getFrom()], m[e->getTo()],
                      octi::graph::CombEdgePL(e));
    }
  }

  for (auto n : *target->getNodes()) {
    size_t routes = 0;
    for (auto e : n->getAdjList()) {
      routes += e->pl().getChilds().front()->pl().getRoutes().size();
    }
    n->pl().setRouteNumber(routes);
  }
}

// _____________________________________________________________________________
void Octilinearizer::writeEdgeOrdering(CombGraph* target) {
  for (auto n : *target->getNodes()) {
    n->pl().setEdgeOrdering(getEdgeOrderingForNode(n));
  }
}

// _____________________________________________________________________________
graph::EdgeOrdering Octilinearizer::getEdgeOrderingForNode(CombNode* n) const {
  return getEdgeOrderingForNode(n, true, std::map<CombNode*, Point>());
}

// _____________________________________________________________________________
graph::EdgeOrdering Octilinearizer::getEdgeOrderingForNode(CombNode* n,
    bool useOrigNextNode, const std::map<CombNode*, Point>& newPos) const {
  graph::EdgeOrdering order;
  for (auto e : n->getAdjList()) {
    auto r = e->pl().getChilds().front();
    Point a = *n->pl().getGeom();
    if (newPos.find(n) != newPos.end()) a = newPos.find(n)->second;

    Point b;
    if (useOrigNextNode) {
      b = *r->getOtherNode(n->pl().getParent())->pl().getGeom();
    } else {
      auto other = e->getOtherNode(n);
      if (e->pl().getGeom()->size() > 2) {
        if (e->getTo() == n) {
          b = e->pl().getGeom()->at(e->pl().getGeom()->size()-2);
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
void Octilinearizer::buildTransitGraph(CombGraph* source,
                                       TransitGraph* target) {
  auto nodes = source->getNodes();

  std::map<TransitNode*, TransitNode*> m;

  for (auto n : *nodes) {
    for (auto f : n->getAdjListOut()) {
      double tot = f->pl().getChilds().size();
      double d = f->pl().getPolyLine().getLength();
      double step = d / tot;

      int i = 0;

      auto pre = n->pl().getParent();

      for (auto e : f->pl().getChilds()) {
        auto from = e->getFrom();
        auto to = e->getTo();

        PolyLine pl = f->pl().getPolyLine().getSegment((step * i) / d,
                                                       (step * (i + 1)) / d);

        if (from == pre) {
          pre = to;
        } else {
          pl.reverse();
          pre = from;
        }

        if (m.find(from) == m.end()) {
          auto payload = from->pl();
          payload.setGeom(pl.getLine().front());
          auto tfrom = new TransitNode(payload);
          m[from] = tfrom;
          target->addNode(tfrom);
        }

        if (m.find(to) == m.end()) {
          auto payload = to->pl();
          payload.setGeom(pl.getLine().back());
          auto tto = new TransitNode(payload);
          m[to] = tto;
          target->addNode(tto);
        }

        auto payload = e->pl();
        payload.setPolyline(pl);
        payload.setGeneration(f->pl().getGeneration());
        target->addEdge(m[from], m[to], payload);

        i++;
      }
    }
  }
}

// _____________________________________________________________________________
void Octilinearizer::combineDeg2(CombGraph* g) {
  std::set<CombNode*> nodes = *g->getNodes();
  for (auto n : nodes) {
    if (n->getAdjList().size() == 2) {
      CombEdge* a = 0;
      CombEdge* b = 0;

      for (auto e : n->getAdjListOut()) {
        if (!a)
          a = g->getEdge(n, e->getOtherNode(n));
        else
          b = g->getEdge(n, e->getOtherNode(n));
      }

      for (auto e : n->getAdjListIn()) {
        if (!a)
          a = g->getEdge(e->getOtherNode(n), n);
        else if (!b)
          b = g->getEdge(e->getOtherNode(n), n);
      }

      // if this combination would turn our graph into a multigraph,
      // dont do it!
      if (g->getEdge(a->getFrom(), b->getOtherNode(n)) ||
          g->getEdge(b->getOtherNode(n), a->getTo()))
        continue;

      auto pl = a->pl();

      // a is the reference
      if (a->getTo() == n) {
        if (b->getTo() != n) {
          pl.getChilds().insert(pl.getChilds().end(),
                                b->pl().getChilds().begin(),
                                b->pl().getChilds().end());
        } else {
          pl.getChilds().insert(pl.getChilds().end(),
                                b->pl().getChilds().rbegin(),
                                b->pl().getChilds().rend());
        }

        pl.setPolyLine(PolyLine(*a->getFrom()->pl().getGeom(),
                                *b->getOtherNode(n)->pl().getGeom()));
        g->addEdge(a->getFrom(), b->getOtherNode(n), pl);
      } else {
        if (b->getTo() == n) {
          pl.getChilds().insert(pl.getChilds().begin(),
                                b->pl().getChilds().begin(),
                                b->pl().getChilds().end());
        } else {
          pl.getChilds().insert(pl.getChilds().begin(),
                                b->pl().getChilds().rbegin(),
                                b->pl().getChilds().rend());
        }

        pl.setPolyLine(PolyLine(*b->getOtherNode(n)->pl().getGeom(),
                                *a->getTo()->pl().getGeom()));
        g->addEdge(b->getOtherNode(n), a->getTo(), pl);
      }

      g->deleteNode(n);
    }
  }
}

// _____________________________________________________________________________
void Octilinearizer::removeEdgesShorterThan(TransitGraph* g, double d) {
start:
  for (auto n1 : *g->getNodes()) {
    for (auto e1 : n1->getAdjList()) {
      if (e1->pl().getPolyline().getLength() < d) {
        if (e1->getOtherNode(n1)->getAdjList().size() > 1 &&
            n1->getAdjList().size() > 1 &&
            (n1->pl().getStops().size() == 0 ||
             e1->getOtherNode(n1)->pl().getStops().size() == 0)) {
          auto otherP = e1->getFrom()->pl().getGeom();

          TransitNode* n = 0;

          if (e1->getTo()->pl().getStops().size() > 0) {
            n = g->mergeNodes(e1->getFrom(), e1->getTo());
          } else {
            n = g->mergeNodes(e1->getTo(), e1->getFrom());
          }

          n->pl().setGeom(util::geo::Point(
              (n->pl().getGeom()->get<0>() + otherP->get<0>()) / 2,
              (n->pl().getGeom()->get<1>() + otherP->get<1>()) / 2));
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
  double gridSize = 50;
  std::cerr << "Removing short edges... ";
  auto begin = std::chrono::steady_clock::now();
  removeEdgesShorterThan(tg, gridSize / 2);
  auto end = std::chrono::steady_clock::now();
  std::cerr << " done ("
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     begin)
                   .count()
            << "ms)" << std::endl;

  std::cerr << "Building combination graph... ";
  begin = std::chrono::steady_clock::now();
  CombGraph cg(true);
  buildCombGraph(tg, &cg);
  end = std::chrono::steady_clock::now();
  std::cerr << " done ("
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     begin)
                   .count()
            << "ms)" << std::endl;

  std::cerr << "Combining deg 2 nodes...";
  begin = std::chrono::steady_clock::now();
  combineDeg2(&cg);
  end = std::chrono::steady_clock::now();
  std::cerr << " done ("
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     begin)
                   .count()
            << "ms)" << std::endl;

  writeEdgeOrdering(&cg);


  auto box = tg->getBBox();

  // pad bbox with 4 grid cell sizes
  box.max_corner().set<0>(box.max_corner().get<0>() + 4 * gridSize);
  box.max_corner().set<1>(box.max_corner().get<1>() + 4 * gridSize);
  box.min_corner().set<0>(box.min_corner().get<0>() - 4 * gridSize);
  box.min_corner().set<1>(box.min_corner().get<1>() - 4 * gridSize);

  Point center;
  bgeo::centroid(box, center);

  double maxDeg = 0;  // 15;
  double step = 5;

  rotate(&cg, center, -maxDeg - step);
  RotatedBox rbox(box, -maxDeg - step, center);

  TransitGraph bestGraph;
  GridGraph* g;
  double bestScore = std::numeric_limits<double>::infinity();

  for (int i = 0; i < 1 + (maxDeg / step) * 2; i++) {
    rotate(&cg, center, step);
    rbox.rotateDeg += step;

    std::cerr << "Rotate degree is " << rbox.rotateDeg << std::endl;

    auto gStart = std::chrono::steady_clock::now();
    g = new GridGraph(getBoundingBox(rbox.getPolygon()), gridSize, pens);
    auto gEnd = std::chrono::steady_clock::now();
    std::cerr << "Build grid graph in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(gEnd -
                                                                       gStart)
                     .count()
              << " ms " << std::endl;
    util::graph::Dijkstra dijkstra;

    // comparator for nodes, based on degree
    struct NodeCompare {
      bool operator()(CombNode* a, CombNode* b) {
        // return a->getAdjList().size() < b->getAdjList().size();
        return a->getAdjList().size() < b->getAdjList().size() ||
               (a->getAdjList().size() == b->getAdjList().size() &&
                a->pl().getRouteNumber() < b->pl().getRouteNumber());
      }
    };

    std::priority_queue<CombNode*, std::vector<CombNode*>, NodeCompare>
        globalPq;
    std::priority_queue<CombNode*, std::vector<CombNode*>, NodeCompare>
        dangling;

    std::set<CombNode*> settled;
    std::map<CombNode*, Point> newPositions;

    for (auto n : *cg.getNodes()) {
      globalPq.push(n);
    }

    std::set<CombEdge*> done;
    size_t gen = 0;
    double totalCost = 0;
    size_t totalIters = 0;
    int dTime = 0;
    int costVecTime = 0;
    int postTime = 0;

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

          dangling.push(e->getOtherNode(n));

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

          auto costVecBegin = std::chrono::steady_clock::now();
          // get surrounding displacement nodes
          double maxDis = getMaxDis(to, e, gridSize);
          std::unordered_set<GridNode*> tos = g->getGridNodesTo(to, maxDis);

          if (tos.size() == 0) {
            LOG(ERROR) << "Could not sort in target node " << to << std::endl;
          }

          std::list<GridEdge*> res;
          GridNode* target = 0;

          // why not distance based? (TODO)
          double movePenPerGrid = tos.size() > 1 ? 10 : 0;


          // open to target node
          for (auto n : tos) {
            double gridDist =
                floor(util::geo::dist(*n->pl().getGeom(), *to->pl().getGeom()) /
                      gridSize);

            double topoPen = 0; //changesTopology(e, *n->pl().getGeom(), newPositions) * 50;

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
          auto costVecEnd = std::chrono::steady_clock::now();
          costVecTime += std::chrono::duration_cast<std::chrono::milliseconds>(
                             costVecEnd - costVecBegin)
                             .count();

          int STOP_GEN = -1;

          auto dBegin = std::chrono::steady_clock::now();
          int iters = dijkstra.shortestPathAStar(
              fromGridNode, tos, GridHeur(g, fromGridNode, tos), &res, &target,
              gen == STOP_GEN);
          auto dEnd = std::chrono::steady_clock::now();
          dTime += std::chrono::duration_cast<std::chrono::milliseconds>(dEnd -
                                                                         dBegin)
                       .count();

          totalIters += iters;


          auto postBegin = std::chrono::steady_clock::now();
          g->removeCostVector(fromGridNode, addCFromInv);
          if (tos.size() == 1) {
            g->removeCostVector(*tos.begin(), addCToInv);
          }

          if (gen == STOP_GEN) {
            util::geo::output::GeoJsonOutput out;
            out.print(*g);
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
          PolyLine pl;
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

          auto postEnd = std::chrono::steady_clock::now();
          postTime += std::chrono::duration_cast<std::chrono::milliseconds>(
                          postEnd - postBegin)
                          .count();

          gen++;
        }

        settled.insert(n);
      }
    }

    std::cerr << " === Total edge cost in graph is " << totalCost
              << " === " << std::endl;

    std::cerr << " === Total dijkstra iterations were " << totalIters
              << " === " << std::endl;

    std::cerr << " === Total dijkstra time was " << dTime
              << " ms === " << std::endl;

    std::cerr << " === Total cost vector calc time was " << costVecTime
              << " ms === " << std::endl;
    std::cerr << " === Total post time was " << postTime
              << " ms === " << std::endl;

    if (totalCost < bestScore) {
      bestScore = totalCost;
      bestGraph = TransitGraph();
      buildTransitGraph(&cg, &bestGraph);
    }
  }

  *gg = g;

  return bestGraph;
}

// _____________________________________________________________________________
void Octilinearizer::buildPolylineFromRes(const std::list<GridEdge*>& res,
                                          PolyLine& pl) {
  for (auto f : res) {
    if (!f->pl().isSecondary()) {
      if (pl.getLine().size() > 0) {
        BezierCurve bc(pl.getLine().back(),
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
double Octilinearizer::getCostFromRes(const std::list<GridEdge*>& res) {
  size_t c = 0;
  double cost = 0;
  for (auto f : res) {
    cost += c > 0 ? f->pl().cost() : 0;
    c++;
  }
  return cost;
}

// _____________________________________________________________________________
double Octilinearizer::addResidentEdges(GridGraph* g, CombEdge* e,
                                        const std::list<GridEdge*>& res) {
  for (auto f : res) {
    assert(f->pl().getResEdges().size() == 0);
    f->pl().addResidentEdge(e);
    assert(g->getEdge(f->getTo(), f->getFrom())->pl().getResEdges().size() ==
           0);
    g->getEdge(f->getTo(), f->getFrom())->pl().addResidentEdge(e);
  }
}

// _____________________________________________________________________________
size_t Octilinearizer::changesTopology(
    CombEdge* eCheck, Point p, const std::map<CombNode*, Point>& newPos) const {
  // collect the affected nodes
  size_t ret = 0;
  std::set<CombNode*> aff;
  CombNode* nOrig= eCheck->getTo();
  auto newPosA = newPos;
  newPosA[nOrig] = p;
  for (auto e : nOrig->getAdjList()) {
    aff.insert(e->getFrom());
    aff.insert(e->getTo());
  }

  for (auto n : aff) {
    if (!getEdgeOrderingForNode(n, false, newPosA).equals(getEdgeOrderingForNode(n, false, newPos))) {
      ret += 1;
    }
  }

  return ret;
}
