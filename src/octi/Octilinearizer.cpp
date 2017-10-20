// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <chrono>
#include "util/geo/output/GeoJsonOutput.h"
#include "octi/Octilinearizer.h"
#include "util/geo/BezierCurve.h"
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
    n->pl().getOrderedEdges().clear();
  }

  for (auto n : *target->getNodes()) {
    for (auto e : n->getAdjListOut()) {
      auto refEdgeA = e->pl().getChilds().front();
      auto refEdgeB = e->pl().getChilds().back();

      // get the angles
      double degA = util::geo::angBetween(
          *n->pl().getParent()->pl().getGeom(),
          *refEdgeA->getOtherNode(n->pl().getParent())->pl().getGeom());
      double degB = util::geo::angBetween(
          *e->getTo()->pl().getParent()->pl().getGeom(),
          *refEdgeB->getOtherNode(e->getTo()->pl().getParent())
               ->pl()
               .getGeom());

      n->pl().addOrderedEdge(e, degA);
      e->getTo()->pl().addOrderedEdge(e, degB);
    }
  }
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
double Octilinearizer::getMaxDis(CombNode* to, CombEdge* e) {
  double tooMuch = 1000;

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

  return 400;
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
TransitGraph Octilinearizer::draw(TransitGraph* tg, const Penalties& pens) {
  std::cerr << "Removing short edges... ";
  auto begin = std::chrono::steady_clock::now();
  removeEdgesShorterThan(tg, 200);
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
  double gridSize = 250;

  Point center;
  bgeo::centroid(box, center);

  double maxDeg = 0;  // 15;
  double step = 5;

  rotate(&cg, center, -maxDeg - step);
  RotatedBox rbox(box, -maxDeg - step, center);

  TransitGraph bestGraph;
  double bestScore = std::numeric_limits<double>::infinity();
  bool verb = false;

  for (int i = 0; i < 1 + (maxDeg / step) * 2; i++) {
    rotate(&cg, center, step);
    rbox.rotateDeg += step;

    std::cerr << "Rotate degree is " << rbox.rotateDeg << std::endl;

    std::map<CombNode*, GridNode*> m;
    std::map<CombNode*, TransitNode*> oldNew;
    std::set<GridNode*> used;

    gridgraph::GridGraph g(getBoundingBox(rbox.getPolygon()), gridSize, pens);

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

    for (auto n : *cg.getNodes()) {
      globalPq.push(n);
    }

    std::set<CombEdge*> done;
    size_t gen = 0;
    double totalCost = 0;

    while (!globalPq.empty()) {
      auto n = globalPq.top();
      globalPq.pop();
      dangling.push(n);

      while (!dangling.empty()) {
        auto n = dangling.top();
        dangling.pop();

        if (settled.find(n) != settled.end()) continue;

        for (auto ee : n->pl().getOrderedEdges()) {
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
            std::cerr << "A" << std::endl;
          }

          if (m.find(from) == m.end() && m.find(to) != m.end()) {
            auto tmp = from;
            from = to;
            to = tmp;
            reversed = !reversed;
          }

          std::cerr << "\n\n ++ Generation " << gen << std::endl;
          std::cerr << "Edge from " << from->pl().toString() << " to " << to->pl().toString() << std::endl;

          if (m.find(from) == m.end()) {
            auto cands = g.getNearestCandidatesFor(*from->pl().getGeom(), 400);

            while (!cands.empty()) {
              if (!cands.top().n->pl().isClosed() && used.find(cands.top().n) == used.end()) {
                m[from] = cands.top().n;
                used.insert(cands.top().n);
                break;
              }
              cands.pop();
            }

            if (m.find(from) == m.end()) {
              LOG(ERROR) << "Could not sort in node " << from << std::endl;
            }
          }

          // get surrounding displacement nodes
          std::unordered_map<GridNode*, bool> tos;

          if (m.find(to) == m.end()) {
            double maxDis = getMaxDis(to, e);

            auto cands = g.getNearestCandidatesFor(*to->pl().getGeom(), maxDis);

            while (!cands.empty()) {
              if (!cands.top().n->pl().isClosed() && used.find(cands.top().n) == used.end()) {
                tos[cands.top().n] = true;
              }
              cands.pop();
            }
          } else {
            tos[m.find(to)->second] = true;
          }

          std::list<GridEdge*> res;
          GridNode* target = 0;

          double movePenPerGrid = tos.size() > 1 ? 10 : 0;

          // open to target node
          for (auto n : tos) {
            double gridDist = floor(util::geo::dist(*n.first->pl().getGeom(), *to->pl().getGeom()) / gridSize);
            g.openNodeSink(n.first, gridDist * movePenPerGrid);
            g.openNode(n.first);
          }

          // open from source node
          g.openNodeSink(m[from], 0);
          g.openNode(m[from]);

          std::cerr << " === Calculating topo penalties... ";
          begin = std::chrono::steady_clock::now();
          double addCFrom[8] = {0, 0, 0, 0, 0, 0, 0, 0};
          double addCTo[8] = {0, 0, 0, 0, 0, 0, 0, 0};
          double addCFromInv[8] = {0, 0, 0, 0, 0, 0, 0, 0};
          double addCToInv[8] = {0, 0, 0, 0, 0, 0, 0, 0};
          g.spacingPenalty(m[from], from, e, addCFrom);
          g.topoBlockPenalty(m[from], from, e, addCFrom);
          normalizeCostVector(addCFrom);
          g.outDegDeviationPenalty(m[from], from, e, addCFrom);

          g.addCostVector(m[from], addCFrom, addCFromInv);

          if (tos.size() == 1) {
            g.spacingPenalty(tos.begin()->first, to, e, addCTo);
            g.topoBlockPenalty(tos.begin()->first, to, e, addCTo);
            normalizeCostVector(addCTo);
            g.outDegDeviationPenalty(tos.begin()->first, to, e, addCTo);

            g.addCostVector(tos.begin()->first, addCTo, addCToInv);
          }

          end = std::chrono::steady_clock::now();
          if (verb)
            std::cerr << " done ("
                      << std::chrono::duration_cast<std::chrono::milliseconds>(
                             end - begin)
                             .count()
                      << "ms)" << std::endl;
          
          std::cerr << tos.size() << std::endl;


          int STOP_GEN = -1;

          if (verb) std::cerr << "Finding shortest path (w heur)... ";
          begin = std::chrono::steady_clock::now();

          std::list<GridEdge*> res2;
          GridNode* target2 = 0;
          // int iters2 = dijkstra.shortestPath(m[from], tos, &res2,
                                     // &target2);
          int iters = dijkstra.shortestPathAStar(m[from], tos, GridHeur(&g, m[from], tos), &res,
                                     &target, gen == STOP_GEN);
          end = std::chrono::steady_clock::now();
          if (verb)
            std::cerr << " done ("
                      << std::chrono::duration_cast<std::chrono::milliseconds>(
                             end - begin)
                             .count()
                      << "ms)" << std::endl;

          std::cerr << "HEUR: took " << iters << " iterations." << std::endl;
          // std::cerr << "BASELINE: took " << iters2 << " iterations." << std::endl;

          double costA = 0;
          double costB = 0;

          for (auto r : res) {
            costA += r->pl().cost();
          }

          for (auto r2 : res2) {
            costB += r2->pl().cost();
          }
          
          if (gen == STOP_GEN) {
            util::geo::output::GeoJsonOutput out;
            out.print(g);
            exit(0);
          }

          //assert(std::round(costA) == std::round(costB));

          g.removeCostVector(m[from], addCFromInv);
          if (tos.size() == 1) {
            g.removeCostVector(tos.begin()->first, addCToInv);
          }

          if (target == 0) {
            LOG(ERROR) << "Could not sort in node " << to << std::endl;
            continue;
          }

          m[to] = target;
          used.insert(target);

          // write everything to the result graph
          if (verb)
            std::cerr << "Writing found path as edge to result graph... ";
          begin = std::chrono::steady_clock::now();
          PolyLine pl;
          double cost = 0;
          size_t c = 0;
          for (auto f : res) {
            cost += c > 0 ? f->pl().cost() : 0;
            c++;

            assert(f->pl().getResEdges().size() == 0);
            f->pl().addResidentEdge(e);
            assert(g.getEdge(f->getTo(), f->getFrom())->pl().getResEdges().size() == 0);
            g.getEdge(f->getTo(), f->getFrom())->pl().addResidentEdge(e);

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

          end = std::chrono::steady_clock::now();
          double estimatedCost = g.heurCost(
              g.getNodeCoords(m[from]).first, g.getNodeCoords(m[from]).second,
              g.getNodeCoords(target).first, g.getNodeCoords(target).second);

          std::cerr << " done ("
                    << std::chrono::duration_cast<std::chrono::milliseconds>(
                           end - begin)
                           .count()
                    << "ms), path cost (excluding take off and touch down) "
                    << cost << ", estimated cost was " << estimatedCost
                    << std::endl;
          assert(std::round(cost) >= std::round(estimatedCost));

          totalCost += cost;

          if (res.size())
            pl << *res.back()->getTo()->pl().getParent()->pl().getGeom();

          // close the target node
          for (auto n : tos) {
            g.closeNodeSink(n.first);
          }
          g.closeNode(target);
          g.closeNodeSink(m[from]);
          g.closeNode(m[from]);

          if (verb) std::cerr << "Balance edge...";
          begin = std::chrono::steady_clock::now();
          for (auto f : res) {
            if (f->pl().isSecondary()) continue;
            g.balanceEdge(f->getFrom()->pl().getParent(),
                          f->getTo()->pl().getParent());
          }
          end = std::chrono::steady_clock::now();
          if (verb)
            std::cerr << " done ("
                      << std::chrono::duration_cast<std::chrono::milliseconds>(
                             end - begin)
                             .count()
                      << "ms)" << std::endl;

          if (reversed) pl.reverse();
          e->pl().setPolyLine(pl);
          e->pl().setGeneration(gen);

          gen++;
        }

        settled.insert(n);
      }
    }

    std::cerr << " === Total edge cost in graph is " << totalCost
              << " === " << std::endl;

    if (totalCost < bestScore) {
      bestScore = totalCost;
      bestGraph = TransitGraph();
      buildTransitGraph(&cg, &bestGraph);
    }
  }

  return bestGraph;
}
