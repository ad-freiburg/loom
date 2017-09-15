// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "octi/config/ConfigReader.h"
#include "octi/graph/CombEdgePL.h"
#include "octi/graph/CombNodePL.h"
#include "octi/graph/Graph.h"
#include "octi/gridgraph/GridGraph.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoJsonOutput.h"
#include "util/graph/Dijkstra.h"
#include "util/log/Log.h"

using std::string;
using namespace octi;

using util::graph::Node;
using octi::graph::NodePL;
using octi::graph::EdgePL;

typedef octi::graph::Graph TransitGraph;
typedef util::graph::Node<octi::graph::NodePL, octi::graph::EdgePL> TransitNode;
typedef util::graph::Edge<octi::graph::NodePL, octi::graph::EdgePL> TransitEdge;

typedef util::graph::Graph<octi::graph::CombNodePL, octi::graph::CombEdgePL>
    CombGraph;
typedef util::graph::Node<octi::graph::CombNodePL, octi::graph::CombEdgePL>
    CombNode;
typedef util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>
    CombEdge;

using octi::gridgraph::GridNode;
using octi::gridgraph::GridEdge;

// _____________________________________________________________________________
void buildCombGraph(TransitGraph* source, CombGraph* target) {
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
}

// _____________________________________________________________________________
void writeEdgeOrdering(CombGraph* target) {
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
          *refEdgeB->getOtherNode(e->getTo()->pl().getParent())->pl().getGeom());

      n->pl().addOrderedEdge(e, degA);
      e->getTo()->pl().addOrderedEdge(e, degB);
    }
  }
}

// _____________________________________________________________________________
void buildTransitGraph(CombGraph* source, TransitGraph* target) {
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
        target->addEdge(m[from], m[to], payload);

        i++;
      }
    }
  }
}

// _____________________________________________________________________________
void combineDeg2(CombGraph* g) {
  std::set<CombNode*> nodes = *g->getNodes();
  for (auto n : nodes) {
    if (n->getAdjList().size() == 2) {
      CombEdge* a = 0;
      CombEdge* b = 0;

      for (auto e : n->getAdjList()) {
        if (!a)
          a = g->getEdge(n, e->getOtherNode(n));
        else
          b = g->getEdge(n, e->getOtherNode(n));
      }

      auto pl = a->pl();

      // a is the referen
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
void removeEdgesShorterThan(TransitGraph* g, double d) {
start:
  for (auto n1 : *g->getNodes()) {
    for (auto e1 : n1->getAdjList()) {
      if (e1->pl().getPolyline().getLength() < d) {
        if (e1->getOtherNode(n1)->getAdjList().size() > 1 &&
            n1->getAdjList().size() > 1) {
          auto otherP = e1->getFrom()->pl().getGeom();
          auto n = g->mergeNodes(e1->getFrom(), e1->getTo());
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
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  config::Config cfg;

  config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  util::geo::output::GeoJsonOutput out;

  TransitGraph tg(&(std::cin));
  CombGraph cg;
  removeEdgesShorterThan(&tg, 150);
  buildCombGraph(&tg, &cg);
  combineDeg2(&cg);
  writeEdgeOrdering(&cg);

  std::map<CombNode*, GridNode*> m;

  std::map<CombNode*, TransitNode*> oldNew;

  std::set<GridNode*> used;

  double gridSize = 150;

  gridgraph::GridGraph g(tg.getBBox(), gridSize);

  g.balance();

  util::graph::Dijkstra dijkstra;

  // comparator for nodes, based on degree
  struct NodeCompare {
    bool operator()(CombNode* a, CombNode* b) {
      return a->getAdjList().size() < b->getAdjList().size();
    }
  };

  std::priority_queue<CombNode*, std::vector<CombNode*>, NodeCompare> pq;

  for (auto n : *cg.getNodes()) {
    pq.push(n);
  }

  std::set<CombEdge*> done;

  while (!pq.empty()) {
    auto n = pq.top();
    pq.pop();

    for (auto e : n->getAdjList()) {
      if (done.find(e) != done.end()) continue;
      done.insert(e);

      auto from = e->getFrom();
      auto to = e->getTo();

      bool reversed = false;

      if (from->getAdjList().size() < to->getAdjList().size()) {
        auto tmp = from;
        from = to;
        to = tmp;
        reversed = true;
      }

      if (m.find(from) == m.end()) {
        auto cands = g.getNearestCandidatesFor(*from->pl().getGeom(), 300);

        while (!cands.empty()) {
          if (used.find(cands.top().n) == used.end()) {
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
      std::set<GridNode*> tos;

      if (m.find(to) == m.end()) {
        double maxDis = to->getAdjList().size() == 1
                            ? util::geo::len(*e->pl().getGeom()) / 2
                            : 300;

        auto cands = g.getNearestCandidatesFor(*to->pl().getGeom(), maxDis);

        while (!cands.empty()) {
          if (used.find(cands.top().n) == used.end()) {
            tos.insert(cands.top().n);
          }
          cands.pop();
        }
      } else {
        tos.insert(m.find(to)->second);
      }

      std::list<GridEdge*> res;
      GridNode* target = 0;

      std::cerr << "Balancing... " << std::endl;
      g.balance();
      g.topoPenalty(m[from], from, e);
      if (tos.size() == 1) g.topoPenalty(*tos.begin(), to, e);
      std::cerr << "... done. " << std::endl;

      dijkstra.shortestPath(m[from], tos, &res, &target);

      if (target == 0) {
        LOG(ERROR) << "Could not sort in node " << to << std::endl;
      }

      m[to] = target;
      used.insert(target);

      // write everything to the result graph
      PolyLine pl;
      bool first = false;
      for (auto f : res) {
        f->pl().addResidentEdge(e);
        g.getEdge(f->getTo(), f->getFrom())->pl().addResidentEdge(e);
        if (!first) {
          pl << *f->getTo()->pl().getParent()->pl().getGeom();
        }
        pl << *f->getFrom()->pl().getParent()->pl().getGeom();
        first = true;
      }

      if (reversed) pl.reverse();
      e->pl().setPolyLine(pl);

    }
  }

  TransitGraph output;
  buildTransitGraph(&cg, &output);
  //out.print(g);
  out.print(output);

  return (0);
}
