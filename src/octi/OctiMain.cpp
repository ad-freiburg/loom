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
#include "octi/graph/Graph.h"
#include "octi/gridgraph/GridGraph.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoJsonOutput.h"
#include "util/graph/Dijkstra.h"
#include "util/log/Log.h"

using std::string;
using namespace octi;

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

  /*
   *   gridgraph::GridGraph g(util::geo::Box(util::geo::Point(0, 0),
   * util::geo::Point(100, 100)), 10);
   *   g.balance();
   *
   *   util::graph::Dijkstra dijkstra;
   *
   *   std::list<util::graph::Edge<gridgraph::NodePL, gridgraph::EdgePL>*> res;
   *
   *   dijkstra.shortestPath(g.getNode(0, 0), g.getNode(5, 8), &res);
   *
   *   for (auto e : res) {
   *     e->pl().addRoute("A");
   *   }
   *
   *   res.clear();
   *   g.balance();
   *
   *   dijkstra.shortestPath(g.getNode(0, 8), g.getNode(9, 0), &res);
   *
   *   for (auto e : res) {
   *     e->pl().addRoute("B");
   *   }
   *
   *   res.clear();
   *   g.balance();
   *
   *   out.print(g);
   */

  graph::Graph tg(&(std::cin));

  std::map<util::graph::Node<octi::graph::NodePL, octi::graph::EdgePL>*,
           util::graph::Node<octi::gridgraph::NodePL, octi::gridgraph::EdgePL>*>
      m;

  std::set<util::graph::Node<octi::gridgraph::NodePL, octi::gridgraph::EdgePL>*>
      used;

  double gridSize = 200;

  gridgraph::GridGraph g(tg.getBBox(), gridSize);

  g.balance();

  util::graph::Dijkstra dijkstra;

  // comparator for nodes, based on degree
  struct NodeCompare {
    bool operator()(
        util::graph::Node<octi::graph::NodePL, octi::graph::EdgePL>* a,
        util::graph::Node<octi::graph::NodePL, octi::graph::EdgePL>* b) {
      return a->getAdjList().size() < b->getAdjList().size();
    }
  };

  std::priority_queue<
      util::graph::Node<octi::graph::NodePL, octi::graph::EdgePL>*,
      std::vector<util::graph::Node<octi::graph::NodePL, octi::graph::EdgePL>*>,
      NodeCompare>
      pq;

  for (auto n : *tg.getNodes()) {
    pq.push(n);
  }

  std::set<
      util::graph::Edge<octi::graph::NodePL, octi::graph::EdgePL>*> done;

  while (!pq.empty()) {
    auto n = pq.top();
    pq.pop();

    for (auto e : n->getAdjList()) {
      if (done.find(e) != done.end()) continue;
      done.insert(e);

      auto from = e->getFrom();
      auto to = e->getTo();

      if (from->getAdjList().size() < to->getAdjList().size()) {
        auto tmp = from;
        from = to;
        to = tmp;
      }

      if (m.find(from) == m.end()) {
        auto cands =
            g.getNearestCandidatesFor(*from->pl().getGeom(), 300);

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
      std::set<
          util::graph::Node<octi::gridgraph::NodePL, octi::gridgraph::EdgePL>*>
          tos;

      if (m.find(to) == m.end()) {
        double maxDis = to->getAdjList().size() == 1 ? 400 : 300;

        auto cands =
            g.getNearestCandidatesFor(*to->pl().getGeom(), maxDis);

        while (!cands.empty()) {
          if (used.find(cands.top().n) == used.end()) {
            tos.insert(cands.top().n);
          }
          cands.pop();
        }
      } else {
        tos.insert(m.find(to)->second);
      }

      std::list<util::graph::Edge<gridgraph::NodePL, gridgraph::EdgePL>*> res;
      util::graph::Node<octi::gridgraph::NodePL, octi::gridgraph::EdgePL>*
          target = 0;

      dijkstra.shortestPath(m[from], tos, &res, &target);

      if (target == 0) {
        LOG(ERROR) << "Could not sort in node " << to << std::endl;
      }

      m[to] = target;
      used.insert(target);

      for (auto e : res) {
        e->pl().addRoute("B");
      }

      g.balance();
    }
  }

  out.print(g);
  // out.print(tg);

  return (0);
}
