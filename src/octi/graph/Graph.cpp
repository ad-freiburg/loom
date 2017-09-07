// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/graph/EdgePL.h"
#include "octi/graph/Graph.h"
#include "octi/graph/NodePL.h"
#include "util/graph/Edge.h"
#include "util/graph/Node.h"
#include "util/log/Log.h"

using namespace octi::graph;
using util::graph::Node;
using util::graph::Edge;

// _____________________________________________________________________________
Graph::Graph(std::istream* s) {
  _bbox = util::geo::minbox();

  json j;
  (*s) >> j;

  std::map<std::string, Node<NodePL, EdgePL>*> idMap;

  if (j["type"] == "FeatureCollection") {
    // first pass, nodes
    for (auto feature : j["features"]) {
      auto props = feature["properties"];
      auto geom = feature["geometry"];
      if (geom["type"] == "Point") {
        std::string id = props["id"];

        std::vector<double> coords = geom["coordinates"];

        Node<NodePL, EdgePL>* n = new Node<NodePL, EdgePL>(
            NodePL(util::geo::Point(coords[0], coords[1])));

        addNode(n);
        idMap[id] = n;
      }
    }

    // second pass, edges
    for (auto feature : j["features"]) {
      auto props = feature["properties"];
      auto geom = feature["geometry"];
      if (geom["type"] == "LineString") {
        if (props["lines"].is_null() || props["lines"].size() == 0) continue;
        std::string from = props["from"];
        std::string to = props["to"];

        std::vector<std::vector<double> > coords = geom["coordinates"];

        PolyLine pl;
        for (auto coord : coords) {
          double x = coord[0], y = coord[1];
          Point p(x, y);
          pl << p;
          expandBBox(p);
        }

        Node<NodePL, EdgePL>* fromN = idMap[from];
        Node<NodePL, EdgePL>* toN = idMap[to];

        if (!fromN) {
          LOG(ERROR) << "Node \"" << from << "\" not found." << std::endl;
          continue;
        }

        if (!toN) {
          LOG(ERROR) << "Node \"" << to << "\" not found." << std::endl;
          continue;
        }

        Edge<NodePL, EdgePL>* e = addEdge(fromN, toN, EdgePL(pl));
      }
    }
  }

  combineDeg2();
}

// _____________________________________________________________________________
void Graph::expandBBox(const Point& p) {
  bgeo::expand(_bbox, boost::geometry::make<bgeo::model::box<Point>>(
                          p.get<0>(), p.get<1>(), p.get<0>(), p.get<1>()));
}

// _____________________________________________________________________________
const util::geo::Box& Graph::getBBox() const {
  return _bbox;
}

// _____________________________________________________________________________
void Graph::combineDeg2() {
  auto nodes = *getNodes();
  for (auto n : nodes) {
    if (n->getAdjList().size() == 2) {
      Node<NodePL, EdgePL>* a = 0;
      Node<NodePL, EdgePL>* b = 0;

      for (auto e : n->getAdjList()) {
        if (!a) a = e->getOtherNode(n);
        else b = e->getOtherNode(n);
      }

      addEdge(a, b, EdgePL(PolyLine(*a->pl().getGeom(), *b->pl().getGeom())));
      deleteNode(n);
    }
  }
}
