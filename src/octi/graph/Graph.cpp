// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/graph/Graph.h"
#include "util/graph/Node.h"

using namespace octi::graph;

// _____________________________________________________________________________
Graph::Graph(std::istream* s) {
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

        Node<NodePL, EdgePL>* n = new Node<NodePL, EdgePL>(new NodePL(coords[0], coords[1]));

        addNode(n);
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
          g->expandBBox(p);
        }

        Node* fromN = g->getNodeById(from);
        Node* toN = g->getNodeById(to);

        if (!fromN) {
          LOG(WARN) << "Node \"" << from << "\" not found." << std::endl;
          continue;
        }

        if (!toN) {
          LOG(WARN) << "Node \"" << to << "\" not found." << std::endl;
          continue;
        }

        if (dist(fromN->getPos(), pl.back()) <
            dist(fromN->getPos(), pl.front())) {
          LOG(WARN) << "Geometry for edge from " << fromN->getId() << " to "
                    << toN->getId() << " seems "
                    << " to have the wrong orientation! This may lead to "
                    << " strange results." << std::endl;
        }

        Edge* e =
            g->addEdge(fromN, toN, pl, _cfg->lineWidth, _cfg->lineSpacing);

      }
    }
}

