// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdint.h>
#include <ostream>
#include "log/Log.h"
#include "transitmap/config/TransitMapConfig.h"
#include "./JsonOutput.h"
#include "transitmap/geo/PolyLine.h"
#include "./../graph/Graph.h"
#include "transitmap/graph/Edge.h"
#include "json/json.hpp"

using json = nlohmann::json;
using namespace skeletonbuilder;

// _____________________________________________________________________________
JsonOutput::JsonOutput(const config::Config* cfg)
: _cfg(cfg) {

}

// _____________________________________________________________________________
void JsonOutput::print(const graph::Graph& outG) {
  json geoj;
  geoj["type"] = "FeatureCollection";
  geoj["features"] = json::array();

  // first pass, nodes
  for (graph::Node* n : outG.getNodes()) {
    json feature;
    feature["type"] = "Feature";

    feature["geometry"]["type"] = "Point";
    std::vector<double> coords;
    coords.push_back(n->getPos().get<0>());
    coords.push_back(n->getPos().get<1>());
    feature["geometry"]["coordinates"] = coords;

    feature["properties"] = json::object();
    feature["properties"]["id"] = boost::lexical_cast<std::string>(n);

    if (n->getStops().size() > 0) {
      feature["properties"]["station_id"] = (*n->getStops().begin())->getId();
      feature["properties"]["station_label"] = (*n->getStops().begin())->getName();
    }

    auto arr = json::array();

    for (graph::Edge* e : n->getAdjList()) {
      for (auto r : *e->getEdgeTripGeoms()->front().getTripsUnordered()) {
        for (graph::Edge* f : n->getAdjList()) {
          if (e == f) continue;
          for (auto rr : *f->getEdgeTripGeoms()->front().getTripsUnordered()) {
            if (r.route == rr.route &&
              (r.direction == 0 || rr.direction == 0 ||
                (r.direction == n && rr.direction != n) ||
                (r.direction != n && rr.direction == n))
              // && false) {
              && !n->isConnOccuring(r.route, e, f)) {
              auto obj = json::object();
              obj["route"] = boost::lexical_cast<std::string>(r.route);
              obj["edge1_node"] = boost::lexical_cast<std::string>(e->getFrom() == n ? e->getTo() : e->getFrom());
              obj["edge2_node"] = boost::lexical_cast<std::string>(f->getFrom() == n ? f->getTo() : f->getFrom());
              arr.push_back(obj);
            } else if (r.route == rr.route && n->getStops().size() == 0) {
              std::cerr << "-: " << r.direction << " " << rr.direction << std::endl;
              std::cerr << "A: " << (r.direction == n && rr.direction != n) << std::endl;
              std::cerr << "B: " << (r.direction != n && rr.direction == n) << std::endl;
              std::cerr << "C: " << (n->isConnOccuring(r.route, e, f)) << std::endl;
            }
          }
        }
      }
    }

    if (arr.size()) feature["properties"]["excluded_line_connections"] = arr;

    geoj["features"].push_back(feature);
  }

  // second pass, edges
  for (graph::Node* n : outG.getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      if (e->getEdgeTripGeoms()->size() > 0) {
        json feature;
        feature["type"] = "Feature";
        feature["properties"]["from"] = boost::lexical_cast<std::string>(e->getFrom());
        feature["properties"]["to"] = boost::lexical_cast<std::string>(e->getTo());

        feature["geometry"]["type"] = "LineString";
        feature["geometry"]["coordinates"] = json::array();

        for (auto p : e->getEdgeTripGeoms()->front().getGeom().getLine()) {
          std::vector<double> coords;
          coords.push_back(p.get<0>());
          coords.push_back(p.get<1>());
          feature["geometry"]["coordinates"].push_back(coords);
        }

        feature["properties"]["lines"] = json::array();

        for (auto r : *e->getEdgeTripGeoms()->front().getTripsUnordered()) {
          json route = json::object();
          route["id"] = boost::lexical_cast<std::string>(r.route);
          route["label"] = r.route->getShortName();
          route["color"] = r.route->getColorString();

          if (r.direction != 0) {
            route["direction"] = boost::lexical_cast<std::string>(r.direction);
          }

          feature["properties"]["lines"].push_back(route);
        }

        geoj["features"].push_back(feature);
      }
    }
  }

  std::cout << geoj.dump(2);
}