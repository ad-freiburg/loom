// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdint.h>
#include <ostream>
#include "gtfs2topo/config/GraphBuilderConfig.h"
#include "gtfs2topo/graph/BuildGraph.h"
#include "gtfs2topo/output/GeoJsonOutput.h"
#include "json/json.hpp"
#include "util/String.h"
#include "util/log/Log.h"
#include "gtfs2topo/graph/NodePL.h"
#include "gtfs2topo/graph/EdgePL.h"

using json = nlohmann::json;
using util::toString;
using namespace gtfs2topo;

// _____________________________________________________________________________
GeoJsonOutput::GeoJsonOutput(const config::Config* cfg) : _cfg(cfg) {}

// _____________________________________________________________________________
void GeoJsonOutput::print(const graph::BuildGraph& outG) {
  json geoj;
  geoj["type"] = "FeatureCollection";
  geoj["features"] = json::array();

  // first pass, nodes
  for (graph::Node* n : outG.getNodes()) {
    json feature;
    feature["type"] = "Feature";

    feature["geometry"]["type"] = "Point";
    std::vector<double> coords;
    coords.push_back(n->pl().getPos().get<0>());
    coords.push_back(n->pl().getPos().get<1>());
    feature["geometry"]["coordinates"] = coords;

    feature["properties"] = json::object();
    feature["properties"]["id"] = toString(n);

    if (n->pl().getStops().size() > 0) {
      feature["properties"]["station_id"] = (*n->pl().getStops().begin())->getId();
      feature["properties"]["station_label"] =
          (*n->pl().getStops().begin())->getName();
    }

    auto arr = json::array();

    for (graph::Edge* e : n->getAdjList()) {
      if (e->pl().getEdgeTripGeoms()->size() == 0) continue;
      for (auto r : *e->pl().getEdgeTripGeoms()->front().getTripsUnordered()) {
        for (graph::Edge* f : n->getAdjList()) {
          if (e == f) continue;
          if (f->pl().getEdgeTripGeoms()->size() == 0) continue;
          for (auto rr : *f->pl().getEdgeTripGeoms()->front().getTripsUnordered()) {
            if (!_cfg->ignoreDirections && r.route == rr.route &&
                (r.direction == 0 || rr.direction == 0 ||
                 (r.direction == n && rr.direction != n) ||
                 (r.direction != n && rr.direction == n)) &&
                !n->pl().isConnOccuring(r.route, e, f)) {
              auto obj = json::object();
              obj["route"] = toString(r.route);
              obj["edge1_node"] =
                  toString(e->getFrom() == n ? e->getTo() : e->getFrom());
              obj["edge2_node"] =
                  toString(f->getFrom() == n ? f->getTo() : f->getFrom());
              arr.push_back(obj);
            }
          }
        }
      }
    }

    if (arr.size()) feature["properties"]["excluded_line_conns"] = arr;
    geoj["features"].push_back(feature);
  }

  // second pass, edges
  for (graph::Node* n : outG.getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      if (e->pl().getEdgeTripGeoms()->size() > 0) {
        json feature;
        feature["type"] = "Feature";
        feature["properties"]["from"] = toString(e->getFrom());
        feature["properties"]["to"] = toString(e->getTo());

        feature["geometry"]["type"] = "LineString";
        feature["geometry"]["coordinates"] = json::array();

        for (auto p : e->pl().getEdgeTripGeoms()->front().getGeom().getLine()) {
          std::vector<double> coords;
          coords.push_back(p.get<0>());
          coords.push_back(p.get<1>());
          feature["geometry"]["coordinates"].push_back(coords);
        }

        feature["properties"]["lines"] = json::array();

        for (auto r : *e->pl().getEdgeTripGeoms()->front().getTripsUnordered()) {
          json route = json::object();
          route["id"] = toString(r.route);
          route["label"] = r.route->getShortName();
          route["color"] = r.route->getColorString();

          if (!_cfg->ignoreDirections && r.direction != 0) {
            route["direction"] = toString(r.direction);
          }

          feature["properties"]["lines"].push_back(route);
        }

        geoj["features"].push_back(feature);
      }
    }
  }

  std::cout << geoj.dump(2);
}
