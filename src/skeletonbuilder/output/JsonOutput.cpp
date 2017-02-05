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

  for (graph::Node* n : outG.getNodes()) {
    json feature;
    feature["type"] = "Feature";
    feature["propertiers"]["id"] = boost::lexical_cast<std::string>(n);

    feature["geometry"]["type"] = "Point";
    std::vector<double> coords;
    coords.push_back(n->getPos().get<0>());
    coords.push_back(n->getPos().get<0>());
    feature["geometry"]["coordinates"] = coords;

    geoj["features"].push_back(feature);
  }

  for (graph::Node* n : outG.getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      if (e->getEdgeTripGeoms()->size() > 0) {
        json feature;
        feature["type"] = "Feature";
        feature["propertiers"]["from"] = boost::lexical_cast<std::string>(e->getFrom());
        feature["propertiers"]["to"] = boost::lexical_cast<std::string>(e->getTo());

        feature["geometry"]["type"] = "LineString";
        feature["geometry"]["coordinates"] = json::array();

        for (auto p : e->getEdgeTripGeoms()->front().getGeom().getLine()) {
          std::vector<double> coords;
          coords.push_back(p.get<0>());
          coords.push_back(p.get<1>());
          feature["geometry"]["coordinates"].push_back(coords);
        }

        geoj["features"].push_back(feature);
      }
    }
  }

  std::cout << geoj;
}