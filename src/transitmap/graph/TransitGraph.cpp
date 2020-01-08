// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <set>
#include <string>
#include "json/json.hpp"
#include "shared/linegraph/Route.h"
#include "transitmap/graph/Edge.h"
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/TransitGraph.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"

using util::geo::Point;
using util::geo::Box;
using util::geo::dist;
using transitmapper::graph::TransitGraph;
using transitmapper::graph::Node;
using transitmapper::graph::Edge;
using shared::linegraph::Route;
using transitmapper::graph::OrderingConfig;

// _____________________________________________________________________________
TransitGraph::TransitGraph() { _bbox = util::geo::minbox<double>(); }

// _____________________________________________________________________________
TransitGraph::~TransitGraph() {
  for (auto n : _nodes) {
    delete n;
  }
}

// _____________________________________________________________________________
const OrderingConfig& TransitGraph::getConfig() const { return _config; }

// _____________________________________________________________________________
void TransitGraph::setConfig(const OrderingConfig& c) { _config = c; }

// _____________________________________________________________________________
void TransitGraph::addNd(Node* n) {
  _nodes.insert(n);
  // expand the bounding box to hold this new node
  expandBBox(n->getPos());
}

// _____________________________________________________________________________
void TransitGraph::expandBBox(const DPoint& p) {
  _bbox = util::geo::extendBox(p, _bbox);
}

// _____________________________________________________________________________
Node* TransitGraph::getNodeById(const std::string& id) const {
  for (auto n : _nodes) {
    if (n->getId() == id) return n;
  }

  return 0;
}

// _____________________________________________________________________________
Edge* TransitGraph::addEdg(Node* from, Node* to, PolyLine<double> pl, double w,
                           double s) {
  if (from == to) return 0;
  Edge* e = getEdg(from, to);
  if (!e) {
    e = new Edge(from, to, pl, w, s);
    from->addEdg(e);
    to->addEdg(e);
    _bbox =
        util::geo::extendBox(util::geo::getBoundingBox(pl.getLine()), _bbox);
  }
  return e;
}

// _____________________________________________________________________________
void TransitGraph::delEdg(Node* from, Node* to) {
  Edge* toDel = getEdg(from, to);
  if (!toDel) return;

  from->removeEdge(toDel);
  to->removeEdge(toDel);

  assert(!getEdg(from, to));

  delete toDel;
}

// _____________________________________________________________________________
void TransitGraph::addRoute(const Route* r) {
  if (!getRoute(r->getId())) {
    _routes[r->getId()] = r;
  }
}

// _____________________________________________________________________________
const Route* TransitGraph::getRoute(const std::string& id) const {
  auto f = _routes.find(id);
  if (f == _routes.end()) return 0;

  return f->second;
}

// _____________________________________________________________________________
Edge* TransitGraph::getEdg(Node* from, Node* to) {
  for (auto e : from->getAdjListOut()) {
    if (e->getTo() == to) return e;
  }

  // also search in the opposite direction, we are handling an undirected
  // graph here
  for (auto e : from->getAdjListIn()) {
    if (e->getFrom() == to) return e;
  }

  return 0;
}

// _____________________________________________________________________________
const std::set<Node*>& TransitGraph::getNds() const { return _nodes; }

// _____________________________________________________________________________
std::set<Node*>* TransitGraph::getNds() { return &_nodes; }

// _____________________________________________________________________________
const DBox& TransitGraph::getBBox() const { return _bbox; }

// _____________________________________________________________________________
size_t TransitGraph::getNumNodes() const {
  return getNumNodes(true) + getNumNodes(false);
}

// _____________________________________________________________________________
size_t TransitGraph::getNumRoutes() const { return _routes.size(); }

// _____________________________________________________________________________
size_t TransitGraph::getMaxCardinality() const {
  size_t ret = 0;
  for (auto n : getNds()) {
    for (auto e : n->getAdjListOut()) {
      if (e->getCardinality() > ret) ret = e->getCardinality();
    }
  }
  return ret;
}

// _____________________________________________________________________________
size_t TransitGraph::maxDeg() const {
  size_t ret = 0;
  for (auto n : getNds()) {
    if (n->getAdjListOut().size() + n->getAdjListIn().size() > ret) {
      ret = n->getAdjListOut().size() + n->getAdjListIn().size();
    }
  }
  return ret;
}

// _____________________________________________________________________________
size_t TransitGraph::getNumEdges() const {
  size_t ret = 0;

  for (auto n : getNds()) {
    ret += n->getAdjListOut().size();
  }

  return ret;
}

// _____________________________________________________________________________
size_t TransitGraph::getNumNodes(bool topo) const {
  size_t ret = 0;
  for (auto n : _nodes) {
    if (n->getAdjListIn().size() + n->getAdjListOut().size() == 0) continue;
    if ((n->getStops().size() == 0) ^ !topo) ret++;
  }

  return ret;
}

// _____________________________________________________________________________
bool TransitGraph::readFromJson(std::istream* s) {
  nlohmann::json j;
  (*s) >> j;

  if (j["type"] == "FeatureCollection") {
    // first pass, nodes
    for (auto feature : j["features"]) {
      auto props = feature["properties"];
      auto geom = feature["geometry"];
      if (geom["type"] == "Point") {
        std::string id = util::toString(props["id"]);

        // check if node already exists
        if (getNodeById(id)) continue;

        std::vector<double> coords = geom["coordinates"];

        Node* n = new Node(id, coords[0], coords[1]);

        shared::linegraph::Station i("", "", n->getPos());
        if (!props["station_id"].is_null() ||
            !props["station_label"].is_null()) {
          if (!props["station_id"].is_null())
            i.id = util::toString(props["station_id"]);
          if (!props["station_label"].is_null())
            i.name = util::toString(props["station_label"]);
          n->addStop(i);
        }

        addNd(n);
      }
    }

    // second pass, edges
    for (auto feature : j["features"]) {
      auto props = feature["properties"];
      auto geom = feature["geometry"];
      if (geom["type"] == "LineString") {
        if (props["lines"].is_null() || props["lines"].size() == 0) continue;
        std::string from = util::toString(props["from"]);
        std::string to = util::toString(props["to"]);

        std::vector<std::vector<double>> coords = geom["coordinates"];

        PolyLine<double> pl;
        for (auto coord : coords) {
          double x = coord[0], y = coord[1];
          DPoint p(x, y);
          pl << p;
          expandBBox(p);
        }

        // TODO
        // pl.applyChaikinSmooth(_cfg->inputSmoothing);

        Node* fromN = getNodeById(from);
        Node* toN = getNodeById(to);

        if (!fromN) {
          continue;
        }

        if (!toN) {
          continue;
        }

        // TODO
        Edge* e = addEdg(fromN, toN, pl, 0, 0);

        assert(e);
        assert(getNodeById(from));
        assert(getNodeById(to));

        for (auto route : props["lines"]) {
          std::string id;
          if (!route["id"].is_null()) {
            id = util::toString(route["id"]);
          } else if (!route["label"].is_null()) {
            id = util::toString(route["label"]);
          } else if (!route["color"].is_null()) {
            id = route["color"];
          } else
            continue;

          const Route* r = getRoute(id);
          if (!r) {
            std::string label = route["label"].is_null() ? "" : route["label"];
            std::string color = route["color"];
            r = new Route(id, label, color);
            addRoute(r);
          }

          Node* dir = 0;

          if (!route["direction"].is_null()) {
            dir = getNodeById(util::toString(route["direction"]));
          }

          if (!route["style"].is_null()) {
            style::LineStyle ls;
            auto style = route["style"];
            std::string dashArray;
            if (!style["dash-array"].is_null()) {
              dashArray = style["dash-array"];
            }

            if (!style["css"].is_null()) {
              ls.setCss(style["css"]);
            }

            ls.setDashArray(dashArray);

            e->addRoute(r, dir, ls);
          } else {
            e->addRoute(r, dir);
          }
        }
      }
    }

    // third pass, exceptions (TODO: do this in the first part, store in some
    // data strcuture,
    //  add here!)
    for (auto feature : j["features"]) {
      auto props = feature["properties"];
      auto geom = feature["geometry"];
      if (geom["type"] == "Point") {
        std::string id = util::toString(props["id"]);

        Node* n = getNodeById(id);

        if (!n) continue;

        if (!props["excluded_line_conns"].is_null()) {
          for (auto excl : props["excluded_line_conns"]) {
            std::string rid = util::toString(excl["route"]);
            std::string nid1 = util::toString(excl["edge1_node"]);
            std::string nid2 = util::toString(excl["edge2_node"]);

            const Route* r = getRoute(rid);

            if (!r) {
              continue;
            }

            Node* n1 = getNodeById(nid1);
            Node* n2 = getNodeById(nid2);

            if (!n1) {
              continue;
            }

            if (!n2) {
              continue;
            }

            Edge* a = n->getEdg(n1);
            Edge* b = n->getEdg(n2);

            if (!a) {
              continue;
            }

            if (!b) {
              continue;
            }

            n->addRouteConnException(r, a, b);
          }
        }
      }
    }

  } else {
    return false;
  }

  return true;
}
