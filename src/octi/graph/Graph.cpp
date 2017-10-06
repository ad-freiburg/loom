// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/graph/EdgePL.h"
#include "octi/graph/Graph.h"
#include "octi/graph/NodePL.h"
#include "util/graph/Edge.h"
#include "util/graph/Node.h"
#include "util/log/Log.h"

#include "transitmap/style/LineStyle.h"
#include "transitmap/style/LineStyle.h"

using namespace octi::graph;
using util::graph::Node;
using util::graph::Edge;

// _____________________________________________________________________________
Graph::Graph() { _bbox = util::geo::minbox(); }

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

        StationInfo i("", "");
        if (!props["station_id"].is_null() ||
            !props["station_label"].is_null()) {
          if (!props["station_id"].is_null()) i.id = props["station_id"];
          if (!props["station_label"].is_null())
            i.name = props["station_label"];
          n->pl().addStop(i);
        }

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

        std::vector<std::vector<double>> coords = geom["coordinates"];

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

        for (auto route : props["lines"]) {
          std::string id;
          if (!route["id"].is_null()) {
            id = route["id"];
          } else if (!route["label"].is_null()) {
            id = route["label"];
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

          Node<NodePL, EdgePL>* dir = 0;

          if (!route["direction"].is_null()) {
            dir = idMap[route["direction"]];
          }

          if (!route["style"].is_null()) {
            transitmapper::style::LineStyle ls;
            auto style = route["style"];
            std::string dashArray;
            if (!style["dash-array"].is_null()) {
              dashArray = style["dash-array"];
            }

            if (!style["css"].is_null()) {
              ls.setCss(style["css"]);
            }

            ls.setDashArray(dashArray);

            e->pl().addRoute(r, dir, ls);
          } else {
            e->pl().addRoute(r, dir);
          }
        }
      }
    }
  }

  buildGrids();
}

// _____________________________________________________________________________
void Graph::buildGrids() {
  _nodeGrid = NodeGrid(200, 200, _bbox);
  _edgeGrid = EdgeGrid(200, 200, _bbox);

  for (auto n : *getNodes()) {
    _nodeGrid.add(*n->pl().getGeom(), n);
    for (auto e : n->getAdjListOut()) {
      _edgeGrid.add(*e->pl().getGeom(), e);
    }
  }
}

// _____________________________________________________________________________
void Graph::expandBBox(const Point& p) {
  bgeo::expand(_bbox, boost::geometry::make<bgeo::model::box<Point>>(
                          p.get<0>(), p.get<1>(), p.get<0>(), p.get<1>()));
}

// _____________________________________________________________________________
const util::geo::Box& Graph::getBBox() const { return _bbox; }


// _____________________________________________________________________________
void Graph::topologizeIsects() {
  proced.clear();
  while (getNextIntersection().a) {
    auto i = getNextIntersection();
    auto x = new util::graph::Node<NodePL, EdgePL>(NodePL(i.bp.p));
    addNode(x);

    double pa = i.a->pl().getPolyline().projectOn(i.bp.p).totalPos;

    auto ba = addEdge(i.b->getFrom(), x,
            EdgePL(i.b->pl().getPolyline().getSegment(0, i.bp.totalPos)));
    auto bb = addEdge(x, i.b->getTo(),
            EdgePL(i.b->pl().getPolyline().getSegment(i.bp.totalPos, 1)));

    _edgeGrid.add(*ba->pl().getGeom(), ba);
    _edgeGrid.add(*bb->pl().getGeom(), bb);

    for (auto r : i.b->pl().getRoutes()) {
      if (r.direction == i.b->getFrom()) {
        ba->pl().addRoute(r.route, i.b->getFrom());
        bb->pl().addRoute(r.route, x);
      } else {
        ba->pl().addRoute(r.route, x);
        bb->pl().addRoute(r.route, i.b->getTo());
      }
    }

    auto aa = addEdge(i.a->getFrom(), x,
            EdgePL(i.a->pl().getPolyline().getSegment(0, pa)));
    auto ab = addEdge(x, i.a->getTo(), EdgePL(i.a->pl().getPolyline().getSegment(pa, 1)));

    _edgeGrid.add(*aa->pl().getGeom(), aa);
    _edgeGrid.add(*ab->pl().getGeom(), ab);

    for (auto r : i.a->pl().getRoutes()) {
      if (r.direction == i.b->getFrom()) {
        aa->pl().addRoute(r.route, i.a->getFrom());
        ab->pl().addRoute(r.route, x);
      } else {
        aa->pl().addRoute(r.route, x);
        ab->pl().addRoute(r.route, i.a->getTo());
      }
    }

    _edgeGrid.remove(i.a);
    _edgeGrid.remove(i.b);
    deleteEdge(i.a->getFrom(), i.a->getTo());
    deleteEdge(i.b->getFrom(), i.b->getTo());
  }
}

// _____________________________________________________________________________
ISect Graph::getNextIntersection() {
  for (auto n1 : *getNodes()) {
    for (auto e1 : n1->getAdjListOut()) {
      if (proced.find(e1) != proced.end()) continue;

      std::set<util::graph::Edge<NodePL, EdgePL>*> neighbors;
      _edgeGrid.getNeighbors(e1, 0, &neighbors);

      for (auto e2 : neighbors) {
        if (proced.find(e2) != proced.end()) continue;
        if (e1 != e2) {
          auto is =
              e1->pl().getPolyline().getIntersections(e2->pl().getPolyline());
          if (is.size()) {
            ISect ret;
            ret.a = e1;
            ret.b = e2;
            ret.bp = *is.begin();
            if (ret.bp.totalPos > 0.001 && 1 - ret.bp.totalPos > 0.001) {
              return ret;
            }
          }
        }
      }
      proced.insert(e1);
    }
  }

  ISect ret;
  ret.a = 0;
  ret.b = 0;
  return ret;
}


// _____________________________________________________________________________
void Graph::addRoute(const Route* r) { _routes[r->getId()] = r; }

// _____________________________________________________________________________
const Route* Graph::getRoute(const std::string& id) const {
  if (_routes.find(id) != _routes.end()) return _routes.find(id)->second;
  return 0;
}
