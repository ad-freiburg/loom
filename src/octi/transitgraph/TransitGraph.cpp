// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "dot/Parser.h"
#include "json/json.hpp"
#include "octi/transitgraph/TransitEdgePL.h"
#include "octi/transitgraph/TransitGraph.h"
#include "octi/transitgraph/TransitNodePL.h"
#include "util/String.h"
#include "util/graph/Edge.h"
#include "util/graph/Node.h"
#include "util/log/Log.h"

#include "transitmap/style/LineStyle.h"

using namespace octi::transitgraph;
using octi::transitgraph::TransitNode;
using octi::transitgraph::TransitEdge;

// _____________________________________________________________________________
TransitGraph::TransitGraph() {}

// _____________________________________________________________________________
void TransitGraph::readFromDot(std::istream* s) {
  _bbox = util::geo::Box<double>();

  dot::parser::Parser dp(s);
  std::map<std::string, TransitNode*> idMap;

  size_t eid = 0;

  while (dp.has()) {
    auto ent = dp.get();

    if (ent.type == dot::parser::EMPTY) {
      continue;
    } else if (ent.type == dot::parser::NODE) {
      // only use nodes with a position atm
      if (ent.attrs.find("pos") == ent.attrs.end()) continue;
      std::string coords = ent.attrs["pos"];
      std::replace(coords.begin(), coords.end(), ',', ' ');

      std::stringstream ss;
      ss << coords;

      double x, y;

      ss >> x;
      ss >> y;

      TransitNode* n = 0;
      if (idMap.find(ent.ids.front()) != idMap.end()) {
        n = idMap[ent.ids.front()];
      }

      if (!n) {
        n = addNd(util::geo::Point<double>(x, y));
        idMap[ent.ids[0]] = n;
      }

      StationInfo i("", "");
      if (ent.attrs.find("station_id") != ent.attrs.end() ||
          ent.attrs.find("label") != ent.attrs.end()) {
        if (ent.attrs.find("station_id") != ent.attrs.end())
          i.id = ent.attrs["station_id"];
        if (ent.attrs.find("label") != ent.attrs.end())
          i.name = ent.attrs["label"];
        n->pl().addStop(i);
      }
    } else if (ent.type == dot::parser::EDGE) {
      eid++;
      std::string prevId = ent.ids.front();
      if (idMap.find(prevId) == idMap.end())
        idMap[prevId] = addNd(util::geo::Point<double>(0, 0));

      for (size_t i = 1; i < ent.ids.size(); ++i) {
        std::string curId = ent.ids[i];

        if (idMap.find(curId) == idMap.end())
          idMap[curId] = addNd(util::geo::Point<double>(0, 0));
        auto e = getEdg(idMap[prevId], idMap[curId]);

        if (!e) {
          PolyLine<double> pl;
          e = addEdg(idMap[curId], idMap[prevId], pl);
        }
        std::string id;
        if (ent.attrs.find("id") != ent.attrs.end()) {
          id = ent.attrs["id"];
        } else if (ent.attrs.find("label") != ent.attrs.end()) {
          id = ent.attrs["label"];
        } else {
          id = util::toString(eid);
        }

        const Route* r = getRoute(id);
        if (!r) {
          std::string label = ent.attrs.find("label") != ent.attrs.end()
                                  ? ""
                                  : ent.attrs["label"];
          std::string color = ent.attrs.find("label") != ent.attrs.end()
                                  ? ""
                                  : ent.attrs["label"];
          r = new Route(id, label, color);
          addRoute(r);
        }

        TransitNode* dir = 0;

        if (ent.graphType == dot::parser::DIGRAPH ||
            ent.graphType == dot::parser::STRICT_DIGRAPH) {
          dir = idMap[curId];
        }

        e->pl().addRoute(r, dir);
      }
    }
  }

  for (auto n : *getNds()) {
    for (auto e : n->getAdjListOut()) {
      PolyLine<double> pl;
      pl << *e->getFrom()->pl().getGeom();
      pl << *e->getTo()->pl().getGeom();

      e->pl().setPolyline(pl);
      expandBBox(pl.front());
      expandBBox(pl.back());
    }
  }

  buildGrids();
}

// _____________________________________________________________________________
void TransitGraph::readFromJson(std::istream* s) {
  _bbox = util::geo::Box<double>();

  nlohmann::json j;
  (*s) >> j;

  std::map<std::string, TransitNode*> idMap;

  if (j["type"] == "FeatureCollection") {
    // first pass, nodes
    for (auto feature : j["features"]) {
      auto props = feature["properties"];
      auto geom = feature["geometry"];
      if (geom["type"] == "Point") {
        std::string id = props["id"];

        std::vector<double> coords = geom["coordinates"];

        TransitNode* n = addNd(util::geo::DPoint(coords[0], coords[1]));

        StationInfo i("", "");
        if (!props["station_id"].is_null() ||
            !props["station_label"].is_null()) {
          if (!props["station_id"].is_null()) i.id = props["station_id"];
          if (!props["station_label"].is_null())
            i.name = props["station_label"];
          n->pl().addStop(i);
        }

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

        PolyLine<double> pl;
        for (auto coord : coords) {
          double x = coord[0], y = coord[1];
          Point<double> p(x, y);
          pl << p;
          expandBBox(p);
        }

        TransitNode* fromN = idMap[from];
        TransitNode* toN = idMap[to];

        if (!fromN) {
          LOG(ERROR) << "Node \"" << from << "\" not found." << std::endl;
          continue;
        }

        if (!toN) {
          LOG(ERROR) << "Node \"" << to << "\" not found." << std::endl;
          continue;
        }

        TransitEdge* e = addEdg(fromN, toN, pl);

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

          TransitNode* dir = 0;

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
void TransitGraph::buildGrids() {
  _nodeGrid = NodeGrid(200, 200, _bbox);
  _edgeGrid = EdgeGrid(200, 200, _bbox);

  for (auto n : *getNds()) {
    _nodeGrid.add(*n->pl().getGeom(), n);
    for (auto e : n->getAdjListOut()) {
      _edgeGrid.add(*e->pl().getGeom(), e);
    }
  }
}

// _____________________________________________________________________________
void TransitGraph::expandBBox(const Point<double>& p) {
  _bbox = util::geo::extendBox(p, _bbox);
}

// _____________________________________________________________________________
const util::geo::DBox& TransitGraph::getBBox() const { return _bbox; }

// _____________________________________________________________________________
void TransitGraph::topologizeIsects() {
  proced.clear();
  while (getNextIntersection().a) {
    auto i = getNextIntersection();
    auto x = addNd(i.bp.p);

    double pa = i.a->pl().getPolyline().projectOn(i.bp.p).totalPos;

    auto ba = addEdg(i.b->getFrom(), x,
                     i.b->pl().getPolyline().getSegment(0, i.bp.totalPos));
    auto bb = addEdg(x, i.b->getTo(),
                     i.b->pl().getPolyline().getSegment(i.bp.totalPos, 1));

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

    auto aa =
        addEdg(i.a->getFrom(), x, i.a->pl().getPolyline().getSegment(0, pa));
    auto ab =
        addEdg(x, i.a->getTo(), i.a->pl().getPolyline().getSegment(pa, 1));

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

    assert(getEdg(i.a->getFrom(), i.a->getTo()));
    assert(getEdg(i.b->getFrom(), i.b->getTo()));
    delEdg(i.a->getFrom(), i.a->getTo());
    delEdg(i.b->getFrom(), i.b->getTo());
  }
}

// _____________________________________________________________________________
ISect TransitGraph::getNextIntersection() {
  for (auto n1 : *getNds()) {
    for (auto e1 : n1->getAdjList()) {
      if (e1->getFrom() != n1) continue;
      if (proced.find(e1) != proced.end()) continue;

      std::set<TransitEdge*> neighbors;
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
void TransitGraph::addRoute(const Route* r) { _routes[r->getId()] = r; }

// _____________________________________________________________________________
const Route* TransitGraph::getRoute(const std::string& id) const {
  if (_routes.find(id) != _routes.end()) return _routes.find(id)->second;
  return 0;
}
