// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "dot/Parser.h"
#include "json/json.hpp"
#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineGraph.h"
#include "shared/linegraph/LineNodePL.h"
#include "util/String.h"
#include "util/graph/Edge.h"
#include "util/graph/Node.h"
#include "util/log/Log.h"

#include "shared/style/LineStyle.h"

using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using shared::linegraph::LineEdge;
using shared::linegraph::Line;
using shared::linegraph::Partner;
using shared::linegraph::LineOcc;
using shared::linegraph::ISect;
using util::geo::Point;
using util::geo::DPoint;

// _____________________________________________________________________________
LineGraph::LineGraph() {}

// _____________________________________________________________________________
void LineGraph::readFromDot(std::istream* s, double smooth) {
  _bbox = util::geo::Box<double>();

  dot::parser::Parser dp(s);
  std::map<std::string, LineNode*> idMap;

  size_t eid = 0;

  while (dp.has()) {
    auto ent = dp.get();

    if (ent.type == dot::parser::EMPTY) {
      continue;
    } else if (ent.type == dot::parser::NODE) {
      // only use nodes with a position
      if (ent.attrs.find("pos") == ent.attrs.end()) continue;
      std::string coords = ent.attrs["pos"];
      std::replace(coords.begin(), coords.end(), ',', ' ');

      std::stringstream ss;
      ss << coords;

      double x, y;

      ss >> x;
      ss >> y;

      LineNode* n = 0;
      if (idMap.find(ent.ids.front()) != idMap.end()) {
        n = idMap[ent.ids.front()];
      }

      if (!n) {
        n = addNd(util::geo::Point<double>(x, y));
        idMap[ent.ids[0]] = n;
      }

      expandBBox(*n->pl().getGeom());

      Station i("", "", *n->pl().getGeom());
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

        const Line* r = getLine(id);
        if (!r) {
          std::string label = ent.attrs.find("label") != ent.attrs.end()
                                  ? ""
                                  : ent.attrs["label"];
          std::string color = ent.attrs.find("label") != ent.attrs.end()
                                  ? ""
                                  : ent.attrs["label"];
          r = new Line(id, label, color);
          addLine(r);
        }

        LineNode* dir = 0;

        if (ent.graphType == dot::parser::DIGRAPH ||
            ent.graphType == dot::parser::STRICT_DIGRAPH) {
          dir = idMap[curId];
        }

        e->pl().addLine(r, dir);
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
void LineGraph::readFromJson(std::istream* s, double smooth) {
  _bbox = util::geo::Box<double>();

  nlohmann::json j;
  (*s) >> j;

  std::map<std::string, LineNode*> idMap;

  if (j["type"] == "FeatureCollection") {
    // first pass, nodes
    for (auto feature : j["features"]) {
      auto props = feature["properties"];
      auto geom = feature["geometry"];
      if (geom["type"] == "Point") {
        std::string id = util::toString(props["id"]);

        std::vector<double> coords = geom["coordinates"];

        LineNode* n = addNd(util::geo::DPoint(coords[0], coords[1]));
        expandBBox(*n->pl().getGeom());

        Station i("", "", *n->pl().getGeom());
        if (!props["station_id"].is_null() ||
            !props["station_label"].is_null()) {
          if (!props["station_id"].is_null())
            i.id = util::toString(props["station_id"]);
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
        std::string from = util::toString(props["from"]);
        std::string to = util::toString(props["to"]);

        std::vector<std::vector<double>> coords = geom["coordinates"];

        PolyLine<double> pl;
        for (auto coord : coords) {
          double x = coord[0], y = coord[1];
          Point<double> p(x, y);
          pl << p;
          expandBBox(p);
        }

        pl.applyChaikinSmooth(smooth);

        LineNode* fromN = idMap[from];
        LineNode* toN = idMap[to];

        if (!fromN) {
          LOG(ERROR) << "Node \"" << from << "\" not found." << std::endl;
          continue;
        }

        if (!toN) {
          LOG(ERROR) << "Node \"" << to << "\" not found." << std::endl;
          continue;
        }

        LineEdge* e = addEdg(fromN, toN, pl);

        for (auto line : props["lines"]) {
          std::string id;
          if (!line["id"].is_null()) {
            id = util::toString(line["id"]);
          } else if (!line["label"].is_null()) {
            id = util::toString(line["label"]);
          } else if (!line["color"].is_null()) {
            id = line["color"];
          } else
            continue;

          const Line* l = getLine(id);
          if (!l) {
            std::string label = line["label"].is_null() ? "" : line["label"];
            std::string color = line["color"];
            l = new Line(id, label, color);
            addLine(l);
          }

          LineNode* dir = 0;

          if (!line["direction"].is_null()) {
            dir = idMap[util::toString(line["direction"])];
          }

          if (!line["style"].is_null()) {
            shared::style::LineStyle ls;
            auto style = line["style"];
            std::string dashArray;
            if (!style["dash-array"].is_null()) {
              dashArray = style["dash-array"];
            }

            if (!style["css"].is_null()) {
              ls.setCss(style["css"]);
            }

            ls.setDashArray(dashArray);

            e->pl().addLine(l, dir, ls);
          } else {
            e->pl().addLine(l, dir);
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

        if (!idMap.count(id)) continue;
        LineNode* n = idMap[id];

        if (!props["excluded_line_conns"].is_null()) {
          for (auto excl : props["excluded_line_conns"]) {
            std::string rid = util::toString(excl["route"]);
            std::string nid1 = util::toString(excl["edge1_node"]);
            std::string nid2 = util::toString(excl["edge2_node"]);

            const Line* r = getLine(rid);

            if (!r) {
              LOG(WARN) << "line connection exclude defined in node " << id
                        << " for line " << rid << ", but no such line exists.";
              continue;
            }

            if (!idMap.count(nid1)) {
              LOG(WARN) << "line connection exclude defined in node " << id
                        << " for edge from " << nid1
                        << ", but no such node exists.";
              continue;
            }

            if (!idMap.count(nid2)) {
              LOG(WARN) << "line connection exclude defined in node " << id
                        << " for edge from " << nid2
                        << ", but no such node exists.";
              continue;
            }

            LineNode* n1 = idMap[nid1];
            LineNode* n2 = idMap[nid2];

            LineEdge* a = getEdg(n, n1);
            LineEdge* b = getEdg(n, n2);

            if (!a) {
              LOG(WARN) << "line connection exclude defined in node " << id
                        << " for edge from " << nid1
                        << ", but no such edge exists.";
              continue;
            }

            if (!b) {
              LOG(WARN) << "line connection exclude defined in node " << id
                        << " for edge from " << nid2
                        << ", but no such edge exists.";
              continue;
            }

            n->pl().addConnExc(r, a, b);
          }
        }
      }
    }
  }

  buildGrids();
}

// _____________________________________________________________________________
void LineGraph::buildGrids() {
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
void LineGraph::expandBBox(const Point<double>& p) {
  _bbox = util::geo::extendBox(p, _bbox);
}

// _____________________________________________________________________________
const util::geo::DBox& LineGraph::getBBox() const { return _bbox; }

// _____________________________________________________________________________
void LineGraph::topologizeIsects() {
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

    for (auto l : i.b->pl().getLines()) {
      if (l.direction == i.b->getFrom()) {
        ba->pl().addLine(l.line, i.b->getFrom());
        bb->pl().addLine(l.line, x);
      } else {
        ba->pl().addLine(l.line, x);
        bb->pl().addLine(l.line, i.b->getTo());
      }
    }

    auto aa =
        addEdg(i.a->getFrom(), x, i.a->pl().getPolyline().getSegment(0, pa));
    auto ab =
        addEdg(x, i.a->getTo(), i.a->pl().getPolyline().getSegment(pa, 1));

    _edgeGrid.add(*aa->pl().getGeom(), aa);
    _edgeGrid.add(*ab->pl().getGeom(), ab);

    for (auto l : i.a->pl().getLines()) {
      if (l.direction == i.b->getFrom()) {
        aa->pl().addLine(l.line, i.a->getFrom());
        ab->pl().addLine(l.line, x);
      } else {
        aa->pl().addLine(l.line, x);
        ab->pl().addLine(l.line, i.a->getTo());
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
ISect LineGraph::getNextIntersection() {
  for (auto n1 : *getNds()) {
    for (auto e1 : n1->getAdjList()) {
      if (e1->getFrom() != n1) continue;
      if (proced.find(e1) != proced.end()) continue;

      std::set<LineEdge*> neighbors;
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
void LineGraph::addLine(const Line* l) { _lines[l->id()] = l; }

// _____________________________________________________________________________
const Line* LineGraph::getLine(const std::string& id) const {
  if (_lines.find(id) != _lines.end()) return _lines.find(id)->second;
  return 0;
}

// _____________________________________________________________________________
LineNode* LineGraph::sharedNode(const LineEdge* a, const LineEdge* b) {
  LineNode* r = 0;
  if (a->getFrom() == b->getFrom() || a->getFrom() == b->getTo())
    r = a->getFrom();
  if (a->getTo() == b->getFrom() || a->getTo() == b->getTo()) r = a->getTo();
  return r;
}

// _____________________________________________________________________________
std::vector<LineOcc> LineGraph::getCtdLinesIn(const Line* r,
                                              const LineNode* dir,
                                              const LineEdge* fromEdge,
                                              const LineEdge* toEdge) {
  std::vector<LineOcc> ret;
  const auto* n = sharedNode(fromEdge, toEdge);
  if (!n || n->getDeg() == 1) return ret;

  for (const auto& to : toEdge->pl().getLines()) {
    if (to.line == r) {
      if (to.direction == 0 || dir == 0 || (to.direction == n && dir != n) ||
          (to.direction != n && dir == n)) {
        if (n->pl().connOccurs(r, fromEdge, toEdge)) ret.push_back(to);
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
std::vector<LineOcc> LineGraph::getCtdLinesIn(const LineEdge* fromEdge,
                                              const LineEdge* toEdge) {
  std::vector<LineOcc> ret;
  const auto* n = sharedNode(fromEdge, toEdge);
  if (!n) return ret;

  for (const auto& to : fromEdge->pl().getLines()) {
    auto r = getCtdLinesIn(to.line, to.direction, fromEdge, toEdge);
    ret.insert(ret.end(), r.begin(), r.end());
  }

  return ret;
}

// _____________________________________________________________________________
size_t LineGraph::getLDeg(const LineNode* nd) {
  size_t ret = 0;
  for (auto e : nd->getAdjList()) ret += e->pl().getLines().size();
  return ret;
}

// _____________________________________________________________________________
size_t LineGraph::getMaxLineNum(const LineNode* nd) {
  size_t ret = 0;
  for (auto e : nd->getAdjList())
    if (e->pl().getLines().size() > ret) ret = e->pl().getLines().size();
  return ret;
}

// _____________________________________________________________________________
size_t LineGraph::getMaxLineNum() {
  size_t ret = 0;
  for (auto nd : *getNds()) {
    size_t lineNum = getMaxLineNum(nd);
    if (lineNum > ret) ret = lineNum;
  }
  return ret;
}

// _____________________________________________________________________________
size_t LineGraph::maxDeg() const {
  size_t ret = 0;
  for (auto nd : getNds())
    if (nd->getDeg() > ret) ret = nd->getDeg();
  return ret;
}

// _____________________________________________________________________________
std::vector<const Line*> LineGraph::getSharedLines(const LineEdge* a,
                                                   const LineEdge* b) {
  std::vector<const Line*> ret;
  for (auto& to : a->pl().getLines()) {
    if (b->pl().hasLine(to.line)) ret.push_back(to.line);
  }

  return ret;
}

// _____________________________________________________________________________
std::vector<Partner> LineGraph::getPartners(const LineNode* n,
                                            const NodeFront* f,
                                            const LineOcc& ro) {
  std::vector<Partner> ret;
  for (const auto& nf : n->pl().fronts()) {
    if (&nf == f) continue;

    // TODO: if that is so, then why do we have the parameter n?
    assert(f->n == n);

    for (const LineOcc& to :
         getCtdLinesIn(ro.line, ro.direction, f->edge, nf.edge)) {
      Partner p(f, nf.edge, to.line);
      p.front = &nf;
      p.edge = nf.edge;
      p.line = to.line;
      ret.push_back(p);
    }
  }
  return ret;
}

// _____________________________________________________________________________
size_t LineGraph::getNumLines() const { return _lines.size(); }

// _____________________________________________________________________________
size_t LineGraph::getNumNds() const { return getNds().size(); }

// _____________________________________________________________________________
size_t LineGraph::getNumNds(bool topo) const { return 0; }
