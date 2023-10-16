// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "3rdparty/json.hpp"
#include "dot/Parser.h"
#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineGraph.h"
#include "shared/linegraph/LineNodePL.h"
#include "shared/style/LineStyle.h"
#include "util/Misc.h"
#include "util/String.h"
#include "util/graph/Algorithm.h"
#include "util/graph/Edge.h"
#include "util/graph/Node.h"
#include "util/log/Log.h"

using shared::linegraph::EdgeGrid;
using shared::linegraph::EdgeOrdering;
using shared::linegraph::ISect;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using shared::linegraph::LineOcc;
using shared::linegraph::NodeGrid;
using shared::linegraph::Partner;
using util::randomHtmlColor;
using util::geo::DPoint;
using util::geo::Point;
using util::graph::Algorithm;

// _____________________________________________________________________________
void LineGraph::readFromDot(std::istream* s) {
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
        n = addNd({util::geo::Point<double>(x, y),
                   std::numeric_limits<uint32_t>::max()});
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
        idMap[prevId] = addNd({util::geo::Point<double>(0, 0),
                               std::numeric_limits<uint32_t>::max()});

      for (size_t i = 1; i < ent.ids.size(); ++i) {
        std::string curId = ent.ids[i];

        if (idMap.find(curId) == idMap.end())
          idMap[curId] = addNd({util::geo::Point<double>(0, 0),
                                std::numeric_limits<uint32_t>::max()});
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
        } else if (ent.attrs.find("color") != ent.attrs.end()) {
          id = ent.attrs["color"];
        } else {
          id = util::toString(eid);
        }

        const Line* r = getLine(id);
        if (!r) {
          std::string label = ent.attrs.find("label") == ent.attrs.end()
                                  ? ""
                                  : ent.attrs["label"];
          std::string color = ent.attrs.find("color") == ent.attrs.end()
                                  ? ""
                                  : ent.attrs["color"];
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

  for (auto n : getNds()) {
    for (auto e : n->getAdjListOut()) {
      PolyLine<double> pl;
      pl << *e->getFrom()->pl().getGeom();
      pl << *e->getTo()->pl().getGeom();

      e->pl().setPolyline(pl);
      expandBBox(pl.front());
      expandBBox(pl.back());
    }
  }

  _bbox = util::geo::pad(_bbox, 100);
  buildGrids();
}
// _____________________________________________________________________________
void LineGraph::readFromTopoJson(nlohmann::json::array_t objects,
                                 nlohmann::json::array_t arcs,
                                 bool webMercCoords) {
  UNUSED(objects);
  UNUSED(arcs);
  UNUSED(webMercCoords);
  throw(std::runtime_error("TopoJSON input not yet implemented."));
}

// _____________________________________________________________________________
void LineGraph::readFromTopoJson(nlohmann::json::array_t objects,
                                 nlohmann::json::array_t arcs) {
  return readFromTopoJson(objects, arcs, false);
}

// _____________________________________________________________________________
void LineGraph::readFromGeoJson(nlohmann::json::array_t features) {
  return readFromGeoJson(features, false);
}

// _____________________________________________________________________________
void LineGraph::readFromGeoJson(nlohmann::json::array_t features,
                                bool webMercCoords) {
  _bbox = util::geo::Box<double>();

  std::map<std::string, LineNode*> idMap;

  for (auto feature : features) {
    auto props = feature["properties"];
    auto geom = feature["geometry"];
    if (geom["type"] == "Point") {
      std::string id;
      if (props.count("id")) id = props["id"].get<std::string>();

      std::vector<double> coords = geom["coordinates"];

      util::geo::DPoint point(coords[0], coords[1]);
      if (!webMercCoords) point = util::geo::latLngToWebMerc(point);

      if (id.empty()) {
        id = std::to_string(static_cast<int>(point.getX())) + "|" +
             std::to_string(static_cast<int>(point.getY()));
      }

      if (idMap.count(id)) continue;

      LineNode* n = addNd({point, std::numeric_limits<uint32_t>::max()});
      expandBBox(*n->pl().getGeom());

      if (props["component"].is_number())
        n->pl().setComponent(props["component"].get<size_t>());

      Station i("", "", *n->pl().getGeom());

      std::string sid = getStationId(props);
      std::string label = getStationLabel(props);
      if (!sid.empty() || !label.empty()) {
        i.id = sid;
        i.name = label;

        n->pl().addStop(i);
      }

      idMap[id] = n;
    }
  }

  // second pass, edges
  for (auto feature : features) {
    auto props = feature["properties"];
    auto geom = feature["geometry"];
    if (geom["type"] == "LineString") {
      std::string from =
          props["from"].is_null() ? "" : props["from"].get<std::string>();
      std::string to =
          props["to"].is_null() ? "" : props["to"].get<std::string>();

      if (geom["coordinates"].is_null()) continue;

      std::vector<std::vector<double>> coords = geom["coordinates"];

      size_t component = std::numeric_limits<uint32_t>::max();

      if (props["component"].is_number())
        component = props["component"].get<size_t>();

      PolyLine<double> pl;
      for (auto coord : coords) {
        double x = coord[0], y = coord[1];
        Point<double> p(x, y);
        if (!webMercCoords) p = util::geo::latLngToWebMerc(p);
        pl << p;
        expandBBox(p);
      }

      if (from.empty()) {
        from = std::to_string(static_cast<int>(pl.front().getX())) + "|" +
               std::to_string(static_cast<int>(pl.front().getY()));
        if (!idMap.count(from))
          idMap[from] = addNd({pl.getLine().front(), component});
      }

      if (to.empty()) {
        to = std::to_string(static_cast<int>(pl.back().getX())) + "|" +
             std::to_string(static_cast<int>(pl.back().getY()));
        if (!idMap.count(to))
          idMap[to] = addNd({pl.getLine().back(), component});
      }

      // pl.applyChaikinSmooth(3);

      LineNode* fromN = 0;
      LineNode* toN = 0;

      if (from.size()) {
        fromN = idMap[from];
        if (!fromN) {
          LOG(ERROR) << "Node \"" << from << "\" not found.";
          continue;
        }
      } else {
        fromN = addNd({pl.getLine().front(), component});
      }

      if (to.size()) {
        toN = idMap[to];
        if (!toN) {
          LOG(ERROR) << "Node \"" << to << "\" not found.";
          continue;
        }
      } else {
        toN = addNd({pl.getLine().back(), component});
      }

      if (fromN == toN) {
        LOGTO(DEBUG, std::cerr) << "Self edges are not supported, dropping...";
        continue;
      }

      LineEdge* e = addEdg(fromN, toN, pl);

      e->pl().setComponent(component);

      if (props["dontcontract"].is_number() && props["dontcontract"].get<int>())
        e->pl().setDontContract(true);

      extractLines(props, e, idMap);

      // if no lines were extracted, completely delete edge
      if (e->pl().getLines().empty()) delEdg(e->getFrom(), e->getTo());
    }
  }

  // third pass, exceptions (TODO: do this in the first part, store in some
  // data structure, add here!)
  for (auto feature : features) {
    auto props = feature["properties"];
    auto geom = feature["geometry"];
    if (geom["type"] == "Point") {
      std::string id;
      if (props.count("id")) id = props["id"].get<std::string>();

      if (id.empty()) {
        std::vector<double> coords = geom["coordinates"];

        util::geo::DPoint point(coords[0], coords[1]);
        if (!webMercCoords) point = util::geo::latLngToWebMerc(point);

        id = std::to_string(static_cast<int>(point.getX())) + "|" +
             std::to_string(static_cast<int>(point.getY()));
      }

      if (!idMap.count(id)) continue;
      LineNode* n = idMap[id];

      if (!props["not_serving"].is_null()) {
        for (auto excl : props["not_serving"]) {
          std::string lid = excl.get<std::string>();

          const Line* r = getLine(lid);

          if (!r) {
            LOG(WARN) << "line " << lid << " marked as not served in in node "
                      << id << ", but no such line exists.";
            continue;
          }

          n->pl().addLineNotServed(r);
        }
      }

      if (!props["excluded_conn"].is_null()) {
        for (auto excl : props["excluded_conn"]) {
          std::string lid = excl["line"].get<std::string>();
          std::string nid1 = excl["node_from"].get<std::string>();
          std::string nid2 = excl["node_to"].get<std::string>();

          const Line* r = getLine(lid);

          if (!r) {
            LOG(WARN) << "line connection exclude defined in node " << id
                      << " for line " << lid << ", but no such line exists.";
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

  _bbox = util::geo::pad(_bbox, 100);

  buildGrids();
}

// _____________________________________________________________________________
void LineGraph::readFromJson(std::istream* s) { return readFromJson(s, false); }

// _____________________________________________________________________________
void LineGraph::readFromJson(std::istream* s, bool useWebMercCoords) {
  nlohmann::json j;
  (*s) >> j;

  if (j["type"] == "FeatureCollection") {
    readFromGeoJson(j["features"], useWebMercCoords);
    if (j.count("properties")) _graphProps = j["properties"];
  }
  if (j["type"] == "Topology")
    readFromTopoJson(j["objects"], j["arcs"], useWebMercCoords);
}

// _____________________________________________________________________________
void LineGraph::buildGrids() {
  _nodeGrid = NodeGrid();
  _edgeGrid = EdgeGrid();

  for (auto n : getNds()) {
    _nodeGrid.add(*n->pl().getGeom(), n);
    for (auto e : n->getAdjListOut()) {
      if (e->getFrom() != n) continue;
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
    auto x = addNd({i.bp.p, i.a->pl().getComponent()});

    double pa = i.a->pl().getPolyline().projectOn(i.bp.p).totalPos;

    auto ba = addEdg(i.b->getFrom(), x, i.b->pl());
    ba->pl().setPolyline(i.b->pl().getPolyline().getSegment(0, i.bp.totalPos));

    auto bb = addEdg(x, i.b->getTo(), i.b->pl());
    bb->pl().setPolyline(i.b->pl().getPolyline().getSegment(i.bp.totalPos, 1));

    edgeRpl(i.b->getFrom(), i.b, ba);
    edgeRpl(i.b->getTo(), i.b, bb);

    nodeRpl(ba, i.b->getTo(), x);
    nodeRpl(bb, i.b->getFrom(), x);

    _edgeGrid.add(*ba->pl().getGeom(), ba);
    _edgeGrid.add(*bb->pl().getGeom(), bb);

    auto aa = addEdg(i.a->getFrom(), x, i.a->pl());
    aa->pl().setPolyline(i.a->pl().getPolyline().getSegment(0, pa));
    auto ab = addEdg(x, i.a->getTo(), i.a->pl());
    ab->pl().setPolyline(i.a->pl().getPolyline().getSegment(pa, 1));

    edgeRpl(i.a->getFrom(), i.a, aa);
    edgeRpl(i.a->getTo(), i.a, ab);

    nodeRpl(aa, i.b->getTo(), x);
    nodeRpl(ab, i.b->getFrom(), x);

    _edgeGrid.remove(aa);
    _edgeGrid.remove(ab);

    _edgeGrid.add(*aa->pl().getGeom(), aa);
    _edgeGrid.add(*ab->pl().getGeom(), ab);

    _edgeGrid.remove(i.a);
    _edgeGrid.remove(i.b);

    assert(getEdg(i.a->getFrom(), i.a->getTo()));
    assert(getEdg(i.b->getFrom(), i.b->getTo()));
    delEdg(i.a->getFrom(), i.a->getTo());
    delEdg(i.b->getFrom(), i.b->getTo());
  }
}

// _____________________________________________________________________________
std::set<LineEdge*> LineGraph::getNeighborEdges(const util::geo::DLine& line,
                                                double d) const {
  std::set<LineEdge*> neighbors;
  _edgeGrid.get(line, d, &neighbors);

  return neighbors;
}

// _____________________________________________________________________________
ISect LineGraph::getNextIntersection() {
  for (auto n1 : getNds()) {
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
            // if the intersection is near a shared node, ignore
            auto shrdNd = sharedNode(e1, e2);
            if (shrdNd &&
                util::geo::dist(*shrdNd->pl().getGeom(), ret.bp.p) < 100) {
              continue;
            }

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
void LineGraph::fillMissingColors() {
  for (auto& l : _lines) {
    auto ll = const_cast<Line*>(l.second);
    if (ll->color().empty()) {
      std::string c = randomHtmlColor();
      ll->setColor(c);
    }
  }
}

// _____________________________________________________________________________
bool LineGraph::lineCtd(const LineEdge* frEdg, const LineEdge* toEdg,
                        const Line* line) {
  if (!frEdg->pl().hasLine(line) || !toEdg->pl().hasLine(line)) return false;
  auto frLn = frEdg->pl().lineOcc(line);
  auto toLn = toEdg->pl().lineOcc(line);
  return lineCtd(frEdg, frLn, toEdg, toLn);
}

// _____________________________________________________________________________
bool LineGraph::lineCtd(const LineEdge* frEdg, const LineOcc& frLn,
                        const LineEdge* toEdg, const LineOcc& toLn) {
  if (frLn.line != toLn.line) return false;
  const auto* n = sharedNode(frEdg, toEdg);
  if (!n || n->getDeg() == 1) return false;
  return (frLn.direction == 0 || toLn.direction == 0 ||
          (frLn.direction == n && toLn.direction != n) ||
          (frLn.direction != n && toLn.direction == n)) &&
         n->pl().connOccurs(frLn.line, frEdg, toEdg);
}

// _____________________________________________________________________________
std::vector<LineOcc> LineGraph::getCtdLinesIn(const LineOcc& frLn,
                                              const LineEdge* frEdge,
                                              const LineEdge* toEdge) {
  std::vector<LineOcc> ret;
  const auto* n = sharedNode(frEdge, toEdge);
  if (!n || n->getDeg() == 1) return ret;

  for (const auto& toLn : toEdge->pl().getLines()) {
    if (lineCtd(frEdge, frLn, toEdge, toLn)) ret.push_back(toLn);
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
    auto r = getCtdLinesIn(to, fromEdge, toEdge);
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
size_t LineGraph::getMaxLineNum() const {
  size_t ret = 0;
  for (auto nd : getNds()) {
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
size_t LineGraph::numLines() const { return _lines.size(); }

// _____________________________________________________________________________
size_t LineGraph::numNds() const { return getNds().size(); }

// _____________________________________________________________________________
size_t LineGraph::numNds(bool topo) const {
  size_t ret = 0;
  for (auto nd : getNds()) {
    if ((nd->pl().stops().size() == 0) ^ !topo) ret++;
  }
  return ret;
}

// _____________________________________________________________________________
size_t LineGraph::numEdgs() const {
  size_t ret = 0;
  for (auto nd : getNds()) {
    for (auto e : nd->getAdjList()) {
      if (e->getFrom() != nd) continue;
      ret++;
    }
  }
  return ret;
}

// _____________________________________________________________________________
NodeGrid* LineGraph::getNdGrid() { return &_nodeGrid; }

// _____________________________________________________________________________
const NodeGrid& LineGraph::getNdGrid() const { return _nodeGrid; }

// _____________________________________________________________________________
EdgeGrid* LineGraph::getEdgGrid() { return &_edgeGrid; }

// _____________________________________________________________________________
const EdgeGrid& LineGraph::getEdgGrid() const { return _edgeGrid; }

// _____________________________________________________________________________
std::set<const shared::linegraph::Line*> LineGraph::servedLines(
    const shared::linegraph::LineNode* n) {
  std::set<const shared::linegraph::Line*> ret;

  for (auto e : n->getAdjList()) {
    for (auto l : e->pl().getLines()) {
      if (n->pl().lineServed(l.line)) {
        ret.insert(l.line);
      }
    }
  }
  return ret;
}

// _____________________________________________________________________________
EdgeOrdering LineGraph::edgeOrdering(LineNode* n, bool useOrigNextNode) {
  EdgeOrdering order;
  util::geo::DPoint a = *n->pl().getGeom();

  for (auto e : n->getAdjList()) {
    util::geo::DPoint b;
    if (useOrigNextNode) {
      b = *e->getOtherNd(n)->pl().getGeom();
    } else {
      auto other = e->getOtherNd(n);
      if (e->pl().getGeom()->size() > 2) {
        if (e->getTo() == n) {
          b = e->pl().getGeom()->at(e->pl().getGeom()->size() - 2);
        } else {
          b = e->pl().getGeom()->at(1);
        }
      } else {
        b = *other->pl().getGeom();
      }
    }

    // get the angles
    double deg = util::geo::angBetween(a, b) - M_PI / 2;
    // make sure the ordering start at 12 o'clock
    if (deg <= 0) deg += M_PI * 2;

    order.add(e, deg);
  }

  return order;
}

// _____________________________________________________________________________
void LineGraph::splitNode(LineNode* n, size_t maxDeg) {
  assert(maxDeg > 2);

  if (n->getAdjList().size() > maxDeg) {
    std::vector<std::pair<LineEdge*, double>> combine;

    const auto& eo = edgeOrdering(n, true);
    const auto& orig = eo.getOrderedSet();

    combine.insert(combine.begin(), orig.begin() + (maxDeg - 1), orig.end());

    // for the new geometry, take the average angle
    double refAngle = 0;
    util::geo::MultiPoint<double> mp;

    for (auto eo : combine) {
      mp.push_back(*eo.first->getOtherNd(n)->pl().getGeom());
    }
    refAngle = util::geo::angBetween(*n->pl().getGeom(), mp);

    auto geom = *n->pl().getGeom();
    geom.setX(geom.getX() + 10 * cos(refAngle));
    geom.setY(geom.getY() + 10 * sin(refAngle));
    // add a new node
    auto cn = addNd({geom, n->pl().getComponent()});

    // add the new trunk edge
    auto ce =
        addEdg(n, cn, LineEdgePL({*n->pl().getGeom(), *cn->pl().getGeom()}));

    for (auto eo : combine) {
      LineEdge* newEdg = 0;
      if (eo.first->getFrom() == n) {
        newEdg = addEdg(cn, eo.first->getOtherNd(n), eo.first->pl());
      } else {
        newEdg = addEdg(eo.first->getOtherNd(n), cn, eo.first->pl());
      }

      // replace direction markers in the new edge
      nodeRpl(newEdg, n, cn);

      // replace exception containing the old edge in the remaining target node
      edgeRpl(eo.first->getOtherNd(n), eo.first, newEdg);

      for (auto lo : eo.first->pl().getLines()) {
        ce->pl().addLine(lo.line, 0);
      }

      ce->pl().setComponent(eo.first->pl().getComponent());

      // in the old node, replace any exception occurence of this node with
      // the new trunk edge
      edgeRpl(n, eo.first, ce);
      _edgeGrid.remove(eo.first);
      delEdg(eo.first->getFrom(), eo.first->getTo());
    }

    // no continuation lines across the new node, only to the trunk edge
    for (auto lo : ce->pl().getLines()) {
      for (auto ea : cn->getAdjList()) {
        if (ea == ce) continue;
        for (auto eb : cn->getAdjList()) {
          if (eb == ce) continue;
          if (ea == eb) continue;
          cn->pl().addConnExc(lo.line, ea, eb);
        }
      }
    }

    assert(n->getDeg() <= maxDeg);

    // recursively split until max deg is satisfied
    if (cn->getDeg() > maxDeg) splitNode(cn, maxDeg);
  }
}

// _____________________________________________________________________________
void LineGraph::splitNodes(size_t maxDeg) {
  std::vector<LineNode*> toSplit;
  for (auto n : getNds()) {
    if (n->getDeg() > maxDeg) toSplit.push_back(n);
  }
  for (auto n : toSplit) splitNode(n, maxDeg);
}

// _____________________________________________________________________________
void LineGraph::edgeDel(LineNode* n, const LineEdge* oldE) {
  // remove from from
  for (auto& r : n->pl().getConnExc()) {
    auto exFr = r.second.begin();
    while (exFr != r.second.end()) {
      if (exFr->first == oldE) {
        exFr = r.second.erase(exFr);
      } else {
        exFr++;
      }
    }
  }

  // remove from to
  for (auto& r : n->pl().getConnExc()) {
    for (auto& exFr : r.second) {
      auto exTo = exFr.second.begin();
      while (exTo != exFr.second.end()) {
        if (*exTo == oldE) {
          exTo = exFr.second.erase(exTo);
        } else {
          exTo++;
        }
      }
    }
  }
}

// _____________________________________________________________________________
void LineGraph::edgeRpl(LineNode* n, const LineEdge* oldE,
                        const LineEdge* newE) {
  if (oldE == newE) return;

  // replace in from
  for (auto& r : n->pl().getConnExc()) {
    auto exFr = r.second.begin();
    while (exFr != r.second.end()) {
      if (exFr->first == oldE) {
        std::swap(r.second[newE], exFr->second);
        exFr = r.second.erase(exFr);
      } else {
        exFr++;
      }
    }
  }

  // replace in to
  for (auto& r : n->pl().getConnExc()) {
    for (auto& exFr : r.second) {
      auto exTo = exFr.second.begin();
      while (exTo != exFr.second.end()) {
        if (*exTo == oldE) {
          exFr.second.insert(newE);
          exTo = exFr.second.erase(exTo);
        } else {
          exTo++;
        }
      }
    }
  }

  // replace in node fronts
  for (auto nf : n->pl().fronts()) {
    if (nf.edge == oldE) {
      n->pl().delFrontFor(nf.edge);
      nf.edge = const_cast<LineEdge*>(newE);
      n->pl().addFront(nf);
    }
  }
}

// _____________________________________________________________________________
void LineGraph::nodeRpl(LineEdge* e, const LineNode* oldN,
                        const LineNode* newN) {
  auto ro = e->pl().getLines().begin();
  while (ro != e->pl().getLines().end()) {
    if (ro->direction == oldN) {
      shared::linegraph::LineOcc newRo = *ro;
      newRo.direction = newN;

      e->pl().updateLineOcc(newRo);
    }
    ro++;
  }
}

// _____________________________________________________________________________
void LineGraph::contractStrayNds() {
  std::vector<LineNode*> toDel;
  for (auto n : getNds()) {
    if (n->pl().stops().size()) continue;
    if (n->getAdjList().size() != 2) continue;

    auto eA = n->getAdjList().front();
    auto eB = n->getAdjList().back();

    // check if all lines continue over this node
    bool lineEqual = true;

    for (const auto& lo : eA->pl().getLines()) {
      if (!lineCtd(eA, eB, lo.line)) {
        lineEqual = false;
        break;
      }
    }

    if (!lineEqual) continue;

    for (const auto& lo : eB->pl().getLines()) {
      if (!lineCtd(eB, eA, lo.line)) {
        lineEqual = false;
        break;
      }
    }

    if (!lineEqual) continue;

    toDel.push_back(n);
  }

  for (auto n : toDel) {
    if (n->getAdjList().size() == 2) {
      auto a = n->getAdjList().front();
      auto b = n->getAdjList().back();

      // if this combination would turn our graph into a multigraph,
      // dont do it!
      if (getEdg(a->getOtherNd(n), b->getOtherNd(n))) continue;

      auto pl = a->pl();

      if (a->getTo() == n) {
        if (a->pl().getGeom() && b->pl().getGeom()) {
          util::geo::Line<double> l;
          l.insert(l.end(), a->pl().getGeom()->begin(),
                   a->pl().getGeom()->end());
          if (b->getFrom() == n)
            l.insert(l.end(), b->pl().getGeom()->begin(),
                     b->pl().getGeom()->end());
          else
            l.insert(l.end(), b->pl().getGeom()->rbegin(),
                     b->pl().getGeom()->rend());
          pl.setPolyline(PolyLine<double>(l));
        } else {
          pl.setPolyline(PolyLine<double>(*a->getFrom()->pl().getGeom(),
                                          *b->getOtherNd(n)->pl().getGeom()));
        }

        assert(!getEdg(a->getFrom(), b->getOtherNd(n)));
        auto newE = addEdg(a->getFrom(), b->getOtherNd(n), pl);
        edgeRpl(a->getFrom(), b, newE);
        edgeRpl(b->getOtherNd(n), b, newE);
        edgeRpl(a->getFrom(), a, newE);
        edgeRpl(b->getOtherNd(n), a, newE);
      } else {
        if (a->pl().getGeom() && b->pl().getGeom()) {
          util::geo::Line<double> l;
          if (b->getTo() == n)
            l.insert(l.end(), b->pl().getGeom()->begin(),
                     b->pl().getGeom()->end());
          else
            l.insert(l.end(), b->pl().getGeom()->rbegin(),
                     b->pl().getGeom()->rend());
          l.insert(l.end(), a->pl().getGeom()->begin(),
                   a->pl().getGeom()->end());
          pl.setPolyline(PolyLine<double>(l));
        } else {
          pl.setPolyline(PolyLine<double>(*a->getTo()->pl().getGeom(),
                                          *b->getOtherNd(n)->pl().getGeom()));
        }
        assert(!getEdg(b->getOtherNd(n), a->getTo()));
        auto newE = addEdg(b->getOtherNd(n), a->getTo(), pl);
        edgeRpl(a->getTo(), b, newE);
        edgeRpl(b->getOtherNd(n), b, newE);
        edgeRpl(a->getTo(), a, newE);
        edgeRpl(b->getOtherNd(n), a, newE);
      }

      for (auto e : n->getAdjList()) {
        _edgeGrid.remove(e);
      }

      delNd(n);
    }
  }
}

// _____________________________________________________________________________
void LineGraph::contractEdge(LineEdge* e) {
  auto n1 = e->getFrom();
  auto n2 = e->getTo();
  auto otherP = n2->pl().getGeom();
  auto newGeom = DPoint((n1->pl().getGeom()->getX() + otherP->getX()) / 2,
                        (n1->pl().getGeom()->getY() + otherP->getY()) / 2);
  LineNode* n = 0;

  // n2 is always the target node
  auto n1Pl = n1->pl();

  if (e->getTo()->pl().stops().size() > 0) {
    auto servedLines = LineGraph::servedLines(e->getTo());
    n = mergeNds(e->getFrom(), e->getTo());
    for (auto l : LineGraph::servedLines(n)) {
      if (!servedLines.count(l)) n->pl().addLineNotServed(l);
    }
  } else if (e->getFrom()->pl().stops().size() > 0) {
    auto servedLines = LineGraph::servedLines(e->getFrom());
    n = mergeNds(e->getTo(), e->getFrom());

    for (auto l : LineGraph::servedLines(n)) {
      if (!servedLines.count(l)) n->pl().addLineNotServed(l);
    }
  } else {
    n = mergeNds(e->getTo(), e->getFrom());
  }

  n->pl().setGeom(newGeom);
}

// _____________________________________________________________________________
LineNode* LineGraph::mergeNds(LineNode* a, LineNode* b) {
  auto eConn = getEdg(a, b);

  std::vector<std::pair<const Line*, std::pair<LineNode*, LineNode*>>> ex;

  auto lineAdjA = servedLines(a);
  auto lineAdjB = servedLines(b);

  NotServedLines notServedNew;

  for (const auto& l : lineAdjA) {
    if (!a->pl().lineServed(l) && (!b->pl().lineServed(l) || !lineAdjB.count(l) )) {
      notServedNew.insert(l);
    }
  }

  for (const auto& l : lineAdjB) {
    if (!b->pl().lineServed(l) && (!a->pl().lineServed(l) || !lineAdjA.count(l) )) {
      notServedNew.insert(l);
    }
  }

  if (eConn) {
    for (auto& ex : a->pl().getConnExc()) {
      auto line = ex.first;
      for (auto& exFr : ex.second) {
        auto fr = exFr.first;

        auto ti = exFr.second.begin();
        while (ti != exFr.second.end()) {
          auto to = *ti;
          if (fr != eConn && to != eConn && lineCtd(fr, eConn, line) &&
              lineCtd(to, eConn, line) && terminatesAt(eConn, b, line)) {
            ti = exFr.second.erase(ti);
          } else {
            ti++;
          }
        }
      }
    }

    for (auto& ex : b->pl().getConnExc()) {
      auto line = ex.first;
      for (auto& exFr : ex.second) {
        auto fr = exFr.first;

        auto ti = exFr.second.begin();
        while (ti != exFr.second.end()) {
          auto to = *ti;
          if (fr != eConn && to != eConn && lineCtd(fr, eConn, line) &&
              lineCtd(to, eConn, line) && terminatesAt(eConn, a, line)) {
            ti = exFr.second.erase(ti);
          } else {
            ti++;
          }
        }
      }
    }

    for (auto fr : a->getAdjList()) {
      if (fr == eConn) continue;
      for (auto lo : fr->pl().getLines()) {
        for (auto to : b->getAdjList()) {
          if (to == eConn) continue;
          if (fr->pl().hasLine(lo.line) && to->pl().hasLine(lo.line) &&
              (!lineCtd(fr, eConn, lo.line) || !lineCtd(eConn, to, lo.line))) {
            auto frNd = fr->getOtherNd(a);
            auto toNd = to->getOtherNd(b);
            auto exA = getEdg(b, frNd);
            auto exB = getEdg(b, toNd);
            if (exA && exB && exA->pl().hasLine(lo.line) &&
                exB->pl().hasLine(lo.line) && lineCtd(exA, exB, lo.line))
              continue;
            exA = getEdg(a, frNd);
            exB = getEdg(a, toNd);
            if (exA && exB && exA->pl().hasLine(lo.line) &&
                exB->pl().hasLine(lo.line) && lineCtd(exA, exB, lo.line))
              continue;
            ex.push_back({lo.line, {frNd, toNd}});
          }
        }
      }
    }
  }

  if (eConn) {
    edgeDel(a, eConn);
    edgeDel(b, eConn);
    _edgeGrid.remove(eConn);
    delEdg(a, b);
  }

  for (auto e : a->getAdjList()) {
    if (e->getFrom() != a) continue;
    if (e->getTo() == b) continue;
    auto ex = getEdg(b, e->getTo());  // check if edge already exists
    auto newE = addEdg(b, e->getTo(), e->pl());
    if (ex) {
      for (auto lo : e->pl().getLines()) {
        if (lo.direction == a) lo.direction = b;
        newE->pl().addLine(lo.line, lo.direction);
      }
    } else {
      // else add to grid
      _edgeGrid.add(*newE->pl().getGeom(), newE);
    }
    edgeRpl(b, e, newE);
    edgeRpl(e->getTo(), e, newE);
    nodeRpl(newE, a, b);
  }

  for (auto e : a->getAdjList()) {
    if (e->getTo() != a) continue;
    if (e->getFrom() == b) continue;
    auto ex = getEdg(e->getFrom(), b);  // check if edge already exists
    auto newE = addEdg(e->getFrom(), b, e->pl());
    if (ex) {
      for (auto lo : e->pl().getLines()) {
        if (lo.direction == a) lo.direction = b;
        newE->pl().addLine(lo.line, lo.direction);
      }
    } else {
      // else add to grid
      _edgeGrid.add(*newE->pl().getGeom(), newE);
    }
    edgeRpl(b, e, newE);
    edgeRpl(e->getFrom(), e, newE);
    nodeRpl(newE, a, b);
  }

  for (auto e : a->getAdjList()) _edgeGrid.remove(e);
  delNd(a);

  for (const auto& newEx : ex) {
    b->pl().addConnExc(newEx.first, getEdg(b, newEx.second.first),
                       getEdg(b, newEx.second.second));
  }

  b->pl().setNotServed(notServedNew);

  return b;
}

// _____________________________________________________________________________
std::vector<Partner> LineGraph::getPartners(const LineNode* nd,
                                            const LineEdge* e,
                                            const LineOcc& lo) {
  std::vector<Partner> ret;
  for (auto toEdg : nd->getAdjList()) {
    if (toEdg == e) continue;

    for (const LineOcc& to : getCtdLinesIn(lo, e, toEdg)) {
      Partner p(toEdg, to.line);
      ret.push_back(p);
    }
  }
  return ret;
}

// _____________________________________________________________________________
void LineGraph::contractEdges(double d) { contractEdges(d, false); }

// _____________________________________________________________________________
void LineGraph::contractEdges(double d, bool onlyNonStatConns) {
  // TODO: the problem here is that contractEdge(e) may delete and replace edges
  // in the graph, which is why we cannot just build a list of edges eligible
  // for contraction and iterate over it - we would have to record the changes
  // made in contractEdge(e) and propagate it back.

  std::vector<std::pair<size_t, LineEdge*>> cands;

breakfor:
  cands.clear();

  for (auto n1 : getNds()) {
    for (auto e : n1->getAdjList()) {
      if (e->getFrom() != n1) continue;
      if (onlyNonStatConns && (e->getFrom()->pl().stops().size() ||
                               e->getTo()->pl().stops().size()))
        continue;
      if (onlyNonStatConns && (e->getFrom()->pl().stops().size() ||
                               e->getTo()->pl().stops().size()))
        continue;

      if (!e->pl().dontContract() && e->pl().getPolyline().shorterThan(d)) {
        if (e->getOtherNd(n1)->getAdjList().size() > 1 &&
            (n1->pl().stops().size() == 0 || n1->getAdjList().size() > 1) &&
            (n1->pl().stops().size() == 0 ||
             e->getOtherNd(n1)->pl().stops().size() == 0 ||
             n1->pl().stops().front().name ==
                 e->getOtherNd(n1)->pl().stops().front().name)) {
          // first contract edges with lower number of adjacent nodes,
          // on ties use shorter edge
          cands.push_back({e->getFrom()->getDeg() + e->getTo()->getDeg() +
                               e->pl().getPolyline().getLength(),
                           e});
        }
      }
    }
  }

  std::sort(cands.begin(), cands.end());
  for (auto e : cands) {
    contractEdge(e.second);
    goto breakfor;
  }
}

// _____________________________________________________________________________
bool LineGraph::isTerminus(const LineNode* nd) {
  for (auto e : nd->getAdjList()) {
    for (const auto& lo : e->pl().getLines()) {
      if (terminatesAt(e, nd, lo.line)) return true;
    }
  }

  return false;
}

// _____________________________________________________________________________
bool LineGraph::terminatesAt(const LineEdge* fromEdge, const LineNode* terminus,
                             const Line* line) {
  for (const auto& toEdg : terminus->getAdjList()) {
    if (toEdg == fromEdge) continue;

    if (lineCtd(fromEdge, toEdg, line)) {
      return false;
    }
  }

  return true;
}

// _____________________________________________________________________________
double LineGraph::searchSpaceSize() const {
  double ret = 1;

  for (auto n : getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      ret *= util::factorial(e->pl().getLines().size());
    }
  }

  return ret;
}

// _____________________________________________________________________________
size_t LineGraph::numConnExcs() const {
  size_t ret = 0;
  for (auto n : getNds()) ret += n->pl().numConnExcs();
  return ret;
}

// _____________________________________________________________________________
void LineGraph::extractLines(const nlohmann::json::object_t& props, LineEdge* e,
                             const std::map<std::string, LineNode*>& idMap) {
  auto i = props.find("lines");

  if (i == props.end()) {
    extractLine(props, e, idMap);
  } else {
    for (auto line : i->second) {
      extractLine(line, e, idMap);
    }
  }
}

// _____________________________________________________________________________
std::string LineGraph::getLineId(const nlohmann::json::object_t& line) {
  std::string id;

  if (line.count("id")) {
    id = line.at("id").get<std::string>();
  } else if (line.count("\"id\"")) {
    id = line.at("\"id\"").get<std::string>();
  } else if (line.count("?id")) {
    id = line.at("?id").get<std::string>();
  } else if (line.count("\"?id\"")) {
    id = line.at("\"?id\"").get<std::string>();
  }

  id = util::trim(id, "\"");

  if (id.empty()) {
    if (!getLineLabel(line).empty()) {
      id = getLineLabel(line);
    } else if (!getLineColor(line).empty()) {
      id = getLineColor(line);
    }
  }

  if (id.empty()) id = "<EMPTY PLACEHOLDER>";

  return util::trim(id, "\"");
}

// _____________________________________________________________________________
std::string LineGraph::getLineColor(const nlohmann::json::object_t& line) {
  std::string color = "";

  if (line.count("color")) {
    color = line.at("color").get<std::string>();
  } else if (line.count("\"color\"")) {
    color = line.at("\"color\"").get<std::string>();
  } else if (line.count("?color")) {
    color = line.at("?color").get<std::string>();
  } else if (line.count("\"?color\"")) {
    color = line.at("\"?color\"").get<std::string>();
  } else if (line.count("?colour")) {
    color = line.at("?colour").get<std::string>();
  } else if (line.count("\"?colour\"")) {
    color = line.at("\"?colour\"").get<std::string>();
  } else if (line.count("colour")) {
    color = line.at("colour").get<std::string>();
  } else if (line.count("\"colour\"")) {
    color = line.at("\"colour\"").get<std::string>();
  } else if (line.count("col")) {
    color = line.at("col").get<std::string>();
  } else if (line.count("\"col\"")) {
    color = line.at("\"col\"").get<std::string>();
  } else if (line.count("?col")) {
    color = line.at("?col").get<std::string>();
  } else if (line.count("\"?col\"")) {
    color = line.at("\"?col\"").get<std::string>();
  }

  return util::normHtmlColor(util::trim(util::trim(color, "\""), "#"));
}

// _____________________________________________________________________________
std::string LineGraph::getLineLabel(const nlohmann::json::object_t& line) {
  std::string label;

  if (line.count("label")) {
    label = line.at("label").get<std::string>();
  } else if (line.count("\"label\"")) {
    label = line.at("\"label\"").get<std::string>();
  } else if (line.count("?label")) {
    label = line.at("?label").get<std::string>();
  } else if (line.count("\"?label\"")) {
    label = line.at("\"?label\"").get<std::string>();
  } else if (line.count("lbl")) {
    label = line.at("lbl").get<std::string>();
  } else if (line.count("\"lbl\"")) {
    label = line.at("\"lbl\"").get<std::string>();
  } else if (line.count("?lbl")) {
    label = line.at("?lbl").get<std::string>();
  } else if (line.count("\"?lbl\"")) {
    label = line.at("\"?lbl\"").get<std::string>();
  }

  return util::trim(label, "\"");
}

// _____________________________________________________________________________
std::string LineGraph::getStationId(const nlohmann::json::object_t& props) {
  std::string id;

  if (props.count("station_id")) {
    id = props.at("station_id").get<std::string>();
  } else if (props.count("\"station_id\"")) {
    id = props.at("\"station_id\"").get<std::string>();
  } else if (props.count("?station_id")) {
    id = props.at("?station_id").get<std::string>();
  } else if (props.count("\"?station_id\"")) {
    id = props.at("\"?station_id\"").get<std::string>();
  }

  return util::trim(id, "\"");
}

// _____________________________________________________________________________
std::string LineGraph::getStationLabel(const nlohmann::json::object_t& line) {
  std::string label;

  if (line.count("station_label")) {
    label = line.at("station_label").get<std::string>();
  } else if (line.count("\"station_label\"")) {
    label = line.at("\"station_label\"").get<std::string>();
  } else if (line.count("?station_label")) {
    label = line.at("?station_label").get<std::string>();
  } else if (line.count("\"?station_label\"")) {
    label = line.at("\"?station_label\"").get<std::string>();
  } else if (line.count("station_lbl")) {
    label = line.at("station_lbl").get<std::string>();
  } else if (line.count("\"station_lbl\"")) {
    label = line.at("\"station_lbl\"").get<std::string>();
  } else if (line.count("?station_lbl")) {
    label = line.at("?station_lbl").get<std::string>();
  } else if (line.count("\"?station_lbl\"")) {
    label = line.at("\"?station_lbl\"").get<std::string>();
  }

  return util::trim(label, "\"");
}

// _____________________________________________________________________________
void LineGraph::extractLine(const nlohmann::json::object_t& line, LineEdge* e,
                            const std::map<std::string, LineNode*>& idMap) {
  std::string id = getLineId(line);
  std::string color = getLineColor(line);
  std::string label = getLineLabel(line);

  const Line* l = getLine(id);
  if (!l) {
    l = new Line(id, label, color);
    addLine(l);
  }

  LineNode* dir = 0;

  if (line.count("direction") && idMap.count(line.at("direction"))) {
    dir = idMap.at(line.at("direction").get<std::string>());
  }

  if (line.count("style") || line.count("outline-style")) {
    shared::style::LineStyle ls;

    if (line.count("style")) ls.setCss(line.at("style"));
    if (line.count("outline-style")) ls.setOutlineCss(line.at("outline-style"));

    e->pl().addLine(l, dir, ls);
  } else {
    e->pl().addLine(l, dir);
  }
}

// _____________________________________________________________________________
std::vector<LineGraph> LineGraph::distConnectedComponents(double d,
                                                          bool write) {
  return distConnectedComponents(d, write, 0);
}

// _____________________________________________________________________________
std::vector<LineGraph> LineGraph::distConnectedComponents(double d, bool write,
                                                          size_t* offset) {
  std::vector<LineGraph> ret;

  size_t idOffset = 0;

  if (offset) idOffset = *offset;

  // first pass, collect components
  const auto& origComps = Algorithm::connectedComponents(*this);
  std::unordered_map<LineNode*, size_t> ndToComp;
  for (size_t comp = 0; comp < origComps.size(); comp++) {
    for (auto nd : origComps[comp]) {
      ndToComp[nd] = comp;
    }
  }

  std::vector<LineEdge*> addedEdgs;

  for (auto nd : getNds()) {
    // connect each node with nodes within distance
    std::set<LineNode*> cands;
    _nodeGrid.get(*nd->pl().getGeom(), d, &cands);

    for (auto cand : cands) {
      if (cand != nd && !getEdg(nd, cand) && !getEdg(cand, nd) &&
          ndToComp[nd] != ndToComp[cand] &&
          util::geo::dist(*nd->pl().getGeom(), *cand->pl().getGeom()) <= d) {
        addedEdgs.push_back(addEdg(nd, cand));
      }
    }
  }

  const auto& geoComps = Algorithm::connectedComponents(*this);

  // delete the newly added edges
  for (auto e : addedEdgs) {
    _edgeGrid.remove(e);
    delEdg(e->getFrom(), e->getTo());
  }

  ret.resize(geoComps.size());

  for (size_t comp = 0; comp < geoComps.size(); comp++) {
    auto* tg = &ret[comp];

    std::unordered_map<LineNode*, LineNode*> nm;
    std::unordered_map<LineEdge*, LineEdge*> em;

    // add nodes
    for (auto nd : geoComps[comp]) {
      if (write) nd->pl().setComponent(idOffset + comp);
      nm[nd] = tg->addNd(nd->pl());
      tg->expandBBox(*nd->pl().getGeom());
    }

    // add edges
    for (auto nd : geoComps[comp]) {
      for (auto edg : nd->getAdjList()) {
        if (edg->getFrom() != nd) continue;
        if (write) edg->pl().setComponent(idOffset + comp);
        if (edg->pl().getLines().size() == 0) {
          continue;
        }

        em[edg] = tg->addEdg(nm[edg->getFrom()], nm[edg->getTo()], edg->pl());
        tg->expandBBox(edg->pl().getGeom()->front());
        tg->expandBBox(edg->pl().getGeom()->back());

        tg->edgeRpl(em[edg]->getFrom(), edg, em[edg]);
        tg->edgeRpl(em[edg]->getTo(), edg, em[edg]);
        tg->nodeRpl(em[edg], edg->getTo(), nm[edg->getTo()]);
        tg->nodeRpl(em[edg], edg->getFrom(), nm[edg->getFrom()]);
      }
    }
  }

  if (offset) *offset = idOffset + geoComps.size() + 1;

  return ret;
}

// _____________________________________________________________________________
void LineGraph::snapOrphanStations() {
  double MAXD = 1;

  for (auto nd : getNds()) {
    if (nd->getDeg() != 0 || nd->pl().stops().size() == 0) continue;

    std::set<LineEdge*> cands;

    _edgeGrid.get(*nd->pl().getGeom(), MAXD * 2, &cands);

    for (auto e : cands) {
      if (util::geo::dist(*nd->pl().getGeom(), *e->pl().getGeom()) <= MAXD) {
        double pa =
            e->pl().getPolyline().projectOn(*nd->pl().getGeom()).totalPos;
        if (getEdg(e->getFrom(), nd) || getEdg(nd, e->getTo())) continue;

        assert(e->getFrom() != nd);
        assert(e->getTo() != nd);

        auto ba = addEdg(e->getFrom(), nd, e->pl());
        ba->pl().setPolyline(e->pl().getPolyline().getSegment(0, pa));

        auto bb = addEdg(nd, e->getTo(), e->pl());
        bb->pl().setPolyline(e->pl().getPolyline().getSegment(pa, 1));

        edgeRpl(e->getFrom(), e, ba);
        edgeRpl(e->getTo(), e, bb);

        nodeRpl(ba, e->getTo(), nd);
        nodeRpl(bb, e->getFrom(), nd);

        assert(ba != bb);

        _edgeGrid.add(*ba->pl().getGeom(), ba);
        _edgeGrid.add(*bb->pl().getGeom(), bb);

        _edgeGrid.remove(e);

        assert(getEdg(e->getFrom(), e->getTo()));
        delEdg(e->getFrom(), e->getTo());
      }
    }
  }
}

// _____________________________________________________________________________
void LineGraph::smooth(double smooth) {
  for (auto n : getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      auto pl = e->pl().getPolyline();
      pl.smoothenOutliers(50);
      pl.simplify(smooth);
      pl.applyChaikinSmooth(3);
      pl.simplify(1);
      e->pl().setPolyline(pl);
    }
  }
}

// _____________________________________________________________________________
void LineGraph::removeDeg1Nodes() {
  std::vector<LineNode*> toDel;
  for (auto n : getNds()) {
    if (n->getDeg() == 0) {
      for (auto e : n->getAdjList()) {
        _edgeGrid.remove(e);
      }

      toDel.push_back(n);
    }
  }

  for (auto n : toDel) {
    _nodeGrid.remove(n);
    delNd(n);
  }
}
