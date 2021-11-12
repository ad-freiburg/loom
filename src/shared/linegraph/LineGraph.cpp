// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "dot/Parser.h"
#include "json/json.hpp"
#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineGraph.h"
#include "shared/linegraph/LineNodePL.h"
#include "shared/style/LineStyle.h"
#include "util/String.h"
#include "util/graph/Edge.h"
#include "util/graph/Node.h"
#include "util/log/Log.h"

using shared::linegraph::EdgeGrid;
using shared::linegraph::EdgeOrdering;
using shared::linegraph::Partner;
using shared::linegraph::ISect;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using shared::linegraph::LineOcc;
using shared::linegraph::NodeGrid;
using util::geo::DPoint;
using util::geo::Point;

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
          if (!props["station_id"].is_null()) {
            if (props["station_id"].type() == nlohmann::json::value_t::string) {
              i.id = props["station_id"].get<std::string>();
            } else {
              i.id = util::toString(props["station_id"]);
            }
          }
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

        if (props["dontcontract"].is_number() &&
            props["dontcontract"].get<int>())
          e->pl().setDontContract(true);

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
    // data structure, add here!)
    for (auto feature : j["features"]) {
      auto props = feature["properties"];
      auto geom = feature["geometry"];
      if (geom["type"] == "Point") {
        std::string id = util::toString(props["id"]);

        if (!idMap.count(id)) continue;
        LineNode* n = idMap[id];

        if (!props["not_serving"].is_null()) {
          for (auto excl : props["not_serving"]) {
            std::string lid = util::toString(excl);

            const Line* r = getLine(lid);

            if (!r) {
              LOG(WARN) << "line " << lid << " marked as not served in in node "
                        << id << ", but no such line exists.";
              continue;
            }

            n->pl().addLineNotServed(r);
          }
        }

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
  size_t gridSize =
      std::max(_bbox.getUpperRight().getX() - _bbox.getLowerLeft().getX(),
               _bbox.getUpperRight().getY() - _bbox.getLowerLeft().getY()) /
      10;

  _nodeGrid = NodeGrid(gridSize, gridSize, _bbox);
  _edgeGrid = EdgeGrid(gridSize, gridSize, _bbox);

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
  // TODO: prevent lines from continuing at the new intersection node!
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
std::set<LineEdge*> LineGraph::getNeighborEdges(const util::geo::DLine& line,
                                                double d) const {
  std::set<LineEdge*> neighbors;
  _edgeGrid.get(line, d, &neighbors);

  return neighbors;
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
const EdgeGrid& LineGraph::getEdgGrid() const { return &_edgeGrid; }

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
    auto cn = addNd(geom);

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
      nodeRpl(newEdg, n, cn);

      for (auto lo : eo.first->pl().getLines()) {
        ce->pl().addLine(lo.line, 0);
      }

      edgeRpl(n, eo.first, ce);
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
  for (auto n : *getNds()) {
    if (n->getDeg() > maxDeg) toSplit.push_back(n);
  }
  for (auto n : toSplit) splitNode(n, maxDeg);
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
    } else {
      ro++;
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
std::vector<Partner> LineGraph::getPartners(const LineNode* nd, const LineEdge* e,
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
void LineGraph::contractEdges(double d) {
  // TODO: the problem here is that contractEdge(e) may delete and replace edges
  // in the graph, which is why we cannot just build a list of edges eligible
  // for contraction and iterate over it - we would have to record the changes
  // made in contractEdge(e) and propagate it back.

breakfor:
  for (auto n1 : *getNds()) {
    for (auto e : n1->getAdjList()) {
      if (e->getFrom() != n1) continue;
      if (!e->pl().dontContract() && e->pl().getPolyline().getLength() < d) {
        if (e->getOtherNd(n1)->getAdjList().size() > 1 &&
            n1->getAdjList().size() > 1 &&
            (n1->pl().stops().size() == 0 ||
             e->getOtherNd(n1)->pl().stops().size() == 0)) {
          contractEdge(e);
          goto breakfor;
        }
      }
    }
  }
}

// _____________________________________________________________________________
double LineGraph::searchSpaceSize() const {
  double ret = 1;

  for (auto n :getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      ret *= util::factorial(e->pl().getLines().size());
    }
  }

  return ret;
}
