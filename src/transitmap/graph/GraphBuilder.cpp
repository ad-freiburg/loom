// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <proj_api.h>
#include <istream>
#include <set>
#include <stack>
#include <vector>
#include "GraphBuilder.h"
#include "json/json.hpp"
#include "shared/linegraph/Route.h"
#include "transitmap/config/TransitMapConfig.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using namespace util::geo;
using transitmapper::graph::GraphBuilder;
using shared::linegraph::NodeFront;
using shared::linegraph::LineNode;
using shared::linegraph::LineEdge;
using shared::linegraph::Route;

const static char* WGS84_PROJ =
    "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

// _____________________________________________________________________________
GraphBuilder::GraphBuilder(const config::Config* cfg) : _cfg(cfg) {
  _mercProj = pj_init_plus(WGS84_PROJ);
}

// _____________________________________________________________________________
bool GraphBuilder::build(std::istream* s, graph::TransitGraph* g) {
  // nlohmann::json j;
  // (*s) >> j;

  // if (j["type"] == "FeatureCollection") {
    // // first pass, nodes
    // for (auto feature : j["features"]) {
      // auto props = feature["properties"];
      // auto geom = feature["geometry"];
      // if (geom["type"] == "Point") {
        // std::string id = util::toString(props["id"]);

        // // check if node already exists
        // if (g->getNodeById(id)) continue;

        // std::vector<double> coords = geom["coordinates"];

        // Node* n = new Node(id, coords[0], coords[1]);

        // shared::linegraph::Station i("", "", n->getPos());
        // if (!props["station_id"].is_null() ||
            // !props["station_label"].is_null()) {
          // if (!props["station_id"].is_null())
            // i.id = util::toString(props["station_id"]);
          // if (!props["station_label"].is_null())
            // i.name = util::toString(props["station_label"]);
          // n->addStop(i);
        // }

        // g->addNd(n);
      // }
    // }

    // // second pass, edges
    // for (auto feature : j["features"]) {
      // auto props = feature["properties"];
      // auto geom = feature["geometry"];
      // if (geom["type"] == "LineString") {
        // if (props["lines"].is_null() || props["lines"].size() == 0) continue;
        // std::string from = util::toString(props["from"]);
        // std::string to = util::toString(props["to"]);

        // std::vector<std::vector<double>> coords = geom["coordinates"];

        // PolyLine<double> pl;
        // for (auto coord : coords) {
          // double x = coord[0], y = coord[1];
          // DPoint p(x, y);
          // pl << p;
          // g->expandBBox(p);
        // }

        // pl.applyChaikinSmooth(_cfg->inputSmoothing);

        // Node* fromN = g->getNodeById(from);
        // Node* toN = g->getNodeById(to);

        // if (!fromN) {
          // LOG(WARN) << "Node \"" << from << "\" not found.";
          // continue;
        // }

        // if (!toN) {
          // LOG(WARN) << "Node \"" << to << "\" not found.";
          // continue;
        // }

        // if (dist(fromN->getPos(), pl.back()) <
            // dist(fromN->getPos(), pl.front())) {
          // LOG(WARN) << "Geometry for edge from " << fromN->getId() << " to "
                    // << toN->getId() << " seems "
                    // << " to have the wrong orientation! This may lead to "
                    // << " strange results.";
        // }

        // Edge* e = g->addEdg(fromN, toN, pl, _cfg->lineWidth, _cfg->lineSpacing);

        // assert(e);
        // assert(g->getNodeById(from));
        // assert(g->getNodeById(to));

        // for (auto route : props["lines"]) {
          // std::string id;
          // if (!route["id"].is_null()) {
            // id = util::toString(route["id"]);
          // } else if (!route["label"].is_null()) {
            // id = util::toString(route["label"]);
          // } else if (!route["color"].is_null()) {
            // id = route["color"];
          // } else
            // continue;

          // const Route* r = g->getRoute(id);
          // if (!r) {
            // std::string label = route["label"].is_null() ? "" : route["label"];
            // std::string color = route["color"];
            // r = new Route(id, label, color);
            // g->addRoute(r);
          // }

          // Node* dir = 0;

          // if (!route["direction"].is_null()) {
            // dir = g->getNodeById(util::toString(route["direction"]));
            // if (!dir)
              // LOG(WARN) << "Direction " << route["direction"]
                        // << " defined in edge " << util::toString(props["id"])
                        // << " is not valid.";
          // }

          // if (!route["style"].is_null()) {
            // style::LineStyle ls;
            // auto style = route["style"];
            // std::string dashArray;
            // if (!style["dash-array"].is_null()) {
              // dashArray = style["dash-array"];
            // }

            // if (!style["css"].is_null()) {
              // ls.setCss(style["css"]);
            // }

            // ls.setDashArray(dashArray);

            // e->addRoute(r, dir, ls);
          // } else {
            // e->addRoute(r, dir);
          // }
        // }
      // }
    // }

    // // third pass, exceptions (TODO: do this in the first part, store in some
    // // data strcuture,
    // //  add here!)
    // for (auto feature : j["features"]) {
      // auto props = feature["properties"];
      // auto geom = feature["geometry"];
      // if (geom["type"] == "Point") {
        // std::string id = util::toString(props["id"]);

        // Node* n = g->getNodeById(id);

        // if (!n) continue;

        // if (!props["excluded_line_conns"].is_null()) {
          // for (auto excl : props["excluded_line_conns"]) {
            // std::string rid = util::toString(excl["route"]);
            // std::string nid1 = util::toString(excl["edge1_node"]);
            // std::string nid2 = util::toString(excl["edge2_node"]);

            // const Route* r = g->getRoute(rid);

            // if (!r) {
              // LOG(WARN) << "line connection exclude defined in node " << id
                        // << " for line " << rid << ", but no such line exists.";
              // continue;
            // }

            // Node* n1 = g->getNodeById(nid1);
            // Node* n2 = g->getNodeById(nid2);

            // if (!n1) {
              // LOG(WARN) << "line connection exclude defined in node " << id
                        // << " for edge from " << nid1
                        // << ", but no such node exists.";
              // continue;
            // }

            // if (!n2) {
              // LOG(WARN) << "line connection exclude defined in node " << id
                        // << " for edge from " << nid2
                        // << ", but no such node exists.";
              // continue;
            // }

            // Edge* a = n->getEdg(n1);
            // Edge* b = n->getEdg(n2);

            // if (!a) {
              // LOG(WARN) << "line connection exclude defined in node " << id
                        // << " for edge from " << nid1
                        // << ", but no such edge exists.";
              // continue;
            // }

            // if (!b) {
              // LOG(WARN) << "line connection exclude defined in node " << id
                        // << " for edge from " << nid2
                        // << ", but no such edge exists.";
              // continue;
            // }

            // n->addRouteConnException(r, a, b);
          // }
        // }
      // }
    // }

  // } else {
    // LOG(ERROR) << "Could not read input.";
    // return false;
  // }

  // return true;
}

// _____________________________________________________________________________
void GraphBuilder::writeMainDirs(TransitGraph* graph) {
  for (auto n : *graph->getNds()) {
    std::set<LineEdge*> eSet;
    eSet.insert(n->getAdjListIn().begin(), n->getAdjListIn().end());
    eSet.insert(n->getAdjListOut().begin(), n->getAdjListOut().end());

    for (LineEdge* e : eSet) {
      NodeFront f(e);
      PolyLine<double> pl;

      f.refEtgLengthBefExp = util::geo::len(*e->pl().getGeom());

      if (e->getTo() == n) {
        pl = PolyLine<double>(*e->pl().getGeom()).getOrthoLineAtDist(util::geo::len(*e->pl().getGeom()),
                                             graph->getTotalWidth(e));
      } else {
        pl = PolyLine<double>(*e->pl().getGeom()).getOrthoLineAtDist(0, graph->getTotalWidth(e));
        pl.reverse();
      }

      f.setInitialGeom(pl);

      n->pl().addMainDir(f);
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::expandOverlappinFronts(TransitGraph* g) {
  // now, look at the nodes entire front geometries and expand them
  // until nothing overlaps
  double step = 1;

  while (true) {
    bool stillFree = false;
    for (auto n : *g->getNds()) {
      if (n->pl().getStops().size() && _cfg->tightStations) continue;
      std::set<NodeFront*> overlaps = nodeGetOverlappingFronts(g, n);
      for (auto f : overlaps) {
        stillFree = true;
        if (f->edge->getTo() == n) {
          f->geom = geo::PolyLine<double>(*f->edge->pl().getGeom()).getOrthoLineAtDist(
              util::geo::len(*f->edge->pl().getGeom()) - step, g->getTotalWidth(f->edge));
        } else {
          f->geom = geo::PolyLine<double>(*f->edge->pl().getGeom()).getOrthoLineAtDist(
              step, g->getTotalWidth(f->edge));
          f->geom.reverse();
        }

        // cut the edges to fit the new front
        freeNodeFront(n, f);
      }
    }
    if (!stillFree) break;
  }
}

// _____________________________________________________________________________
void GraphBuilder::createMetaNodes(TransitGraph* g) {
  // std::vector<NodeFront> cands;
  // while ((cands = getNextMetaNodeCand(g)).size() > 0) {
    // // remove all edges completely contained
    // for (auto nf : cands) {
      // auto onfs = getClosedNodeFronts(nf.n);

      // for (auto onf : onfs) {
        // const LineEdge* e = onf.edge;

        // bool found = false;

        // for (auto nf : cands) {
          // if (nf.n == e->getFrom()) {
            // found = true;
            // break;
          // }
        // }

        // if (!found) continue;

        // for (auto nf : cands) {
          // if (nf.n == e->getTo()) {
            // found = true;
            // break;
          // }
        // }

        // g->delEdg(e->getTo(), e->getFrom());
      // }
    // }

    // // first node has new ref node id
    // LineNode* ref = g->addNd(cands[0].n->getId(), cands[0].n->getPos());

    // for (auto nf : cands) {
      // for (auto onf : getOpenNodeFronts(nf.n)) {
        // if (onf.edge->getTo() == nf.n) {
          // onf.edge->setTo(ref);
        // } else {
          // onf.edge->setFrom(ref);
        // }

        // for (auto& to : *onf.edge->getRoutes()) {
          // if (to.direction == nf.n) {
            // to.direction = ref;
          // }
        // }

        // g->getNds()->erase(onf.n);
        // onf.n = ref;
        // ref->addMainDir(onf);
        // ref->addEdg(onf.edge);
      // }
    // }
  // }
}

// _____________________________________________________________________________
std::vector<NodeFront> GraphBuilder::getNextMetaNodeCand(
    TransitGraph* g) const {
  for (auto n : *g->getNds()) {
    if (n->pl().getStops().size()) continue;
    if (getOpenNodeFronts(g, n).size() != 1) continue;

    std::set<const LineNode*> potClique;

    std::stack<const LineNode*> nodeStack;
    nodeStack.push(n);

    while (!nodeStack.empty()) {
      const LineNode* n = nodeStack.top();
      nodeStack.pop();

      if (n->pl().getStops().size() == 0) {
        potClique.insert(n);
        for (auto nff : getClosedNodeFronts(g, n)) {
          const LineNode* m;

          if (nff.edge->getTo() == n) {
            m = nff.edge->getFrom();
          } else {
            m = nff.edge->getTo();
          }

          if (potClique.find(m) == potClique.end()) {
            nodeStack.push(m);
          }
        }
      }
    }

    if (isClique(g, potClique)) {
      std::vector<NodeFront> ret;

      for (auto n : potClique) {
        if (getOpenNodeFronts(g, n).size() > 0) {
          ret.push_back(getOpenNodeFronts(g, n)[0]);
        } else {
          for (auto nf : getClosedNodeFronts(g, n)) {
            ret.push_back(nf);
          }
        }
      }

      return ret;
    }
  }

  return std::vector<NodeFront>();
}

// _____________________________________________________________________________
bool GraphBuilder::isClique(const TransitGraph* g, std::set<const LineNode*> potClique) const {
  if (potClique.size() < 2) return false;

  for (const LineNode* a : potClique) {
    for (const LineNode* b : potClique) {
      if (util::geo::dist(*a->pl().getGeom(), *b->pl().getGeom()) >
          (_cfg->lineWidth + _cfg->lineSpacing) * 10) {
        return false;
      }
    }
  }

  std::set<const LineNode*> periphery;

  for (const LineNode* n : potClique) {
    for (auto nf : getClosedNodeFronts(g, n)) {
      if (nf.edge->getTo() == n) {
        if (potClique.find(nf.edge->getFrom()) == potClique.end()) {
          return false;
        }
      } else {
        if (potClique.find(nf.edge->getTo()) == potClique.end()) {
          return false;
        }
      }
    }

    // catch cases where two clique nodes share the same node at teir
    // open front
    for (auto nf : getOpenNodeFronts(g, n)) {
      if (nf.edge->getTo() == n) {
        if (periphery.count(nf.edge->getFrom())) return false;
        periphery.insert(nf.edge->getFrom());
      } else {
        if (periphery.count(nf.edge->getTo())) return false;
        periphery.insert(nf.edge->getTo());
      }
    }

    for (auto nf : getOpenNodeFronts(g, n)) {
      if (nf.edge->getTo() == n) {
        if (potClique.find(nf.edge->getFrom()) != potClique.end()) {
          return false;
        }
      } else {
        if (potClique.find(nf.edge->getTo()) != potClique.end()) {
          return false;
        }
      }
    }
  }

  return true;
}

// _____________________________________________________________________________
std::vector<NodeFront> GraphBuilder::getOpenNodeFronts(const graph::TransitGraph* g, const LineNode* n) const {
  std::vector<NodeFront> res;
  for (auto nf : n->pl().getMainDirs()) {
    if (util::geo::len(*nf.edge->pl().getGeom()) >
            (g->getWidth(nf.edge) + g->getSpacing(nf.edge)) ||
        (nf.edge->getOtherNd(n)->pl().getNodeFrontFor(nf.edge)->geom.distTo(
             *nf.edge->getOtherNd(n)->pl().getGeom()) >
         6 * (g->getWidth(nf.edge) + g->getSpacing(nf.edge))) ||
        (nf.edge->getTo()->pl().getStops().size() > 0) ||
        (nf.edge->getFrom()->pl().getStops().size() > 0)) {
      res.push_back(nf);
    }
  }

  return res;
}

// _____________________________________________________________________________
std::vector<NodeFront> GraphBuilder::getClosedNodeFronts(const graph::TransitGraph* g, const LineNode* n) const {
  std::vector<NodeFront> res;
  for (auto nf : n->pl().getMainDirs()) {
    if (!(util::geo::len(*nf.edge->pl().getGeom()) >
          (g->getWidth(nf.edge) + g->getSpacing(nf.edge))) &&
        !(nf.edge->getOtherNd(n)->pl().getNodeFrontFor(nf.edge)->geom.distTo(
              *nf.edge->getOtherNd(n)->pl().getGeom()) >
          6 * (g->getWidth(nf.edge) + g->getSpacing(nf.edge))) &&
        (nf.edge->getTo()->pl().getStops().size() == 0) &&
        (nf.edge->getFrom()->pl().getStops().size() == 0)) {
      res.push_back(nf);
    }
  }

  return res;
}

// _____________________________________________________________________________
std::set<NodeFront*> GraphBuilder::nodeGetOverlappingFronts(const TransitGraph* g, 
    const LineNode* n) const {
  std::set<NodeFront*> ret;
  double minLength = 6;

  for (size_t i = 0; i < n->pl().getMainDirs().size(); ++i) {
    const NodeFront& fa = n->pl().getMainDirs()[i];

    for (size_t j = 0; j < n->pl().getMainDirs().size(); ++j) {
      const NodeFront& fb = n->pl().getMainDirs()[j];

      if (fa.geom.equals(fb.geom, 5) || j == i) continue;

      if (nodeFrontsOverlap(g, fa, fb)) {
        if (util::geo::len(*fa.edge->pl().getGeom()) > minLength &&
            fa.geom.distTo(*n->pl().getGeom()) < 2.5 * g->getMaxNodeFrontWidth(n)) {
          ret.insert(const_cast<NodeFront*>(&fa));
        }
        if (util::geo::len(*fb.edge->pl().getGeom()) > minLength &&
            fb.geom.distTo(*n->pl().getGeom()) < 2.5 * g->getMaxNodeFrontWidth(n)) {
          ret.insert(const_cast<NodeFront*>(&fb));
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
bool GraphBuilder::nodeFrontsOverlap(const TransitGraph* g, const NodeFront& a,
                                     const NodeFront& b) const {
  size_t numShr = g->getSharedRoutes(a.edge, b.edge).size();

  DPoint aa = a.geom.front();
  DPoint ab = a.geom.back();
  DPoint ba = b.geom.front();
  DPoint bb = b.geom.back();

  bool intersects = lineIntersects(aa, ab, ba, bb);
  if (!intersects) return false;

  DPoint i = intersection(aa, ab, ba, bb);

  if (numShr && a.geom.distTo(i) < (g->getWidth(a.edge) + g->getSpacing(a.edge)))
    return true;
  if (b.geom.distTo(a.geom) <
      fmax((g->getWidth(b.edge) + g->getSpacing(b.edge)) * 3,
           (g->getTotalWidth(b.edge) + g->getTotalWidth(a.edge)) / 4))
    return true;

  return false;
}

// _____________________________________________________________________________
void GraphBuilder::freeNodeFront(const LineNode* n, NodeFront* f) {
  PolyLine<double> cutLine = f->geom;

  std::set<LinePoint<double>, LinePointCmp<double>> iSects =
      cutLine.getIntersections(*f->edge->pl().getGeom());
  if (iSects.size() > 0) {
    if (f->edge->getTo() != n) {
      // cut at beginning
      f->edge->pl().setGeom(
          geo::PolyLine<double>(*f->edge->pl().getGeom()).getSegment(iSects.begin()->totalPos, 1).getLine());

      assert(cutLine.distTo(f->edge->pl().getGeom()->front()) < 0.1);

    } else {
      // cut at end
      f->edge->pl().setGeom(
          geo::PolyLine<double>(*f->edge->pl().getGeom()).getSegment(0, (--iSects.end())->totalPos).getLine());

      assert(cutLine.distTo(f->edge->pl().getGeom()->back()) < 0.1);
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::writeInitialConfig(TransitGraph* g) {
  OrderingConfig c;
  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjListOut()) {
      Ordering order(e->pl().getRoutes().size());
      for (size_t i = 0; i < e->pl().getRoutes().size(); i++) order[i] = i;
      c[e] = order;
    }
  }

  g->setConfig(c);
}
