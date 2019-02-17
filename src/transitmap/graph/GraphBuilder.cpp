// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <proj_api.h>
#include <istream>
#include <set>
#include <stack>
#include <vector>
#include "transitmap/config/TransitMapConfig.h"
#include "GraphBuilder.h"
#include "json/json.hpp"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using namespace transitmapper;
using namespace graph;
using namespace util::geo;
using json = nlohmann::json;

const static char* WGS84_PROJ = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

// _____________________________________________________________________________
GraphBuilder::GraphBuilder(const config::Config* cfg) : _cfg(cfg) {
  _mercProj = pj_init_plus(WGS84_PROJ);
}

// _____________________________________________________________________________
bool GraphBuilder::build(std::istream* s, graph::TransitGraph* g) {
  json j;
  (*s) >> j;

  if (j["type"] == "FeatureCollection") {
    // first pass, nodes
    for (auto feature : j["features"]) {
      auto props = feature["properties"];
      auto geom = feature["geometry"];
      if (geom["type"] == "Point") {
        std::string id = props["id"];

        // check if node already exists
        if (g->getNodeById(id)) continue;

        std::vector<double> coords = geom["coordinates"];

        Node* n = new Node(id, coords[0], coords[1]);

        StationInfo i("", "");
        if (!props["station_id"].is_null() ||
            !props["station_label"].is_null()) {
          if (!props["station_id"].is_null()) i.id = props["station_id"];
          if (!props["station_label"].is_null())
            i.name = props["station_label"];
          n->addStop(i);
        }

        g->addNode(n);
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

        std::vector<std::vector<double> > coords = geom["coordinates"];

        PolyLine<double> pl;
        for (auto coord : coords) {
          double x = coord[0], y = coord[1];
          DPoint p(x, y);
          pl << p;
          g->expandBBox(p);
        }

        pl.applyChaikinSmooth(_cfg->inputSmoothing);

        Node* fromN = g->getNodeById(from);
        Node* toN = g->getNodeById(to);

        if (!fromN) {
          LOG(WARN) << "Node \"" << from << "\" not found.";
          continue;
        }

        if (!toN) {
          LOG(WARN) << "Node \"" << to << "\" not found.";
          continue;
        }

        if (dist(fromN->getPos(), pl.back()) <
            dist(fromN->getPos(), pl.front())) {
          LOG(WARN) << "Geometry for edge from " << fromN->getId() << " to "
                    << toN->getId() << " seems "
                    << " to have the wrong orientation! This may lead to "
                    << " strange results.";
        }

        Edge* e =
            g->addEdge(fromN, toN, pl, _cfg->lineWidth, _cfg->lineSpacing);

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

          const Route* r = g->getRoute(id);
          if (!r) {
            std::string label = route["label"].is_null() ? "" : route["label"];
            std::string color = route["color"];
            r = new Route(id, label, color);
            g->addRoute(r);
          }

          Node* dir = 0;

          if (!route["direction"].is_null()) {
            dir = g->getNodeById(route["direction"]);
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
        std::string id = props["id"];

        Node* n = g->getNodeById(id);

        if (!n) continue;

        if (!props["excluded_line_conns"].is_null()) {
          for (auto excl : props["excluded_line_conns"]) {
            std::string rid = excl["route"];
            std::string nid1 = excl["edge1_node"];
            std::string nid2 = excl["edge2_node"];

            const Route* r = g->getRoute(rid);

            if (!r) {
              LOG(WARN) << "line connection exclude defined in node " << id
                        << " for line " << rid << ", but no such line exists.";
              continue;
            }

            Node* n1 = g->getNodeById(nid1);
            Node* n2 = g->getNodeById(nid2);

            if (!n1) {
              LOG(WARN) << "line connection exclude defined in node " << id
                        << " for edge from " << nid1
                        << ", but no such node exists.";
              continue;
            }

            if (!n2) {
              LOG(WARN) << "line connection exclude defined in node " << id
                        << " for edge from " << nid2
                        << ", but no such node exists.";
              continue;
            }

            Edge* a = n->getEdge(n1);
            Edge* b = n->getEdge(n2);

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

            n->addRouteConnException(r, a, b);
          }
        }
      }
    }

  } else {
    LOG(ERROR) << "Could not read input.";
    return false;
  }

  return true;
}

// _____________________________________________________________________________
void GraphBuilder::writeMainDirs(TransitGraph* graph) {
  for (auto n : *graph->getNodes()) {
    std::set<Edge*> eSet;
    eSet.insert(n->getAdjListIn().begin(), n->getAdjListIn().end());
    eSet.insert(n->getAdjListOut().begin(), n->getAdjListOut().end());

    for (Edge* e : eSet) {
      NodeFront f(e, n);
      PolyLine<double> pl;

      f.refEtgLengthBefExp = e->getGeom().getLength();

      if (e->getTo() == n) {
        pl = e->getGeom().getOrthoLineAtDist(e->getGeom().getLength(),
                                             e->getTotalWidth());
      } else {
        pl = e->getGeom().getOrthoLineAtDist(0, e->getTotalWidth());
        pl.reverse();
      }

      f.setInitialGeom(pl);

      n->addMainDir(f);
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::writeStationGeoms(TransitGraph* graph) {
  for (auto n : *graph->getNodes()) {
    n->generateStationHull((_cfg->lineSpacing + _cfg->lineWidth) * 0.8,
                           _cfg->simpleRenderForTwoEdgeNodes);
  }
}

// _____________________________________________________________________________
void GraphBuilder::expandOverlappinFronts(TransitGraph* g) {
  // now, look at the nodes entire front geometries and expand them
  // until nothing overlaps
  double step = 1;

  while (true) {
    bool stillFree = false;
    for (auto n : *g->getNodes()) {
      if (n->getStops().size() && _cfg->tightStations) continue;
      std::set<NodeFront*> overlaps = nodeGetOverlappingFronts(n);
      for (auto f : overlaps) {
        stillFree = true;
        if (f->edge->getTo() == n) {
          f->geom = f->edge->getGeom().getOrthoLineAtDist(
              f->edge->getGeom().getLength() - step, f->edge->getTotalWidth());
        } else {
          f->geom = f->edge->getGeom().getOrthoLineAtDist(
              step, f->edge->getTotalWidth());
          f->geom.reverse();
        }

        // cut the edges to fit the new front
        freeNodeFront(f);
      }
    }
    if (!stillFree) break;
  }
}

// _____________________________________________________________________________
void GraphBuilder::createMetaNodes(TransitGraph* g) {
  std::vector<NodeFront> cands;
  while ((cands = getNextMetaNodeCand(g)).size() > 0) {
    // remove all edges completely contained
    for (auto nf : cands) {
      auto onfs = getClosedNodeFronts(nf.n);

      for (auto onf : onfs) {
        const Edge* e = onf.edge;

        bool found = false;

        for (auto nf : cands) {
          if (nf.n == e->getFrom()) {
            found = true;
            break;
          }
        }

        if (!found) continue;

        for (auto nf : cands) {
          if (nf.n == e->getTo()) {
            found = true;
            break;
          }
        }

        g->deleteEdge(e->getTo(), e->getFrom());
      }
    }

    // first node has new ref node id
    Node* ref = new Node(cands[0].n->getId(), cands[0].n->getPos());
    g->addNode(ref);

    for (auto nf : cands) {
      for (auto onf : getOpenNodeFronts(nf.n)) {
        if (onf.edge->getTo() == nf.n) {
          onf.edge->setTo(ref);
        } else {
          onf.edge->setFrom(ref);
        }

        for (auto& to : *onf.edge->getRoutes()) {
          if (to.direction == nf.n) {
            to.direction = ref;
          }
        }

        g->getNodes()->erase(onf.n);
        onf.n = ref;
        ref->addMainDir(onf);
        ref->addEdge(onf.edge);
      }
    }
  }
}

// _____________________________________________________________________________
std::vector<NodeFront> GraphBuilder::getNextMetaNodeCand(
    TransitGraph* g) const {
  for (auto n : *g->getNodes()) {
    if (n->getStops().size()) continue;
    if (getOpenNodeFronts(n).size() != 1) continue;

    std::set<const Node*> potClique;

    std::stack<const Node*> nodeStack;
    nodeStack.push(n);

    while (!nodeStack.empty()) {
      const Node* n = nodeStack.top();
      nodeStack.pop();

      if (n->getStops().size() == 0) {
        potClique.insert(n);
        for (auto nff : getClosedNodeFronts(n)) {
          const Node* m;

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

    if (isClique(potClique)) {
      std::vector<NodeFront> ret;

      for (auto n : potClique) {
        if (getOpenNodeFronts(n).size() > 0) {
          ret.push_back(getOpenNodeFronts(n)[0]);
        } else {
          for (auto nf : getClosedNodeFronts(n)) {
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
bool GraphBuilder::isClique(std::set<const Node*> potClique) const {
  if (potClique.size() < 2) return false;

  for (const Node* a : potClique) {
    for (const Node* b : potClique) {
      if (util::geo::dist(a->getPos(), b->getPos()) >
          (_cfg->lineWidth + _cfg->lineSpacing) * 10) {
        return false;
      }
    }
  }

  std::set<const Node*> periphery;

  for (const Node* n : potClique) {
    for (auto nf : getClosedNodeFronts(n)) {
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
    for (auto nf : getOpenNodeFronts(n)) {
      if (nf.edge->getTo() == n) {
        if (periphery.count(nf.edge->getFrom())) return false;
        periphery.insert(nf.edge->getFrom());
      } else {
        if (periphery.count(nf.edge->getTo())) return false;
        periphery.insert(nf.edge->getTo());
      }
    }


    for (auto nf : getOpenNodeFronts(n)) {
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
std::vector<NodeFront> GraphBuilder::getOpenNodeFronts(const Node* n) const {
  std::vector<NodeFront> res;
  for (auto nf : n->getMainDirs()) {
    if (nf.edge->getGeom().getLength() >
            (nf.edge->getWidth() + nf.edge->getSpacing()) ||
        (nf.edge->getOther(n)->getNodeFrontFor(nf.edge)->geom.distTo(
             nf.edge->getOther(n)->getPos()) >
         6 * (nf.edge->getWidth() + nf.edge->getSpacing())) ||
        (nf.edge->getTo()->getStops().size() > 0) ||
        (nf.edge->getFrom()->getStops().size() > 0)) {
      res.push_back(nf);
    }
  }

  return res;
}

// _____________________________________________________________________________
std::vector<NodeFront> GraphBuilder::getClosedNodeFronts(const Node* n) const {
  std::vector<NodeFront> res;
  for (auto nf : n->getMainDirs()) {
    if (!(nf.edge->getGeom().getLength() >
          (nf.edge->getWidth() + nf.edge->getSpacing())) &&
        !(nf.edge->getOther(n)->getNodeFrontFor(nf.edge)->geom.distTo(
              nf.edge->getOther(n)->getPos()) >
          6 * (nf.edge->getWidth() + nf.edge->getSpacing())) &&
        (nf.edge->getTo()->getStops().size() == 0) &&
        (nf.edge->getFrom()->getStops().size() == 0)) {
      res.push_back(nf);
    }
  }

  return res;
}

// _____________________________________________________________________________
std::set<NodeFront*> GraphBuilder::nodeGetOverlappingFronts(
    const Node* n) const {
  std::set<NodeFront*> ret;
  double minLength = 6;

  for (size_t i = 0; i < n->getMainDirs().size(); ++i) {
    const NodeFront& fa = n->getMainDirs()[i];

    for (size_t j = 0; j < n->getMainDirs().size(); ++j) {
      const NodeFront& fb = n->getMainDirs()[j];

      if (fa.geom.equals(fb.geom, 5) || j == i) continue;

      if (nodeFrontsOverlap(fa, fb)) {
        if (fa.edge->getGeom().getLength() > minLength &&
            fa.geom.distTo(n->getPos()) < 2.5 * n->getMaxNodeFrontWidth()) {
          ret.insert(const_cast<NodeFront*>(&fa));
        }
        if (fb.edge->getGeom().getLength() > minLength &&
            fb.geom.distTo(n->getPos()) < 2.5 * n->getMaxNodeFrontWidth()) {
          ret.insert(const_cast<NodeFront*>(&fb));
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
bool GraphBuilder::nodeFrontsOverlap(const NodeFront& a,
                                     const NodeFront& b) const {
  size_t numShr = a.edge->getShrdRoutes(*b.edge).size();

  DPoint aa = a.geom.front();
  DPoint ab = a.geom.back();
  DPoint ba = b.geom.front();
  DPoint bb = b.geom.back();

  bool intersects = lineIntersects(aa, ab, ba, bb);
  if (!intersects) return false;

  DPoint i = intersection(aa, ab, ba, bb);

  if (numShr && a.geom.distTo(i) < (a.edge->getWidth() + a.edge->getSpacing()))
    return true;
  if (b.geom.distTo(a.geom) <
      fmax((b.edge->getWidth() + b.edge->getSpacing()) * 3,
           (b.edge->getTotalWidth() + a.edge->getTotalWidth()) / 4))
    return true;

  return false;
}

// _____________________________________________________________________________
void GraphBuilder::freeNodeFront(NodeFront* f) {
  PolyLine<double> cutLine = f->geom;

  std::set<LinePoint<double>, LinePointCmp<double>> iSects =
      cutLine.getIntersections(f->edge->getGeom());
  if (iSects.size() > 0) {
    if (f->edge->getTo() != f->n) {
      // cut at beginning
      f->edge->setGeom(
          f->edge->getGeom().getSegment(iSects.begin()->totalPos, 1));

      assert(cutLine.distTo(f->edge->getGeom().front()) < 0.1);

    } else {
      // cut at end
      f->edge->setGeom(
          f->edge->getGeom().getSegment(0, (--iSects.end())->totalPos));

      assert(cutLine.distTo(f->edge->getGeom().back()) < 0.1);
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::combinePartnerRoutes(graph::TransitGraph* g) {
  auto partners = getPartnerRoutes(g);
  std::set<const Route*> processed;

  for (const auto& p : partners) {
    if (processed.find(p.first) != processed.end()) continue;
    if (p.second.size() == 0) continue;

    const_cast<Route*>(p.first)->setNumCollapsedPartners(p.second.size());

    for (const auto pp : p.second) {
      LOG(INFO) << "Combining " << pp << "(" << pp->getLabel() << ","
                << pp->getColor() << ")"
                << " and " << p.first << " (" << p.first->getLabel() << ","
                << p.first->getColor() << ").";
      processed.insert(pp);
      const_cast<Route*>(pp)->setRelativeTo(p.first);
    }
  }
}

// _____________________________________________________________________________
std::map<const Route*, std::set<const Route*> > GraphBuilder::getPartnerRoutes(
    graph::TransitGraph* g) const {
  std::map<const Route*, std::set<const Route*> > partners;

  for (graph::Node* n : *g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (const auto& ra : *(e->getRoutes())) {
        std::set<const Route*> partnersHere;

        for (const auto& rb : *(e->getRoutes())) {
          if (ra.direction == rb.direction) {
            partnersHere.insert(rb.route);
          }
        }

        if (partners.find(ra.route) == partners.end()) {
          partners[ra.route] = partnersHere;
        } else {
          std::set<const Route*> interSect;
          std::set_intersection(partners[ra.route].begin(),
                                partners[ra.route].end(), partnersHere.begin(),
                                partnersHere.end(),
                                std::inserter(interSect, interSect.begin()));
          partners[ra.route] = interSect;
        }
      }
    }
  }

  for (const auto& p : partners) {
    std::set<const Route*> realPartners;
    for (const auto& pp : p.second) {
      if (p.first != pp && partners[pp].find(p.first) != partners[pp].end()) {
        realPartners.insert(pp);
      }
    }

    partners[p.first] = realPartners;
  }

  return partners;
}

// _____________________________________________________________________________
void GraphBuilder::writeInitialConfig(TransitGraph* g) {
  OrderingConfig c;
  for (graph::Node* n : *g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      Ordering order(e->getCardinality());
      for (size_t i = 0; i < e->getCardinality(); i++) order[i] = i;
      c[e] = order;
    }
  }

  g->setConfig(c);
}
