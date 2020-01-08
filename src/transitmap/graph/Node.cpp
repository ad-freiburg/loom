// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include <algorithm>
#include "transitmap/graph/Edge.h"
#include "transitmap/graph/Node.h"
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/Penalties.h"
#include "shared/linegraph/Route.h"
#include "transitmap/graph/TransitGraph.h"
#include "util/geo/BezierCurve.h"
#include "util/geo/Geo.h"

using transitmapper::graph::NodeFront;
using transitmapper::graph::Node;
using transitmapper::graph::Edge;
using transitmapper::graph::Partner;
using transitmapper::graph::InnerGeometry;
using shared::linegraph::Route;

using util::geo::Point;
using util::geo::Polygon;
using util::geo::BezierCurve;
using util::geo::PolyLine;
using util::geo::Line;

// _____________________________________________________________________________
DPoint NodeFront::getTripOccPos(const Route* r, const OrderingConfig& c) const {
  return getTripOccPos(r, c, false);
}

// _____________________________________________________________________________
DPoint NodeFront::getTripOccPos(const Route* r, const OrderingConfig& c,
                                bool origGeom) const {
  assert(c.find(edge) != c.end());

  size_t p = edge->getRoutePosUnder(r, c.find(edge)->second);
  return getTripPos(edge, p, n == edge->getTo(), origGeom);
}

// _____________________________________________________________________________
DPoint NodeFront::getTripPos(const Edge* e, size_t pos, bool inv) const {
  return getTripPos(e, pos, inv, false);
}

// _____________________________________________________________________________
double NodeFront::getOutAngle() const {
  double checkDist = 10;
  if (edge->getFrom() == n) {
    return angBetween(n->getPos(), edge->getGeom().getPointAtDist(checkDist).p);
  } else {
    return angBetween(
        n->getPos(),
        edge->getGeom()
            .getPointAtDist(edge->getGeom().getLength() - checkDist)
            .p);
  }
}

// _____________________________________________________________________________
DPoint NodeFront::getTripPos(const Edge* e, size_t pos, bool inv,
                             bool origG) const {
  double p;
  if (!inv) {
    p = (e->getWidth() + e->getSpacing()) * pos + e->getWidth() / 2;
  } else {
    p = (e->getWidth() + e->getSpacing()) * (e->getCardinality() - 1 - pos) +
        e->getWidth() / 2;
  }
  // use interpolate here directly for speed
  if (origG) {
    return origGeom.interpolate(origGeom.getLine().front(),
                                origGeom.getLine().back(), p);
  } else {
    return geom.interpolate(geom.getLine().front(), geom.getLine().back(), p);
  }
}

// _____________________________________________________________________________
double Node::getMaxNodeFrontWidth() const {
  double ret = 0;
  for (const NodeFront& g : _mainDirs) {
    if (g.edge->getTotalWidth() > ret) ret = g.edge->getTotalWidth();
  }
  return ret;
}

// _____________________________________________________________________________
size_t Node::getMaxNodeFrontCardinality() const {
  size_t ret = 0;
  for (const NodeFront& g : _mainDirs) {
    if (g.edge->getCardinality() > ret) ret = g.edge->getCardinality();
  }
  return ret;
}

// _____________________________________________________________________________
Node::Node(const std::string& id, DPoint pos) : _id(id), _pos(pos) {}

// _____________________________________________________________________________
Node::Node(const std::string& id, double x, double y) : _id(id), _pos(x, y) {}

// _____________________________________________________________________________
Node::Node(const std::string& id, DPoint pos, shared::linegraph::Station s)
    : _id(id), _pos(pos) {
  addStop(s);
}

// _____________________________________________________________________________
Node::Node(const std::string& id, double x, double y, shared::linegraph::Station s)
    : _id(id), _pos(x, y) {
  addStop(s);
}

// _____________________________________________________________________________
Node::~Node() {
  for (auto e = _adjListOut.begin(); e != _adjListOut.end();) {
    Edge* eP = *e;

    if (eP->getFrom() == this) {
      // careful with invalidating iterators
      e = _adjListOut.erase(e);
    } else {
      eP->getFrom()->removeEdge(eP);
      e++;
    }

    eP->getTo()->removeEdge(eP);

    delete eP;
  }

  for (auto e = _adjListIn.begin(); e != _adjListIn.end();) {
    Edge* eP = *e;

    if (eP->getTo() == this) {
      // careful with invalidating iterators
      e = _adjListIn.erase(e);
    } else {
      eP->getTo()->removeEdge(eP);
      e++;
    }

    eP->getFrom()->removeEdge(eP);

    delete eP;
  }
}

// _____________________________________________________________________________
void Node::addStop(shared::linegraph::Station s) { _stops.push_back(s); }

// _____________________________________________________________________________
const std::vector<shared::linegraph::Station>& Node::getStops() const { return _stops; }

// _____________________________________________________________________________
void Node::addEdg(Edge* e) {
  if (e->getFrom() == this) _adjListOut.insert(e);
  if (e->getTo() == this) _adjListIn.insert(e);
}

// _____________________________________________________________________________
void Node::removeEdge(Edge* e) {
  if (e->getFrom() == this) _adjListOut.erase(e);
  if (e->getTo() == this) _adjListIn.erase(e);

  for (size_t i = 0; i < _mainDirs.size(); i++) {
    if (_mainDirs[i].edge == e) {
      _mainDirs.erase(_mainDirs.begin() + i);
    }
  }

  // TODO: remove from _routeConnExceptions
}

// _____________________________________________________________________________
const DPoint& Node::getPos() const { return _pos; }

// _____________________________________________________________________________
void Node::setPos(const DPoint& p) { _pos = p; }

// _____________________________________________________________________________
const std::string& Node::getId() const { return _id; }

// _____________________________________________________________________________
void Node::addMainDir(NodeFront f) { _mainDirs.push_back(f); }

// _____________________________________________________________________________
const NodeFront* Node::getNodeFrontFor(const Edge* e) const {
  for (auto& nf : getMainDirs()) {
    if (nf.edge == e) {
      return &nf;
    }
  }

  return 0;
}

// _____________________________________________________________________________
std::set<Edge*> Node::getAdjList() const {
  std::set<Edge*> ret;
  ret.insert(getAdjListIn().begin(), getAdjListIn().end());
  ret.insert(getAdjListOut().begin(), getAdjListOut().end());

  return ret;
}

// _____________________________________________________________________________
std::vector<Partner> Node::getPartners(const NodeFront* f,
                                       const RouteOccurance& ro) const {
  std::vector<Partner> ret;
  for (const auto& nf : getMainDirs()) {
    if (&nf == f) continue;

    for (const RouteOccurance& to :
         nf.edge->getCtdRoutesIn(this, ro.route, ro.direction, f->edge)) {
      Partner p(f, nf.edge, to.route);
      p.front = &nf;
      p.edge = nf.edge;
      p.route = to.route;
      ret.push_back(p);
    }
  }
  return ret;
}


// _____________________________________________________________________________
void Node::addRouteConnException(const Route* r, const Edge* edgeA,
                                 const Edge* edgeB) {
  _routeConnExceptions[r][edgeA].insert(edgeB);
  // index the other direction also, will lead to faster lookups later on
  _routeConnExceptions[r][edgeB].insert(edgeA);
}

// _____________________________________________________________________________
bool Node::connOccurs(const Route* r, const Edge* edgeA,
                      const Edge* edgeB) const {
  const auto& i = _routeConnExceptions.find(r);
  if (i == _routeConnExceptions.end()) return true;

  const auto& ii = i->second.find(edgeA);
  if (ii == i->second.end()) return true;

  return ii->second.find(edgeB) == ii->second.end();
}

// _____________________________________________________________________________
Edge* Node::getEdg(const Node* other) const {
  for (auto e = _adjListOut.begin(); e != _adjListOut.end(); ++e) {
    Edge* eP = *e;

    if (eP->getTo() == other) {
      return eP;
    }
  }

  for (auto e = _adjListIn.begin(); e != _adjListIn.end(); ++e) {
    Edge* eP = *e;

    if (eP->getFrom() == other) {
      return eP;
    }
  }

  return 0;
}

// _____________________________________________________________________________
size_t Node::getConnCardinality() const {
  size_t ret = 0;
  std::map<const Route*, std::set<const NodeFront*>> processed;

  for (size_t i = 0; i < getMainDirs().size(); ++i) {
    const NodeFront& nf = getMainDirs()[i];

    for (size_t j = 0; j < nf.edge->getCardinality(); j++) {
      const RouteOccurance& routeOcc = (*nf.edge->getRoutes())[j];

      std::vector<Partner> partners = getPartners(&nf, routeOcc);

      for (const Partner& p : partners) {
        if (processed[routeOcc.route].find(p.front) !=
            processed[routeOcc.route].end()) {
          continue;
        }
        ret++;
      }

      processed[routeOcc.route].insert(&nf);
    }
  }

  return ret;
}

