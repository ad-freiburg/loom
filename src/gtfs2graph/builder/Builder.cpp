// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <proj_api.h>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "gtfs2graph/builder/Builder.h"
#include "gtfs2graph/graph/BuildGraph.h"
#include "gtfs2graph/graph/EdgePL.h"
#include "gtfs2graph/graph/NodePL.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/log/Log.h"

using namespace gtfs2graph;
using namespace graph;

using util::geo::Box;
using util::geo::DBox;
using util::geo::DPoint;
using util::geo::extendBox;
using util::geo::Grid;
using util::geo::Point;
using util::geo::PolyLine;
using util::geo::SharedSegments;

using graph::Edge;
using graph::Node;

using ad::cppgtfs::gtfs::Feed;
using ad::cppgtfs::gtfs::Shape;
using ad::cppgtfs::gtfs::Stop;
using ad::cppgtfs::gtfs::StopTime;
using ad::cppgtfs::gtfs::Trip;

const static char* WGS84_PROJ =
    "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

// _____________________________________________________________________________
Builder::Builder(const config::Config* cfg) : _cfg(cfg) {
  _mercProj = pj_init_plus(WGS84_PROJ);
  _graphProj = pj_init_plus(cfg->projectionString.c_str());
}

// _____________________________________________________________________________
void Builder::consume(const Feed& f, BuildGraph* g) {
  DBox graphBox(getProjectedPoint(f.getMinLat(), f.getMinLon(), _graphProj),
                getProjectedPoint(f.getMaxLat(), f.getMaxLon(), _graphProj));

  NodeGrid ngrid(2000, 2000, graphBox);

  size_t i = 0;

  for (auto t = f.getTrips().begin(); t != f.getTrips().end(); ++t) {
    i++;
    // ignore trips with only one stop
    if (t->second->getStopTimes().size() < 2) continue;
    if (!_cfg->useMots.count(t->second->getRoute()->getType())) continue;

    auto st = t->second->getStopTimes().begin();

    auto prev = *st;
    const Edge* prevEdge = 0;
    addStop(prev.getStop(), g, &ngrid);
    ++st;

    if (i % 100 == 0)
      LOGTO(INFO, std::cerr) << "@ trip " << i << "/" << f.getTrips().size();

    for (; st != t->second->getStopTimes().end(); ++st) {
      const auto& cur = *st;

      Node* fromNode = getNodeByStop(g, prev.getStop());
      Node* toNode = addStop(cur.getStop(), g, &ngrid);

      // TODO: we should also allow this, for round-trips
      if (fromNode == toNode) continue;

      Edge* exE = g->getEdg(fromNode, toNode);

      if (!exE) {
        exE = g->addEdg(fromNode, toNode, EdgePL());
        exE->pl().setEdge(exE);
      }

      Node* directionNode = toNode;

      std::pair<bool, PolyLine<double>> edgeGeom;
      edgeGeom = getSubPolyLine(prev.getStop(), cur.getStop(), t->second,
                                prev.getShapeDistanceTravelled(),
                                cur.getShapeDistanceTravelled());

      if (prevEdge) {
        fromNode->pl().connOccurs(t->second->getRoute(), prevEdge, exE);
      }

      exE->pl().addTrip(t->second, edgeGeom.second, directionNode);

      prev = cur;
      prevEdge = exE;
    }
  }
}

// _____________________________________________________________________________
DPoint Builder::getProjectedPoint(double lat, double lng, projPJ p) const {
  double x = lng * DEG_TO_RAD, y = lat * DEG_TO_RAD;
  pj_transform(_mercProj, p, 1, 1, &x, &y, 0);
  return DPoint(x, y);
}

// _____________________________________________________________________________
void Builder::simplify(BuildGraph* g) {
  // delete edges without a reference ETG
  std::vector<Edge*> toDel;
  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (!e->pl().getRefETG()) toDel.push_back(e);
    }
  }

  for (auto e : toDel) g->delEdg(e->getFrom(), e->getTo());

  // try to merge both-direction edges into a single one
  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      e->pl().simplify();
    }
  }
}

// _____________________________________________________________________________
std::pair<bool, PolyLine<double>> Builder::getSubPolyLine(const Stop* a,
                                                          const Stop* b,
                                                          Trip* t, double distA,
                                                          double distB) {
  DPoint ap = getProjectedPoint(a->getLat(), a->getLng(), _graphProj);
  DPoint bp = getProjectedPoint(b->getLat(), b->getLng(), _graphProj);

  if (!t->getShape()) {
    return std::pair<bool, PolyLine<double>>(false, PolyLine<double>(ap, bp));
  }

  auto pl = _polyLines.find(t->getShape());
  if (pl == _polyLines.end()) {
    // generate polyline for this shape
    pl = _polyLines
             .insert(std::pair<Shape*, PolyLine<double>>(t->getShape(),
                                                         PolyLine<double>()))
             .first;

    for (const auto& sp : t->getShape()->getPoints()) {
      pl->second << getProjectedPoint(sp.lat, sp.lng, _graphProj);
    }
  }

  PolyLine<double> p;

  p = pl->second.getSegment(ap, bp);

  return std::pair<bool, PolyLine<double>>(true, p);
}

// _____________________________________________________________________________
Node* Builder::addStop(const Stop* curStop, BuildGraph* g, NodeGrid* grid) {
  Node* n = getNodeByStop(g, curStop);
  if (n) return n;

  DPoint p =
      getProjectedPoint(curStop->getLat(), curStop->getLng(), _graphProj);

  if (n) {
    n->pl().addStop(curStop);
    _stopNodes[curStop] = n;
  } else {
    n = g->addNd(NodePL(p, curStop));
    n->pl().setNode(n);
    grid->add(n->pl().getPos(), n);
    _stopNodes[curStop] = n;
  }

  return n;
}

// _____________________________________________________________________________
Node* Builder::getNodeByStop(const BuildGraph* g, const gtfs::Stop* s) const {
  if (_stopNodes.find(s) != _stopNodes.end()) return _stopNodes.find(s)->second;

  for (const auto n : g->getNds()) {
    if (n->pl().getStops().find(const_cast<gtfs::Stop*>(s)) !=
        n->pl().getStops().end()) {
      return n;
    }
  }
  return 0;
}
