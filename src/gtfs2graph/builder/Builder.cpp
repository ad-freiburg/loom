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

using util::geo::Point;
using util::geo::DPoint;
using util::geo::Grid;
using util::geo::Box;
using util::geo::DBox;
using util::geo::extendBox;
using util::geo::PolyLine;
using util::geo::SharedSegments;

using graph::Node;
using graph::Edge;

using ad::cppgtfs::gtfs::Stop;
using ad::cppgtfs::gtfs::StopTime;
using ad::cppgtfs::gtfs::Trip;
using ad::cppgtfs::gtfs::Shape;
using ad::cppgtfs::gtfs::Feed;

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

  bool SANITY_CHECK = true;

  for (auto t = f.getTrips().begin(); t != f.getTrips().end(); ++t) {
    if (t->second->getStopTimes().size() < 2) continue;
    if (!_cfg->useMots.count(t->second->getRoute()->getType())) continue;

    if (SANITY_CHECK && !checkTripSanity(t->second)) continue;

    auto st = t->second->getStopTimes().begin();

    auto prev = *st;
    const Edge* prevEdge = 0;
    addStop(prev.getStop(), _cfg->stationAggrLevel, g, &ngrid);
    ++st;

    for (; st != t->second->getStopTimes().end(); ++st) {
      const auto& cur = *st;

      Node* fromNode = getNodeByStop(g, prev.getStop(), _cfg->stationAggrLevel);
      Node* toNode = addStop(cur.getStop(), _cfg->stationAggrLevel, g, &ngrid);

      // TODO: we should also allow this, for round-trips
      if (fromNode == toNode) continue;

      Edge* exE = g->getEdg(fromNode, toNode);

      if (!exE) {
        exE = g->addEdg(fromNode, toNode, EdgePL());
        exE->pl().setEdge(exE);
      }

      Node* directionNode = toNode;

      if (!exE->pl().addTrip(t->second, directionNode)) {
        std::pair<bool, PolyLine<double>> edgeGeom;
        if (_cfg->stationAggrLevel) {
          const Stop* frs = prev.getStop()->getParentStation()
                                ? prev.getStop()->getParentStation()
                                : prev.getStop();
          const Stop* tos = cur.getStop()->getParentStation()
                                ? cur.getStop()->getParentStation()
                                : cur.getStop();
          edgeGeom = getSubPolyLine(frs, tos, t->second,
                                    prev.getShapeDistanceTravelled(),
                                    cur.getShapeDistanceTravelled());
        } else {
          edgeGeom = getSubPolyLine(prev.getStop(), cur.getStop(), t->second,
                                    prev.getShapeDistanceTravelled(),
                                    cur.getShapeDistanceTravelled());
        }

        // only take geometries that could be found using the
        // shape geometry, not some fallback
        if (!SANITY_CHECK || edgeGeom.first) {
          if (prevEdge) {
            fromNode->pl().connOccurs(t->second->getRoute(), prevEdge, exE);
          }

          exE->pl().addTrip(t->second, edgeGeom.second, directionNode);
        }
      }
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

  double totalTripDist = t->getShape()->getPoints().rbegin()->travelDist -
                         t->getShape()->getPoints().begin()->travelDist;

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

    pl->second.simplify(20);
    pl->second.smoothenOutliers(50);
    pl->second.fixTopology(50);
  }

  // if ((pl->second.distTo(ap) > 200) || (pl->second.distTo(bp) > 200)) {
    /**
     * something is not right, the distance from the station to its geometry
     * is excessive. fall back to straight line connection
     */
    // PolyLine<double> p = PolyLine<double>(ap, bp);
    // p.smoothenOutliers(50);
    // return std::pair<bool, PolyLine<double>>(false, p);
  // }

  PolyLine<double> p;

  if (!_cfg->ignoreGtfsDistances && distA > -1 && distA > -1 &&
      totalTripDist > 0) {
    p = pl->second.getSegment(
        (distA - t->getShape()->getPoints().begin()->travelDist) /
            totalTripDist,
        (distB - t->getShape()->getPoints().begin()->travelDist) /
            totalTripDist);
  } else {
    p = pl->second.getSegment(ap, bp);
  }

  return std::pair<bool, PolyLine<double>>(true, p);
}

// _____________________________________________________________________________
Node* Builder::addStop(const Stop* curStop, uint8_t aggrLevel, BuildGraph* g,
                       NodeGrid* grid) {
  if (aggrLevel && curStop->getParentStation() != 0) {
    Node* n = addStop(curStop->getParentStation(), aggrLevel, g, grid);
    _stopNodes[curStop] = n;
    return n;
  }

  Node* n = getNodeByStop(g, curStop, aggrLevel);
  if (n) return n;

  DPoint p =
      getProjectedPoint(curStop->getLat(), curStop->getLng(), _graphProj);

  if (aggrLevel > 1) {
    n = getNearestStop(p, _cfg->stationAggrDistance, grid);
  }

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
bool Builder::checkTripSanity(Trip* t) const {
  return checkShapeSanity(t->getShape());
}

// _____________________________________________________________________________
bool Builder::checkShapeSanity(Shape* s) const {
  if (!s || s->getPoints().size() < 2) return false;
  return true;
}

// _____________________________________________________________________________
Node* Builder::getNodeByStop(const BuildGraph* g, const gtfs::Stop* s,
                             bool getParent) const {
  if (getParent && s->getParentStation())
    return getNodeByStop(g, s->getParentStation());

  return getNodeByStop(g, s);
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

// _____________________________________________________________________________
Node* Builder::getNearestStop(const DPoint& p, double maxD,
                              const NodeGrid* grid) const {
  double curD = DBL_MAX;

  Node* curN = 0;
  std::set<Node*> neighbors;
  grid->get(p, maxD, &neighbors);

  for (auto n : neighbors) {
    double d = util::geo::dist(n->pl().getPos(), p);
    if (d < maxD && d < curD) {
      curN = n;
      curD = d;
    }
  }

  return curN;
}
