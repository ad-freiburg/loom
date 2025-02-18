// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

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

// _____________________________________________________________________________
Builder::Builder(const config::Config* cfg) : _cfg(cfg) {}

// _____________________________________________________________________________
void Builder::consume(const Feed& f, BuildGraph* g) {
  DBox graphBox(getProjP(f.getMinLat(), f.getMinLon()),
                getProjP(f.getMaxLat(), f.getMaxLon()));

  NodeGrid ngrid(2000, 2000, graphBox);

  size_t i = 0;

  for (auto t = f.getTrips().begin(); t != f.getTrips().end(); ++t) {
    i++;
    // ignore trips with only one stop
    if (t->second->getStopTimes().size() < 2) continue;
    if (!_cfg->useMots.count(t->second->getRoute()->getType())) continue;
    if (_cfg->useRoutes.size() > 0 && !_cfg->useRoutes.count(t->second->getRoute()->getId())) continue;

    auto st = t->second->getStopTimes().begin();

    auto prev = *st;
    const Edge* prevEdge = 0;
    addStop(prev.getStop(), g, &ngrid);
    ++st;

    if (i % 100 == 0)
      LOGTO(DEBUG, std::cerr) << "@ trip " << i << "/" << f.getTrips().size();

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
DPoint Builder::getProjP(double lat, double lng) const {
  return util::geo::latLngToWebMerc<double>(lat, lng);
}

// _____________________________________________________________________________
void Builder::simplify(BuildGraph* g) {
  // calculate average number of trip occurences per line
  double avg = 0;
  int c = 0;
  for (auto n : g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      for (auto& etg : *e->pl().getEdgeTripGeoms()) {
        for (auto& r : *etg.getTripsUnordered()) {
          avg += r.trips.size();
          c++;
        }
      }
    }
  }

  avg /= c;

  // try to merge both-direction edges into a single one
  // also prune edges with few trips
  for (auto n : g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      e->pl().simplify(avg * _cfg->pruneThreshold);
    }
  }

  // delete edges without a reference ETG
  std::vector<Edge*> toDel;
  for (auto n : g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (!e->pl().getRefETG()) toDel.push_back(e);
    }
  }

  for (auto e : toDel) g->delEdg(e->getFrom(), e->getTo());
}

// _____________________________________________________________________________
std::pair<bool, PolyLine<double>> Builder::getSubPolyLine(const Stop* a,
                                                          const Stop* b,
                                                          Trip* t, double distA,
                                                          double distB) {
  UNUSED(distA);
  UNUSED(distB);
  DPoint ap = getProjP(a->getLat(), a->getLng());
  DPoint bp = getProjP(b->getLat(), b->getLng());

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
      pl->second << getProjP(sp.lat, sp.lng);
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

  DPoint p = getProjP(curStop->getLat(), curStop->getLng());

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
