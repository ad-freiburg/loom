// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "util/geo/Geo.h"
#include "shared/linegraph/NodeFront.h"

using shared::linegraph::NodeFront;
using transitmapper::graph::OrderingConfig;
using  util::geo::DPoint;

// _____________________________________________________________________________
DPoint NodeFront::getTripOccPos(const Route* r, const OrderingConfig& c,
                                bool origGeom) const {
  assert(c.find(edge) != c.end());

  size_t p = edge->pl().getRoutePosUnder(r, c.find(edge)->second);
  return getTripPos(edge, p, n == edge->getTo(), origGeom);
}

// _____________________________________________________________________________
double NodeFront::getOutAngle() const {
  double checkDist = 10;
  if (edge->getFrom() == n) {
    return angBetween(*n->pl().getGeom(), PolyLine<double>(*edge->pl().getGeom()).getPointAtDist(checkDist).p);
  } else {
    return angBetween(
        *n->pl().getGeom(),
        PolyLine<double>(*edge->pl().getGeom())
            .getPointAtDist(util::geo::len(*edge->pl().getGeom()) - checkDist)
            .p);
  }
}

// _____________________________________________________________________________
DPoint NodeFront::getTripPos(const LineEdge* e, size_t pos, bool inv,
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

