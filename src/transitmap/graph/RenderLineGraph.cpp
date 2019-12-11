// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/graph/RenderLineGraph.h"
#include "util/geo/Geo.h"

// we still need boost geometry for the polygon buffering, which is no
// yet implemented in util/Geo.h
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>
namespace bgeo = boost::geometry;

using transitmapper::graph::RenderLineGraph;

// _____________________________________________________________________________
void RenderLineGraph::addFront(const TransitNode* nd, const Front& fr) {
  _fronts[nd].push_back(fr);

  struct Cmp {
    Cmp(const TransitNode* nd) : nd(nd){};
    bool operator()(const Front& a, const Front& b) {
      double ang1 = util::geo::angBetween(
          *nd->pl().getGeom(),
          util::geo::pointAtDist(*a.edg->pl().getGeom(), 5));
      double ang2 = util::geo::angBetween(
          *nd->pl().getGeom(),
          util::geo::pointAtDist(*b.edg->pl().getGeom(), 5));
      return (ang1 < ang2);
    }
    const TransitNode* nd;
  };

  // re-sort
  std::sort(_fronts[nd].begin(), _fronts[nd].end(), Cmp(nd));

  // TODO: re-sort after node front expansion
}

// _____________________________________________________________________________
void RenderLineGraph::writeStatGeoms() {
  for (auto nd : *getNds()) {
    _hulls[nd] = stationGeom(nd, (_cfg->lineSpacing + _cfg->lineWidth) * 0.8,
                             true, _cfg->simpleRenderForTwoEdgeNodes);
  }
}

// _____________________________________________________________________________
Polygon<double> RenderLineGraph::stationGeom(const TransitNode* node, double d,
                                             bool rectangulize,
                                             bool simple) const {
  auto fronts = _fronts.find(node)->second;
  double cd = d;

  typedef bgeo::model::point<double, 2, bgeo::cs::cartesian> BoostPoint;
  typedef bgeo::model::linestring<BoostPoint> BoostLine;
  typedef bgeo::model::polygon<BoostPoint> BoostPoly;
  typedef bgeo::model::multi_polygon<BoostPoly> BoostMultiPoly;

  BoostMultiPoly ret;
  double pointsPerCircle = 36;
  bgeo::strategy::buffer::distance_symmetric<double> distanceStrat(d);
  bgeo::strategy::buffer::join_round joinStrat(pointsPerCircle);
  bgeo::strategy::buffer::end_round endStrat(pointsPerCircle);
  bgeo::strategy::buffer::point_circle circleStrat(pointsPerCircle);
  bgeo::strategy::buffer::side_straight sideStrat;

  if (!simple || node->getDeg() == 2) {
    MultiLine<double> l;
    for (auto& nf : fronts) {
      l.push_back(util::geo::segment(nf.geom,
                      (cd / 2) / util::geo::len(nf.geom),
                                  (util::geo::len(nf.geom) - cd / 2) /
                                      util::geo::len(nf.geom)));
    }

    Polygon<double> hull = util::geo::convexHull(l);

    if (rectangulize && getMaxNodeFrontCardinality() > 1) {
      MultiLine<double> ll;
      for (auto& nf : fronts) {
        ll.push_back(nf.geom);
      }
      Polygon<double> env = util::geo::convexHull(
          util::geo::shrink(util::geo::getOrientedEnvelopeAvg(ll), cd / 2));

      double incr = (util::geo::area(env) / util::geo::area(hull)) - 1;
      if (ll.size() < 5 || incr < 0.5) {
        hull = env;
      }
    }

    BoostPoly hullBgeo;
    for (const auto& p : hull.getOuter())
      hullBgeo.outer().push_back({p.getX(), p.getY()});
    hullBgeo.outer().push_back(
        {hull.getOuter().front().getX(), hull.getOuter().front().getY()});

    // boost geometry expects polygons in clockwise fashion
    bgeo::correct(hullBgeo);

    bgeo::buffer(hullBgeo, ret, distanceStrat, sideStrat, joinStrat, endStrat,
                 circleStrat);

  } else {
    // for two main dirs, take average
    std::vector<const Line<double>*> pols;

    Line<double> a =
        util::geo::segment(fronts[0].geom, (cd / 2) / util::geo::len(fronts[0].geom),
                                  (util::geo::len(fronts[0].geom) - cd / 2) /
                                      util::geo::len(fronts[0].geom));

    Line<double> b =
        util::geo::segment(fronts[1].geom, (cd / 2) / util::geo::len(fronts[1].geom),
                                  (util::geo::len(fronts[1].geom) - cd / 2) /
                                      util::geo::len(fronts[1].geom));

    assert(a.size() > 1);
    assert(b.size() > 1);

    if (dist(a[0], b[0]) >
        dist(a[1], b[0])) {
      std::reverse(a.begin(), a.end());
    }

    pols.push_back(&a);
    pols.push_back(&b);

    auto avg = util::geo::average(pols);

    BoostLine lineBgeo;
    for (const auto& p : avg) {
      lineBgeo.push_back({p.getX(), p.getY()});
    }

    bgeo::buffer(lineBgeo, ret, distanceStrat, sideStrat, joinStrat, endStrat,
                 circleStrat);
  }

  assert(ret.size() > 0);

  Polygon<double> retPoly;
  for (const auto& p : ret[0].outer()) {
    retPoly.getOuter().push_back({p.get<0>(), p.get<1>()});
  }

  return retPoly;
}
