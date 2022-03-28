// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <set>
#include <string>
#include "3rdparty/json.hpp"
#include "shared/linegraph/Line.h"
#include "shared/rendergraph/OrderCfg.h"
#include "shared/rendergraph/RenderGraph.h"
#include "util/Misc.h"
#include "util/geo/BezierCurve.h"
#include "util/geo/Geo.h"

using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineNode;
using shared::linegraph::LineOcc;
using shared::linegraph::NodeFront;
using shared::linegraph::Partner;
using shared::rendergraph::InnerGeom;
using shared::rendergraph::OrderCfg;
using shared::rendergraph::RenderGraph;
using util::geo::BezierCurve;
using util::geo::dist;
using util::geo::DPoint;
using util::geo::MultiLine;
using util::geo::Polygon;
using util::geo::PolyLine;

// _____________________________________________________________________________
void RenderGraph::smooth() {
  for (auto n : getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      auto pl = e->pl().getPolyline();
      pl.smoothenOutliers(50);
      pl.simplify(1);
      pl.applyChaikinSmooth(1);
      pl.simplify(1);
      e->pl().setPolyline(pl);
    }
  }
}

// _____________________________________________________________________________
void RenderGraph::writePermutation(const OrderCfg& c) {
  for (auto n : getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      e->pl().writePermutation(c.at(e));
    }
  }
}

// _____________________________________________________________________________
bool RenderGraph::isTerminus(const LineNode* n) {
  if (n->getDeg() == 1) return true;
  for (size_t i = 0; i < n->pl().fronts().size(); ++i) {
    const NodeFront& nf = n->pl().fronts()[i];

    for (size_t p = 0; p < nf.edge->pl().getLines().size(); p++) {
      const LineOcc& lineOcc = nf.edge->pl().lineOccAtPos(p);
      std::vector<Partner> partners = getPartners(n, nf.edge, lineOcc);
      if (partners.size() == 0) return true;
    }
  }
  return false;
}

// _____________________________________________________________________________
std::vector<InnerGeom> RenderGraph::innerGeoms(const LineNode* n,
                                               double prec) const {
  std::vector<InnerGeom> ret;
  std::map<const Line*, std::set<const LineEdge*>> processed;

  for (size_t i = 0; i < n->pl().fronts().size(); ++i) {
    const NodeFront& nf = n->pl().fronts()[i];

    for (size_t j = 0; j < nf.edge->pl().getLines().size(); j++) {
      const LineOcc& lineOcc = nf.edge->pl().lineOccAtPos(j);
      Partner o(nf.edge, lineOcc.line);

      std::vector<Partner> partners = getPartners(n, nf.edge, lineOcc);

      for (const Partner& p : partners) {
        if (processed[lineOcc.line].find(p.edge) !=
            processed[lineOcc.line].end()) {
          continue;
        }

        auto is = getInnerLine(n, o, p);

        auto nfo = n->pl().frontFor(p.edge);

        double dmax = std::max(
            util::geo::dist(nfo->geom.getLine(), nf.geom.getLine().front()),
            util::geo::dist(nfo->geom.getLine(), nf.geom.getLine().back()));

        if (prec > 0 && is.geom.getLength() > 0 && dmax > 5) {
          ret.push_back(getInnerBezier(n, o, p, prec));
        } else {
          ret.push_back(is);
        }
      }

      // handle lines where this node is the terminus
      if (partners.size() == 0 && !notCompletelyServed(n)) {
        auto is = getTerminusLine(n, o);
        if (is.geom.getLength() > 0) {
          if (prec > 0) {
            ret.push_back(getTerminusBezier(n, o, prec));
          } else {
            ret.push_back(is);
          }
        }
      }

      processed[lineOcc.line].insert(nf.edge);
    }
  }

  return ret;
}

// _____________________________________________________________________________
InnerGeom RenderGraph::getInnerBezier(const LineNode* n,
                                      const Partner& partnerFrom,
                                      const Partner& partnerTo,
                                      double prec) const {
  double EPSI = 0.001;
  InnerGeom ret = getInnerLine(n, partnerFrom, partnerTo);
  DPoint p = ret.geom.getLine().front();
  DPoint q = ret.geom.getLine().back();
  double d = util::geo::dist(p, q);

  std::pair<double, double> slpA, slpB;

  if (util::geo::len(*partnerFrom.edge->pl().getGeom()) <= 5) return ret;
  if (util::geo::len(*partnerTo.edge->pl().getGeom()) <= 5) return ret;

  if (d <= 5) return ret;

  if (partnerFrom.edge->getTo() == n) {
    slpA = PolyLine<double>(*partnerFrom.edge->pl().getGeom())
               .getSlopeBetweenDists(
                   util::geo::len(*partnerFrom.edge->pl().getGeom()) - 5,
                   util::geo::len(*partnerFrom.edge->pl().getGeom()));
  } else {
    slpA = PolyLine<double>(*partnerFrom.edge->pl().getGeom())
               .getSlopeBetweenDists(5, 0);
  }

  if (partnerTo.edge->getTo() == n) {
    slpB = PolyLine<double>(*partnerTo.edge->pl().getGeom())
               .getSlopeBetweenDists(
                   util::geo::len(*partnerTo.edge->pl().getGeom()) - 5,
                   util::geo::len(*partnerTo.edge->pl().getGeom()));
  } else {
    slpB = PolyLine<double>(*partnerTo.edge->pl().getGeom())
               .getSlopeBetweenDists(5, 0);
  }

  double da = d * 0.55;

  DPoint pa =
      DPoint(p.getX() + slpA.first * d * 2, p.getY() + slpA.second * d * 2);
  DPoint ppa =
      DPoint(q.getX() + slpB.first * d * 2, q.getY() + slpB.second * d * 2);

  if (d > EPSI && util::geo::intersects(p, pa, q, ppa)) {
    DPoint isect = util::geo::intersection(p, pa, q, ppa);

    if (!std::isnan(isect.getX()) && !std::isnan(isect.getY())) {
      double dar = util::geo::dist(isect, p);
      double dbr = util::geo::dist(isect, q);

      double avg = (dar + dbr) / 2;

      da = avg * 0.55;
    }
  }

  if ((2 * da) < EPSI || (2 * da) > d * 2) {
    da = d * 0.55;
  }

  auto pp = DPoint(p.getX() + slpA.first * (da), p.getY() + slpA.second * (da));
  auto qq = DPoint(q.getX() + slpB.first * (da), q.getY() + slpB.second * (da));

  BezierCurve<double> bc(p, pp, qq, q);
  ret.geom = bc.render(prec);

  return ret;
}

// _____________________________________________________________________________
InnerGeom RenderGraph::getTerminusLine(const LineNode* n,
                                       const Partner& partnerFrom) const {
  auto nf = n->pl().frontFor(partnerFrom.edge);
  DPoint p = linePosOn(*nf, partnerFrom.line, false);
  DPoint pp = linePosOn(*nf, partnerFrom.line, true);

  size_t s = partnerFrom.edge->pl().linePos(partnerFrom.line);
  size_t ss = partnerFrom.edge->pl().linePos(partnerFrom.line);

  return InnerGeom(PolyLine<double>(p, pp), partnerFrom, Partner(), s, ss);
}

// _____________________________________________________________________________
InnerGeom RenderGraph::getInnerLine(const LineNode* n,
                                    const Partner& partnerFrom,
                                    const Partner& partnerTo) const {
  auto nfFr = n->pl().frontFor(partnerFrom.edge);
  auto nfTo = n->pl().frontFor(partnerTo.edge);
  DPoint p = linePosOn(*nfFr, partnerFrom.line, false);
  DPoint pp = linePosOn(*nfTo, partnerTo.line, false);

  size_t s = partnerFrom.edge->pl().linePos(partnerFrom.line);
  size_t ss = partnerTo.edge->pl().linePos(partnerTo.line);

  return InnerGeom(PolyLine<double>(p, pp), partnerFrom, partnerTo, s, ss);
}

// _____________________________________________________________________________
InnerGeom RenderGraph::getTerminusBezier(const LineNode* n,
                                         const Partner& partnerFrom,
                                         double prec) const {
  InnerGeom ret = getTerminusLine(n, partnerFrom);
  DPoint p = ret.geom.getLine().front();
  DPoint pp = ret.geom.getLine().back();
  double d = util::geo::dist(p, pp) / 2;

  DPoint b = p;
  DPoint c = pp;
  std::pair<double, double> slopeA;

  if (util::geo::len(*partnerFrom.edge->pl().getGeom()) <= 5) return ret;

  if (partnerFrom.edge->getTo() == n) {
    slopeA = PolyLine<double>(*partnerFrom.edge->pl().getGeom())
                 .getSlopeBetweenDists(
                     util::geo::len(*partnerFrom.edge->pl().getGeom()) - 5,
                     util::geo::len(*partnerFrom.edge->pl().getGeom()));
  } else {
    slopeA = PolyLine<double>(*partnerFrom.edge->pl().getGeom())
                 .getSlopeBetweenDists(5, 0);
  }

  b = DPoint(p.getX() + slopeA.first * d, p.getY() + slopeA.second * d);
  c = DPoint(pp.getX(), pp.getY());

  BezierCurve<double> bc(p, b, c, pp);
  ret.geom = bc.render(prec);

  return ret;
}

// _____________________________________________________________________________
Polygon<double> RenderGraph::getConvexFrontHull(const LineNode* n, double d,
                                                bool rectangulize, bool tight,
                                                size_t pointsPerCircle) const {
  double cd = d;

  if (n->pl().fronts().size() != 2) {
    MultiLine<double> l;
    for (auto& nf : n->pl().fronts()) {
      l.push_back(
          nf.geom
              .getSegment((cd / 2) / nf.geom.getLength(),
                          (nf.geom.getLength() - cd / 2) / nf.geom.getLength())
              .getLine());
    }

    Polygon<double> hull = util::geo::convexHull(l);

    if (rectangulize) {
      double step = tight ? 45 : 1;
      MultiLine<double> ll;
      for (auto& nf : n->pl().fronts()) ll.push_back(nf.geom.getLine());
      Polygon<double> env = util::geo::convexHull(
          util::geo::shrink(util::geo::getOrientedEnvelope(ll, step), cd / 2));

      double incr = (util::geo::area(env) / util::geo::area(hull)) - 1;
      if (ll.size() < 5 || incr < 0.5) {
        hull = env;
      }
    }

    return util::geo::buffer(hull, d, pointsPerCircle);

  } else {
    // for two node fronts, take average
    std::vector<const PolyLine<double>*> pols;

    PolyLine<double> a = n->pl().fronts()[0].origGeom.getSegment(
        (cd / 2) / n->pl().fronts()[0].origGeom.getLength(),
        (n->pl().fronts()[0].origGeom.getLength() - cd / 2) /
            n->pl().fronts()[0].origGeom.getLength());

    PolyLine<double> b = n->pl().fronts()[1].origGeom.getSegment(
        (cd / 2) / n->pl().fronts()[1].origGeom.getLength(),
        (n->pl().fronts()[1].origGeom.getLength() - cd / 2) /
            n->pl().fronts()[1].origGeom.getLength());

    assert(a.getLine().size() > 1);
    assert(b.getLine().size() > 1);

    if (dist(a.getLine()[0], b.getLine()[0]) >
        dist(a.getLine()[1], b.getLine()[0])) {
      a.reverse();
    }

    pols.push_back(&a);
    pols.push_back(&b);

    auto avg = PolyLine<double>::average(pols).getLine();

    return util::geo::buffer(avg, d, pointsPerCircle);
  }
}

// _____________________________________________________________________________
std::vector<util::geo::Polygon<double>> RenderGraph::getIndStopPolys(
    const std::set<const shared::linegraph::Line*>& served, const LineNode* n,
    double d) const {
  double pointsPerCircle = 32;
  std::vector<util::geo::Polygon<double>> retPolys;

  std::vector<std::pair<size_t, const LineEdge*>> orderedEdgs;
  for (auto e : n->getAdjList()) {
    size_t unserved = 0;
    for (auto lo : e->pl().getLines())
      if (!served.count(lo.line)) unserved++;
    orderedEdgs.push_back({unserved, e});
  }

  std::sort(orderedEdgs.begin(), orderedEdgs.end());

  std::set<const Line*> proced;

  for (auto ep : orderedEdgs) {
    util::geo::Line<double> lineBgeo;
    auto e = ep.second;

    for (size_t p = 0; p < e->pl().getLines().size(); p++) {
      auto lineAtPos = e->pl().lineOccAtPos(p).line;
      if (proced.count(lineAtPos)) continue;
      if (!served.count(lineAtPos)) {
        if (lineBgeo.size()) {
          retPolys.push_back(util::geo::buffer(lineBgeo, d, pointsPerCircle));
          lineBgeo.clear();
        }
        continue;
      }
      proced.insert(lineAtPos);

      auto pos = linePosOn(*n->pl().frontFor(e), lineAtPos, false);

      lineBgeo.push_back({pos.getX(), pos.getY()});
    }

    // last entry
    if (lineBgeo.size()) {
      retPolys.push_back(util::geo::buffer(lineBgeo, d, pointsPerCircle));
    }
  }

  return retPolys;
}

// _____________________________________________________________________________
bool RenderGraph::notCompletelyServed(const LineNode* n) {
  for (auto e : n->getAdjList()) {
    for (auto l : e->pl().getLines()) {
      if (!n->pl().lineServed(l.line)) {
        return true;
      }
    }
  }

  return false;
}

// _____________________________________________________________________________
std::vector<Polygon<double>> RenderGraph::getStopGeoms(
    const LineNode* n, double d, bool tight, size_t pointsPerCircle) const {
  if (notCompletelyServed(n)) {
    // render each stop individually
    auto served = servedLines(n);
    return getIndStopPolys(served, n, d);
  }

  if (n->pl().fronts().size() == 0) return {};
  auto h = getConvexFrontHull(n, d, true, tight, pointsPerCircle);

  // auto nfWidth = getMaxNdFrontWidth(n);

  // prevent too large stations
  // if (util::geo::area(h) > 1.25 * (nfWidth * nfWidth)) {
    // // render each stop individually
    // auto served = servedLines(n);
    // return getIndStopPolys(served, n, d);
  // }

  return {h};
}

// _____________________________________________________________________________
double RenderGraph::getTotalWidth(const LineEdge* e) const {
  return _defWidth * e->pl().getLines().size() +
         _defSpacing * (e->pl().getLines().size() - 1);
}

// _____________________________________________________________________________
size_t RenderGraph::getConnCardinality(const LineNode* n) {
  size_t ret = 0;
  std::map<const Line*, std::set<const LineEdge*>> processed;

  for (const auto& nf : n->pl().fronts()) {
    for (size_t j = 0; j < nf.edge->pl().getLines().size(); j++) {
      const auto& lineOcc = nf.edge->pl().lineOccAtPos(j);

      std::vector<Partner> partners = getPartners(n, nf.edge, lineOcc);

      for (const Partner& p : partners) {
        if (processed[lineOcc.line].find(p.edge) !=
            processed[lineOcc.line].end()) {
          continue;
        }
        ret++;
      }

      processed[lineOcc.line].insert(nf.edge);
    }
  }

  return ret;
}

// _____________________________________________________________________________
double RenderGraph::getWidth(const shared::linegraph::LineEdge* e) const {
  UNUSED(e);
  return _defWidth;
}

// _____________________________________________________________________________
double RenderGraph::getSpacing(const shared::linegraph::LineEdge* e) const {
  UNUSED(e);
  return _defSpacing;
}

// _____________________________________________________________________________
double RenderGraph::getMaxNdFrontWidth() const {
  double ret = 0;
  for (auto n : getNds()) {
    for (const auto& g : n->pl().fronts()) {
      if (getTotalWidth(g.edge) > ret) ret = getTotalWidth(g.edge);
    }
  }
  return ret;
}

// _____________________________________________________________________________
double RenderGraph::getMaxNdFrontWidth(const LineNode* n) const {
  double ret = 0;
  for (const auto& g : n->pl().fronts()) {
    if (getTotalWidth(g.edge) > ret) ret = getTotalWidth(g.edge);
  }
  return ret;
}

// _____________________________________________________________________________
DPoint RenderGraph::linePosOn(const NodeFront& nf, const Line* r,
                              bool origGeom) const {
  size_t p = nf.edge->pl().linePos(r);
  return linePosOn(nf, nf.edge, p, nf.n == nf.edge->getTo(), origGeom);
}

// _____________________________________________________________________________
DPoint RenderGraph::linePosOn(const NodeFront& nf, const LineEdge* e,
                              size_t pos, bool inv, bool origG) const {
  double p;
  if (!inv) {
    p = (getWidth(e) + getSpacing(e)) * pos + getWidth(e) / 2;
  } else {
    p = (getWidth(e) + getSpacing(e)) * (e->pl().getLines().size() - 1 - pos) +
        getWidth(e) / 2;
  }
  // use interpolate here directly for speed
  if (origG) {
    return nf.origGeom.interpolate(nf.origGeom.getLine().front(),
                                   nf.origGeom.getLine().back(), p);
  } else {
    return nf.geom.interpolate(nf.geom.getLine().front(),
                               nf.geom.getLine().back(), p);
  }
}

// _____________________________________________________________________________
double RenderGraph::getOutAngle(const LineNode* n, const LineEdge* e) {
  double checkDist = 10;
  assert(e->getFrom() == n || e->getTo() == n);
  if (e->getFrom() == n) {
    return angBetween(
        *n->pl().getGeom(),
        PolyLine<double>(*e->pl().getGeom()).getPointAtDist(checkDist).p);
  } else {
    return angBetween(
        *n->pl().getGeom(),
        PolyLine<double>(*e->pl().getGeom())
            .getPointAtDist(util::geo::len(*e->pl().getGeom()) - checkDist)
            .p);
  }
}

// _____________________________________________________________________________
void RenderGraph::createMetaNodes() {
  std::vector<NodeFront> cands;
  while ((cands = getNextMetaNodeCand()).size() > 0) {
    // remove all edges completely contained
    for (auto nf : cands) {
      auto onfs = getClosedNodeFronts(nf.n);

      for (auto onf : onfs) {
        LineEdge* e = onf.edge;

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

        if (!found) continue;

        e->getTo()->pl().delFrontFor(e);
        e->getFrom()->pl().delFrontFor(e);
        delEdg(e->getTo(), e->getFrom());
        getEdgGrid()->remove(e);
      }
    }

    // first node has new ref node id
    LineNode* ref = addNd(*cands[0].n->pl().getGeom());

    std::set<LineNode*> toDel;

    for (auto nf : cands) {
      for (auto onf : getOpenNodeFronts(nf.n)) {
        LineEdge* e;
        LineNode* other;
        if (onf.edge->getTo() == nf.n) {
          e = addEdg(onf.edge->getFrom(), ref, onf.edge->pl());
          other = onf.edge->getFrom();
        } else {
          e = addEdg(ref, onf.edge->getTo(), onf.edge->pl());
          other = onf.edge->getTo();
        }

        // update the directions, if necessary
        std::set<shared::linegraph::LineOcc> del;
        for (auto& to : e->pl().getLines()) {
          if (to.direction == nf.n) del.insert(to);
        }

        // also update the edge of the other node front
        NodeFront otherFr = *other->pl().frontFor(onf.edge);
        other->pl().delFrontFor(onf.edge);
        otherFr.edge = e;
        other->pl().addFront(otherFr);

        // remove the original edge
        delEdg(onf.edge->getFrom(), onf.edge->getTo());
        getEdgGrid()->remove(onf.edge);

        // update the original edge in the checked node front
        onf.edge = e;

        // update the lines
        for (auto ro : del) {
          e->pl().delLine(ro.line);
          e->pl().addLine(ro.line, ref);
        }

        toDel.insert(onf.n);

        // update the node front node
        onf.n = ref;

        // add the new node front to the new node
        ref->pl().addFront(onf);
      }
    }

    // delete the nodes marked for deletion
    for (auto toDelNd : toDel) {
      getNdGrid()->remove(toDelNd);
      delNd(toDelNd);
    }
  }
}

// _____________________________________________________________________________
std::vector<NodeFront> RenderGraph::getNextMetaNodeCand() const {
  for (auto n : getNds()) {
    if (n->pl().stops().size()) continue;
    if (getOpenNodeFronts(n).size() != 1) continue;

    std::set<const LineNode*> potClique;

    std::stack<const LineNode*> nodeStack;
    nodeStack.push(n);

    while (!nodeStack.empty()) {
      const LineNode* n = nodeStack.top();
      nodeStack.pop();

      if (n->pl().stops().size() == 0) {
        potClique.insert(n);
        for (auto nff : getClosedNodeFronts(n)) {
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
bool RenderGraph::isClique(std::set<const LineNode*> potClique) const {
  if (potClique.size() < 2) return false;

  for (const LineNode* a : potClique) {
    for (const LineNode* b : potClique) {
      if (util::geo::dist(*a->pl().getGeom(), *b->pl().getGeom()) >
          (_defWidth + _defSpacing) * 10) {
        return false;
      }
    }
  }

  std::set<const LineNode*> periphery;

  for (const LineNode* n : potClique) {
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

    // catch cases where two clique nodes share the same node at their
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
std::vector<NodeFront> RenderGraph::getOpenNodeFronts(const LineNode* n) const {
  std::vector<NodeFront> res;
  for (auto nf : n->pl().fronts()) {
    if (util::geo::len(*nf.edge->pl().getGeom()) >
            (getWidth(nf.edge) + getSpacing(nf.edge)) ||
        (nf.edge->getOtherNd(n)->pl().frontFor(nf.edge)->geom.distTo(
             *nf.edge->getOtherNd(n)->pl().getGeom()) >
         6 * (getWidth(nf.edge) + getSpacing(nf.edge))) ||
        (nf.edge->getTo()->pl().stops().size() > 0) ||
        (nf.edge->getFrom()->pl().stops().size() > 0)) {
      res.push_back(nf);
    }
  }

  return res;
}

// _____________________________________________________________________________
std::vector<NodeFront> RenderGraph::getClosedNodeFronts(
    const LineNode* n) const {
  std::vector<NodeFront> res;
  for (auto nf : n->pl().fronts()) {
    if (!(util::geo::len(*nf.edge->pl().getGeom()) >
          (getWidth(nf.edge) + getSpacing(nf.edge))) &&
        !(nf.edge->getOtherNd(n)->pl().frontFor(nf.edge)->geom.distTo(
              *nf.edge->getOtherNd(n)->pl().getGeom()) >
          6 * (getWidth(nf.edge) + getSpacing(nf.edge))) &&
        (nf.edge->getTo()->pl().stops().size() == 0) &&
        (nf.edge->getFrom()->pl().stops().size() == 0)) {
      res.push_back(nf);
    }
  }

  return res;
}
