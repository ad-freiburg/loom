// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <proj_api.h>
#include "shared/transitgraph/TransitGraph.h"
#include "topo/builder/Builder.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/log/Log.h"

using namespace topo;

using util::geo::Point;
using util::geo::DPoint;
using util::geo::Grid;
using util::geo::Box;
using util::geo::DBox;
using util::geo::extendBox;
using util::geo::PolyLine;
using util::geo::SharedSegments;

// _____________________________________________________________________________
Builder::Builder(const config::TopoConfig* cfg) : _cfg(cfg) {}

// _____________________________________________________________________________
ShrdSegWrap Builder::getNextSharedSegment(TransitGraph* g, bool final,
                                          EdgeGrid* grid) const {
  int i = 0;
  double dmin = _cfg->maxAggrDistance;

  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      i++;
      if (_indEdges.find(e) != _indEdges.end()) continue;

      std::set<TransitEdge*> neighbors;
      grid->getNeighbors(e, fmax(5, dmin * 10), &neighbors);

      for (auto toTest : neighbors) {
        if (_indEdgesPairs.find({e, toTest}) != _indEdgesPairs.end() ||
            _indEdges.find(toTest) != _indEdges.end()) {
          continue;
        }

        if (e != toTest) {
          double dmax = 5;

          if (final) {
            dmax = fmax(dmin, e->pl().getRoutes().size() * (dmin / 2) +
                                  toTest->pl().getRoutes().size() * (dmin / 2));
          }

          const auto& s = e->pl().getPolyline().getSharedSegments(
              toTest->pl().getPolyline(), dmax);

          if (s.segments.size() > 0) {
            _pEdges[{e, toTest}]++;
            _pEdges[{toTest, e}]++;

            if (_pEdges[{e, toTest}] > 20) {
              _indEdgesPairs.insert({e, toTest});
              _indEdgesPairs.insert({toTest, e});
            }

            return ShrdSegWrap(e, toTest, s.segments.front());
          }
        }
      }

      // nothing overlapping found for this edge anymore, mark as processed
      _indEdges.insert(e);
    }
  }

  return ShrdSegWrap();
}

// _____________________________________________________________________________
bool Builder::crossesAt(const TransitNode* a, const TransitEdge* ea,
                                const TransitEdge* eb) const {
  if (!(ea->getTo() == a || ea->getFrom() == a) ||
      !(eb->getTo() == a || eb->getFrom() == a))
    return false;

  double lookDist = 100, minInAngle = 0.2;

  const auto e = ea;
  const auto f = eb;

  double angleE = 0, angleF = 0, angleEPartner = 0, angleFPartner = 0;

  if (e->getTo() == a) {
    // look at end of geom
    DPoint a = e->pl().getPolyline().back();
    DPoint b = e->pl()
                   .getPolyline()
                   .getPointAtDist(e->pl().getPolyline().getLength() - lookDist)
                   .p;
    angleE = util::geo::angBetween(a, b);
  } else {
    // look at front of geom
    DPoint a = e->pl().getPolyline().front();
    DPoint b = e->pl().getPolyline().getPointAtDist(lookDist).p;
    angleE = util::geo::angBetween(a, b);
  }

  if (f->getTo() == a) {
    // look at end of geom
    DPoint a = f->pl().getPolyline().back();
    DPoint b = f->pl()
                   .getPolyline()
                   .getPointAtDist(f->pl().getPolyline().getLength() - lookDist)
                   .p;
    angleF = util::geo::angBetween(a, b);
  } else {
    // look at front of geom
    DPoint a = f->pl().getPolyline().front();
    DPoint b = f->pl().getPolyline().getPointAtDist(lookDist).p;
    angleF = util::geo::angBetween(a, b);
  }

  double inDiff = fabs(angleE - angleF);

  if (inDiff < minInAngle) return false;

  bool eFound = false, fFound = false;

  // now lets get some edge that is not e or f and has all of e's routes

  for (auto edge : a->getAdjList()) {
    if (edge == ea || edge == eb) continue;

    if (routeEq(edge, e)) {
      if (edge->getTo() != a) {
        // look at end of geom
        DPoint a = edge->pl().getPolyline().back();
        DPoint b =
            edge->pl()
                .getPolyline()
                .getPointAtDist(edge->pl().getPolyline().getLength() - lookDist)
                .p;
        angleEPartner = util::geo::angBetween(a, b);
      } else {
        // look at front of geom
        DPoint a = edge->pl().getPolyline().front();
        DPoint b = edge->pl().getPolyline().getPointAtDist(lookDist).p;
        angleEPartner = util::geo::angBetween(a, b);
      }
      eFound = true;
      break;
    }
  }

  for (auto edge : a->getAdjList()) {
    if (edge == ea || edge == eb) continue;

    if (routeEq(edge, f)) {
      if (edge->getTo() != a) {
        // look at end of geom
        DPoint a = edge->pl().getPolyline().back();
        DPoint b =
            edge->pl()
                .getPolyline()
                .getPointAtDist(edge->pl().getPolyline().getLength() - lookDist)
                .p;
        angleFPartner = util::geo::angBetween(a, b);
      } else {
        // look at front of geom
        DPoint a = edge->pl().getPolyline().front();
        DPoint b = edge->pl().getPolyline().getPointAtDist(lookDist).p;
        angleFPartner = util::geo::angBetween(a, b);
      }
      fFound = true;
      break;
    }
  }

  if (!eFound || !fFound) return false;

  double outDiff = fabs(angleEPartner - angleFPartner);

  if (outDiff < minInAngle) return false;

  if (fabs(angleE - angleEPartner) < 0.2 &&
      fabs(angleF - angleFPartner) < 0.2) {
    return true;
  }

  return false;
}

// _____________________________________________________________________________
bool Builder::routeEq(const TransitEdge* a, const TransitEdge* b) const {
  // shortcut
  if (a->pl().getRoutes().size() != a->pl().getRoutes().size()) return false;

  const auto shrNd = TransitGraph::sharedNode(a, b);

  // TODO: remove quadratic code
  for (const auto& ra : a->pl().getRoutes()) {
    bool found = false;
    for (const auto& rb : b->pl().getRoutes()) {
      if (ra.route == rb.route && ra.style == rb.style) {
        // if the route does not continue over the shared node, the routes
        // are not equivalent!
        if (!shrNd->pl().connOccurs(ra.route, a, b)) return false;

        if (ra.direction == 0 && rb.direction == 0) {
          found = true;
          break;
        }
        if (ra.direction == shrNd && rb.direction != 0 &&
            rb.direction != shrNd) {
          found = true;
          break;
        }
        if (ra.direction != shrNd && ra.direction != 0 &&
            rb.direction == shrNd) {
          found = true;
          break;
        }
      }
    }
    if (!found) return false;
  }
  return true;
}

// _____________________________________________________________________________
bool Builder::createTopologicalNodes(TransitGraph* g, bool final) {
  ShrdSegWrap w;
  _indEdges.clear();
  _pEdges.clear();
  bool found = false;

  EdgeGrid grid = getGeoIndex(g);

  while ((w = getNextSharedSegment(g, final, &grid)).e) {
    const auto curEdgeGeom = w.e;
    const auto cmpEdgeGeom = w.f;

    const auto& eap = w.s.first.first;
    const auto& ebp = w.s.second.first;
    const auto& fap = w.s.first.second;
    const auto& fbp = w.s.second.second;

    auto ea =
        curEdgeGeom->pl().getPolyline().getSegment(0, eap.totalPos);
    auto ec =
        curEdgeGeom->pl().getPolyline().getSegment(ebp.totalPos, 1);

    PolyLine<double> fa, fc;

    if (fap.totalPos > fbp.totalPos) {
      fa = cmpEdgeGeom->pl().getPolyline().getSegment(fap.totalPos, 1);
      fc = cmpEdgeGeom->pl().getPolyline().getSegment(0, fbp.totalPos);
    } else {
      fa = cmpEdgeGeom->pl().getPolyline().getSegment(0, fap.totalPos);
      fc = cmpEdgeGeom->pl().getPolyline().getSegment(fbp.totalPos, 1);
    }

    auto ab = getAveragedFromSharedSeg(w);

    // new nodes at the start and end of the shared segment
    TransitNode *a = 0, *b = 0;

    double maxSnapDist = 20;

    if (ea.getLength() < maxSnapDist) {
      a = w.e->getFrom();
    }

    if (ec.getLength() < maxSnapDist) {
      b = w.e->getTo();
    }

    if (fa.getLength() < maxSnapDist) {
      a = !(fap.totalPos > fbp.totalPos) ? w.f->getFrom() : w.f->getTo();
    }

    if (fc.getLength() < maxSnapDist) {
      b = !(fap.totalPos > fbp.totalPos) ? w.f->getTo() : w.f->getFrom();
    }

    if (!a) {
      a = g->addNd({w.s.first.first.p});
    }

    if (!b) {
      b = g->addNd({w.s.second.first.p});
    }

    if (a == b) {
      continue;
    }

    // TODO: what about the nodes we added above??
    if (crossesAt(a, w.e, w.f)) continue;
    if (crossesAt(b, w.e, w.f)) continue;

    TransitEdgePL eaEdgeGeom(ea);
    TransitEdgePL abEdgeGeom(ab);
    TransitEdgePL ecEdgeGeom(ec);

    const TransitNode *faDir = 0, *fcDir = 0;

    if (fap.totalPos > fbp.totalPos) {
      faDir = w.f->getTo();
      fcDir = b;
    } else {
      faDir = a;
      fcDir = w.f->getTo();
    }

    TransitEdgePL faEdgeGeom(fa);
    TransitEdgePL fcEdgeGeom(fc);

    auto wefrom = w.e->getFrom();
    auto weto = w.e->getTo();
    auto wffrom = w.f->getFrom();
    auto wfto = w.f->getTo();

    for (const auto& r : curEdgeGeom->pl().getRoutes()) {
      if (!r.direction) {
        eaEdgeGeom.addRoute(r.route, 0, r.style);
        abEdgeGeom.addRoute(r.route, 0, r.style);
        ecEdgeGeom.addRoute(r.route, 0, r.style);
      } else if (r.direction == weto) {
        eaEdgeGeom.addRoute(r.route, a, r.style);
        abEdgeGeom.addRoute(r.route, b, r.style);
        ecEdgeGeom.addRoute(r.route, weto, r.style);
      } else {
        eaEdgeGeom.addRoute(r.route, wefrom, r.style);
        abEdgeGeom.addRoute(r.route, a, r.style);
        ecEdgeGeom.addRoute(r.route, b, r.style);
      }
    }

    for (const auto& r : cmpEdgeGeom->pl().getRoutes()) {
      if (!r.direction) {
        faEdgeGeom.addRoute(r.route, 0, r.style);
        abEdgeGeom.addRoute(r.route, 0, r.style);
        fcEdgeGeom.addRoute(r.route, 0, r.style);
      } else if ((r.direction == wfto)) {
        if (fap.totalPos > fbp.totalPos) {
          faEdgeGeom.addRoute(r.route, wfto, r.style);
          abEdgeGeom.addRoute(r.route, a, r.style);
          fcEdgeGeom.addRoute(r.route, b, r.style);
        } else {
          faEdgeGeom.addRoute(r.route, a, r.style);
          abEdgeGeom.addRoute(r.route, b, r.style);
          fcEdgeGeom.addRoute(r.route, wfto, r.style);
        }
      } else {
        if (fap.totalPos > fbp.totalPos) {
          faEdgeGeom.addRoute(r.route, a, r.style);
          abEdgeGeom.addRoute(r.route, b, r.style);
          fcEdgeGeom.addRoute(r.route, b, r.style);
        } else {
          faEdgeGeom.addRoute(r.route, wffrom, r.style);
          abEdgeGeom.addRoute(r.route, a, r.style);
          fcEdgeGeom.addRoute(r.route, wfto, r.style);
        }
      }
    }

    // delete old edges
    grid.remove(g->getEdg(w.e->getFrom(), w.e->getTo()));
    grid.remove(g->getEdg(w.f->getFrom(), w.f->getTo()));
    g->delEdg(w.e->getFrom(), w.e->getTo());
    g->delEdg(w.f->getFrom(), w.f->getTo());

    TransitEdge *eaE = 0, *abE = 0, *ebE = 0, *faE = 0, *fbE = 0;

    // add new edges
    assert(a != b);
    abE = g->addEdg(a, b, TransitEdgePL());
    // if (abE) abE->pl().setEdge(abE);

    if (a != wefrom) {
      eaE = g->addEdg(wefrom, a);
      edgeRpl(wefrom, w.e, eaE);
    } else {
      edgeRpl(a, w.e, abE);
    }

    if (b != weto) {
      ebE = g->addEdg(b, weto);
      edgeRpl(weto, w.e, ebE);
    } else {
      assert(b == weto);
      edgeRpl(b, w.e, abE);
    }

    if (fap.totalPos > fbp.totalPos) {
      if (a != wfto) {
        faE = g->addEdg(a, wfto);
        edgeRpl(wfto, w.f, faE);
      } else {
        assert(a == wfto);
        edgeRpl(a, w.f, abE);
      }

      if (b != wffrom) {
        fbE = g->addEdg(wffrom, b);
        edgeRpl(wffrom, w.f, fbE);
      } else {
        edgeRpl(b, w.f, abE);
      }
    } else {
      if (a != wffrom) {
        faE = g->addEdg(wffrom, a);
        edgeRpl(wffrom, w.f, faE);
      } else {
        assert(a == wffrom);
        edgeRpl(a, w.f, abE);
      }

      if (b != wfto) {
        fbE = g->addEdg(b, wfto);
        edgeRpl(wfto, w.f, fbE);
      } else {
        edgeRpl(b, w.f, abE);
      }
    }

    if (abE) {
      abE->pl() = abEdgeGeom;
      grid.add(*abE->pl().getGeom(), abE);
    } else {
      // we use abE below without checking!
      assert(false);
    }

    if (eaE) {
      // eaE->pl().addEdgeTripGeom(eaEdgeGeom);
      eaE->pl() = eaEdgeGeom;
      // eaE->pl().simplify();
      grid.add(*eaE->pl().getGeom(), eaE);
      // a->pl().sewConnectionsTogether(eaE, abE);
    }

    if (ebE) {
      ebE->pl() = ecEdgeGeom;
      // ebE->pl().simplify();
      grid.add(*ebE->pl().getGeom(), ebE);
      // b->pl().sewConnectionsTogether(abE, ebE);
    }

    if (faE) {
      faE->pl() = faEdgeGeom;
      // faE->pl().simplify();
      grid.add(*faE->pl().getGeom(), faE);
      // a->pl().sewConnectionsTogether(faE, abE);
    }

    if (fbE) {
      fbE->pl() = fcEdgeGeom;
      // fbE->pl().simplify();
      grid.add(*fbE->pl().getGeom(), fbE);
      // b->pl().sewConnectionsTogether(abE, fbE);
    }

    found = true;
  }

  return found;
}

// _____________________________________________________________________________
void Builder::averageNodePositions(TransitGraph* g) {
  for (auto n : *g->getNds()) {
    double x = 0, y = 0;
    size_t c = 0;

    for (auto e : n->getAdjList()) {
      if (e->getTo() != n) {
        x += e->pl().getPolyline().front().getX();
        y += e->pl().getPolyline().front().getY();
      } else {
        x += e->pl().getPolyline().back().getX();
        y += e->pl().getPolyline().back().getY();
      }
      c++;
    }

    if (c > 0)
      n->pl().setGeom(
          DPoint(x / static_cast<double>(c), y / static_cast<double>(c)));
  }
}

// _____________________________________________________________________________
void Builder::removeEdgeArtifacts(TransitGraph* g) {
  while (contractNodes(g)) {}
    ;
}

// _____________________________________________________________________________
void Builder::removeNodeArtifacts(TransitGraph* g) {
  while (contractEdges(g)) {}
    ;
}

// _____________________________________________________________________________
bool Builder::contractNodes(TransitGraph* g) {
  double MIN_SEG_LENGTH = 20;

  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (e->pl().getPolyline().getLength() < MIN_SEG_LENGTH) {
        auto from = e->getFrom();
        auto to = e->getTo();
        if (from->pl().getStops().size() == 0 ||
            to->pl().getStops().size() == 0) {
          combineNodes(from, to, g);
          return true;
        }
      }
    }
  }
  return false;
}

// _____________________________________________________________________________
bool Builder::contractEdges(TransitGraph* g) {
  for (auto n : *g->getNds()) {
    std::vector<TransitEdge*> edges;
    edges.insert(edges.end(), n->getAdjList().begin(), n->getAdjList().end());
    if (edges.size() == 2 && n->pl().getStops().size() == 0) {
      if (routeEq(edges[0], edges[1])) {
        combineEdges(edges[0], edges[1], n, g);
        return true;
      }
    }
  }
  return false;
}

// _____________________________________________________________________________
bool Builder::combineEdges(TransitEdge* a, TransitEdge* b, TransitNode* n,
                           TransitGraph* g) {
  assert((a->getTo() == n || a->getFrom() == n) &&
         (b->getTo() == n || b->getFrom() == n));

  TransitEdge* newEdge = 0;
  util::geo::PolyLine<double> newPl;

  // TODO: there is some copying going on below, which is not always necessary.
  // insert a non-const getLine to polyline and re-use existing polylines where
  // possible

  if (a->getTo() == n && b->getTo() != n) {
    //   a       b
    // ----> n ---->
    auto lineA = a->pl().getPolyline().getLine();
    const auto& lineB = b->pl().getPolyline().getLine();
    lineA.insert(lineA.end(), lineB.begin(), lineB.end());
    newPl = util::geo::PolyLine<double>(lineA);

    newEdge = g->addEdg(a->getFrom(), b->getTo(), a->pl());
    routeDirRepl(n, newEdge->getTo(), newEdge);

    edgeRpl(newEdge->getFrom(), a, newEdge);
    edgeRpl(newEdge->getTo(), b, newEdge);
  }

  if (a->getTo() != n && b->getTo() == n) {
    //   a       b
    // <---- n <----
    auto lineB = b->pl().getPolyline().getLine();
    const auto& lineA = a->pl().getPolyline().getLine();
    lineB.insert(lineB.end(), lineA.begin(), lineA.end());
    newPl = util::geo::PolyLine<double>(lineB);

    newEdge = g->addEdg(b->getFrom(), a->getTo(), b->pl());
    routeDirRepl(n, newEdge->getTo(), newEdge);

    edgeRpl(newEdge->getFrom(), b, newEdge);
    edgeRpl(newEdge->getTo(), a, newEdge);
  }

  if (a->getFrom() == n && b->getFrom() == n) {
    //   a       b
    // <---- n ---->
    auto lineA = a->pl().getPolyline().getLine();
    std::reverse(lineA.begin(), lineA.end());
    const auto& lineB = b->pl().getPolyline().getLine();
    lineA.insert(lineA.end(), lineB.begin(), lineB.end());
    newPl = util::geo::PolyLine<double>(lineA);

    newEdge = g->addEdg(a->getTo(), b->getTo(), b->pl());
    routeDirRepl(n, newEdge->getFrom(), newEdge);

    edgeRpl(newEdge->getFrom(), a, newEdge);
    edgeRpl(newEdge->getTo(), b, newEdge);
  }

  if (a->getTo() == n && b->getTo() == n) {
    //   a       b
    // ----> n <----
    auto lineA = a->pl().getPolyline().getLine();
    const auto& lineB = b->pl().getPolyline().getLine();
    lineA.insert(lineA.end(), lineB.rbegin(), lineB.rend());
    newPl = util::geo::PolyLine<double>(lineA);

    newEdge = g->addEdg(a->getFrom(), b->getFrom(), a->pl());
    routeDirRepl(n, newEdge->getTo(), newEdge);

    edgeRpl(newEdge->getFrom(), a, newEdge);
    edgeRpl(newEdge->getTo(), b, newEdge);
  }

  // set new polyline and smoothen a bit
  newPl.smoothenOutliers(50);
  newPl.applyChaikinSmooth(3);
  newPl.simplify(3);
  newEdge->pl().setPolyline(newPl);

  g->delEdg(a->getFrom(), a->getTo());
  g->delEdg(b->getFrom(), b->getTo());
  g->delNd(n);

  return true;
}

// _____________________________________________________________________________
bool Builder::combineNodes(TransitNode* a, TransitNode* b, TransitGraph* g) {
  assert(a->pl().getStops().size() == 0 || b->pl().getStops().size() == 0);

  if (a->pl().getStops().size() != 0) {
    TransitNode* c = a;
    a = b;
    b = c;
  }

  TransitEdge* connecting = g->getEdg(a, b);
  assert(connecting);

  // we will delete a and the connecting edge {a, b}.
  // b will be the contracted node

  // go through all route exceptions in a
  for (const auto& ro : a->pl().getConnExc()) {
    for (const auto& exFr : ro.second) {
      for (const auto& exTo : exFr.second) {
        if (exFr.first == connecting) {
          for (auto e : b->getAdjList()) {
            if (e == connecting) continue;
            if (e == exTo) continue;
            if (!e->pl().hasRoute(ro.first)) continue;
            b->pl().addConnExc(ro.first, e, exTo);
          }
        } else if (exTo == connecting) {
          for (auto e : b->getAdjList()) {
            if (e == connecting) continue;
            if (e == exFr.first) continue;
            if (!e->pl().hasRoute(ro.first)) continue;
            b->pl().addConnExc(ro.first, e, exFr.first);
          }
        } else {
          b->pl().addConnExc(ro.first, exFr.first, exTo);
        }
      }
    }
  }

  // go through all route exceptions in b
  for (const auto& ro : b->pl().getConnExc()) {
    for (const auto& exFr : ro.second) {
      for (const auto& exTo : exFr.second) {
        if (exFr.first == connecting) {
          for (auto e : a->getAdjList()) {
            if (e == connecting) continue;
            if (e == exTo) continue;
            if (!e->pl().hasRoute(ro.first)) continue;
            b->pl().addConnExc(ro.first, e, exTo);
          }
        } else if (exTo == connecting) {
          for (auto e : a->getAdjList()) {
            if (e == connecting) continue;
            if (e == exFr.first) continue;
            if (!e->pl().hasRoute(ro.first)) continue;
            b->pl().addConnExc(ro.first, e, exFr.first);
          }
        } else {
          b->pl().addConnExc(ro.first, exFr.first, exTo);
        }
      }
    }
  }

  for (auto* oldE : a->getAdjList()) {
    if (oldE->getFrom() != a) continue;
    if (connecting == oldE) continue;

    assert(b != oldE->getTo());

    // add a new edge going from b to the non-a node
    auto* newE = g->addEdg(b, oldE->getTo(), oldE->pl());

    // update route dirs
    routeDirRepl(a, b, newE);

    // replace each occurance of oldE in the non-a node with the newE
    edgeRpl(newE->getTo(), oldE, newE);
    edgeRpl(newE->getFrom(), oldE, newE);
  }

  for (auto* oldE : a->getAdjList()) {
    if (oldE->getTo() != a) continue;
    if (connecting == oldE) continue;

    assert(b != oldE->getFrom());

    auto* newE = g->addEdg(oldE->getFrom(), b, oldE->pl());

    // update route dirs
    routeDirRepl(a, b, newE);

    // replace each occurance of oldE in the non-a node with the newE
    edgeRpl(newE->getTo(), oldE, newE);
    edgeRpl(newE->getFrom(), oldE, newE);
  }

  g->delEdg(a, b);
  g->delNd(a);

  return true;
}

// _____________________________________________________________________________
PolyLine<double> Builder::getAveragedFromSharedSeg(const ShrdSegWrap& w) const {
  const auto& geomA = w.e;
  const auto& geomB = w.f;

  PolyLine<double> a = geomA->pl().getPolyline().getSegment(
      w.s.first.first.totalPos, w.s.second.first.totalPos);

  PolyLine<double> b;

  if (w.s.first.second.totalPos > w.s.second.second.totalPos) {
    b = geomB->pl().getPolyline().getSegment(w.s.second.second.totalPos,
                                             w.s.first.second.totalPos);
    b.reverse();
  } else {
    b = geomB->pl().getPolyline().getSegment(w.s.first.second.totalPos,
                                             w.s.second.second.totalPos);
  }

  bool dominationHeuristic = false;

  if (dominationHeuristic) {
    // test if one geom is the dominant line here
    bool aDomin = lineDominatesSharedSeg(w, w.e);
    bool bDomin = lineDominatesSharedSeg(w, w.f);
    if (aDomin && !bDomin) {
      return a;
    }

    if (bDomin && !aDomin) {
      return b;
    }
  }

  std::vector<const PolyLine<double>*> avg{&a, &b};
  std::vector<double> weights{
      1.0 * geomA->pl().getRoutes().size() * geomA->pl().getRoutes().size(),
      1.0 * geomB->pl().getRoutes().size() * geomB->pl().getRoutes().size()};

  PolyLine<double> ret = PolyLine<double>::average(avg, weights);
  ret.simplify(5);
  return ret;
}

// _____________________________________________________________________________
bool Builder::lineDominatesSharedSeg(const ShrdSegWrap& w,
                                     TransitEdge* e) const {
  if (e != w.e && e != w.f) return false;

  double LOOKAHEAD = 50, DELTA = .5;

  DPoint a, b, c, d;

  if (e == w.e) {
    const auto& geom = w.e;
    double la = LOOKAHEAD / geom->pl().getPolyline().getLength();

    a = geom->pl().getPolyline().getPointAt(w.s.first.first.totalPos - la).p;
    b = geom->pl().getPolyline().getPointAt(w.s.first.first.totalPos).p;
    c = geom->pl().getPolyline().getPointAt(w.s.second.first.totalPos).p;
    d = geom->pl().getPolyline().getPointAt(w.s.second.first.totalPos + la).p;
  } else if (e == w.f) {
    const auto& geom = w.f;
    double la = LOOKAHEAD / geom->pl().getPolyline().getLength();

    if (w.s.first.second.totalPos > w.s.second.second.totalPos) {
      a = geom->pl().getPolyline().getPointAt(w.s.first.second.totalPos + la).p;
      d = geom->pl()
              .getPolyline()
              .getPointAt(w.s.second.second.totalPos - la)
              .p;
    } else {
      a = geom->pl().getPolyline().getPointAt(w.s.first.second.totalPos - la).p;
      d = geom->pl()
              .getPolyline()
              .getPointAt(w.s.second.second.totalPos + la)
              .p;
    }
    b = geom->pl().getPolyline().getPointAt(w.s.first.second.totalPos).p;
    c = geom->pl().getPolyline().getPointAt(w.s.second.second.totalPos).p;
  }

  double ang = util::geo::angBetween(a, b) - util::geo::angBetween(c, d);
  double tang = fabs(tan(ang));

  return tang < DELTA;
}

// _____________________________________________________________________________
EdgeGrid Builder::getGeoIndex(const TransitGraph* g) const {
  EdgeGrid grid(120, 120, getGraphBoundingBox(g));

  for (auto n : g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      grid.add(e->pl().getPolyline().getLine(), e);
    }
  }

  return grid;
}

// _____________________________________________________________________________
DBox Builder::getGraphBoundingBox(const TransitGraph* g) const {
  DBox b;

  for (auto n : g->getNds()) {
    b = extendBox(*n->pl().getGeom(), b);
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      b = extendBox(e->pl().getPolyline().getLine(), b);
    }
  }

  return b;
}

// _____________________________________________________________________________
void Builder::edgeRpl(TransitNode* n, const TransitEdge* oldE,
                      const TransitEdge* newE) const {
  // replace in from
  for (auto& r : n->pl().getConnExc()) {
    for (auto& exFr : r.second) {
      if (exFr.first == oldE) {
        std::swap(r.second[newE], exFr.second);
        r.second.erase(exFr.first);
      }
    }
  }

  // replace in to
  for (auto& r : n->pl().getConnExc()) {
    for (auto& exFr : r.second) {
      for (auto& exTo : exFr.second) {
        if (exTo == oldE) {
          exFr.second.erase(exTo);
          exFr.second.insert(newE);
        }
      }
    }
  }
}

// _____________________________________________________________________________
void Builder::routeDirRepl(TransitNode* oldN, TransitNode* newN,
                           TransitEdge* e) const {
  for (auto& ro : e->pl().getRoutes()) {
    if (ro.direction == oldN) {
      shared::transitgraph::RouteOcc newRo = ro;
      newRo.direction = newN;

      // delete old
      e->pl().getRoutes().erase(ro);

      // add new
      e->pl().getRoutes().insert(newRo);
    }
  }
}
