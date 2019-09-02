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

double MAX_SNAP_DIST = 20;
double MIN_SEG_LENGTH = 20;

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

      std::cerr << neighbors.size() << std::endl;
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

  std::cerr << "oops..." << std::endl;
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
  return createTopologicalNodes(g, final, 10000000);
}

// _____________________________________________________________________________
bool Builder::createTopologicalNodes(TransitGraph* g, bool final,
                                     size_t steps) {
  ShrdSegWrap w;
  _indEdges.clear();
  _pEdges.clear();
  bool found = false;

  EdgeGrid grid = getGeoIndex(g);

  while (steps > 0 && (w = getNextSharedSegment(g, final, &grid)).e) {
    std::cerr << "---------------------" << std::endl;
    steps--;
    const auto curEdgeGeom = w.e;
    const auto cmpEdgeGeom = w.f;

    const auto& eap = w.s.first.first;
    const auto& ebp = w.s.second.first;
    const auto& fap = w.s.first.second;
    const auto& fbp = w.s.second.second;

    auto ea = curEdgeGeom->pl().getPolyline().getSegment(0, eap.totalPos);
    auto ec = curEdgeGeom->pl().getPolyline().getSegment(ebp.totalPos, 1);

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

    // if (ea.getLength() < MAX_SNAP_DIST && snapElig(w.e->getFrom(), w.f, g)) {
      // std::cerr << "SNAP A" << std::endl;
      // a = w.e->getFrom();
    // }

    // if (ec.getLength() < MAX_SNAP_DIST && snapElig(w.e->getTo(), w.f, g)) {
      // std::cerr << "SNAP B" << std::endl;
      // b = w.e->getTo();
    // }

    // if (fap.totalPos <= fbp.totalPos) {
      // if (fa.getLength() < MAX_SNAP_DIST && snapElig(w.f->getFrom(), w.e, g)) {
        // std::cerr << "SNAP A" << std::endl;
        // a = w.f->getFrom();
      // }

      // if (fc.getLength() < MAX_SNAP_DIST && snapElig(w.f->getTo(), w.e, g)) {
        // std::cerr << "SNAP B" << std::endl;
        // b = w.f->getTo();
      // }
    // } else {
      // if (fa.getLength() < MAX_SNAP_DIST && snapElig(w.f->getTo(), w.e, g)) {
        // std::cerr << "SNAP A" << std::endl;
        // a = w.f->getTo();
      // }

      // if (fc.getLength() < MAX_SNAP_DIST && snapElig(w.f->getFrom(), w.e, g)) {
        // std::cerr << "SNAP B" << std::endl;
        // b = w.f->getFrom();
      // }
    // }

    if (!a) {
      std::cerr << "NOSNAP A" << std::endl;
      a = g->addNd({w.s.first.first.p});
    }

    if (!b) {
      std::cerr << "NOSNAP B" << std::endl;
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

    bool eaShort = false;
    bool ecShort = false;
    bool faShort = false;
    bool fcShort = false;
    bool abShort = false;

    // bool eaShort = eaEdgeGeom.getPolyline().getLength() < 5000;
    // bool ecShort = ecEdgeGeom.getPolyline().getLength() < 5000;
    // bool faShort = faEdgeGeom.getPolyline().getLength() < 5000;
    // bool fcShort = fcEdgeGeom.getPolyline().getLength() < 5000;
    // bool abShort = abEdgeGeom.getPolyline().getLength() < 5000;

    for (const auto& r : curEdgeGeom->pl().getRoutes()) {
      if (!r.direction) {
        eaEdgeGeom.addRoute(r.route, 0, r.style);
        abEdgeGeom.addRoute(r.route, 0, r.style);
        ecEdgeGeom.addRoute(r.route, 0, r.style);
      } else if (r.direction == weto) {
        eaEdgeGeom.addRoute(r.route, eaShort ? 0 : a, r.style);
        abEdgeGeom.addRoute(r.route, abShort ? 0 : b, r.style);
        ecEdgeGeom.addRoute(r.route, ecShort ? 0 : weto, r.style);
      } else {
        eaEdgeGeom.addRoute(r.route, eaShort ? 0 : wefrom, r.style);
        abEdgeGeom.addRoute(r.route, abShort ? 0 : a, r.style);
        ecEdgeGeom.addRoute(r.route, ecShort ? 0 : b, r.style);
      }
    }

    for (const auto& r : cmpEdgeGeom->pl().getRoutes()) {
      if (!r.direction) {
        faEdgeGeom.addRoute(r.route, 0, r.style);
        abEdgeGeom.addRoute(r.route, 0, r.style);
        fcEdgeGeom.addRoute(r.route, 0, r.style);
      } else if ((r.direction == wfto)) {
        if (fap.totalPos > fbp.totalPos) {
          std::cerr << "A " << r.route->getId() << std::endl;
          faEdgeGeom.addRoute(r.route, faShort ? 0 : wfto, r.style);
          abEdgeGeom.addRoute(r.route, abShort ? 0 : a, r.style);
          fcEdgeGeom.addRoute(r.route, fcShort ? 0 : b, r.style);
        } else {
          std::cerr << "B " << r.route->getId() << std::endl;
          faEdgeGeom.addRoute(r.route, faShort ? 0 : a, r.style);
          abEdgeGeom.addRoute(r.route, abShort ? 0 : b, r.style);
          fcEdgeGeom.addRoute(r.route, fcShort ? 0 : wfto, r.style);
        }
      } else {
        if (fap.totalPos > fbp.totalPos) {
          std::cerr << fap.totalPos << " vs " << fbp.totalPos << std::endl;
          std::cerr << "C " << r.route->getId() << std::endl;
          faEdgeGeom.addRoute(r.route, faShort ? 0 : a, r.style);
          abEdgeGeom.addRoute(r.route, b, r.style);
          fcEdgeGeom.addRoute(r.route, fcShort ? 0 : wffrom, r.style);
          // fcEdgeGeom.addRoute(r.route, b, r.style);
        } else {
          std::cerr << "D " << r.route->getId() << std::endl;
          faEdgeGeom.addRoute(r.route, faShort ? 0 : wffrom, r.style);
          abEdgeGeom.addRoute(r.route, a, r.style);
          fcEdgeGeom.addRoute(r.route, fcShort ? 0 : b, r.style);
          // fcEdgeGeom.addRoute(r.route, wfto, r.style);
        }
      }
    }

    std::cerr << "---" << std::endl;

    for (const auto& r : abEdgeGeom.getRoutes()) {
      std::cerr << r.route->getLabel() << std::endl;
    }

    std::cerr << "====" << std::endl;

    // delete old edges
    grid.remove(g->getEdg(w.e->getFrom(), w.e->getTo()));
    grid.remove(g->getEdg(w.f->getFrom(), w.f->getTo()));
    g->delEdg(w.e->getFrom(), w.e->getTo());
    g->delEdg(w.f->getFrom(), w.f->getTo());

    TransitEdge *eaE = 0, *abE = 0, *ebE = 0, *faE = 0, *fbE = 0, *helper = 0;

    // add new edges
    assert(a != b);
    abE = g->addEdg(a, b, TransitEdgePL());

    assert(eap.totalPos < ebp.totalPos);

    if (a != wefrom) {
      if (g->getEdg(wefrom, a)) {
        std::tie(eaE, helper) = split(eaEdgeGeom, wefrom, a, g);
        edgeRpl(wefrom, w.e, helper);
      } else {
        eaE = g->addEdg(wefrom, a);
        assert(eaE != abE);
        edgeRpl(wefrom, w.e, eaE);
      }
    } else {
      // this is the case when e.from->a was below the merge threshold
      edgeRpl(a, w.e, abE);
    }

    if (b != weto) {
      if (g->getEdg(b, weto)) {
        std::tie(ebE, helper) = split(ecEdgeGeom, b, weto, g);
        edgeRpl(weto, w.e, helper);
      } else {
        ebE = g->addEdg(b, weto);
        edgeRpl(weto, w.e, ebE);
      }
    } else {
      assert(b == weto);
      edgeRpl(b, w.e, abE);
    }

    if (fap.totalPos > fbp.totalPos) {
      if (a != wfto) {
        if (g->getEdg(a, wfto)) {
          std::tie(faE, helper) = split(faEdgeGeom, a, wfto, g);
          edgeRpl(wfto, w.f, helper);
        } else {
          faE = g->addEdg(a, wfto);
          assert(faE != abE);
          edgeRpl(wfto, w.f, faE);
        }
      } else {
        assert(a == wfto);
        edgeRpl(a, w.f, abE);
      }

      if (b != wffrom) {
        if (g->getEdg(wffrom, b)) {
          std::tie(fbE, helper) = split(fcEdgeGeom, wffrom, b, g);
          edgeRpl(wffrom, w.f, helper);
        } else {
          fbE = g->addEdg(wffrom, b);
          assert(fbE != abE);
          edgeRpl(wffrom, w.f, fbE);
        }
      } else {
        assert(b == wffrom);
        edgeRpl(b, w.f, abE);
      }
    } else {
      if (a != wffrom) {
        if (g->getEdg(wffrom, a)) {
          std::tie(faE, helper) = split(faEdgeGeom, wffrom, a, g);
          edgeRpl(wffrom, w.f, helper);
        } else {
          faE = g->addEdg(wffrom, a);
          assert(faE != abE);
          edgeRpl(wffrom, w.f, faE);
        }
      } else {
        assert(a == wffrom);
        edgeRpl(a, w.f, abE);
      }

      if (b != wfto) {
        if (g->getEdg(b, wfto)) {
          std::tie(fbE, helper) = split(fcEdgeGeom, b, wfto, g);
          edgeRpl(wfto, w.f, helper);
        } else {
          fbE = g->addEdg(b, wfto);
          assert(fbE != abE);
          edgeRpl(wfto, w.f, fbE);
        }
      } else {
        assert(b == wfto);
        edgeRpl(b, w.f, abE);
      }
    }

    if (abE) {
      abE->pl() = abEdgeGeom;  // no merge needed here
      grid.add(*abE->pl().getGeom(), abE);
    } else {
      assert(false);
    }

    // TODO: plmerge is not necessary here

    if (eaE) {
      plMerge(eaE->pl(), eaEdgeGeom);
      grid.add(*eaE->pl().getGeom(), eaE);
    }

    if (ebE) {
      plMerge(ebE->pl(), ecEdgeGeom);
      grid.add(*ebE->pl().getGeom(), ebE);
    }

    if (faE) {
      plMerge(faE->pl(), faEdgeGeom);
      grid.add(*faE->pl().getGeom(), faE);
    }

    if (fbE) {
      plMerge(fbE->pl(), fcEdgeGeom);
      grid.add(*fbE->pl().getGeom(), fbE);
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
  while (contractNodes(g)) {
    std::cerr << g->getNds()->size() << std::endl;
  };
}

// _____________________________________________________________________________
void Builder::removeNodeArtifacts(TransitGraph* g) {
  while (contractEdges(g)) {
    std::cerr << g->getNds()->size() << std::endl;
  };
}

// _____________________________________________________________________________
bool Builder::contractNodes(TransitGraph* g) {
  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (e->pl().getPolyline().getLength() < MIN_SEG_LENGTH) {
        auto from = e->getFrom();
        auto to = e->getTo();
        if ((from->pl().getStops().size() == 0 ||
             to->pl().getStops().size() == 0) &&
            !isTriFace(e, g)) {
          if (combineNodes(from, to, g)) return true;
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

  std::cerr << "Combining " << a << " and " << b << " at node " << n
            << std::endl;

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

  // look for lines on {a, b} which serve as a continuation blocker,
  // for example
  //     1->      <-1
  // a ------> b -----> c
  // these have to be integrated as explicit connection exceptions to b
  // to prevent their loss later on
  //
  // this will lead to some superfluous connection exceptions - it may be
  // good to filter them out at some point :)

  explicateNonCons(connecting, b);
  explicateNonCons(connecting, a);

  // remove restrictions in a which are circumvented because b is a terminus
  terminusPass(a, connecting);

  // remove restrictions in a which are circumvented because b is a terminus
  terminusPass(b, connecting);

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
            b->pl().addConnExc(ro.first, e, exTo);
          }
        } else if (exTo == connecting) {
          for (auto e : b->getAdjList()) {
            if (e == connecting) continue;
            if (e == exFr.first) continue;
            b->pl().addConnExc(ro.first, e, exFr.first);
          }
        } else {
          b->pl().addConnExc(ro.first, exFr.first, exTo);
        }
      }
    }
  }

  std::vector<std::pair<const Route*,
                        std::pair<const TransitEdge*, const TransitEdge*>>>
      toDel;

  // go through all route exceptions in b
  for (const auto& ro : b->pl().getConnExc()) {
    for (const auto& exFr : ro.second) {
      for (const auto& exTo : exFr.second) {
        if (exFr.first == connecting) {
          for (auto e : a->getAdjList()) {
            if (e == connecting) continue;
            if (e == exTo) continue;
            b->pl().addConnExc(ro.first, e, exTo);
          }
          toDel.push_back({ro.first, {exFr.first, exTo}});
        } else if (exTo == connecting) {
          for (auto e : a->getAdjList()) {
            if (e == connecting) continue;
            if (e == exFr.first) continue;
            b->pl().addConnExc(ro.first, e, exFr.first);
            toDel.push_back({ro.first, {exFr.first, exTo}});
          }
          toDel.push_back({ro.first, {exFr.first, exTo}});
        } else {
          // not needed!
          // b->pl().addConnExc(ro.first, exFr.first, exTo);
        }
      }
    }
  }

  for (auto* oldE : a->getAdjList()) {
    if (oldE->getFrom() != a) continue;
    if (connecting == oldE) continue;

    std::cerr << a->getAdjList().size() << std::endl;

    assert(b != oldE->getTo());
    auto* newE = g->getEdg(b, oldE->getTo());

    if (!newE) {
      // add a new edge going from b to the non-a node
      newE = g->addEdg(b, oldE->getTo(), oldE->pl());

      // update route dirs
      routeDirRepl(a, b, newE);

      // replace each occurance of oldE in the non-a node with the newE
      edgeRpl(newE->getTo(), oldE, newE);
      edgeRpl(newE->getFrom(), oldE, newE);
    } else {
      // edge is already existing
      foldEdges(oldE, newE);

      // update route dirs
      routeDirRepl(a, b, newE);

      // replace each occurance of oldE in the non-a node with the newE
      edgeRpl(newE->getTo(), oldE, newE);
      edgeRpl(newE->getFrom(), oldE, newE);
    }
  }

  for (auto* oldE : a->getAdjList()) {
    if (oldE->getTo() != a) continue;
    if (connecting == oldE) continue;

    assert(b != oldE->getFrom());
    auto* newE = g->getEdg(oldE->getFrom(), b);

    if (!newE) {
      newE = g->addEdg(oldE->getFrom(), b, oldE->pl());

      // update route dirs
      routeDirRepl(a, b, newE);

      // replace each occurance of oldE in the non-a node with the newE
      edgeRpl(newE->getTo(), oldE, newE);
      edgeRpl(newE->getFrom(), oldE, newE);
    } else {
      // edge is already existing
      foldEdges(oldE, newE);

      // update route dirs
      routeDirRepl(a, b, newE);

      // replace each occurance of oldE in the non-a node with the newE
      edgeRpl(newE->getTo(), oldE, newE);
      edgeRpl(newE->getFrom(), oldE, newE);
    }
  }

  // remove the original exceptions
  for (const auto& ro : toDel) {
    b->pl().delConnExc(ro.first, ro.second.first, ro.second.second);
  }

  g->delEdg(a, b);
  g->delNd(a);

  cleanEx(b);

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
  if (oldE == newE) return;
  // replace in from
  for (auto& r : n->pl().getConnExc()) {
    auto exFr = r.second.begin();
    while (exFr != r.second.end()) {
      if (exFr->first == oldE) {
        std::swap(r.second[newE], exFr->second);
        exFr = r.second.erase(exFr);
      } else {
        exFr++;
      }
    }
  }

  // replace in to
  for (auto& r : n->pl().getConnExc()) {
    for (auto& exFr : r.second) {
      auto exTo = exFr.second.begin();
      while (exTo != exFr.second.end()) {
        if (*exTo == oldE) {
          exFr.second.insert(newE);
          exTo = exFr.second.erase(exTo);
        } else {
          exTo++;
        }
      }
    }
  }
}

// _____________________________________________________________________________
void Builder::routeDirRepl(TransitNode* oldN, TransitNode* newN,
                           TransitEdge* e) const {
  auto ro = e->pl().getRoutes().begin();
  while (ro != e->pl().getRoutes().end()) {
    if (ro->direction == oldN) {
      shared::transitgraph::RouteOcc newRo = *ro;
      newRo.direction = newN;

      // delete old
      ro = e->pl().getRoutes().erase(ro);

      // add new
      e->pl().getRoutes().insert(newRo);
    } else {
      ro++;
    }
  }
}

// _____________________________________________________________________________
bool Builder::isTriFace(const TransitEdge* a, const TransitGraph* g) const {
  // checks whether an edge is part of something like this:
  //         1
  //     a -----> b
  //  1  |        ^
  //     c -------|
  //          2
  // where the contraction of any two nodes would lead to a multigraph.

  std::set<TransitNode *> frNds, toNds;

  for (auto e : a->getFrom()->getAdjList()) {
    frNds.insert(e->getOtherNd(a->getFrom()));
  }

  for (auto e : a->getTo()->getAdjList()) {
    toNds.insert(e->getOtherNd(a->getTo()));
  }

  std::vector<TransitNode*> iSect;
  std::set_intersection(frNds.begin(), frNds.end(), toNds.begin(), toNds.end(),
                   std::back_inserter(iSect));

  if (iSect.size() == 0) return false;

  for (auto nd : iSect) {
    if (g->getEdg(nd, a->getFrom())
                 ->pl()
                 .getPolyline()
                 .getLength() > MIN_SEG_LENGTH) return true;
    if (g->getEdg(nd, a->getTo())
                 ->pl()
                 .getPolyline()
                 .getLength() > MIN_SEG_LENGTH) return true;
  }
  return false;
}

// _____________________________________________________________________________
void Builder::explicateNonCons(const TransitEdge* m, TransitNode* hub) const {
  assert(m->getFrom() == hub || m->getTo() == hub);

  for (const auto& ro : m->pl().getRoutes()) {
    if (ro.direction == 0) continue;
    for (auto ea : hub->getAdjList()) {
      if (ea == m) continue;
      if (!ea->pl().hasRoute(ro.route)) continue;
      const auto& ero = ea->pl().getRouteOcc(ro.route);

      assert(ro.direction == m->getTo() || ro.direction == m->getFrom());
      assert(ero.direction == 0 || ero.direction == ea->getTo() || ero.direction == ea->getFrom());

      if (ro.direction == ero.direction ||
          (ero.direction != 0 && m->getOtherNd(ro.direction) == ea->getOtherNd(ero.direction))) {
        hub->pl().addConnExc(ro.route, m, ea);
      }
    }
  }

  for (auto ea : hub->getAdjList()) {
    for (const auto& ro : ea->pl().getRoutes()) {
      if (ea == m) continue;
      if (!m->pl().hasRoute(ro.route)) {
        hub->pl().addConnExc(ro.route, m, ea);
      }
    }
  }
}

// _____________________________________________________________________________
bool Builder::lostIn(const Route* r, const TransitEdge* a) const {
  auto fr = a->getFrom();
  auto to = a->getTo();

  auto aro = a->pl().getRouteOcc(r);

  bool lost = true;

  for (auto e : fr->getAdjList()) {
    if (e == a) continue;
    if (!e->pl().hasRoute(r)) continue;
    if (!fr->pl().connOccurs(r, e, a)) continue;

    const auto& ero = e->pl().getRouteOcc(r);
    if (ero.direction == 0) {
      lost = false;
      break;
    }
    if (aro.direction == 0) {
      lost = false;
      break;
    }

    if (ero.direction == fr && aro.direction != fr) {
      lost = false;
      break;
    }
    if (aro.direction == fr && ero.direction != fr) {
      lost = false;
      break;
    }
  }

  if (lost) return true;

  for (auto e : to->getAdjList()) {
    if (e == a) continue;
    if (!e->pl().hasRoute(r)) continue;
    if (!to->pl().connOccurs(r, e, a)) continue;

    const auto& ero = e->pl().getRouteOcc(r);
    if (ero.direction == 0) return false;
    if (aro.direction == 0) return false;

    if (ero.direction == to && aro.direction != to) return false;
    if (aro.direction == to && ero.direction != to) return false;
  }

  return true;
}

// _____________________________________________________________________________
bool Builder::foldEdges(TransitEdge* a, TransitEdge* b) {
  const auto shrNd = TransitGraph::sharedNode(a, b);
  const auto majNonShrNd = b->getOtherNd(shrNd);
  assert(shrNd);

  //
  //                  b
  //         shrNd --------> v
  //          \             /
  //           \           /
  //            \         /
  //           a \       /
  //              \     /
  //               \   /
  //                \ /
  //            majNonShrNd
  //


  std::cerr << "Folding " << a << " and " << b << std::endl;

  // b is the new edge

  explicateNonCons(a, shrNd);
  explicateNonCons(b, shrNd);

  // remove "lost lines" which are on one of the fold edges, but cannot
  // leave it
  // std::vector<const Route *> lostA, lostB;
  // for (auto ro : a->pl().getRoutes()) {
    // if (lostIn(ro.route, a)) {
      // lostA.push_back(ro.route);
    // }
  // }
  // for (auto ro : b->pl().getRoutes()) {
    // if (lostIn(ro.route, b)) {
      // lostB.push_back(ro.route);
    // }
  // }

  // for (auto r : lostA) a->pl().delRoute(r);
  // for (auto r : lostB) b->pl().delRoute(r);

  auto keptExcsShrd = keptExcs(shrNd, a, b);
  auto keptExcsMajNonShrd = keptExcs(majNonShrNd, a, b);

  shrNd->pl().getConnExc().clear();
  for (auto ex : keptExcsShrd) {
    shrNd->pl().addConnExc(ex.first, ex.second.first, ex.second.second);
  }

  majNonShrNd->pl().getConnExc().clear();
  for (auto ex : keptExcsMajNonShrd) {
    majNonShrNd->pl().addConnExc(ex.first, ex.second.first, ex.second.second);
  }

  for (auto ro : a->pl().getRoutes()) {
    if (!b->pl().hasRoute(ro.route)) {
      // simply add the route
      if (ro.direction == 0)
        b->pl().addRoute(ro.route, 0);
      else if (ro.direction == shrNd)
        b->pl().addRoute(ro.route, shrNd);
      else if (ro.direction != shrNd)
        b->pl().addRoute(ro.route, b->getOtherNd(shrNd));
    } else {
      auto old = b->pl().getRouteOcc(ro.route);
      if (ro.direction == 0 && old.direction != 0) {
        // now goes in both directions
        b->pl().delRoute(ro.route);
        b->pl().addRoute(ro.route, 0);
      }

      if (ro.direction == shrNd && old.direction != shrNd) {
        // now goes in both directions
        b->pl().delRoute(ro.route);
        b->pl().addRoute(ro.route, 0);
      }

      if (ro.direction != shrNd && old.direction == shrNd) {
        // now goes in both directions
        b->pl().delRoute(ro.route);
        b->pl().addRoute(ro.route, 0);
      }
    }
  }

  return true;
}

// _____________________________________________________________________________
std::vector<
    std::pair<const Route*, std::pair<const TransitEdge*, const TransitEdge*>>>
Builder::keptExcs(TransitNode* nd, const TransitEdge* a, const TransitEdge* b) {
  std::vector<std::pair<const Route*,
                        std::pair<const TransitEdge*, const TransitEdge*>>>
      ret;
  for (auto ex : nd->pl().getConnExc()) {
    for (auto fr : ex.second) {
      for (auto to : fr.second) {
        if (fr.first == to) continue;
        if (fr.first == a && to == b) continue;
        if (fr.first == b && to == a) continue;

        // the exception does not concern the folded edges, so we
        // keep it
        if (a != fr.first && a != to && b != fr.first && b != to) {
          ret.push_back({ex.first, {fr.first, to}});
          continue;
        }

        // if the route is not contained in any edge, continue
        if (!a->pl().hasRoute(ex.first) && !b->pl().hasRoute(ex.first))
          continue;

        // if the route is only contained in one edge
        if (a->pl().hasRoute(ex.first) ^ b->pl().hasRoute(ex.first)) {
          if (a->pl().hasRoute(ex.first) && a != fr.first && a != to) {
            // but the edge doesnt occur in the exception, so drop it
            continue;
          }
          if (b->pl().hasRoute(ex.first) && b != fr.first && b != to) {
            // but the edge doesnt occur in the exceptions, so drop it
            continue;
          }
          auto efr = fr.first;
          auto eto = to;
          if (efr == a) efr = b;
          if (eto == a) eto = b;

          ret.push_back({ex.first, {efr, eto}});
          continue;
        }

        // if the route is contained in both edges, we must check if the
        // exception is overridden by the other edge
        if (a->pl().hasRoute(ex.first) && b->pl().hasRoute(ex.first)) {
          if (fr.first == a || to == a) {
            // concerns a
            auto otherEdge = fr.first == a ? to : fr.first;
            if (nd->pl().connOccurs(ex.first, b, otherEdge)) {
              // okay, there is no exception.
              // Important: we can be sure that a connection is possible then,
              // because we explicated direction-induced restrictions before!

              // -> drop this exception!
            } else {
              // else, keep it
              auto nfr = fr.first;
              auto nto = to;
              if (nfr == a) nfr = b;
              if (nto == a) nto = b;
              ret.push_back({ex.first, {nfr, nto}});
            }

          } else if (fr.first == b || to == b) {
            // concerns b

            auto otherEdge = fr.first == b ? to : fr.first;
            if (nd->pl().connOccurs(ex.first, a, otherEdge)) {
              // okay, there is no exception.
              // Important: we can be sure that a connection is possible then,
              // because we explicated direction-induced restrictions before!

              // -> drop this exception!
            } else {
              // else, keep it
              auto nfr = fr.first;
              auto nto = to;
              if (nfr == a) nfr = b;
              if (nto == a) nto = b;
              ret.push_back({ex.first, {nfr, nto}});
            }
          }
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
void Builder::terminusPass(TransitNode* nd, const TransitEdge* edg) {
  assert(nd == edg->getFrom() || nd == edg->getTo());

  auto otherNd = edg->getOtherNd(nd);

  std::vector<std::pair<const Route*,
                        std::pair<const TransitEdge*, const TransitEdge*>>>
      toDel;

  for (const auto& ro : nd->pl().getConnExc()) {
    for (const auto& exFr : ro.second) {
      for (const auto& exTo : exFr.second) {
        // don't handle exceptions to/from the connecting edge
        if (exFr.first == edg || exTo == edg) continue;

        if (edg->pl().hasRoute(ro.first)) {
          auto edgRo = edg->pl().getRouteOcc(ro.first);

          // only handle lines which go in both directions on edg
          if (edgRo.direction != 0) continue;

          // check if line has a terminus in nd
          bool term = true;
          for (auto eEdg : otherNd->getAdjList()) {
            if (eEdg == edg) continue;
            if (eEdg->pl().hasRoute(edgRo.route) &&
                otherNd->pl().connOccurs(edgRo.route, eEdg, edg)) {
              term = false;
              break;
            }
          }

          if (term) {
            toDel.push_back({ro.first, {exFr.first, exTo}});
          }
        }
      }
    }
  }

  for (const auto& ro : toDel) {
    nd->pl().delConnExc(ro.first, ro.second.first, ro.second.second);
  }
}

// _____________________________________________________________________________
void Builder::plMerge(TransitEdgePL& a, const TransitEdgePL& b) const {
  // merge routes
  for (auto ro : b.getRoutes()) {
    if (!a.hasRoute(ro.route)) {
      a.addRoute(ro.route, ro.direction);
    } else {
      auto old = a.getRouteOcc(ro.route);
      if (ro.direction != old.direction) {
        // now goes in both directions
        a.delRoute(ro.route);
        a.addRoute(ro.route, 0);
      }
    }
  }

  // merge geom
  if (a.getPolyline().getLine().size() == 0) {
    // a did not yet have a polyline
    a.setPolyline(b.getPolyline());
  } else {
    assert(false);
    // average the two
    std::vector<const PolyLine<double>*> lines{&a.getPolyline(),
                                               &b.getPolyline()};
    a.setPolyline(PolyLine<double>::average(lines));
  }
}

// _____________________________________________________________________________
 std::pair<TransitEdge*, TransitEdge*> Builder::split(TransitEdgePL& a, TransitNode* fr, TransitNode* to,
                            TransitGraph* g) const {
  TransitEdge* ret;
  auto right = a.getPolyline().getSegment(0.5, 1);
  a.setPolyline(a.getPolyline().getSegment(0, 0.5));
  auto helper = g->addNd(a.getPolyline().back());
  ret = g->addEdg(fr, helper);
  auto ro = a.getRoutes().begin();
  auto helperEdg = g->addEdg(helper, to, right);

  while (ro != a.getRoutes().end()) {
    if (ro->direction == to) {
      auto* route = ro->route;  // store because of deletion below
      ro = a.getRoutes().erase(ro);
      a.addRoute(route, helper);
      helperEdg->pl().addRoute(route, to);
    } else if (ro->direction == fr) {
      helperEdg->pl().addRoute(ro->route, helper);
      ro++;
    } else {
      helperEdg->pl().addRoute(ro->route, 0);
      ro++;
    }
  }

  return std::pair<TransitEdge*, TransitEdge*>(ret,helperEdg);
}

// _____________________________________________________________________________
bool Builder::snapElig(TransitNode* nat, const TransitEdge* eOther,
                       const TransitGraph* g) const {
  return (!g->getEdg(nat, eOther->getTo()) ||
          !g->getEdg(nat, eOther->getFrom()) ||
          g->getEdg(nat, eOther->getTo()) == g->getEdg(nat, eOther->getFrom()));
}

// _____________________________________________________________________________
void Builder::cleanEx(TransitNode* nd) const {
  std::vector<std::pair<const Route*,
                        std::pair<const TransitEdge*, const TransitEdge*>>>
      toDel;

  // go through all route exceptions in b
  for (const auto& ro : nd->pl().getConnExc()) {
    for (const auto& exFr : ro.second) {
      for (const auto& exTo : exFr.second) {
        std::cerr << nd << ", exFr " << exFr.first->getFrom() << ", " << exFr.first->getTo() << ", exTo "  << exTo->getFrom() << ", " << exTo->getTo() << std::endl;
        if (!(exFr.first->getFrom() == nd || exFr.first->getTo() == nd)) {
          std::cerr << "ALERT!" << std::endl;
          assert(false);
        }
        if (!(exTo->getFrom() == nd || exTo->getTo() == nd)) {
          std::cerr << "ALERT!" << std::endl;
          assert(false);
        }

        // exception from edge to itself
        if (exFr.first == exTo) toDel.push_back({ro.first, {exFr.first, exTo}});

        // exceptions for routes that don't occur in one of the edges
        else if (!exTo->pl().hasRoute(ro.first)) toDel.push_back({ro.first, {exFr.first, exTo}});
        else if (!exFr.first->pl().hasRoute(ro.first)) toDel.push_back({ro.first, {exFr.first, exTo}});
        else {
          auto roFr = exFr.first->pl().getRouteOcc(ro.first);
          auto roTo = exTo->pl().getRouteOcc(ro.first);

          if (roFr.direction == 0 || roTo.direction == 0) continue;

          if (roFr.direction == roTo.direction) toDel.push_back({ro.first, {exFr.first, exTo}});
          else if (exFr.first->getOtherNd(roFr.direction) == exTo->getOtherNd(roTo.direction)) toDel.push_back({ro.first, {exFr.first, exTo}});
        }
      }
    }
  }

  for (const auto& ro : toDel) {
    nd->pl().delConnExc(ro.first, ro.second.first, ro.second.second);
  }
}
