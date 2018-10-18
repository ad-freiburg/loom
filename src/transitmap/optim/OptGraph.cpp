// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/OptGraph.h"
#include "util/String.h"
#include "util/graph/Algorithm.h"

using namespace transitmapper;
using namespace optim;
using transitmapper::graph::RouteOccurance;

// _____________________________________________________________________________
EtgPart OptGraph::getFirstEdg(const OptEdge* optEdg) {
  for (const auto e : optEdg->pl().etgs) {
    if (e.etg->getFrom() == optEdg->getFrom()->pl().node ||
        e.etg->getTo() == optEdg->getFrom()->pl().node) {
      return e;
    }
  }
  assert(false);
}

// _____________________________________________________________________________
EtgPart OptGraph::getLastEdg(const OptEdge* optEdg) {
  for (const auto e : optEdg->pl().etgs) {
    if (e.etg->getFrom() == optEdg->getTo()->pl().node ||
        e.etg->getTo() == optEdg->getTo()->pl().node) {
      return e;
    }
  }
  assert(false);
}

// _____________________________________________________________________________
const std::vector<transitmapper::graph::RouteOccurance>& OptEdgePL::getRoutes()
    const {
  if (!partialRoutes.size()) return *etgs[0].etg->getRoutes();

  return partialRoutes;
}

// _____________________________________________________________________________
size_t OptEdgePL::getCardinality() const {
  if (!partialRoutes.size()) {
    assert(etgs.size());
    return etgs[0].etg->getCardinality(true);
  }

  size_t ret = 0;

  for (const auto& ro : getRoutes()) {
    if (ro.route->relativeTo() == 0) ret++;
  }

  return ret;
}

// _____________________________________________________________________________
std::string OptEdgePL::getStrRepr() const {
  const void* address = static_cast<const void*>(this);
  std::stringstream ss;
  ss << address;

  return ss.str();
}

// _____________________________________________________________________________
Edge* OptGraph::getAdjEdg(const OptEdge* e, const OptNode* n) {
  if (e->getFrom() == n) {
    return getFirstEdg(e).etg;
  } else if (e->getTo() == n) {
    return getLastEdg(e).etg;
  }

  return 0;
}

// _____________________________________________________________________________
OptGraph::OptGraph(TransitGraph* toOptim) : _g(toOptim) { build(); }

// _____________________________________________________________________________
void OptGraph::build() {
  for (graph::Node* n : *_g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      graph::Node* fromTn = e->getFrom();
      graph::Node* toTn = e->getTo();

      OptNode* from = getNodeForTransitNode(fromTn);
      OptNode* to = getNodeForTransitNode(toTn);

      if (!from) {
        from = addNd(fromTn);
      }

      if (!to) {
        to = addNd(toTn);
      }

      OptEdge* edge = addEdg(from, to);

      edge->pl().etgs.push_back(EtgPart(e, e->getTo() == toTn));
    }
  }

  writeEdgeOrder();
}

// _____________________________________________________________________________
OptNode* OptGraph::getNodeForTransitNode(const Node* tn) const {
  for (auto n : getNds()) {
    if (n->pl().node == tn) return n;
  }

  return 0;
}

// _____________________________________________________________________________
void OptGraph::split() {
  std::vector<OptEdge*> toCut;

  // collect edges to cut
  for (OptNode* n : *getNds()) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;

      // we only cut if both nodes have a deg > 1 - otherwise, we cut
      // forever...
      if (e->getFrom()->getDeg() < 2 || e->getTo()->getDeg() < 2) continue;

      if (e->pl().getCardinality() == 1) toCut.push_back(e);
    }
  }

  for (OptEdge* e : toCut) {
    // cut this edge

    // add two new opt nodes, with slight offsets at the center of this
    // edge (this is for debug visualization purposes only)

    util::geo::PolyLine<double> edgeGeom(*e->getFrom()->pl().getGeom(),
                                         *e->getTo()->pl().getGeom());

    OptNode* leftN = addNd(edgeGeom.getPointAt(0.45).p);
    OptNode* rightN = addNd(edgeGeom.getPointAt(0.55).p);
    OptNode* eFrom = e->getFrom();
    OptNode* eTo = e->getTo();

    addEdg(eFrom, leftN, e->pl());
    auto newE = addEdg(rightN, eTo, e->pl());
    for (auto& etg : newE->pl().etgs) etg.wasCut = true;

    delEdg(eFrom, eTo);

    updateEdgeOrder(leftN);
    updateEdgeOrder(rightN);

    updateEdgeOrder(eFrom);
    updateEdgeOrder(eTo);
  }
}

// _____________________________________________________________________________
void OptGraph::simplify() {
  while (simplifyStep()) {
  }
}

// _____________________________________________________________________________
void OptGraph::untangle() {
  while (untangleYStep()) {
  }

  while (untangleDogBoneStep()) {
  }
}

// _____________________________________________________________________________
bool OptGraph::simplifyStep() {
  for (OptNode* n : *getNds()) {
    if (n->getDeg() == 2) {
      OptEdge* first = n->getAdjList().front();
      OptEdge* second = n->getAdjList().back();

      assert(n->pl().node);

      if (dirRouteEqualIn(first, n, second)) {
        OptNode* newFrom = 0;
        OptNode* newTo = 0;

        bool firstReverted;
        bool secondReverted;

        // add new edge
        if (first->getTo() != n) {
          newFrom = first->getTo();
          firstReverted = true;
        } else {
          newFrom = first->getFrom();
          firstReverted = false;
        }

        if (second->getTo() != n) {
          newTo = second->getTo();
          secondReverted = false;
        } else {
          newTo = second->getFrom();
          secondReverted = true;
        }

        if (newFrom == newTo) continue;

        OptEdge* newEdge = addEdg(newFrom, newTo);

        // add etgs...
        for (EtgPart& etgp : first->pl().etgs) {
          newEdge->pl().etgs.push_back(
              EtgPart(etgp.etg, (etgp.dir ^ firstReverted), etgp.order, etgp.wasCut));
        }

        for (EtgPart& etgp : second->pl().etgs) {
          newEdge->pl().etgs.push_back(
              EtgPart(etgp.etg, (etgp.dir ^ secondReverted), etgp.order, etgp.wasCut));
        }

        newEdge->pl().depth = std::max(first->pl().depth, second->pl().depth);

        assert(newFrom != n);
        assert(newTo != n);

        delNd(n);

        updateEdgeOrder(newFrom);
        updateEdgeOrder(newTo);

        return true;
      }
    }
  }
  return false;
}

// _____________________________________________________________________________
TransitGraph* OptGraph::getGraph() const { return _g; }

// _____________________________________________________________________________
size_t OptGraph::getNumNodes() const { return getNds().size(); }

// _____________________________________________________________________________
size_t OptGraph::getNumNodes(bool topo) const {
  size_t ret = 0;
  for (auto n : getNds()) {
    if ((n->pl().node->getStops().size() == 0) ^ !topo) ret++;
  }

  return ret;
}

// _____________________________________________________________________________
size_t OptGraph::getNumEdges() const {
  size_t ret = 0;

  for (auto n : getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      ret += 1;
    }
  }

  return ret;
}

// _____________________________________________________________________________
size_t OptGraph::getNumRoutes() const {
  std::set<const graph::Route*> routes;

  for (auto n : getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      for (const auto& to : *getFirstEdg(e).etg->getRoutes()) {
        if (to.route->relativeTo()) continue;
        routes.insert(to.route);
      }
    }
  }
  return routes.size();
}

// _____________________________________________________________________________
size_t OptGraph::getMaxCardinality() const {
  size_t ret = 0;
  for (auto n : getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (getFirstEdg(e).etg->getCardinality(true) > ret) {
        ret = getFirstEdg(e).etg->getCardinality(true);
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
const util::geo::Line<double>* OptEdgePL::getGeom() { return 0; }

// _____________________________________________________________________________
util::json::Dict OptEdgePL::getAttrs() {
  util::json::Dict ret;
  std::string lines;
  for (const auto& r : getRoutes()) {
    if (r.route->relativeTo()) lines += "(" + r.route->relativeTo()->getLabel() + "+" + r.route->getLabel() + "), ";
    else lines += r.route->getLabel() + ", ";
  }
  ret["lines"] = lines;
  ret["depth"] = util::toString(depth);
  return ret;
}

// _____________________________________________________________________________
const util::geo::Point<double>* OptNodePL::getGeom() { return &p; }

// _____________________________________________________________________________
util::json::Dict OptNodePL::getAttrs() {
  util::json::Dict ret;
  if (!node) {
    ret["stat_name"] = "<dummy>";
  } else {
    if (node->getStops().size()) {
      ret["stat_name"] = node->getStops()[0].name;
    } else {
      ret["stat_name"] = "";
    }
  }
  std::string edges;
  for (auto e : orderedEdges) {
    edges += util::toString(e) + ",";
  }

  ret["edge_order"] = edges;
  return ret;
}

// _____________________________________________________________________________
bool OptGraph::untangleYStep() {
  double DO = 400;  // only relevant for debug output
  for (OptNode* na : *getNds()) {
    if (na->getDeg() != 1) continue;  // only look at terminus nodes

    // the only outgoing edge
    OptEdge* ea = na->getAdjList().front();
    OptNode* nb = ea->getOtherNd(na);

    if (isYAt(ea, nb)) {
      assert(nb->pl().node);
      assert(na->pl().node);

      // the geometry of the main leg
      util::geo::PolyLine<double> pl(*nb->pl().getGeom(), *na->pl().getGeom());
      double bandW = (nb->getDeg() - 1) * (DO / (ea->pl().depth + 1));
      auto orthoPl = pl.getOrthoLineAtDist(0, bandW);
      auto orthoPlOrig = pl.getOrthoLineAtDist(pl.getLength(), bandW);

      std::vector<OptNode*> centerNds(nb->getDeg() - 1);
      std::vector<OptNode*> origNds(nb->getDeg() - 1);

      // for each minor leg of the Y, create a new node at the Y center
      for (size_t i = 0; i < centerNds.size(); i++) {
        double p = (centerNds.size() - 1 - i) / (double)(centerNds.size());
        centerNds[i] = addNd(orthoPl.getPointAt(p).p);
        centerNds[i]->pl().node = nb->pl().node;
      }

      // for each minor leg of the Y, create a new node at the Y center
      for (size_t i = 0; i < origNds.size(); i++) {
        double p = (origNds.size() - 1 - i) / (double)(origNds.size());
        origNds[i] = addNd(orthoPlOrig.getPointAt(p).p);
        origNds[i]->pl().node = na->pl().node;
      }

      // each leg, in clockwise fashion
      auto minLgs = clockwEdges(ea, nb);
      assert(minLgs.size() == nb->pl().orderedEdges.size() - 1);

      for (size_t i = 0; i < minLgs.size(); i++) {
        OptEdge* newLeg = 0;
        if (minLgs[i]->getFrom() == nb)
          newLeg = addEdg(centerNds[i], minLgs[i]->getTo(), minLgs[i]->pl());
        else
          newLeg = addEdg(minLgs[i]->getFrom(), centerNds[i], minLgs[i]->pl());
        delEdg(minLgs[i]->getFrom(), minLgs[i]->getTo());
        minLgs[i] = newLeg;
      }

      size_t offset = 0;
      for (size_t i = 0; i < minLgs.size(); i++) {
        size_t j = i;
        if (ea->getFrom() == nb) {
          if (ea->pl().etgs[0].dir) j = minLgs.size() - 1 - i;
          addEdg(centerNds[j], origNds[j], getView(ea, nb, minLgs[j], offset));
        } else {
          if (!ea->pl().etgs[0].dir) j = minLgs.size() - 1 - i;
          addEdg(origNds[j], centerNds[j], getView(ea, nb, minLgs[j], offset));
        }
        offset += minLgs[j]->pl().getRoutes().size();
      }

      // delete remaining stuff
      delNd(nb);
      delNd(na);

      // update orderings
      for (auto n : centerNds) updateEdgeOrder(n);
      for (auto n : origNds) updateEdgeOrder(n);

      return true;
    }
  }
  return false;
}

// _____________________________________________________________________________
OptEdgePL OptGraph::getView(OptEdge* parent, OptNode* origin, OptEdge* leg,
                            size_t offset) {
  OptEdgePL ret(parent->pl());
  ret.depth++;

  for (auto& etg : ret.etgs) {
    if (parent->pl().etgs[0].dir ^ etg.dir)
      etg.order += parent->pl().getRoutes().size() -
                   leg->pl().getRoutes().size() - offset;
    else
      etg.order += offset;
  }

  for (auto ro : leg->pl().getRoutes()) {
    auto e = getAdjEdg(parent, origin);
    auto routes = e->getCtdRoutesIn(origin->pl().node, ro.route, ro.direction,
                                    getAdjEdg(leg, origin));
    ret.partialRoutes.insert(ret.partialRoutes.begin(), routes.begin(),
                             routes.end());
  }

  return ret;
}

// _____________________________________________________________________________
bool OptGraph::untangleDogBoneStep() {
  double DO = 400;  // only relevant for debug output
  for (OptNode* na : *getNds()) {
    if (na->getDeg() < 3) continue;  // only look at terminus nodes

    for (OptEdge* mainLeg : na->getAdjList()) {
      if (mainLeg->getFrom() != na) continue;
      if (isDogBone(mainLeg)) {
        OptNode* nb = mainLeg->getOtherNd(na);
        // the geometry of the main leg
        util::geo::PolyLine<double> pl(*na->pl().getGeom(),
                                       *nb->pl().getGeom());
        double bandW = (nb->getDeg() - 1) * (DO / (mainLeg->pl().depth + 1));
        auto orthoPlA = pl.getOrthoLineAtDist(0, bandW);
        auto orthoPlB = pl.getOrthoLineAtDist(pl.getLength(), bandW);

        std::vector<OptNode*> aNds(na->getDeg() - 1);
        std::vector<OptNode*> bNds(nb->getDeg() - 1);

        // map positions in b to positions in a
        std::vector<size_t> bToA(nb->getDeg() - 1);

        // for each minor leg at node a, create a new node
        for (size_t i = 0; i < aNds.size(); i++) {
          double p = (aNds.size() - 1 - i) / (double)(aNds.size());
          aNds[i] = addNd(orthoPlA.getPointAt(p).p);
          aNds[i]->pl().node = na->pl().node;
        }

        // each leg in a, in clockwise fashion
        auto minLgsA = clockwEdges(mainLeg, na);
        assert(minLgsA.size() == na->pl().orderedEdges.size() - 1);

        // for each minor leg at node b, create a new node
        for (size_t i = 0; i < bNds.size(); i++) {
          double p = (bNds.size() - 1 - i) / (double)(bNds.size());
          bNds[i] = addNd(orthoPlB.getPointAt(p).p);
          bNds[i]->pl().node = nb->pl().node;
        }

        // each leg in b, in counter-clockwise fashion
        auto minLgsB = clockwEdges(mainLeg, nb);
        assert(minLgsB.size() == nb->pl().orderedEdges.size() - 1);
        std::reverse(minLgsB.begin(), minLgsB.end());

        if (!dirContinuedOver(minLgsA.front(), mainLeg, minLgsB.front())) {
          std::cout << "NOT: "
                    << minLgsA.front()->pl().getRoutes().front().route << " -> "
                    << minLgsB.front()->pl().getRoutes().front().route
                    << std::endl;
          for (size_t i = 0; i < bNds.size(); i++)
            bToA[i] = bNds.size() - 1 - i;
        } else {
          for (size_t i = 0; i < bNds.size(); i++) bToA[i] = i;
        }

        for (size_t i = 0; i < minLgsA.size(); i++) {
          OptEdge* newLeg = 0;
          if (minLgsA[i]->getFrom() == na)
            newLeg = addEdg(aNds[i], minLgsA[i]->getTo(), minLgsA[i]->pl());
          else
            newLeg = addEdg(minLgsA[i]->getFrom(), aNds[i], minLgsA[i]->pl());
          delEdg(minLgsA[i]->getFrom(), minLgsA[i]->getTo());
          minLgsA[i] = newLeg;
        }

        for (size_t i = 0; i < minLgsB.size(); i++) {
          OptEdge* newLeg = 0;
          if (minLgsB[i]->getFrom() == nb)
            newLeg = addEdg(bNds[i], minLgsB[i]->getTo(), minLgsB[i]->pl());
          else
            newLeg = addEdg(minLgsB[i]->getFrom(), bNds[i], minLgsB[i]->pl());
          delEdg(minLgsB[i]->getFrom(), minLgsB[i]->getTo());
          minLgsB[i] = newLeg;
        }

        size_t offset = 0;
        for (size_t i = 0; i < minLgsA.size(); i++) {
          size_t j = i;
          if (mainLeg->getFrom() == na) {
            if (mainLeg->pl().etgs[0].dir) j = minLgsA.size() - 1 - i;
            addEdg(aNds[j], bNds[bToA[j]],
                   getView(mainLeg, na, minLgsA[j], offset));
          } else {
            if (!mainLeg->pl().etgs[0].dir) j = minLgsA.size() - 1 - i;
            addEdg(bNds[bToA[j]], aNds[j],
                   getView(mainLeg, na, minLgsA[j], offset));
          }
          offset += minLgsA[j]->pl().getRoutes().size();
        }

        // delete remaining stuff
        delNd(nb);
        delNd(na);

        // update orderings
        for (auto n : aNds) updateEdgeOrder(n);
        for (auto n : bNds) updateEdgeOrder(n);

        return true;
      }
    }
  }

  return false;
}

// _____________________________________________________________________________
void OptGraph::writeEdgeOrder() {
  for (auto nd : *getNds()) {
    updateEdgeOrder(nd);
  }
}

// _____________________________________________________________________________
void OptGraph::updateEdgeOrder(OptNode* n) {
  n->pl().orderedEdges.clear();

  if (n->getDeg() == 1) {
    n->pl().orderedEdges.push_back(n->getAdjList().front());
    return;
  }

  for (auto e : n->getAdjList()) {
    n->pl().orderedEdges.push_back(e);
  }
  std::sort(n->pl().orderedEdges.begin(), n->pl().orderedEdges.end(), cmpEdge);
}

// _____________________________________________________________________________
bool OptGraph::dirRouteContains(const OptEdge* a, const OptNode* dirN,
                                const OptEdge* b) {
  for (auto& to : b->pl().getRoutes()) {
    if (!getAdjEdg(a, dirN)
             ->getSameDirRoutesIn(dirN->pl().node, to.route, to.direction,
                                  getAdjEdg(b, dirN))
             .size()) {
      return false;
    }
  }
  return true;
}

// _____________________________________________________________________________
bool OptGraph::dirRouteEqualIn(const OptEdge* a, const OptNode* dirN,
                               const OptEdge* b) {
  if (a->pl().getCardinality() != b->pl().getCardinality()) return false;

  for (auto& to : b->pl().getRoutes()) {
    if (!getAdjEdg(a, dirN)
             ->getSameDirRoutesIn(dirN->pl().node, to.route, to.direction,
                                  getAdjEdg(b, dirN))
             .size()) {
      return false;
    }
  }
  return true;
}

// _____________________________________________________________________________
bool OptGraph::isYAt(OptEdge* eLeg, OptNode* n) const {
  if (n->getDeg() < 3) return false;
  if (eLeg->getOtherNd(n)->getDeg() != 1) return false;
  // if (eLeg->pl().getCardinality() < 2) return false;

  for (OptEdge* e : n->getAdjList()) {
    if (e == eLeg) continue;

    if (!dirRouteContains(eLeg, n, e)) return false;
  }
  return true;
}

// _____________________________________________________________________________
bool OptGraph::isDogBone(OptEdge* leg) const {
  if (leg->getFrom()->getDeg() != leg->getTo()->getDeg()) return false;
  if (leg->getFrom()->getDeg() < 3) return false;

  for (OptEdge* ea : leg->getFrom()->getAdjList()) {
    if (ea == leg) continue;
    bool found = false;
    for (OptEdge* eb : leg->getTo()->getAdjList()) {
      if (eb == leg) continue;
      if (dirContinuedOver(ea, leg, eb)) {
        found = true;
        break;
      }
    }
    if (!found) return false;
  }

  for (OptEdge* ea : leg->getTo()->getAdjList()) {
    if (ea == leg) continue;
    bool found = false;
    for (OptEdge* eb : leg->getFrom()->getAdjList()) {
      if (eb == leg) continue;
      if (dirContinuedOver(ea, leg, eb)) {
        found = true;
        break;
      }
    }
    if (!found) return false;
  }

  return true;
}

// _____________________________________________________________________________
std::vector<OptEdge*> OptGraph::clockwEdges(OptEdge* noon, OptNode* n) {
  std::vector<OptEdge*> clockwise;

  size_t passed = n->getDeg();
  for (size_t i = 0; i < n->pl().orderedEdges.size(); i++) {
    auto e = n->pl().orderedEdges[i];
    if (e == noon) {
      passed = i;
      continue;
    }
    if (passed > n->getDeg() - 1)
      clockwise.push_back(e);
    else {
      clockwise.insert(clockwise.begin() + (i - 1 - passed), e);
    }
  }

  return clockwise;
}

// _____________________________________________________________________________
bool OptGraph::dirContinuedOver(const OptEdge* a, const OptEdge* b,
                                const OptEdge* c) {
  OptNode* ab = sharedNode(a, b);
  OptNode* bc = sharedNode(b, c);
  assert(ab);
  assert(bc);

  for (auto& to : a->pl().getRoutes()) {
    for (auto& too : getAdjEdg(a, ab)->getSameDirRoutesIn(
             ab->pl().node, to.route, to.direction, getAdjEdg(b, ab))) {
      if (!getAdjEdg(c, bc)
               ->getSameDirRoutesIn(bc->pl().node, too.route, too.direction,
                                    getAdjEdg(c, bc))
               .size()) {
        return false;
      }
    }
  }
  return true;
}

// _____________________________________________________________________________
OptNode* OptGraph::sharedNode(const OptEdge* a, const OptEdge* b) {
  OptNode* r = 0;
  if (a->getFrom() == b->getFrom() || a->getFrom() == b->getTo())
    r = a->getFrom();
  if (a->getTo() == b->getFrom() || a->getTo() == b->getTo()) r = a->getTo();
  return r;
}
