// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/OptGraph.h"
#include "util/graph/Algorithm.h"
#include "util/String.h"

using namespace transitmapper;
using namespace optim;

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
EtgPart OptGraph::getLastEdg(const OptEdge* optEdg) { for (const auto e : optEdg->pl().etgs) {
    if (e.etg->getFrom() == optEdg->getTo()->pl().node ||
        e.etg->getTo() == optEdg->getTo()->pl().node) {
      return e;
    }
  }
  assert(false);
}

// _____________________________________________________________________________
const std::vector<transitmapper::graph::RouteOccurance>& OptEdgePL::getRoutes() const {
  if (!partialRoutes.size())
    return *etgs[0].etg->getRoutes();

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

      // we only cut if both nodes have a deg > 1
      if (e->getFrom()->getDeg() < 2 || e->getTo()->getDeg() < 2) continue;

      if (e->pl().etgs.front().etg->getCardinality(true) == 1)
        toCut.push_back(e);
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

    OptEdge* sibling = addEdg(eFrom, leftN, e->pl());
    addEdg(rightN, eTo, e->pl())->pl().siameseSibl = sibling;

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
}

// _____________________________________________________________________________
bool OptGraph::simplifyStep() {
  for (OptNode* n : *getNds()) {
    if (n->getDeg() == 2) {
      OptEdge* first = n->getAdjList().front();
      OptEdge* second = n->getAdjList().back();

      assert(n->pl().node);

      if (getAdjEdg(first, n)->dirRouteEqualIn(getAdjEdg(second, n),
                                               n->pl().node)) {
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
              EtgPart(etgp.etg, (etgp.dir ^ firstReverted)));
        }

        for (EtgPart& etgp : second->pl().etgs) {
          newEdge->pl().etgs.push_back(
              EtgPart(etgp.etg, (etgp.dir ^ secondReverted)));
        }

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
    lines += r.route->getLabel() + ", ";
  }
  ret["lines"] = lines;
  ret["order"] = util::toString(order);
  return ret;
}

// _____________________________________________________________________________
const util::geo::Point<double>* OptNodePL::getGeom() {
  return &p;
}

// _____________________________________________________________________________
util::json::Dict OptNodePL::getAttrs() {
  util::json::Dict ret;
  if (!node) {
    ret["stat_name"] = "<dummy>";
  } else {
    if (node->getStops().size()) {
      ret["stat_name"] = node->getStops().front().name;
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
  // for now, just untangle deg 3 nodes!
  for (OptNode* na : *getNds()) {
    if (na->getDeg() != 1) continue;  // only look at terminus nodes

    // the only outgoing edge
    OptEdge* ea = na->getAdjList().front();
    OptNode* nb = ea->getOtherNd(na);

    if (nb->getDeg() < 3) continue;

    assert(nb->pl().node);
    assert(na->pl().node);

    bool isY = true;
    for (OptEdge* e : nb->getAdjList()) {
      if (e == ea) continue;

      if (!getAdjEdg(ea, nb)->dirRouteContains(getAdjEdg(e, nb),
                                               nb->pl().node)) {
        isY = false;
        break;
      }
    }

    if (isY) {
      // the geometry of the main leg
      util::geo::PolyLine<double> pl(*nb->pl().getGeom(), *na->pl().getGeom());
      auto orthoPl = pl.getOrthoLineAtDist(0, (nb->getDeg() - 1) * 200);
      auto orthoPlOrig = pl.getOrthoLineAtDist(pl.getLength(), (nb->getDeg() - 1) * 200);

      std::vector<OptNode*> centerNds(nb->getDeg() - 1);
      std::vector<OptNode*> origNds(nb->getDeg() - 1);

      // for each minor leg of the Y, create a new node at the Y center
      for (size_t i = 0; i < centerNds.size(); i++) {
        centerNds[i] =
            addNd(orthoPl.getPointAt(((double)(centerNds.size() - 1 - i)) / (double)(centerNds.size())).p);
        centerNds[i]->pl().node = nb->pl().node;
      }

      // for each minor leg of the Y, create a new node at the Y center
      for (size_t i = 0; i < origNds.size(); i++) {
        origNds[i] =
            addNd(orthoPlOrig.getPointAt(((double)(origNds.size() - 1 - i)) / (double)(origNds.size())).p);
        origNds[i]->pl().node = na->pl().node;
      }

      // each leg, in clockwise fashion
      std::vector<OptEdge*> minLegs;

      size_t passed = nb->getDeg();
      for (size_t i = 0; i < nb->pl().orderedEdges.size(); i++) {
        auto e = nb->pl().orderedEdges[i];
        if (e == ea) {
          passed = i;
          continue;
        }
        if (passed > nb->getDeg() - 1) minLegs.push_back(e);
        else {
          std::cout << (passed + 1 - i) << std::endl; 
          minLegs.insert(minLegs.begin() + (i - 1 - passed), e);
        }
      }

      assert(minLegs.size() == nb->pl().orderedEdges.size() - 1);

      std::cout << std::endl;
      for (auto e : minLegs) std::cout << e->pl().getRoutes().front().route->getLabel() << std::endl;

      for (size_t i = 0; i < minLegs.size(); i++) {
        OptEdge* newLeg = 0;
        if (minLegs[i]->getFrom() == nb) newLeg = addEdg(centerNds[i], minLegs[i]->getTo(), minLegs[i]->pl());
        else newLeg = addEdg(minLegs[i]->getFrom(), centerNds[i], minLegs[i]->pl());
        delEdg(minLegs[i]->getFrom(), minLegs[i]->getTo());
        minLegs[i] = newLeg;
      }

      size_t offset = 0;
      for (size_t i = 0; i < minLegs.size(); i++) {
        size_t j = i;
        if (ea->getFrom() == nb) {
          if (ea->pl().etgs.front().dir) j = minLegs.size() - 1 - i;
          addEdg(centerNds[j], origNds[j], OptGraph::getOptEdgePLView(ea, nb, minLegs[j], offset));
        } else {
          if (!ea->pl().etgs.front().dir) j = minLegs.size() - 1 - i;
          addEdg(origNds[j], centerNds[j], OptGraph::getOptEdgePLView(ea, nb, minLegs[j], offset));
        }
        offset += minLegs[j]->pl().getRoutes().size();
      }

      delEdg(na, nb);
      delNd(nb);
      delNd(na);

      for (auto n : centerNds) {
        updateEdgeOrder(n);
      }

      for (auto n : origNds) {
        updateEdgeOrder(n);
      }

      return true;
    }
  }
  return false;
}

// _____________________________________________________________________________
OptEdgePL OptGraph::getOptEdgePLView(OptEdge* parent, OptNode* origin, OptEdge* leg, size_t offset) {
  OptEdgePL ret(parent->pl());

  ret.order = offset;

  for (auto ro : leg->pl().getRoutes()) {
    auto e = getAdjEdg(parent, origin);
    auto routes = e->getCtdRoutesIn(origin->pl().node, ro.route, ro.direction, getAdjEdg(leg, origin));
    ret.partialRoutes.insert(ret.partialRoutes.begin(), routes.begin(), routes.end());
  }


  return ret;
}

// _____________________________________________________________________________
bool OptGraph::untangleDogBoneStep() {
  // for now, just untangle deg 3 nodes!
  for (OptNode* n : *getNds()) {
    if (n->getDeg() != 3) continue;
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (e->getTo()->getDeg() != e->getFrom()->getDeg()) continue;
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
