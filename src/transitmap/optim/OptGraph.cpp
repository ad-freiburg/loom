// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/OptGraph.h"
#include "util/graph/Algorithm.h"

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
  } else {
    return getLastEdg(e).etg;
  }
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

      if (e->pl().etgs.front().etg->getCardinality(true) == 1) toCut.push_back(e);
    }
  }

  for (OptEdge* e : toCut) {
    // cut this edge

    // add two new opt nodes, with slight offsets at the center of this
    // edge (this is for debug visualization purposes only)

    util::geo::PolyLine<double> edgeGeom(*e->getFrom()->pl().getGeom(), *e->getTo()->pl().getGeom());

    OptNode* leftN = addNd(edgeGeom.getPointAt(0.45).p);
    OptNode* rightN = addNd(edgeGeom.getPointAt(0.55).p);

    OptEdge* sibling = addEdg(e->getFrom(), leftN, e->pl());
    addEdg(rightN, e->getTo(), e->pl())->pl().siameseSibl = sibling;

    delEdg(e->getFrom(), e->getTo());
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

      if (getAdjEdg(first, n)->dirRouteEqualIn(getAdjEdg(second, n), n->pl().node)) {
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

        newFrom->addEdge(newEdge);
        newTo->addEdge(newEdge);
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
const util::geo::Line<double>* OptEdgePL::getGeom() {
  return 0;
}

// _____________________________________________________________________________
util::json::Dict OptEdgePL::getAttrs() {
  util::json::Dict ret;
  std::string lines;
  for (const auto& r : *etgs.front().etg->getRoutes()) {
    lines += r.route->getLabel() + ", ";
  }
  ret["lines"] = lines;
  return ret;
}

// _____________________________________________________________________________
const util::geo::Point<double>* OptNodePL::getGeom() {
  if (node) return &node->getPos();
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

    if (nb->getDeg() < 2) continue;

    assert(nb->pl().node);

    bool isY = true;
    for (OptEdge* e : nb->getAdjList()) {
      if (e == ea) continue;

      // std::cout << "--" << std::endl;
      // std::cout << getAdjEdg(ea, nb)->toString() << std::endl;
      // std::cout << getAdjEdg(e, nb)->toString() << std::endl;

      if (!getAdjEdg(ea, nb)->dirRouteContains(getAdjEdg(e, nb), nb->pl().node)) {
        isY = false;
        break;
      }
    }

    if (isY) {
      std::cout << "Found Y at node " << nb << std::endl;

      // the geometry of the main leg
      util::geo::PolyLine<double> pl(*nb->pl().getGeom(), *na->pl().getGeom());
      auto orthoPl = pl.getOrthoLineAtDist(0, (nb->getDeg() - 1) * 10);

      std::vector<OptNode*> nds(nb->getDeg() - 1);

      // for each minor leg of the Y, create a new node
      for (size_t i = 0; i < nds.size(); i++) {
        nds[i] = addNd(orthoPl.getPointAt(((double)i) / (double)(nds.size())).p);
      }


      // delEdg(na, nb);

      return false;
    }

  }
  return false;
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
