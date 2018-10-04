// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/OptGraph.h"

using namespace transitmapper;
using namespace optim;

// _____________________________________________________________________________
EtgPart OptGraph::getFirstEdge(const OptEdge* optEdg) {
  for (const auto e : optEdg->pl().etgs) {
    if (e.etg->getFrom() == optEdg->getFrom()->pl().node || e.etg->getTo() == optEdg->getFrom()->pl().node) {
      return e;
    }
  }
}

// _____________________________________________________________________________
EtgPart OptGraph::getLastEdge(const OptEdge* optEdg) {
  for (const auto e : optEdg->pl().etgs) {
    if (e.etg->getFrom() == optEdg->getTo()->pl().node || e.etg->getTo() == optEdg->getTo()->pl().node) {
      return e;
    }
  }
}

// _____________________________________________________________________________
std::string OptEdgePL::getStrRepr() const {
  const void* address = static_cast<const void*>(this);
  std::stringstream ss;
  ss << address;

  return ss.str();
}

// _____________________________________________________________________________
Edge* OptGraph::getAdjacentEdge(const OptEdge* e, const OptNode* n) {
  if (e->getFrom() == n) {
    return getFirstEdge(e).etg;
  } else {
    return getLastEdge(e).etg;
  }
}

// _____________________________________________________________________________
OptGraph::OptGraph(TransitGraph* toOptim)
 : _g(toOptim) { build(); }

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
void OptGraph::simplify() {
  while (simplifyStep()) {
  }
}

// _____________________________________________________________________________
bool OptGraph::simplifyStep() {
  for (OptNode* n : *getNds()) {
    if (n->getDeg() == 2) {
      OptEdge* first = 0;
      OptEdge* second = 0;

      for (OptEdge* e : n->getAdjList()) {
        if (!first)
          first = e;
        else
          second = e;
      }

      bool equal = first->pl().etgs.front().etg->getCardinality() ==
                   second->pl().etgs.front().etg->getCardinality();

      for (auto& to : *getAdjacentEdge(first, n)->getTripsUnordered()) {
        if (!getAdjacentEdge(second, n)
                 ->getSameDirRoutesIn(n->pl().node, to.route, to.direction,
                                      getAdjacentEdge(first, n))
                 .size()) {
          equal = false;
          break;
        }
      }

      if (equal) {
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
        // delEdg(first->getFrom(), first->getTo());
        // delEdg(second->getFrom(), second->getTo());

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
      ret +=1;
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
      for (const auto& to : *getFirstEdge(e).etg->getTripsUnordered()) {
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
      if (getFirstEdge(e).etg->getCardinality(true) > ret) {
        ret = getFirstEdge(e).etg->getCardinality(true);
      }
    }
  }

  return ret;
}
