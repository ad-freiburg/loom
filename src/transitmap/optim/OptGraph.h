// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_OPTIM_OPTGRAPH_H_
#define TRANSITMAP_GRAPH_OPTIM_OPTGRAPH_H_

#include <set>
#include <string>

#include "transitmap/graph/Edge.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/optim/Scorer.h"
#include "util/graph/UndirGraph.h"
#include "util/json/Writer.h"

using transitmapper::graph::TransitGraph;
using transitmapper::graph::Node;
using transitmapper::graph::Edge;
using util::graph::UndirGraph;

namespace transitmapper {
namespace optim {

struct OptNodePL;
struct OptEdgePL;

typedef util::graph::Node<OptNodePL, OptEdgePL> OptNode;
typedef util::graph::Edge<OptNodePL, OptEdgePL> OptEdge;

typedef std::map<const transitmapper::optim::OptEdge*, std::vector<const graph::Route*>>
    OptOrderingConfig;

struct OptRO {
  OptRO(const graph::Route* r, const Node* dir) : route(r), direction(dir) {}
  const graph::Route* route;
  const Node* direction;  // 0 if in both directions

  bool operator==(const OptRO& b) const {
    return b.route == route;
  }
  bool operator==(const graph::RouteOccurance& b) const {
    return b.route == route;
  }
};

struct EtgPart {
  Edge* etg;
  bool dir;
  size_t order;

  // there is another edge with determines the
  // ordering in this edge - important to prevent double
  // writing of ordering later on
  bool wasCut;

  EtgPart(Edge* etg, bool dir) : etg(etg), dir(dir), order(0), wasCut(false){};
  EtgPart(Edge* etg, bool dir, size_t order, bool wasCut)
      : etg(etg), dir(dir), order(order), wasCut(wasCut){};
};

struct OptEdgePL {
  OptEdgePL() : depth(0), firstEtg(0), lastEtg(0) {};

  // all original ETGs from the transit graph contained in this edge
  // Guarantee: they are all equal in terms of (directed) routes
  std::vector<EtgPart> etgs;

  size_t depth;

  size_t firstEtg;
  size_t lastEtg;

  size_t getCardinality() const;
  std::string toStr() const;
  const std::vector<OptRO>& getRoutes() const;

  // partial routes
  // For the ETGs contained in .etgs, only these route occurances are
  // actually contained in this edge. Their relative ordering is defined by
  // .order
  std::vector<OptRO> routes;

  std::string getStrRepr() const;

  const util::geo::Line<double>* getGeom();
  util::json::Dict getAttrs();
};

struct OptNodePL {
  const Node* node;
  util::geo::Point<double> p;

  // the edges arriving at this node, in clockwise fashion, based
  // on the geometry in the original graph
  std::vector<OptEdge*> orderedEdges;

  OptNodePL(util::geo::Point<double> p) : node(0), p(p){};
  OptNodePL(const Node* node) : node(node), p(node->getPos()){};
  OptNodePL() : node(0){};

  const util::geo::Point<double>* getGeom();
  util::json::Dict getAttrs();
};

class OptGraph : public UndirGraph<OptNodePL, OptEdgePL> {
 public:
  OptGraph(TransitGraph* toOptim, const Scorer* scorer);

  TransitGraph* getGraph() const;

  size_t getNumNodes() const;
  size_t getNumNodes(bool topo) const;
  size_t getNumEdges() const;
  size_t getNumRoutes() const;
  size_t getMaxCardinality() const;

  double getMaxCrossPen() const;
  double getMaxSplitPen() const;

  void simplify();
  void untangle();

  static graph::Edge* getAdjEdg(const OptEdge* e, const OptNode* n);
  static EtgPart getAdjEtgp(const OptEdge* e, const OptNode* n);
  static bool hasCtdRoutesIn(const graph::Route* r,
                                              const Node* dir,
                                              const OptEdge* fromEdge,
                                              const OptEdge* toEdge);
  static std::vector<OptRO> getCtdRoutesIn(
      const graph::Route* r, const Node* dir, const OptEdge* fromEdge,
      const OptEdge* toEdge);
  static std::vector<OptRO> getSameDirRoutesIn(
      const graph::Route* r, const Node* dir, const OptEdge* fromEdge,
      const OptEdge* toEdge);
  static EtgPart getFirstEdg(const OptEdge*);
  static EtgPart getLastEdg(const OptEdge*);

  // apply splitting rules
  void split();

 private:
  TransitGraph* _g;
  const Scorer* _scorer;

  OptNode* getNodeForTransitNode(const Node* tn) const;

  void build();
  void writeEdgeOrder();
  void updateEdgeOrder(OptNode* n);
  bool simplifyStep();

  bool untangleFullCross();
  bool untangleYStep();
  bool untanglePartialYStep();
  bool untangleDogBoneStep();
  bool untanglePartialDogBoneStep();
  bool untangleStumpStep();

  std::vector<OptNode*> explodeNodeAlong(OptNode* nd, const PolyLine<double>& pl, size_t n);

  std::vector<OptEdge*> branchesAt(OptEdge* e, OptNode* n) const;
  bool branchesAtInto(OptEdge* e, OptNode* n,
                      std::vector<OptEdge*> branchesA) const;
  bool partiallyBranchesAtInto(OptEdge* e, OptNode* n,
                               std::vector<OptEdge*> branchesA) const;
  std::vector<OptEdge*> partiallyBranchesAt(OptEdge* e, OptNode* n) const;

  std::pair<OptEdge*, OptEdge*> isFullCross(OptNode* n) const;
  bool isYAt(OptEdge* e, OptNode* n) const;
  bool isPartialYAt(OptEdge* e, OptNode* n) const;
  OptEdge* isStump(OptEdge* e) const;
  OptEdge* isStumpAt(OptEdge* e, OptNode* n) const;

  bool isDogBone(OptEdge* e) const;
  OptNode* isPartialDogBone(OptEdge* e) const;


  static void upFirstLastEdg(OptEdge*);

  static OptEdgePL getView(OptEdge* parent, OptEdge* leg, size_t offset);
  static OptEdgePL getPartialView(OptEdge* parent, OptEdge* leg, size_t offset);

  std::vector<size_t> mapPositions(std::vector<OptEdge*> a, OptEdge* leg,
                                   std::vector<OptEdge*> b) const;

  static bool dirRouteEndsIn(const OptEdge* a, const OptEdge* b);
  static bool dirRouteContains(const OptEdge* a, const OptEdge* b);

  static bool dirRouteEqualIn(const OptEdge* a, const OptEdge* b);
  static bool dirContinuedOver(const OptEdge* a, const OptEdge* b,
                               const OptEdge* c);
  static bool dirPartialContinuedOver(const OptEdge* a, const OptEdge* b);
  static bool dirContinuedOver(const OptRO& ro,
                               const OptEdge* a, const OptEdge* b);

  static std::vector<OptRO> getCtdRoutesIn(
      const OptEdge* fromEdge, const OptEdge* toEdge);

  static OptNode* sharedNode(const OptEdge* a, const OptEdge* b);

  static std::vector<OptEdge*> clockwEdges(OptEdge* noon, OptNode* n);
  static std::vector<OptEdge*> partialClockwEdges(OptEdge* noon, OptNode* n);
};

// compare the orientation of two edges adjacent to some shared node
inline bool cmpEdge(const OptEdge* a, const OptEdge* b) {
  double angA, angB;

  // n is the shared node
  OptNode* n = 0;
  if (a->getFrom() == b->getFrom() || a->getFrom() == b->getTo())
    n = a->getFrom();
  else
    n = a->getTo();
  assert(n->pl().node);

  auto tgEdgeA = OptGraph::getAdjEdg(a, n);
  assert(tgEdgeA);
  assert(n->pl().node->getNodeFrontFor(tgEdgeA));

  angA = n->pl().node->getNodeFrontFor(tgEdgeA)->getOutAngle();

  auto tgEdgeB = OptGraph::getAdjEdg(b, n);
  assert(tgEdgeB);
  assert(n->pl().node->getNodeFrontFor(tgEdgeB));

  angB = n->pl().node->getNodeFrontFor(tgEdgeB)->getOutAngle();

  if (tgEdgeA == tgEdgeB) {
    // if these edges originally came from the same node front, use their
    // internal ordering
    // TODO: what if the edge is reversed?
    if ((a->getFrom() == n && b->getFrom() == n)) {
      if (OptGraph::getAdjEtgp(a, n).dir) {
        assert(OptGraph::getAdjEtgp(b, n).dir);
        return OptGraph::getAdjEtgp(a, n).order <
               OptGraph::getAdjEtgp(b, n).order;
      } else {
        return OptGraph::getAdjEtgp(a, n).order >
               OptGraph::getAdjEtgp(b, n).order;
      }
    } else if ((a->getTo() == n && b->getTo() == n)) {
      if (OptGraph::getAdjEtgp(a, n).dir) {
        assert(OptGraph::getAdjEtgp(b, n).dir);
        return OptGraph::getAdjEtgp(a, n).order >
               OptGraph::getAdjEtgp(b, n).order;
      } else {
        return OptGraph::getAdjEtgp(a, n).order <
               OptGraph::getAdjEtgp(b, n).order;
      }
    }
  }
  return fmod(angA + M_PI * 1.5, 2 * M_PI) > fmod(angB + M_PI * 1.5, 2 * M_PI);
}
}
}

#endif  // TRANSITMAP_GRAPH_OPTIM_OPTGRAPH_H_
