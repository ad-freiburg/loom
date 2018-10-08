// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_OPTIM_OPTGRAPH_H_
#define TRANSITMAP_GRAPH_OPTIM_OPTGRAPH_H_

#include <set>
#include <string>

#include "util/graph/UndirGraph.h"
#include "util/json/Writer.h"
#include "transitmap/graph/Edge.h"
#include "transitmap/graph/TransitGraph.h"

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

struct EtgPart {
  Edge* etg;
  bool dir;

  EtgPart(Edge* etg, bool dir) : etg(etg), dir(dir){};
};

struct OptEdgePL {
  OptEdgePL() : siameseSibl(0) {};
  std::vector<EtgPart> etgs;

  size_t getCardinality() const;

  // there is another edge with determines the
  // ordering in this edge - important to prevent double
  // writing of ordering later on
  OptEdge* siameseSibl;

  std::string getStrRepr() const;

  const util::geo::Line<double>* getGeom();
  util::json::Dict getAttrs();
};

struct OptNodePL {
  const Node* node;
  util::geo::Point<double> p;

  OptNodePL(util::geo::Point<double> p) : node(0), p(p) {};
  OptNodePL(const Node* node) : node(node){};
  OptNodePL() : node(0){};

  const util::geo::Point<double>* getGeom();
  util::json::Dict getAttrs();
};


class OptGraph : public UndirGraph<OptNodePL, OptEdgePL> {
 public:
  explicit OptGraph(TransitGraph* toOptim);

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

  // apply splitting rules
  void split();

 private:
  TransitGraph* _g;

  OptNode* getNodeForTransitNode(const Node* tn) const;

  void build();
  bool simplifyStep();

  bool untangleYStep();
  bool untangleDogBoneStep();


  static EtgPart getFirstEdg(const OptEdge*);
  static EtgPart getLastEdg(const OptEdge*);
};
}
}

#endif  // TRANSITMAP_GRAPH_OPTIM_OPTGRAPH_H_
