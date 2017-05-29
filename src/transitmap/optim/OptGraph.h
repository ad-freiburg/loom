// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_OPTIM_OPTGRAPH_H_
#define TRANSITMAP_GRAPH_OPTIM_OPTGRAPH_H_

#include <set>
#include <string>

#include "./../graph/Edge.h"
#include "./../graph/TransitGraph.h"

using transitmapper::graph::TransitGraph;
using transitmapper::graph::Node;
using transitmapper::graph::Edge;

namespace transitmapper {
namespace optim {

struct OptNode;

struct EtgPart {
  Edge* etg;
  bool dir;

  EtgPart(Edge* etg, bool dir) : etg(etg), dir(dir){};
};

struct OptEdge {
  std::vector<EtgPart> etgs;

  OptNode* from;
  OptNode* to;
  OptEdge(OptNode* from, OptNode* to) : from(from), to(to){};

  EtgPart getFirstEdge() const;
  EtgPart getLastEdge() const;

  std::string getStrRepr() const;

  graph::Edge* getAdjacentEdge(const OptNode* n) const;
};

struct OptNode {
  const Node* node;

  std::set<OptEdge*> adjListIn;
  std::set<OptEdge*> adjListOut;
  std::set<OptEdge*> adjList;

  OptNode(const Node* node) : node(node){};

  void addEdge(OptEdge* e) {
    adjList.insert(e);
    if (e->from == this) adjListOut.insert(e);
    if (e->to == this) adjListIn.insert(e);
  }

  void deleteEdge(OptEdge* e) {
    adjList.erase(e);
    adjListOut.erase(e);
    adjListIn.erase(e);
  }
};

class OptGraph {
 public:
  explicit OptGraph(TransitGraph* toOptim);

  const std::set<OptNode*>& getNodes() const;
  TransitGraph* getGraph() const;

  size_t getNumNodes() const;
  size_t getNumNodes(bool topo) const;
  size_t getNumEdges() const;
  size_t getNumRoutes() const;
  size_t getMaxCardinality() const;

  void simplify();

 private:
  TransitGraph* _g;
  std::set<OptNode*> _nodes;

  void addNode(OptNode* n);

  OptNode* getNodeForTransitNode(const Node* tn) const;

  void build();

  bool simplifyStep();
};
}
}

#endif  // TRANSITMAP_GRAPH_TRANSITGRAPH_H_
