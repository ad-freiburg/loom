// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_RESTR_RESTRINFERRER_H_
#define TOPO_RESTR_RESTRINFERRER_H_

#include <unordered_map>
#include "shared/transitgraph/TransitGraph.h"
#include "topo/config/TopoConfig.h"
#include "topo/restr/RestrGraph.h"
#include "util/graph/EDijkstra.h"

using topo::config::TopoConfig;
using shared::transitgraph::TransitGraph;
using shared::transitgraph::TransitNode;
using shared::transitgraph::TransitEdge;
using shared::transitgraph::TransitNodePL;
using shared::transitgraph::TransitEdgePL;
using util::graph::EDijkstra;
using topo::restr::RestrEdge;
using topo::restr::RestrNode;
using topo::restr::RestrEdgePL;
using topo::restr::RestrNodePL;
using topo::restr::RestrGraph;

namespace topo {
namespace restr {

typedef std::map<const TransitEdge*, std::set<const TransitEdge*>> OrigEdgs;

struct CostFunc : public EDijkstra::CostFunc<RestrNodePL, RestrEdgePL, double> {
  CostFunc(std::map<const RestrEdge*, double> sourcePos,
           std::map<const RestrEdge*, double> targetPos)
      : _sourcePos(sourcePos), _targetPos(targetPos) {}
  double operator()(const RestrEdge* from, const RestrNode* n,
                    const RestrEdge* to) const {
    if (n) {
      // for edges in _targetPos, add cost to arrive at edge
      return to->pl().geom.getLength();
    } else {
      // for edges in _sourcePos, add cost to leave the edge
    }
    return 0;
  };
  double inf() const { return std::numeric_limits<double>::infinity(); };

  std::map<const RestrEdge *, double> _sourcePos, _targetPos;
};

class RestrInferrer {
 public:
  RestrInferrer(const TopoConfig* cfg);

  void init(TransitGraph* g);
  void infer(TransitGraph* g, const OrigEdgs& origEdgs);

 private:
  const TopoConfig* _cfg;

  // internal copy of the original, unmodified graph
  RestrGraph _g;

  // mapping from the edges in the original graph to our internal
  // graph representation
  std::unordered_map<TransitEdge*, std::pair<RestrEdge*, RestrEdge*>> _eMap;

  // mapping from the nodes in the original graph to our internal
  // graph representation
  std::unordered_map<TransitNode*, RestrNode*> _nMap;

  // check whether a connection ocurred in the original graph
  bool check(const Route* r, const TransitEdge* edg1, const TransitEdge* edg2,
             const OrigEdgs& origEdgs) const;
};
}
}

#endif  // TOPO_RESTR_RESTRINFERRER
