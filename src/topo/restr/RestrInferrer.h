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
  CostFunc(const Route* r, double max) : _max(max), _route(r) {}
  double inf() const { return _max; };
  double operator()(const RestrEdge* from, const RestrNode* n,
                    const RestrEdge* to) const {
    if (!from) return 0;
    if (!from->pl().routes.count(_route)) return inf();
    if (!to->pl().routes.count(_route)) return inf();

    double c = 0;

    if (n) {
      c += 0;  // turn restrictions

      // dont allow going back the same edge
      if (from->getOtherNd(n) == to->getOtherNd(n)) return inf();
    }

    return c + to->pl().geom.getLength();
  };

  std::map<const RestrEdge *, double> _sourcePos, _targetPos;
  double _max;
  const Route* _route;
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
  std::unordered_map<TransitEdge*, std::vector<RestrEdge*>> _eMap;

  // mapping from the nodes in the original graph to our internal
  // graph representation
  std::unordered_map<TransitNode*, RestrNode*> _nMap;

  // check whether a connection ocurred in the original graph
  bool check(const Route* r, const TransitEdge* edg1, const TransitEdge* edg2,
             const OrigEdgs& origEdgs) const;

  void insertHandles(TransitGraph* g, const OrigEdgs& origEdgs);
  void insertHandles(
      const TransitEdge* e, const OrigEdgs& origEdgs,
      std::map<RestrEdge*, std::vector<std::pair<RestrNode*, double>>>*
          handles);

  std::map<const TransitEdge *, std::set<RestrNode *>> _handlesA, _handlesB;
};
}
}

#endif  // TOPO_RESTR_RESTRINFERRER
