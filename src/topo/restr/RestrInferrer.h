// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_RESTR_RESTRINFERRER_H_
#define TOPO_RESTR_RESTRINFERRER_H_

#include <unordered_map>
#include "shared/linegraph/LineGraph.h"
#include "shared/linegraph/Line.h"
#include "topo/config/TopoConfig.h"
#include "topo/restr/RestrGraph.h"
#include "util/graph/EDijkstra.h"

using topo::config::TopoConfig;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using shared::linegraph::LineEdge;
using shared::linegraph::LineNodePL;
using shared::linegraph::LineEdgePL;
using util::graph::EDijkstra;
using topo::restr::RestrEdge;
using topo::restr::RestrNode;
using topo::restr::RestrEdgePL;
using topo::restr::RestrNodePL;
using topo::restr::RestrGraph;

namespace topo {
namespace restr {

typedef std::map<const LineEdge*, std::set<const LineEdge*>> OrigEdgs;
typedef std::pair<RestrNode*, double> Hndl;
typedef std::vector<Hndl> HndlLst;

using shared::linegraph::Line;

struct CostFunc : public EDijkstra::CostFunc<RestrNodePL, RestrEdgePL, double> {
  CostFunc(const Line* r, double max) : _max(max), _line(r) {}
  double inf() const { return _max; };
  double operator()(const RestrEdge* from, const RestrNode* n,
                    const RestrEdge* to) const {
    // don't count start edge
    if (!from) return 0;

    // if an edge does not contain the line we are routing for, set
    // cost to inf
    if (!from->pl().lines.count(_line)) return inf();
    if (!to->pl().lines.count(_line)) return inf();

    double c = 0;

    if (n) {
      auto lRestrs = n->pl().restrs.find(_line);
      if (lRestrs != n->pl().restrs.end() && lRestrs->second.count(from)) {
        if (lRestrs->second.find(from)->second.count(to))  return inf();
      }

      // dont allow going back the same edge
      if (from->getOtherNd(n) == to->getOtherNd(n)) return inf();
    }

    // final cost are turn restriction costs plus length of the (to) edge
    return c + to->pl().geom.getLength();
  };

  double _max;
  const Line* _line;
};

struct HndlCmp {
  bool operator()(const Hndl& a, const Hndl& b) const {
    return a.second < b.second;
  }
};

class RestrInferrer {
 public:
  RestrInferrer(const TopoConfig* cfg, LineGraph* g);

  void init();
  size_t infer(const OrigEdgs& origEdgs);

 private:
  const TopoConfig* _cfg;

  // internal copy of the original, unmodified graph we initialized with
  RestrGraph _rg;

  // the transit graph we operate on, may be modified externally
  LineGraph* _tg;

  // mapping from the edges in the original graph to our internal
  // graph representation
  std::unordered_map<const LineEdge*, std::vector<RestrEdge*>> _eMap;

  // mapping from the nodes in the original graph to our internal
  // graph representation
  std::unordered_map<const LineNode*, RestrNode*> _nMap;

  // check whether a connection ocurred in the original graph
  bool check(const Line* r, const LineEdge* edg1, const LineEdge* edg2) const;

  void addHndls(const OrigEdgs& origEdgs);
  void addHndls(const LineEdge* e, const OrigEdgs& origEdgs,
                std::map<RestrEdge*, HndlLst>* handles);

  void edgeRpl(RestrNode* n, const RestrEdge* oldE,
                             const RestrEdge* newE);

  std::map<const LineEdge *, std::set<RestrNode *>> _handlesA, _handlesB;
};
}
}

#endif  // TOPO_RESTR_RESTRINFERRER
