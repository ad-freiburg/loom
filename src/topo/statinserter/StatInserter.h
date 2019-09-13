// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_STATINSERTER_STATINSERTER_H_
#define TOPO_STATINSERTER_STATINSERTER_H_

#include <proj_api.h>
#include <algorithm>
#include <unordered_map>
#include "shared/transitgraph/TransitGraph.h"
#include "topo/config/TopoConfig.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/geo/PolyLine.h"
#include "util/graph/Graph.h"

using util::geo::Grid;
using util::geo::Box;
using util::geo::Line;
using util::geo::PolyLine;
using util::geo::Point;
using util::geo::Line;
using util::geo::DPoint;
using util::geo::DLine;
using util::geo::DBox;
using util::geo::SharedSegment;

using topo::config::TopoConfig;

using shared::transitgraph::TransitGraph;
using shared::transitgraph::TransitNode;
using shared::transitgraph::Station;
using shared::transitgraph::TransitEdge;
using shared::transitgraph::TransitEdgePair;
using shared::transitgraph::TransitNodePL;
using shared::transitgraph::TransitEdgePL;

typedef Grid<TransitNode*, Point, double> NodeGrid;
typedef Grid<TransitEdge*, Line, double> EdgeGrid;

typedef std::map<const TransitEdge*, std::set<const TransitEdge*>> OrigEdgs;

namespace topo {

struct StationOcc {
  Station station;
  std::set<const TransitEdge*> edges;
};

struct StationCand {
  // either an edge
  TransitEdge* edg;
  double pos;

  // or a node
  TransitNode* nd;

  double dist;

  size_t shouldServ;

  size_t truelyServ;
  size_t falselyServ;
};

class StatInserter {
 public:
  StatInserter(const TopoConfig* cfg, TransitGraph* g);

  void init();
  bool insertStations(const OrigEdgs& origEdgs);

 private:
  const config::TopoConfig* _cfg;
  TransitGraph* _g;

  std::vector<StationCand> candidates(const StationOcc& occ,
                                      const EdgeGrid& idx,
                                      const OrigEdgs& origEdgs);

  DBox bbox() const;
  EdgeGrid geoIndex();

  static double candScore(const StationCand& c);

  std::pair<size_t, size_t> served(const std::vector<TransitEdge*>& adj,
                                   const std::set<const TransitEdge*>& toServe,
                                   const OrigEdgs& origEdgs);

  TransitEdgePair split(TransitEdgePL& a, TransitNode* fr, TransitNode* to,
                        double p);

  void edgeRpl(TransitNode* n, const TransitEdge* oldE,
               const TransitEdge* newE);

  std::vector<std::vector<StationOcc>> _statClusters;
};

}  // namespace topo

#endif  // TOPO_STATINSERTER_STATINSERTER_H_
