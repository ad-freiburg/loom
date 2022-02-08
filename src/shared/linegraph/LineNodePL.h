// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_LINEGRAPH_LINENODEPL_H_
#define SHARED_LINEGRAPH_LINENODEPL_H_

#include "shared/linegraph/Line.h"
#include "shared/linegraph/LineEdgePL.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Edge.h"
#include "util/graph/Node.h"

namespace shared {
namespace linegraph {

typedef util::graph::Edge<LineNodePL, LineEdgePL> LineEdge;
typedef util::graph::Node<LineNodePL, LineEdgePL> LineNode;
typedef std::map<const Line*,
                 std::map<const LineEdge*, std::set<const LineEdge*>>>
    ConnEx;

typedef std::set<const Line*> NotServedLines;

struct NodeFront {
  NodeFront(LineNode* n, LineEdge* e) : n(n), edge(e) {}

  LineNode* n;
  LineEdge* edge;

  // geometry after expansion
  PolyLine<double> geom;

  // geometry before expansion
  PolyLine<double> origGeom;

  void setInitialGeom(const PolyLine<double>& g) {
    geom = g;
    origGeom = g;
  };
  void setGeom(const PolyLine<double>& g) { geom = g; };

  // TODO
  double refEtgLengthBefExp;
};

struct Station {
  Station(const std::string& id, const std::string& name,
          const util::geo::DPoint& pos)
      : id(id), name(name), pos(pos) {}
  std::string id, name;
  util::geo::DPoint pos;
};

struct ConnException {
  ConnException(const LineEdge* from, const LineEdge* to) : fr(from), to(to) {}
  const LineEdge* fr;
  const LineEdge* to;
};

class LineNodePL : util::geograph::GeoNodePL<double> {
 public:
  LineNodePL(){};
  LineNodePL(util::geo::Point<double> pos);

  const util::geo::Point<double>* getGeom() const;
  void setGeom(const util::geo::Point<double>& p);
  util::json::Dict getAttrs() const;

  void addStop(const Station& i);
  const std::vector<Station>& stops() const;
  void clearStops();

  // TODO refactor, all front related stuff should go into rendergraph
  const std::vector<NodeFront>& fronts() const;
  std::vector<NodeFront>& fronts();
  void delFrontFor(const LineEdge* e);
  const NodeFront* frontFor(const LineEdge* e) const;

  void addFront(const NodeFront& f);

  void addConnExc(const Line* r, const LineEdge* edgeA, const LineEdge* edgeB);

  void delConnExc(const Line* r, const LineEdge* edgeA, const LineEdge* edgeB);

  bool connOccurs(const Line* r, const LineEdge* edgeA,
                  const LineEdge* edgeB) const;

  void addLineNotServed(const Line* r);
  bool lineServed(const Line* r) const;

  void clearConnExc();

  size_t numConnExcs() const;

  ConnEx& getConnExc() { return _connEx; }
  const ConnEx& getConnExc() const { return _connEx; }

  std::string toString() const;

 private:
  util::geo::Point<double> _pos;
  std::vector<Station> _is;

  std::map<const LineEdge*, size_t> _edgToNf;
  std::vector<NodeFront> _nodeFronts;

  ConnEx _connEx;
  NotServedLines _notServed;
};
}  // namespace linegraph
}  // namespace shared

#endif  // SHARED_LINEGRAPH_LINENODEPL_H_
