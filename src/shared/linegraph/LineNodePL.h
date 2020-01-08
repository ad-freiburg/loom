// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_LINEGRAPH_LINENODEPL_H_
#define SHARED_LINEGRAPH_LINENODEPL_H_

#include "shared/linegraph/Route.h"
#include "shared/linegraph/LineEdgePL.h"
#include "transitmap/graph/OrderingConfig.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"

namespace shared {
namespace linegraph {

using util::geo::Point;
using util::geo::DPoint;
using util::graph::Node;

// TODO: move to shared
using transitmapper::graph::OrderingConfig;

class LineNodePL;

typedef util::graph::Edge<LineNodePL, LineEdgePL> LineEdge;
typedef std::map<const Route*,
                 std::map<const LineEdge*, std::set<const LineEdge*>>>
    ConnEx;

struct NodeFront {
  NodeFront(LineEdge* e) : edge(e) {}

  DPoint getTripOccPos(const shared::linegraph::Route* r,
                       const OrderingConfig& c, bool origGeom) const;
  DPoint getTripPos(const LineEdge* e, size_t pos, bool inv, bool originGeom) const;

  double getOutAngle() const;

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

struct Partner {
  Partner() : front(0), edge(0), route(0){};
  Partner(const NodeFront* f, const LineEdge* e, const shared::linegraph::Route* r)
      : front(f), edge(e), route(r){};
  const NodeFront* front;
  const LineEdge* edge;
  const shared::linegraph::Route* route;
};

struct InnerGeometry {
  InnerGeometry(PolyLine<double> g, Partner a, Partner b, size_t slotF,
                size_t slotT)
      : geom(g), from(a), to(b), slotFrom(slotF), slotTo(slotT){};
  PolyLine<double> geom;
  Partner from, to;
  size_t slotFrom, slotTo;
};

struct Station {
  Station(const std::string& id, const std::string& name,
          const util::geo::DPoint& pos)
      : id(id), name(name), pos(pos) {}
  std::string id, name;
  util::geo::DPoint pos;
};

struct ConnException {
  ConnException(const LineEdge* from, const LineEdge* to)
      : fr(from), to(to) {}
  const LineEdge* fr;
  const LineEdge* to;
};

class LineNodePL : util::geograph::GeoNodePL<double> {
 public:
  LineNodePL(){};
  LineNodePL(Point<double> pos);

  const Point<double>* getGeom() const;
  void setGeom(const Point<double>& p);
  util::json::Dict getAttrs() const;

  void addStop(const Station& i);
  const std::vector<Station>& getStops() const;
  void clearStops();

  size_t getLineDeg() const;

  // TODO refactor
  const std::vector<NodeFront>& getMainDirs() const { return _mainDirs; }
  std::vector<NodeFront>& getMainDirs() { return _mainDirs; }

  void addMainDir(NodeFront f);

  void addConnExc(const Route* r, const LineEdge* edgeA,
                  const LineEdge* edgeB);

  void delConnExc(const Route* r, const LineEdge* edgeA,
                  const LineEdge* edgeB);

  bool connOccurs(const Route* r, const LineEdge* edgeA,
                  const LineEdge* edgeB) const;

  ConnEx& getConnExc() { return _connEx; }
  const ConnEx& getConnExc() const { return _connEx; }

 private:
  Point<double> _pos;
  std::vector<Station> _is;

  std::vector<NodeFront> _mainDirs;

  ConnEx _connEx;
};
}
}

#endif  // SHARED_LINEGRAPH_LINENODEPL_H_
