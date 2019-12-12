// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_TRANSITGRAPH_TRANSITNODEPL_H_
#define SHARED_TRANSITGRAPH_TRANSITNODEPL_H_

#include "shared/transitgraph/TransitEdgePL.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"

using util::geo::Point;
using util::graph::Node;

namespace shared {
namespace transitgraph {

typedef util::graph::Node<TransitNodePL, TransitEdgePL> TransitNode;
typedef util::graph::Edge<TransitNodePL, TransitEdgePL> TransitEdge;
typedef std::map<const Route*,
                 std::map<const TransitEdge*, std::set<const TransitEdge*>>>
    ConnEx;

struct Station {
  Station(const std::string& id, const std::string& name,
          const util::geo::DPoint& pos)
      : id(id), name(name), pos(pos) {}
  std::string id, name;
  util::geo::DPoint pos;
};

struct ConnException {
  ConnException(const TransitEdge* from, const TransitEdge* to)
      : fr(from), to(to) {}
  const TransitEdge* fr;
  const TransitEdge* to;
};

class TransitNodePL : util::geograph::GeoNodePL<double> {
 public:
  TransitNodePL(){};
  TransitNodePL(Point<double> pos);

  const Point<double>* getGeom() const;
  void setGeom(const Point<double>& p);
  util::json::Dict getAttrs() const;

  void addStop(const Station& i);
  const std::vector<Station>& getStops() const;
  void clearStops();

  void addConnExc(const Route* r, const TransitEdge* edgeA,
                  const TransitEdge* edgeB);

  void delConnExc(const Route* r, const TransitEdge* edgeA,
                  const TransitEdge* edgeB);

  bool connOccurs(const Route* r, const TransitEdge* edgeA,
                  const TransitEdge* edgeB) const;

  ConnEx& getConnExc() { return _connEx; }
  const ConnEx& getConnExc() const { return _connEx; }

 private:
  Point<double> _pos;
  std::vector<Station> _is;

  ConnEx _connEx;
};
}
}

#endif  // SHARED_TRANSITGRAPH_TRANSITNODEPL_H_
