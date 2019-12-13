// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_LINEGRAPH_LINENODEPL_H_
#define SHARED_LINEGRAPH_LINENODEPL_H_

#include "shared/linegraph/Route.h"
#include "shared/linegraph/LineEdgePL.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"

namespace shared {
namespace linegraph {

using util::geo::Point;
using util::graph::Node;

class LineNodePL;

typedef util::graph::Edge<LineNodePL, LineEdgePL> LineEdge;
typedef std::map<const Route*,
                 std::map<const LineEdge*, std::set<const LineEdge*>>>
    ConnEx;

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

  ConnEx _connEx;
};
}
}

#endif  // SHARED_LINEGRAPH_LINENODEPL_H_
