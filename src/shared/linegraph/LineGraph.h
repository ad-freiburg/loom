// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_LINEGRAPH_LINEGRAPH_H_
#define SHARED_LINEGRAPH_LINEGRAPH_H_

#include "3rdparty/json.hpp"
#include "shared/linegraph/EdgeOrdering.h"
#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineNodePL.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/UndirGraph.h"

namespace shared {
namespace linegraph {

typedef util::graph::Node<LineNodePL, LineEdgePL> LineNode;
typedef util::graph::Edge<LineNodePL, LineEdgePL> LineEdge;

typedef std::pair<LineEdge*, LineEdge*> LineEdgePair;

typedef util::geo::Grid<LineNode*, util::geo::Point, double> NodeGrid;
typedef util::geo::Grid<LineEdge*, util::geo::Line, double> EdgeGrid;

struct ISect {
  LineEdge *a, *b;
  util::geo::LinePoint<double> bp;
};

struct Partner {
  Partner() : edge(0), line(0){};
  Partner(const LineEdge* e, const Line* r) : edge(e), line(r){};
  const LineEdge* edge;
  const Line* line;
};

class LineGraph : public util::graph::UndirGraph<LineNodePL, LineEdgePL> {
 public:
  LineGraph() = default;
  LineGraph(const LineGraph& other) = delete;
  void operator=(const LineGraph& other) = delete;

  LineGraph(LineGraph&& other) {
    _bbox = other._bbox;
    proced = other.proced;
    _lines = other._lines;
    _nodeGrid = other._nodeGrid;
    _edgeGrid = other._edgeGrid;

    _nodes = other._nodes;
    other._nodes.clear();
  }

  LineGraph& operator=(LineGraph&& other) {
    _bbox = other._bbox;
    proced = other.proced;
    _lines = other._lines;
    _nodeGrid = other._nodeGrid;
    _edgeGrid = other._edgeGrid;

    _nodes = other._nodes;
    other._nodes.clear();
    return *this;
  }

  virtual void readFromJson(std::istream* s, double smooth);
  virtual void readFromGeoJson(nlohmann::json::array_t, double smooth);
  virtual void readFromTopoJson(nlohmann::json::array_t objects,
                                nlohmann::json::array_t arc, double smooth);
  virtual void readFromDot(std::istream* s, double smooth);

  const util::geo::Box<double>& getBBox() const;
  void topologizeIsects();

  size_t maxDeg() const;

  // TODO: make the following functions private
  void addLine(const Line* r);
  const Line* getLine(const std::string& id) const;
  void expandBBox(const util::geo::Point<double>& p);

  size_t numNds() const;
  size_t numNds(bool topo) const;
  size_t numEdgs() const;
  size_t numLines() const;
  size_t numConnExcs() const;

  static std::vector<LineOcc> getCtdLinesIn(const LineOcc& line,
                                            const LineEdge* fromEdge,
                                            const LineEdge* toEdge);

  static std::vector<LineOcc> getCtdLinesIn(const LineEdge* fromEdge,
                                            const LineEdge* toEdge);

  static bool lineCtd(const LineEdge* fromEdge, const LineOcc& fromLine,
                      const LineEdge* toEdge, const LineOcc& toLine);

  static bool lineCtd(const LineEdge* fromEdge, const LineEdge* toEdge,
                      const Line* line);

  static bool terminatesAt(const LineEdge* fromEdge, const LineNode* terminus,
                           const Line* line);

  static std::vector<const Line*> getSharedLines(const LineEdge* a,
                                                 const LineEdge* b);

  static size_t getLDeg(const LineNode* nd);
  static size_t getMaxLineNum(const LineNode* nd);
  size_t getMaxLineNum() const;

  static std::set<const shared::linegraph::Line*> servedLines(
      const shared::linegraph::LineNode* n);

  static EdgeOrdering edgeOrdering(LineNode* n, bool useOrigNextNode);

  static void edgeDel(LineNode* n, const LineEdge* oldE);
  static void edgeRpl(LineNode* n, const LineEdge* oldE, const LineEdge* newE);
  static void nodeRpl(LineEdge* e, const LineNode* oldN, const LineNode* newN);

  std::set<LineEdge*> getNeighborEdges(const util::geo::DLine& line,
                                       double d) const;

  static std::vector<Partner> getPartners(const LineNode* nd, const LineEdge* e,
                                          const LineOcc& lo);

  NodeGrid* getNdGrid();
  const NodeGrid& getNdGrid() const;

  EdgeGrid* getEdgGrid();
  const EdgeGrid& getEdgGrid() const;

  void splitNode(LineNode* n, size_t maxDeg);
  void splitNodes(size_t maxDeg);

  void contractStrayNds();

  void contractEdges(double d);
  void contractEdges(double d, bool onlyNonStatConns);
  void contractEdge(LineEdge* e);

  double searchSpaceSize() const;

  virtual LineNode* mergeNds(LineNode* a, LineNode* b);

 private:
  util::geo::Box<double> _bbox;

  ISect getNextIntersection();

  void buildGrids();

  // TODO: remove this
  std::set<LineEdge*> proced;
  std::map<std::string, const Line*> _lines;

  NodeGrid _nodeGrid;
  EdgeGrid _edgeGrid;
};

}  // namespace linegraph
}  // namespace shared

#endif  // SHARED_LINEGRAPH_LINEGRAPH_H_
