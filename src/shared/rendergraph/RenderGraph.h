// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_RENDERGRAPH_RENDERGRAPH_H_
#define SHARED_RENDERGRAPH_RENDERGRAPH_H_

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <set>
#include <string>

#include "shared/linegraph/Line.h"
#include "shared/linegraph/LineGraph.h"
#include "shared/rendergraph/OrderCfg.h"
#include "shared/rendergraph/Penalties.h"
#include "util/geo/Geo.h"

namespace shared {
namespace rendergraph {

struct InnerGeom {
  InnerGeom(util::geo::PolyLine<double> g, shared::linegraph::Partner a,
            shared::linegraph::Partner b, size_t slotF, size_t slotT)
      : geom(g), from(a), to(b), slotFrom(slotF), slotTo(slotT){};
  util::geo::PolyLine<double> geom;
  shared::linegraph::Partner from, to;
  size_t slotFrom, slotTo;
};

class RenderGraph : public shared::linegraph::LineGraph {
 public:
  RenderGraph() : _defWidth(5), _defSpacing(5){};
  RenderGraph(double defLineWidth, double defLineSpace)
      : _defWidth(defLineWidth), _defSpacing(defLineSpace){};

  virtual void readFromJson(std::istream* s, double smooth);
  virtual void readFromDot(std::istream* s, double smooth);

  void writePermutation(const OrderCfg&);

  size_t numEdgs() const;

  std::vector<shared::rendergraph::InnerGeom> innerGeoms(
      const shared::linegraph::LineNode* n, double prec) const;

  std::vector<util::geo::Polygon<double>> getStopGeoms(
      const shared::linegraph::LineNode* n, double d, bool simple) const;

  // TODO: maybe move this to LineGraph?
  static size_t getConnCardinality(const shared::linegraph::LineNode* n);

  double getTotalWidth(const shared::linegraph::LineEdge* e) const;

  double getWidth(const shared::linegraph::LineEdge* e) const;
  double getSpacing(const shared::linegraph::LineEdge* e) const;

  double getMaxNdFrontWidth(const shared::linegraph::LineNode* n) const;
  size_t getMaxNdFrontCard(const shared::linegraph::LineNode* n) const;

  util::geo::DPoint linePosOn(const shared::linegraph::NodeFront& nf,
                              const shared::linegraph::Line* r,
                              bool origGeom) const;

  util::geo::DPoint linePosOn(const shared::linegraph::NodeFront& nf,
                              const shared::linegraph::LineEdge* e, size_t pos,
                              bool inv, bool origG) const;

  static bool notCompletelyServed(const shared::linegraph::LineNode* n);

  std::vector<util::geo::Polygon<double>> getIndStopPolys(
      const std::set<const shared::linegraph::Line*>& served,
      const shared::linegraph::LineNode* n, double d) const;

  void createMetaNodes();

  static bool isTerminus(const shared::linegraph::LineNode* n);

 private:
  double _defWidth, _defSpacing;

  shared::rendergraph::InnerGeom getInnerBezier(
      const shared::linegraph::LineNode* n,
      const shared::linegraph::Partner& partnerFrom,
      const shared::linegraph::Partner& partnerTo, double prec) const;

  shared::rendergraph::InnerGeom getInnerStraightLine(
      const shared::linegraph::LineNode* n,
      const shared::linegraph::Partner& partnerFrom,
      const shared::linegraph::Partner& partnerTo) const;

  shared::rendergraph::InnerGeom getTerminusStraightLine(
      const shared::linegraph::LineNode* n,
      const shared::linegraph::Partner& partnerFrom) const;

  shared::rendergraph::InnerGeom getTerminusBezier(
      const shared::linegraph::LineNode* n,
      const shared::linegraph::Partner& partnerFrom, double prec) const;

  util::geo::Polygon<double> getConvexFrontHull(
      const shared::linegraph::LineNode* n, double d, bool rectangulize,
      bool simpleRenderForTwoEdgeNodes) const;

  std::vector<shared::linegraph::NodeFront> getNextMetaNodeCand() const;
  bool isClique(std::set<const shared::linegraph::LineNode*> potClique) const;

  std::vector<shared::linegraph::NodeFront> getOpenNodeFronts(
      const shared::linegraph::LineNode* n) const;

  std::vector<shared::linegraph::NodeFront> getClosedNodeFronts(
      const shared::linegraph::LineNode* n) const;
};
}  // namespace rendergraph
}  // namespace shared

#endif  // SHARED_RENDERGRAPH_RENDERGRAPH_H_
