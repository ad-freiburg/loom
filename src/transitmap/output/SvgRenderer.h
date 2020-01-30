// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OUTPUT_SVGRENDERER_H_
#define TRANSITMAP_OUTPUT_SVGRENDERER_H_

#include <ostream>
#include <set>
#include <string>
#include <vector>
#include "Renderer.h"
#include "shared/linegraph/Line.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/RenderGraph.h"
#include "transitmap/label/Labeller.h"
#include "transitmap/optim/Scorer.h"
#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"
#include "util/xml/XmlWriter.h"

using util::Nullable;

namespace transitmapper {
namespace output {

class SvgRendererException : public std::exception {
 public:
  SvgRendererException(std::string msg) : _msg(msg) {}
  ~SvgRendererException() throw() {}

  virtual const char* what() const throw() { return _msg.c_str(); };

 private:
  std::string _msg;
};

struct InnerClique {
  InnerClique(const shared::linegraph::LineNode* n,
              shared::linegraph::InnerGeom geom)
      : n(n) {
    geoms.push_back(geom);
  };
  std::vector<shared::linegraph::InnerGeom> geoms;

  double getZWeight() const;
  size_t getNumBranchesIn(const shared::linegraph::NodeFront* front) const;
  bool operator<(const InnerClique& rhs) const;

  const shared::linegraph::LineNode* n;
};

struct RenderParams {
  double width;
  double height;
  int64_t xOff;
  int64_t yOff;
};

struct EndMarker {
  EndMarker(const std::string& name, const std::string& color,
            const std::string& path, double width, double height)
      : name(name), color(color), path(path), width(width), height(height) {}
  std::string name;
  std::string color;
  std::string path;
  double width, height;
};

typedef std::map<std::string, std::string> Params;
typedef std::pair<Params, util::geo::PolyLine<double>> PrintDelegate;

struct OutlinePrintPair {
  OutlinePrintPair(PrintDelegate front, PrintDelegate back)
      : front(front), back(back) {}

  PrintDelegate front;
  PrintDelegate back;
};

class SvgRenderer : public Renderer {
 public:
  SvgRenderer(std::ostream* o, const config::Config* cfg,
              const optim::Scorer* scorer);
  virtual ~SvgRenderer(){};

  virtual void print(const graph::RenderGraph& outputGraph);

  void printLine(const util::geo::PolyLine<double>& l,
                 const std::map<std::string, std::string>& ps,
                 const RenderParams& params);
  void printLine(const util::geo::PolyLine<double>& l, const std::string& style,
                 const RenderParams& params);
  void printPoint(const util::geo::DPoint& p, const std::string& style,
                  const RenderParams& params);
  void printPolygon(const util::geo::Polygon<double>& g,
                    const std::string& style, const RenderParams& params);
  void printPolygon(const util::geo::Polygon<double>& g,
                    const std::map<std::string, std::string>& ps,
                    const RenderParams& params);
  void printCircle(const util::geo::DPoint& center, double rad,
                   const std::map<std::string, std::string>& ps,
                   const RenderParams& rparams);
  void printCircle(const util::geo::DPoint& center, double rad,
                   const std::string& style, const RenderParams& rparams);

 private:
  std::ostream* _o;
  util::xml::XmlWriter _w;

  const config::Config* _cfg;
  const optim::Scorer* _scorer;

  std::map<uintptr_t, std::vector<OutlinePrintPair>> _delegates;
  std::vector<std::map<uintptr_t, std::vector<OutlinePrintPair>>>
      _innerDelegates;
  std::vector<EndMarker> _markers;

  void outputNodes(const graph::RenderGraph& outputGraph,
                   const RenderParams& params);
  void outputEdges(const graph::RenderGraph& outputGraph,
                   const RenderParams& params);

  void renderEdgeTripGeom(const graph::RenderGraph& outG,
                          const shared::linegraph::LineEdge* e,
                          const RenderParams& params);

  void renderNodeConnections(const graph::RenderGraph& outG,
                             const shared::linegraph::LineNode* n,
                             const RenderParams& params);

  void renderLinePart(const util::geo::PolyLine<double> p, double width,
                      const shared::linegraph::Line& line,
                      const shared::linegraph::LineEdge* e);

  void renderLinePart(const util::geo::PolyLine<double> p, double width,
                      const shared::linegraph::Line& line,
                      const shared::linegraph::LineEdge* edge,
                      const std::string& endMarker);

  void renderDelegates(const graph::RenderGraph& outG,
                       const RenderParams& params);

  void renderNodeFronts(const graph::RenderGraph& outG,
                        const RenderParams& params);

  void renderLineLabels(const label::Labeller& lbler,
                        const RenderParams& params);

  void renderStationLabels(const label::Labeller& lbler,
                        const RenderParams& params);

  std::multiset<InnerClique> getInnerCliques(
      const shared::linegraph::LineNode* n,
      std::vector<shared::linegraph::InnerGeom> geoms, size_t level) const;

  void renderClique(const InnerClique& c,
                    const shared::linegraph::LineNode* node);

  bool isNextTo(const shared::linegraph::InnerGeom& a,
                const shared::linegraph::InnerGeom b) const;
  bool hasSameOrigin(const shared::linegraph::InnerGeom& a,
                     const shared::linegraph::InnerGeom b) const;

  size_t getNextPartner(const InnerClique& forGeom,
                        const std::vector<shared::linegraph::InnerGeom>& pool,
                        size_t level) const;

  std::string getMarkerPathMale(double w) const;
  std::string getMarkerPathFemale(double w) const;
};
}
}

#endif  // TRANSITMAP_OUTPUT_SVGRENDERER_H_
