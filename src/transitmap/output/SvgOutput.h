// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OUTPUT_SVGOUTPUT_H_
#define TRANSITMAP_OUTPUT_SVGOUTPUT_H_

#include <ostream>
#include <string>
#include "./../config/TransitMapConfig.h"
#include "./../graph/Route.h"
#include "./../graph/TransitGraph.h"
#include "Output.h"
#include "pbutil/geo/Geo.h"
#include "pbutil/geo/PolyLine.h"
#include "pbutil/xml/XmlWriter.h"

using pbutil::Nullable;
using namespace pbutil::geo;

namespace transitmapper {
namespace output {

class SvgOutputException : public std::exception {
 public:
  SvgOutputException(std::string msg) : _msg(msg) {}
  ~SvgOutputException() throw() {}

  virtual const char* what() const throw() { return _msg.c_str(); };

 private:
  std::string _msg;
};

struct InnerClique {
  InnerClique(graph::InnerGeometry geom) {geoms.push_back(geom);};
  std::vector<graph::InnerGeometry> geoms;
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
typedef std::pair<Params, PolyLine> PrintDelegate;
typedef std::pair<PrintDelegate, PrintDelegate> OutlinePrintPair;

class SvgOutput : public Output {
 public:
  SvgOutput(std::ostream* o, const config::Config* cfg);
  virtual ~SvgOutput(){};

  virtual void print(const graph::TransitGraph& outputGraph);

  void printLine(const PolyLine& l,
                 const std::map<std::string, std::string>& ps,
                 const RenderParams& params);
  void printLine(const PolyLine& l, const std::string& style,
                 const RenderParams& params);
  void printPoint(const Point& p, const std::string& style,
                  const RenderParams& params);
  void printPolygon(const Polygon& g, const std::string& style,
                    const RenderParams& params);
  void printPolygon(const Polygon& g,
                    const std::map<std::string, std::string>& ps,
                    const RenderParams& params);

 private:
  std::ostream* _o;
  util::XmlWriter _w;

  const config::Config* _cfg;

  std::map<uintptr_t, std::vector<OutlinePrintPair> > _delegates;
  std::map<uintptr_t, std::vector<OutlinePrintPair> > _innerDelegates;
  std::vector<EndMarker> _markers;

  void outputNodes(const graph::TransitGraph& outputGraph,
                   const RenderParams& params);
  void outputEdges(const graph::TransitGraph& outputGraph,
                   const RenderParams& params);

  void renderEdgeTripGeom(const graph::TransitGraph& outG, const graph::Edge* e,
                          const RenderParams& params);

  void renderNodeConnections(const graph::TransitGraph& outG,
                             const graph::Node* n, const RenderParams& params);

  void renderNodeScore(const graph::TransitGraph& outG, const graph::Node* n,
                       const RenderParams& params);

  void renderLinePart(const PolyLine p, double width, const graph::Route& route,
                      const Nullable<style::LineStyle> style);

  void renderLinePart(const PolyLine p, double width, const graph::Route& route,
                      const std::string& endMarker,
                      const Nullable<style::LineStyle> style);

  void renderDelegates(const graph::TransitGraph& outG,
                       const RenderParams& params);

  void renderNodeFronts(const graph::TransitGraph& outG,
                        const RenderParams& params);
  std::vector<InnerClique> getInnerCliques(
      std::vector<graph::InnerGeometry> geoms, size_t level) const;

  void renderClique(const InnerClique& c, const graph::Node* node);

  bool isNextTo(const graph::InnerGeometry& a, const graph::InnerGeometry b) const;
  bool hasSameOrigin(const graph::InnerGeometry& a, const graph::InnerGeometry b) const;

  size_t getNextPartner(
      const InnerClique& forGeom,
      const std::vector<graph::InnerGeometry>& pool, size_t level) const;

  std::string getMarkerPathMale(double w) const;
  std::string getMarkerPathFemale(double w) const;
};
}
}

#endif  // TRANSITMAP_OUTPUT_SVGOUTPUT_H_
