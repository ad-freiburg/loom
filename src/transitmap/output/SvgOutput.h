// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OUTPUT_SVGOUTPUT_H_
#define TRANSITMAP_OUTPUT_SVGOUTPUT_H_

#include <string>
#include <ostream>
#include "Output.h"
#include "pbutil/geo/Geo.h"
#include "pbutil/geo/PolyLine.h"
#include "./../config/TransitMapConfig.h"
#include "./../util/XmlWriter.h"
#include "./../graph/TransitGraph.h"
#include "./../graph/Route.h"

using pbutil::Nullable;
using namespace pbutil::geo;

namespace transitmapper {
namespace output {

class SvgOutputException : public std::exception {
 public:
  SvgOutputException(std::string msg)
   : _msg(msg) {}
  ~SvgOutputException() throw() {}

  virtual const char* what() const throw() {
    return _msg.c_str();
  };

 private:
  std::string _msg;
};

struct EndMarker {
  EndMarker(const std::string& name, const std::string& color, const std::string& path, double width, double height)
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
  virtual ~SvgOutput() {};

  virtual void print(const graph::TransitGraph& outputGraph);

	void printLine(const PolyLine& l,
								const std::map<std::string, std::string>& ps,
                double w, double h, int64_t xOffs, int64_t yOffs);
	void printLine(const PolyLine& l,
								const std::string& style,
                double w, double h, int64_t xOffs, int64_t yOffs);
  void printPoint(const Point& p, const std::string& style,
                          double w, double h, int64_t xOffs, int64_t yOffs);
  void printPolygon(const Polygon& g,
										const std::string& style,
                    double w, double h, int64_t xOffs, int64_t yOffs);
  void printPolygon(const Polygon& g,
								    const std::map<std::string, std::string>& ps,
                    double w, double h, int64_t xOffs, int64_t yOffs);
 private:
  std::ostream* _o;
  util::XmlWriter _w;

  const config::Config* _cfg;

  std::map<uintptr_t, std::vector<OutlinePrintPair> > _delegates;
  std::vector<EndMarker> _markers;

  void outputNodes(const graph::TransitGraph& outputGraph, double w, double h);
  void outputEdges(const graph::TransitGraph& outputGraph, double w, double h);

  void renderEdgeTripGeom(const graph::TransitGraph& outG,
    const graph::Edge* e, double w, double h);

  void renderNodeConnections(const graph::TransitGraph& outG,
    const graph::Node* n, double w, double h);

  void renderNodeScore(const graph::TransitGraph& outG,
      const graph::Node* n, double w, double h);

  void renderLinePart(const PolyLine p, double width,
    const graph::Route& route, const Nullable<style::LineStyle> style);

  void renderLinePart(const PolyLine p, double width,
    const graph::Route& route, const std::string& endMarker,
    const Nullable<style::LineStyle> style);

  void renderDelegates(const graph::TransitGraph& outG, double w, double h);

  void renderNodeFronts(const graph::TransitGraph& outG, double w, double h);

  std::string getMarkerPathMale(double w) const;
  std::string getMarkerPathFemale(double w) const;
};

}}

#endif  // TRANSITMAP_OUTPUT_SVGOUTPUT_H_
