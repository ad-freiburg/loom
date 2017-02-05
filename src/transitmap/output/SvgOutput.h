// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OUTPUT_SVGOUTPUT_H_
#define TRANSITMAP_OUTPUT_SVGOUTPUT_H_

#include <string>
#include <ostream>
#include "Output.h"
#include "./../config/TransitMapConfig.h"
#include "./../util/XmlWriter.h"
#include "./../util/Geo.h"
#include "./../graph/TransitGraph.h"
#include "./../graph/Route.h"
#include "./../geo/PolyLine.h"

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

typedef std::map<std::string, std::string> Params;
typedef std::pair<Params, geo::PolyLine> PrintDelegate;
typedef std::pair<PrintDelegate, PrintDelegate> OutlinePrintPair;

class SvgOutput : public Output {

 public:
  SvgOutput(std::ostream* o, const config::Config* cfg);
  virtual ~SvgOutput() {};

  virtual void print(const graph::TransitGraph& outputGraph);

	void printLine(const transitmapper::geo::PolyLine& l,
								const std::map<std::string, std::string>& ps,
                double w, double h, int64_t xOffs, int64_t yOffs);
	void printLine(const transitmapper::geo::PolyLine& l,
								const std::string& style,
                double w, double h, int64_t xOffs, int64_t yOffs);
  void printPoint(const util::geo::Point& p, const std::string& style,
                          double w, double h, int64_t xOffs, int64_t yOffs);
  void printPolygon(const util::geo::Polygon& g,
										const std::string& style,
                    double w, double h, int64_t xOffs, int64_t yOffs);
  void printPolygon(const util::geo::Polygon& g,
								    const std::map<std::string, std::string>& ps,
                    double w, double h, int64_t xOffs, int64_t yOffs);
 private:
  std::ostream* _o;
  util::XmlWriter _w;

  const config::Config* _cfg;

  std::map<uintptr_t, std::vector<OutlinePrintPair> > _delegates;

  void outputNodes(const graph::TransitGraph& outputGraph, double w, double h);
  void outputEdges(const graph::TransitGraph& outputGraph, double w, double h);

  void renderEdgeTripGeom(const graph::TransitGraph& outG,
    const graph::Edge* e, double w, double h);

  void renderNodeConnections(const graph::TransitGraph& outG,
    const graph::Node* n, double w, double h);

  void renderNodeScore(const graph::TransitGraph& outG,
      const graph::Node* n, double w, double h);

  void renderLinePart(const geo::PolyLine p, double width,
    const graph::Route& route);

  void renderDelegates(const graph::TransitGraph& outG, double w, double h);

  void renderNodeFronts(const graph::TransitGraph& outG, double w, double h);
};

}}

#endif  // TRANSITMAP_OUTPUT_SVGOUTPUT_H_
