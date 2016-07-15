// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OUTPUT_SVGOUTPUT_H_
#define TRANSITMAP_OUTPUT_SVGOUTPUT_H_

#include <string>
#include <ostream>
#include "output.h"
#include "../util/XmlWriter.h"
#include "../util/Geo.h"
#include "../graph/transitgraph.h"
#include "../geo/PolyLine.h"

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

class SvgOutput : public Output {

 public:
  SvgOutput(std::ostream* o);
  virtual ~SvgOutput() {};

  virtual void print(const graph::TransitGraph& outputGraph);

	void printLine(const transitmapper::geo::PolyLine& l, const std::string& style);
 private:
  std::ostream* _o;
  util::XmlWriter _w;

  void outputNodes(const graph::TransitGraph& outputGraph);
  void outputEdges(const graph::TransitGraph& outputGraph);
};

}}

#endif  // TRANSITMAP_OUTPUT_SVGOUTPUT_H_