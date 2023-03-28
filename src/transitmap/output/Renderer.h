// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OUTPUT_RENDERER_H_
#define TRANSITMAP_OUTPUT_RENDERER_H_

#include "shared/rendergraph/RenderGraph.h"

namespace transitmapper {
namespace output {

class RendererException : public std::exception {
 public:
  RendererException(std::string msg) : _msg(msg) {}
  ~RendererException() throw() {}

  virtual const char* what() const throw() { return _msg.c_str(); };

 private:
  std::string _msg;
};

struct InnerClique {
  InnerClique(const shared::linegraph::LineNode* n,
              shared::rendergraph::InnerGeom geom)
      : n(n) {
    geoms.push_back(geom);
  };
  std::vector<shared::rendergraph::InnerGeom> geoms;

  double getZWeight() const;
  size_t getNumBranchesIn(const shared::linegraph::LineEdge* front) const;
  bool operator<(const InnerClique& rhs) const;

  const shared::linegraph::LineNode* n;
};

struct RenderParams {
  double width;
  double height;
  int64_t xOff;
  int64_t yOff;
};

typedef std::map<std::string, std::string> Params;
typedef std::pair<Params, util::geo::PolyLine<double>> PrintDelegate;

struct OutlinePrintPair {
  OutlinePrintPair(PrintDelegate front, PrintDelegate back)
      : front(front), back(back) {}

  PrintDelegate front;
  PrintDelegate back;
};


class Renderer {
 public:
  virtual ~Renderer() {};

  // print the outputGraph
  virtual void print(const shared::rendergraph::RenderGraph& outG) = 0;
};

}}

#endif  // TRANSITMAP_OUTPUT_RENDERER_H_
