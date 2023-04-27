// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_LINEGRAPH_LINEEDGEPL_H_
#define SHARED_LINEGRAPH_LINEEDGEPL_H_

#include <set>

#include "shared/linegraph/Line.h"
#include "shared/style/LineStyle.h"
#include "util/Nullable.h"
#include "util/geo/GeoGraph.h"
#include "util/geo/PolyLine.h"
#include "util/graph/Node.h"

namespace shared {
namespace linegraph {

using util::Nullable;
using util::geo::PolyLine;
using util::graph::Node;

class LineEdgePL;
class LineNodePL;

struct LineOcc {
  LineOcc() {}
  LineOcc(const Line* r, const Node<LineNodePL, LineEdgePL>* dir)
      : line(r), direction(dir) {}
  LineOcc(const Line* r, const Node<LineNodePL, LineEdgePL>* dir,
          const util::Nullable<shared::style::LineStyle>& ls)
      : line(r), direction(dir), style(ls) {}
  const Line* line;
  const Node<LineNodePL, LineEdgePL>* direction;  // 0 if in both directions

  util::Nullable<shared::style::LineStyle> style;
};

inline bool operator<(const LineOcc& x, const LineOcc& y) {
  return x.line < y.line;
}

class LineEdgePL : util::geograph::GeoEdgePL<double> {
 public:
  LineEdgePL();
  LineEdgePL(const PolyLine<double>& p);

  void addLine(const Line* r, const Node<LineNodePL, LineEdgePL>* dir,
               util::Nullable<shared::style::LineStyle> ls);
  void addLine(const Line* r, const Node<LineNodePL, LineEdgePL>* dir);

  const std::vector<LineOcc>& getLines() const;

  bool hasLine(const Line* r) const;
  void delLine(const Line* r);

  const LineOcc& lineOcc(const Line* r) const;
  const LineOcc& lineOccAtPos(size_t i) const;

  void updateLineOcc(const LineOcc& occ);

  size_t linePos(const Line* r) const;

  const util::geo::Line<double>* getGeom() const;
  void setGeom(const util::geo::Line<double>& l);
  util::json::Dict getAttrs() const;

  const PolyLine<double>& getPolyline() const;
  void setPolyline(const PolyLine<double>& p);

  uint32_t getComponent() const { return _comp; }
  void setComponent(uint32_t id) { _comp = id; }

  void writePermutation(const std::vector<size_t> order);

  void setDontContract(bool dontContract) { _dontContract = dontContract; }
  bool dontContract() { return _dontContract; }

 private:
  std::map<const Line*, size_t> _lineToIdx;
  std::vector<LineOcc> _lines;
  bool _dontContract;
  uint32_t _comp = std::numeric_limits<uint32_t>::max();

  PolyLine<double> _p;
};
}  // namespace linegraph
}  // namespace shared

#endif  // SHARED_LINEGRAPH_LINEEDGEPL_H_
