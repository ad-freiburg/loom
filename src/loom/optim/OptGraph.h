// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_GRAPH_OPTIM_OPTGRAPH_H_
#define LOOM_GRAPH_OPTIM_OPTGRAPH_H_

#include <set>
#include <string>

#include "shared/linegraph/LineGraph.h"
#include "shared/rendergraph/RenderGraph.h"
#include "util/Misc.h"
#include "util/graph/UndirGraph.h"
#include "util/json/Writer.h"

namespace loom {
namespace optim {

struct OptNodePL;
struct OptEdgePL;

class OptGraphScorer;

typedef util::graph::Node<OptNodePL, OptEdgePL> OptNode;
typedef util::graph::Edge<OptNodePL, OptEdgePL> OptEdge;

typedef std::map<const loom::optim::OptEdge*,
                 std::vector<const shared::linegraph::Line*>>
    OptOrderCfg;

struct OptLO {
  OptLO() : line(0), dir(0) {}
  OptLO(const shared::linegraph::Line* r,
        const shared::linegraph::LineNode* dir)
      : line(r), dir(dir) {
    relatives.push_back(r);
  }
  const shared::linegraph::Line* line;
  const shared::linegraph::LineNode* dir;  // 0 if in both directions

  std::vector<const shared::linegraph::Line*> relatives;

  bool operator==(const shared::linegraph::Line* b) const { return b == line; }
  bool operator<(const shared::linegraph::Line* b) const { return b < line; }
  bool operator>(const shared::linegraph::Line* b) const { return b > line; }

  bool operator==(const OptLO& b) const { return b.line == line; }
  bool operator<(const OptLO& b) const { return b.line < line; }
  bool operator>(const OptLO& b) const { return b.line > line; }

  bool operator==(const shared::linegraph::LineOcc& b) const {
    return b.line == line;
  }
};

struct PartnerPath {
  // Important: OptLOs with the same route are r equivalent to each other and
  // to the original route, see above
  std::set<OptLO> partners;
  std::vector<OptEdge*> path;
  std::vector<bool> inv;
};

struct LnEdgPart {
  shared::linegraph::LineEdge* lnEdg;
  bool dir;
  size_t order;

  // there is another edge which determines the
  // ordering in this edge - important to prevent double
  // writing of ordering later on
  bool wasCut;

  LnEdgPart(shared::linegraph::LineEdge* lnEdg, bool dir)
      : lnEdg(lnEdg), dir(dir), order(0), wasCut(false){};
  LnEdgPart(shared::linegraph::LineEdge* lnEdg, bool dir, size_t order,
            bool wasCut)
      : lnEdg(lnEdg), dir(dir), order(order), wasCut(wasCut){};
};

struct OptEdgePL {
  OptEdgePL() : depth(0), firstLnEdg(0), lastLnEdg(0){};

  // all original line edges from the transit graph contained in this edge
  // Guarantee: they are all equal in terms of (directed) routes
  std::vector<LnEdgPart> lnEdgParts;

  size_t depth;

  size_t firstLnEdg;
  size_t lastLnEdg;

  size_t getCardinality() const;
  std::string toStr() const;
  std::vector<OptLO>& getLines();
  const std::vector<OptLO>& getLines() const;

  const OptLO* getLineOcc(const shared::linegraph::Line* l) const;

  // partial routes
  // For the line edge parts contained in lnEdgParts, only these route
  // occurances are actually contained in this edge. Their relative ordering is
  // defined by .order
  std::vector<OptLO> lines;

  std::string getStrRepr() const;

  const util::geo::Line<double>* getGeom();
  util::json::Dict getAttrs();
};

struct OptNodePL {
  OptNodePL(util::geo::Point<double> p) : node(0), p(p){};
  OptNodePL(const shared::linegraph::LineNode* node)
      : node(node), p(*node->pl().getGeom()){};
  OptNodePL() : node(0){};

  size_t circOrder(OptEdge*) const;

  const util::geo::Point<double>* getGeom();
  util::json::Dict getAttrs();

  const shared::linegraph::LineNode* node;
  util::geo::Point<double> p;

  // the edges adjacent to this node, in clockwise fashion, based
  // on the geometry in the original graph
  std::vector<OptEdge*> circOrdering;
  std::map<OptEdge*, size_t> circOrderMap;
};

class OptGraph : public util::graph::UndirGraph<OptNodePL, OptEdgePL> {
 public:
  OptGraph(const OptGraphScorer* scorer) : _scorer(scorer){};

  std::map<const shared::linegraph::LineNode*, OptNode*> build(
      shared::rendergraph::RenderGraph* rg);

  size_t getNumNodes() const;
  size_t getNumNodes(bool topo) const;
  size_t getNumEdges() const;
  size_t getNumLines() const;
  size_t getMaxCardinality() const;

  double getMaxCrossPen() const;
  double getMaxSplitPen() const;

  void contractDeg2Nds();
  void untangle();
  void partnerLines();

  std::vector<PartnerPath> getPartnerLines() const;
  PartnerPath pathFromComp(const std::set<OptNode*>& comp) const;

  static shared::linegraph::LineEdge* getAdjEdg(const OptEdge* e,
                                                const OptNode* n);
  static LnEdgPart getAdjLnEdgPart(const OptEdge* e, const OptNode* n);

  static bool hasCtdLineIn(const shared::linegraph::Line* r,
                           const shared::linegraph::LineNode* dir,
                           const OptEdge* fromEdge, const OptEdge* toEdge);

  static const OptLO* getCtdLineIn(const shared::linegraph::Line* r,
                                   const shared::linegraph::LineNode* dir,
                                   const OptEdge* fromEdge,
                                   const OptEdge* toEdge);

  static const OptLO* getSameDirLineIn(const shared::linegraph::Line* r,
                                       const shared::linegraph::LineNode* dir,
                                       const OptEdge* fromEdge,
                                       const OptEdge* toEdge);

  static LnEdgPart getFirstLnEdgPart(const OptEdge*);
  static LnEdgPart getLastLnEdgPart(const OptEdge*);

  // apply splitting rules
  void splitSingleLineEdgs();
  void terminusDetach();
  void deleteSingleEdgeComponents();

 private:
  const OptGraphScorer* _scorer;
  void writeEdgeOrder();
  void updateEdgeOrder(OptNode* n);
  bool contractDeg2Step();

  bool untangleFullCross();
  bool untangleYStep();
  bool untanglePartialYStep();
  bool untangleDogBoneStep();
  bool untanglePartialDogBoneStep();
  bool untangleStumpStep();

  std::vector<OptNode*> explodeNodeAlong(OptNode* nd,
                                         const util::geo::PolyLine<double>& pl,
                                         size_t n);

  std::vector<OptEdge*> branchesAt(OptEdge* e, OptNode* n) const;
  bool branchesAtInto(OptEdge* e, OptNode* n,
                      std::vector<OptEdge*> branchesA) const;
  bool partiallyBranchesAtInto(OptEdge* e, OptNode* n,
                               std::vector<OptEdge*> branchesA) const;
  std::vector<OptEdge*> partiallyBranchesAt(OptEdge* e, OptNode* n) const;

  std::pair<OptEdge*, OptEdge*> isFullCross(OptNode* n) const;
  bool isYAt(OptEdge* e, OptNode* n) const;
  bool isPartialYAt(OptEdge* e, OptNode* n) const;
  std::pair<OptEdge*, OptEdge*> isStump(OptEdge* e) const;
  std::pair<OptEdge*, OptEdge*> isStumpAt(OptEdge* e, OptNode* n) const;

  std::set<const shared::linegraph::Line*> getLines() const;

  bool isDogBone(OptEdge* e) const;
  OptNode* isPartialDogBone(OptEdge* e) const;

  static void upFirstLastEdg(OptEdge*);

  static OptEdgePL getView(OptEdge* parent, OptEdge* leg, size_t offset);
  static OptEdgePL getPartialView(OptEdge* parent, OptEdge* leg, size_t offset);

  std::vector<size_t> mapPositions(std::vector<OptEdge*> a, OptEdge* leg,
                                   std::vector<OptEdge*> b) const;

  static bool dirLineContains(const OptEdge* a, const OptEdge* b);

  static bool dirLineEqualIn(const OptEdge* a, const OptEdge* b);
  static bool dirContinuedOver(const OptEdge* a, const OptEdge* b,
                               const OptEdge* c);
  static bool dirPartialContinuedOver(const OptEdge* a, const OptEdge* b);
  static bool dirContinuedOver(const OptLO& ro, const OptEdge* a,
                               const OptEdge* b);
  static bool dirContinuedOver(const OptLO& ro, const OptEdge* a,
                               const OptNode* n);

  static bool linesBranchAt(const OptLO& roA, const OptLO& roB,
                            const OptEdge* a, const OptNode* n);

  static bool linesCtnOver(const OptLO& roA, const OptLO& roB, const OptEdge* a,
                           const OptNode* n);

  static bool uniquelyExtendsOver(const OptLO& a, const OptEdge* e, const OptNode* n);

  bool contractCheaper(const OptNode* cont, const OptNode* cheaper,
                       const std::vector<OptLO>& lines) const;

  static bool lineDisjunct(const std::vector<const OptEdge*>& edges);

  static std::vector<OptLO> getCtdLinesIn(const OptEdge* fromEdge,
                                          const OptEdge* toEdge);

  static OptNode* sharedNode(const OptEdge* a, const OptEdge* b);

  static util::Nullable<const OptLO> getLO(const OptEdge* a,
                                           const shared::linegraph::Line*);

  static std::vector<OptEdge*> clockwEdges(OptEdge* noon, OptNode* n);
  static std::vector<OptEdge*> partialClockwEdges(OptEdge* noon, OptNode* n);

  static bool terminatesAt(const OptLO& lo, const OptEdge* e, const OptNode* nd);
  static bool terminatesAt(const OptEdge* e, const OptNode* nd);
};

// compare the orientation of two edges adjacent to some shared node
inline bool cmpEdge(const OptEdge* a, const OptEdge* b) {
  double angA, angB;

  // n is the shared node
  OptNode* n = 0;
  if (a->getFrom() == b->getFrom() || a->getFrom() == b->getTo())
    n = a->getFrom();
  else
    n = a->getTo();
  assert(n->pl().node);

  auto tgEdgeA = OptGraph::getAdjEdg(a, n);
  assert(tgEdgeA);

  angA = shared::rendergraph::RenderGraph::getOutAngle(n->pl().node, tgEdgeA);

  auto tgEdgeB = OptGraph::getAdjEdg(b, n);
  assert(tgEdgeB);

  angB = shared::rendergraph::RenderGraph::getOutAngle(n->pl().node, tgEdgeB);

  if (tgEdgeA == tgEdgeB) {
    // if these edges originally came from the same node front, use their
    // internal ordering
    // TODO: what if the edge is reversed?
    if ((a->getFrom() == n && b->getFrom() == n)) {
      if (OptGraph::getAdjLnEdgPart(a, n).dir) {
        assert(OptGraph::getAdjLnEdgPart(b, n).dir);
        return OptGraph::getAdjLnEdgPart(a, n).order <
               OptGraph::getAdjLnEdgPart(b, n).order;
      } else {
        return OptGraph::getAdjLnEdgPart(a, n).order >
               OptGraph::getAdjLnEdgPart(b, n).order;
      }
    } else if ((a->getTo() == n && b->getTo() == n)) {
      if (OptGraph::getAdjLnEdgPart(a, n).dir) {
        assert(OptGraph::getAdjLnEdgPart(b, n).dir);
        return OptGraph::getAdjLnEdgPart(a, n).order >
               OptGraph::getAdjLnEdgPart(b, n).order;
      } else {
        return OptGraph::getAdjLnEdgPart(a, n).order <
               OptGraph::getAdjLnEdgPart(b, n).order;
      }
    }
  }
  return fmod(angA + M_PI * 1.5, 2 * M_PI) > fmod(angB + M_PI * 1.5, 2 * M_PI);
}
}  // namespace optim
}  // namespace loom

#endif  // LOOM_GRAPH_OPTIM_OPTGRAPH_H_
