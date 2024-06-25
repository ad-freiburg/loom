// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <set>
#include "loom/optim/OptGraph.h"
#include "loom/optim/OptGraphScorer.h"
#include "shared/linegraph/Line.h"
#include "shared/linegraph/LineGraph.h"
#include "shared/rendergraph/RenderGraph.h"
#include "util/String.h"
#include "util/graph/Algorithm.h"
#include "util/log/Log.h"

using loom::optim::LnEdgPart;
using loom::optim::OptEdge;
using loom::optim::OptEdgePL;
using loom::optim::OptGraph;
using loom::optim::OptGraphScorer;
using loom::optim::OptLO;
using loom::optim::OptNode;
using loom::optim::OptNodePL;
using loom::optim::PartnerPath;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using shared::linegraph::LineOcc;
using shared::rendergraph::RenderGraph;
using util::Nullable;
using util::graph::Algorithm;

const static double DO = 100;

// _____________________________________________________________________________
void OptGraph::upFirstLastEdg(OptEdge* optEdg) {
  size_t i = 0;
  for (const auto e : optEdg->pl().lnEdgParts) {
    if (e.lnEdg->getFrom() == optEdg->getFrom()->pl().node ||
        e.lnEdg->getTo() == optEdg->getFrom()->pl().node) {
      break;
    }
    i++;
  }

  size_t j = 0;
  for (const auto e : optEdg->pl().lnEdgParts) {
    if (e.lnEdg->getFrom() == optEdg->getTo()->pl().node ||
        e.lnEdg->getTo() == optEdg->getTo()->pl().node) {
      break;
    }
    j++;
  }

  optEdg->pl().firstLnEdg = i;
  optEdg->pl().lastLnEdg = j;
}

// _____________________________________________________________________________
LnEdgPart OptGraph::getFirstLnEdgPart(const OptEdge* optEdg) {
  return optEdg->pl().lnEdgParts[optEdg->pl().firstLnEdg];
}

// _____________________________________________________________________________
LnEdgPart OptGraph::getLastLnEdgPart(const OptEdge* optEdg) {
  return optEdg->pl().lnEdgParts[optEdg->pl().lastLnEdg];
}

// _____________________________________________________________________________
const std::vector<OptLO>& OptEdgePL::getLines() const { return lines; }

// _____________________________________________________________________________
std::vector<OptLO>& OptEdgePL::getLines() { return lines; }

// _____________________________________________________________________________
size_t OptEdgePL::getCardinality() const { return getLines().size(); }

// _____________________________________________________________________________
const OptLO* OptEdgePL::getLineOcc(const Line* l) const {
  if (lines.size() && lines.back().line > l) return 0;
  if (lines.size() > 0 && lines[0].line == l) return &lines[0];
  if (lines.size() > 1 && lines[1].line == l) return &lines[1];
  if (lines.size() > 2 && lines[2].line == l) return &lines[2];

  if (lines.size() < 4) return 0;

  // we make use of the fact here that the lines are always sorted
  auto it = std::lower_bound(lines.begin(), lines.end(), l);

  if (it->line != l) return 0;

  return &*it;
}

// _____________________________________________________________________________
std::string OptEdgePL::getStrRepr() const {
  const void* address = static_cast<const void*>(this);
  std::stringstream ss;
  ss << address;

  return ss.str();
}

// _____________________________________________________________________________
LnEdgPart OptGraph::getAdjLnEdgPart(const OptEdge* e, const OptNode* n) {
  if (e->getFrom() == n) {
    return getFirstLnEdgPart(e);
  } else if (e->getTo() == n) {
    return getLastLnEdgPart(e);
  }

  // TODO: throw exception
  assert(false);
}

// _____________________________________________________________________________
LineEdge* OptGraph::getAdjEdg(const OptEdge* e, const OptNode* n) {
  if (e->getFrom() == n) {
    return getFirstLnEdgPart(e).lnEdg;
  } else if (e->getTo() == n) {
    return getLastLnEdgPart(e).lnEdg;
  }

  return 0;
}

// _____________________________________________________________________________
std::map<const LineNode*, OptNode*> OptGraph::build(RenderGraph* rg) {
  std::map<const LineNode*, OptNode*> lnNdToOptNd;
  for (auto n : rg->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;

      auto frLnNd = e->getFrom();
      auto toLnNd = e->getTo();

      OptNode* from;
      OptNode* to;

      auto fromI = lnNdToOptNd.find(frLnNd);
      auto toI = lnNdToOptNd.find(toLnNd);

      if (fromI == lnNdToOptNd.end()) {
        fromI = lnNdToOptNd.insert({frLnNd, addNd(frLnNd)}).first;
      }

      if (toI == lnNdToOptNd.end()) {
        toI = lnNdToOptNd.insert({toLnNd, addNd(toLnNd)}).first;
      }

      from = fromI->second;
      to = toI->second;

      OptEdge* edge = addEdg(from, to);

      edge->pl().lnEdgParts.push_back(LnEdgPart(e, true));
      for (auto roOld : e->pl().getLines()) {
        edge->pl().lines.push_back(OptLO(roOld.line, roOld.direction));
      }

      std::sort(edge->pl().lines.begin(), edge->pl().lines.end());
    }
  }

  writeEdgeOrder();

  return lnNdToOptNd;
}

// _____________________________________________________________________________
void OptGraph::terminusDetach() {
  std::vector<std::pair<OptEdge*, OptNode*>> toDetach;

  // collect edges to cut
  for (OptNode* n : getNds()) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;

      // we only cut if the terminating node has a degree greater than 1
      if (e->getFrom()->getDeg() > 1 && terminatesAt(e, e->getFrom())) {
        toDetach.push_back({e, e->getFrom()});
        continue;
      }

      if (e->getTo()->getDeg() > 1 && terminatesAt(e, e->getTo())) {
        toDetach.push_back({e, e->getTo()});
        continue;
      }
    }
  }

  for (auto ePair : toDetach) {
    OptEdge* e = ePair.first;
    OptNode* n = ePair.second;

    if (n->getDeg() < 2)
      continue;  // may happen if we have detached an edge
                 // from the other side

    OptNode* eFrom = e->getFrom();
    OptNode* eTo = e->getTo();

    if (n == eFrom) {
      OptNode* newFrom = addNd(e->getFrom()->pl().node);
      addEdg(newFrom, eTo, e->pl());

      delEdg(eFrom, eTo);
      updateEdgeOrder(newFrom);
      updateEdgeOrder(eFrom);
      updateEdgeOrder(eTo);
    } else {
      OptNode* newTo = addNd(e->getTo()->pl().node);
      addEdg(eFrom, newTo, e->pl());

      delEdg(eFrom, eTo);
      updateEdgeOrder(newTo);
      updateEdgeOrder(eFrom);
      updateEdgeOrder(eTo);
    }
  }
}

// _____________________________________________________________________________
void OptGraph::splitSingleLineEdgs() {
  std::vector<OptEdge*> toCut;

  // collect edges to cut
  for (OptNode* n : getNds()) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;

      // we only cut if both nodes have a deg > 1 - otherwise, we would cut
      // forever...
      if (e->getFrom()->getDeg() < 2 || e->getTo()->getDeg() < 2) continue;

      if (e->pl().getCardinality() == 1) toCut.push_back(e);
    }
  }

  for (OptEdge* e : toCut) {
    // cut this edge

    // add two new opt nodes, with slight offsets at the center of this
    // edge (this is for debug visualization purposes only)

    util::geo::PolyLine<double> edgeGeom(*e->getFrom()->pl().getGeom(),
                                         *e->getTo()->pl().getGeom());

    OptNode* leftN = addNd(edgeGeom.getPointAt(0.45).p);
    OptNode* rightN = addNd(edgeGeom.getPointAt(0.55).p);
    OptNode* eFrom = e->getFrom();
    OptNode* eTo = e->getTo();

    addEdg(eFrom, leftN, e->pl());
    auto newE = addEdg(rightN, eTo, e->pl());
    for (auto& lnEdgPart : newE->pl().lnEdgParts) lnEdgPart.wasCut = true;

    delEdg(eFrom, eTo);

    updateEdgeOrder(leftN);
    updateEdgeOrder(rightN);

    updateEdgeOrder(eFrom);
    updateEdgeOrder(eTo);
  }
}

// _____________________________________________________________________________
void OptGraph::partnerLines() {
  const auto& partners = getPartnerLines();

  for (const auto& p : partners) {
    if (p.partners.size() == 1) continue;

    for (size_t i = 0; i < p.path.size(); i++) {
      // TODO: why isnt there a getLine function f or OptEdgePL?
      auto e = p.path[i];
      auto it = e->pl().getLines().begin();

      while (it != e->pl().getLines().end()) {
        auto& ro = *it;
        if (ro == *p.partners.begin()) {
          ro.relatives.clear();
          for (auto partner : p.partners) ro.relatives.push_back(partner.line);
          if (p.inv[i]) std::reverse(ro.relatives.begin(), ro.relatives.end());
        } else if (p.partners.count(ro)) {
          it = e->pl().getLines().erase(it);
          continue;
        }
        it++;
      }
    }
  }
}

// _____________________________________________________________________________
std::set<const Line*> OptGraph::getLines() const {
  std::set<const Line*> lines;
  // TODO: there has to be a faster way
  for (auto n : getNds()) {
    for (auto e : n->getAdjList()) {
      for (auto l : e->pl().getLines()) {
        lines.insert(l.line);
      }
    }
  }

  return lines;
}

// _____________________________________________________________________________
Nullable<const OptLO> OptGraph::getLO(const OptEdge* a, const Line* r) {
  auto lo = a->pl().getLineOcc(r);
  if (lo) return *lo;
  return Nullable<const OptLO>();
}

// _____________________________________________________________________________
std::vector<PartnerPath> OptGraph::getPartnerLines() const {
  std::vector<PartnerPath> ret;

  for (auto rt : getLines()) {
    // create connected components w.r.t. route (each component consists only of
    // edges containing route rt)

    struct Check : public Algorithm::EdgeCheckFunc<OptNodePL, OptEdgePL> {
      Check(const Line* l) : rt(l){};
      const Line* rt;
      virtual bool operator()(const OptNode* frNd, const OptEdge* edg) const {
        UNUSED(frNd);
        // TODO: this will put lines not continueing over a particular edge
        // in one component, as only the lines' presence on the edge is checked.
        // Later on, we then cannot follow the edge through the entire path,
        // which may lead to collapsing opportunities being missed.
        auto ro = getLO(edg, rt);
        return !ro.isNull();
      };
    };

    const auto& comps = Algorithm::connectedComponents(*this, Check(rt));

    for (const auto& comp : comps) {
      if (comp.size() < 2) continue;
      auto p = pathFromComp(comp);
      if (p.partners.size() && p.path.size()) ret.push_back(p);
    }
  }

  return ret;
}

// _____________________________________________________________________________
PartnerPath OptGraph::pathFromComp(const std::set<OptNode*>& comp) const {
  std::vector<OptNode*> nPath;
  nPath.push_back(*comp.begin());
  PartnerPath pp;

  // search for an entry node or abort if a path splits
  OptNode* entry = 0;
  for (auto n : comp) {
    size_t compDeg = 0;
    for (auto e : n->getAdjList()) {
      if (!comp.count(e->getOtherNd(n))) continue;
      compDeg++;
      pp.partners.insert(e->pl().getLines().begin(), e->pl().getLines().end());
    }
    if (compDeg > 2) return PartnerPath();  // not a simple path!
    if (compDeg == 1 && !entry) entry = n;
    assert(compDeg != 0);
  }

  if (!entry) entry = *comp.begin();  // we have a cycle, take any node as start
  auto cur = entry;

  std::set<OptLO> notPartner;

  // backtrack
  while (true) {
    OptNode* toNd = 0;

    // get the continuing edge cntE
    OptEdge* cntE = 0;
    for (auto e : cur->getAdjList()) {
      // dont follow edge we came from
      if (pp.path.size() && e == pp.path.back()) continue;
      auto toNdCand = e->getOtherNd(cur);
      if (comp.count(toNdCand)) {
        cntE = e;
        break;
      }
    }

    // check if either the new edge, or the old edge, have any lines
    // continueing into any other edge and add them to the notPartner set
    for (auto e : cur->getAdjList()) {
      // dont use the edge we came from
      if (pp.path.size() && e == pp.path.back()) continue;
      if (e == cntE) continue;

      if (pp.path.size()) {
        const auto& ctdLines = getCtdLinesIn(pp.path.back(), e);
        notPartner.insert(ctdLines.begin(), ctdLines.end());
      }

      if (cntE) {
        const auto& ctdLines = getCtdLinesIn(e, cntE);
        notPartner.insert(ctdLines.begin(), ctdLines.end());
      }
    }

    // no continuation edge found, break
    if (!cntE) break;

    if (pp.path.size()) {
      const auto& ctd = getCtdLinesIn(pp.path.back(), cntE);
      std::set<OptLO> newPartners;
      for (auto ctdRt : ctd) {
        if (pp.partners.count(ctdRt)) newPartners.insert(ctdRt);
      }
      pp.partners = newPartners;
    } else {
      const auto& ctd = cntE->pl().getLines();
      pp.partners.insert(ctd.begin(), ctd.end());
    }

    pp.path.push_back(cntE);
    pp.inv.push_back(cntE->getOtherNd(cur) == cntE->getFrom());
    toNd = cntE->getOtherNd(cur);

    cur = toNd;

    // we had a circle and are back at the beginning
    if (cur == entry) break;
  }

  // remove routes continuing outside the component from the partners set
  for (const auto& r : notPartner) pp.partners.erase(r);

  return pp;
}

// _____________________________________________________________________________
void OptGraph::contractDeg2Nds() {
  while (contractDeg2Step()) {
  }
}

// _____________________________________________________________________________
void OptGraph::untangle() {
  untangleDoubleStump();

  untangleOuterStump();

  while (untangleFullX()) {
  }

  untangleY();

  untanglePartialY();

  untangleDogBone();

  untanglePartialDogBone();

  untangleInnerStump();
}

// _____________________________________________________________________________
bool OptGraph::uniquelyExtendsOver(const OptLO& lo, const OptEdge* eA,
                                   const OptNode* nd) {
  size_t c = 0;
  for (const auto eB : nd->getAdjList()) {
    if (eB == eA) continue;
    if (dirContinuedOver(lo, eA, eB)) c++;
    if (c > 1) return false;
  }

  return c == 1;
}

// _____________________________________________________________________________
bool OptGraph::dogBoneCheaper(const OptNode* a, const OptNode* b,
                              const std::vector<OptLO>& lines) const {
  // this currently only checks whether the crossing cost is cheaper, as we do
  // not yet have line-dependend costs
  for (const auto& loA : lines) {
    for (const auto& loB : lines) {
      if (loA == loB) continue;

      if (_scorer->getCrossingPenDiffSeg(a) >
          _scorer->getCrossingPenDiffSeg(b)) {
        return false;
      }
    }
  }

  return true;
}

// _____________________________________________________________________________
bool OptGraph::contractCheaper(const OptNode* cont, const OptNode* cheaper,
                               const std::vector<OptLO>& lines) const {
  if (cheaper->getDeg() == 1) return true;  // trivially a terminus node

  // check if the lines *uniquely* extend over cheaper, otherwise don't contract
  for (const auto& lo : lines) {
    if ((!terminatesAt(lo, getEdg(cont, cheaper), cheaper) &&
         !uniquelyExtendsOver(lo, getEdg(cont, cheaper), cheaper))) {
      return false;
    }
  }

  // for degree > 2, check if true for each line pair - it might be that they
  // cross into different segments there

  for (const auto& loA : lines) {
    for (const auto& loB : lines) {
      if (loA == loB) continue;

      if (!contractCheaper(cont, cheaper, loA, loB)) {
        return false;
      }
    }
  }

  return true;
}

// _____________________________________________________________________________
bool OptGraph::contractCheaper(const OptNode* cont, const OptNode* cheaper,
                               const OptLO& loA, const OptLO& loB) const {
  auto edg = getEdg(cont, cheaper);
  if (linesCtnOver(loA, loB, edg, cheaper)) {
    // we have to consider same segment crossings and separations
    if (_scorer->getSeparationPen(cont) < _scorer->getSeparationPen(cheaper) ||
        _scorer->getCrossingPenSameSeg(cont) <
            _scorer->getCrossingPenSameSeg(cheaper)) {
      return false;
    }
  }

  if (linesBranchAt(loA, loB, edg, cheaper)) {
    // we have to consider diff segment crossings in the cheaper node,
    // but same segment crossings in the cont node (no diff crossings are
    // possible there)
    if (_scorer->getCrossingPenSameSeg(cont) <
        _scorer->getCrossingPenDiffSeg(cheaper)) {
      return false;
    }
  }

  return true;
}

// _____________________________________________________________________________
bool OptGraph::contractDeg2Step() {
  for (OptNode* n : getNds()) {
    if (n->getDeg() == 2) {
      OptEdge* first = n->getAdjList().front();
      OptEdge* second = n->getAdjList().back();

      assert(n->pl().node);

      if (dirLineEqualIn(first, second)) {
        // if both edges have more than 2 lines, only contract if we can move
        // potential crossings to a cheaper location
        if (first->pl().getCardinality() > 1) {
          if (!contractCheaper(n, first->getOtherNd(n),
                               first->pl().getLines()) &&
              !contractCheaper(n, second->getOtherNd(n),
                               first->pl().getLines()))
            continue;
        }

        OptNode* newFrom = 0;
        OptNode* newTo = 0;

        bool firstReverted;
        bool secondReverted;

        // add new edge
        if (first->getTo() != n) {
          newFrom = first->getTo();
          firstReverted = true;
        } else {
          newFrom = first->getFrom();
          firstReverted = false;
        }

        if (second->getTo() != n) {
          newTo = second->getTo();
          secondReverted = false;
        } else {
          newTo = second->getFrom();
          secondReverted = true;
        }

        // Important: dont create a multigraph, dont add self-edges
        if (newFrom == newTo || getEdg(newFrom, newTo)) continue;

        OptEdge* newEdge = addEdg(newFrom, newTo);

        // add lnEdgParts...
        for (LnEdgPart& lnEdgPart : first->pl().lnEdgParts) {
          newEdge->pl().lnEdgParts.push_back(
              LnEdgPart(lnEdgPart.lnEdg, (lnEdgPart.dir ^ firstReverted),
                        lnEdgPart.order, lnEdgPart.wasCut));
        }

        for (LnEdgPart& lnEdgPart : second->pl().lnEdgParts) {
          newEdge->pl().lnEdgParts.push_back(
              LnEdgPart(lnEdgPart.lnEdg, (lnEdgPart.dir ^ secondReverted),
                        lnEdgPart.order, lnEdgPart.wasCut));
        }

        upFirstLastEdg(newEdge);

        newEdge->pl().depth = std::max(first->pl().depth, second->pl().depth);

        newEdge->pl().lines = first->pl().lines;

        // update direction markers
        for (auto& ro : newEdge->pl().lines) {
          if (ro.dir == n->pl().node) ro.dir = newTo->pl().node;
        }

        assert(newFrom != n);
        assert(newTo != n);

        delNd(n);

        updateEdgeOrder(newFrom);
        updateEdgeOrder(newTo);

        return true;
      }
    }
  }

  return false;
}

// _____________________________________________________________________________
size_t OptGraph::getNumNodes() const { return getNds().size(); }

// _____________________________________________________________________________
size_t OptGraph::getNumNodes(bool topo) const {
  size_t ret = 0;
  for (auto n : getNds()) {
    if ((n->pl().node->pl().stops().size() == 0) ^ !topo) ret++;
  }

  return ret;
}

// _____________________________________________________________________________
size_t OptGraph::getNumEdges() const {
  size_t ret = 0;

  for (auto n : getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      ret += 1;
    }
  }

  return ret;
}

// _____________________________________________________________________________
size_t OptGraph::getNumLines() const {
  std::set<const Line*> lines;

  for (auto n : getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      for (auto to : getFirstLnEdgPart(e).lnEdg->pl().getLines()) {
        lines.insert(to.line);
      }
    }
  }

  return lines.size();
}

// _____________________________________________________________________________
size_t OptGraph::getMaxCardinality() const {
  size_t ret = 0;
  for (auto n : getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (getFirstLnEdgPart(e).lnEdg->pl().getLines().size() > ret) {
        ret = getFirstLnEdgPart(e).lnEdg->pl().getLines().size();
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
const util::geo::Line<double>* OptEdgePL::getGeom() { return 0; }

// _____________________________________________________________________________
util::json::Dict OptEdgePL::getAttrs() {
  util::json::Dict ret;
  std::string lines;

  for (const auto& r : getLines()) {
    if (r.relatives.size() > 1)
      lines += r.line->label() + "(x" + util::toString(r.relatives.size()) +
               ")" + "[" + r.line->color() + ", -> " + util::toString(r.dir) +
               "], ";
    else
      lines += r.line->label() + "[" + r.line->color() + ", -> " +
               util::toString(r.dir) + "], ";
  }

  ret["lines"] = lines;
  ret["num_line_edge_parts"] = util::toString(lnEdgParts.size());
  ret["depth"] = util::toString(depth);

  return ret;
}

// _____________________________________________________________________________
std::string OptEdgePL::toStr() const {
  std::string lines;
  for (const auto& r : getLines()) {
    if (r.relatives.size() > 1)
      lines += r.line->label() + "(x" + util::toString(r.relatives.size()) +
               ")" + "[" + r.line->color() + ", -> " + util::toString(r.dir) +
               "], ";
    else
      lines += r.line->label() + "[" + r.line->color() + ", -> " +
               util::toString(r.dir) + "], ";
  }

  return "[" + lines + "]";
}

// _____________________________________________________________________________
const util::geo::Point<double>* OptNodePL::getGeom() { return &p; }

// _____________________________________________________________________________
util::json::Dict OptNodePL::getAttrs() {
  util::json::Dict ret;
  ret["orig_node"] = util::toString(node);
  if (!node) {
    ret["stat_name"] = "<dummy>";
  } else {
    if (node->pl().stops().size()) {
      ret["stat_name"] = node->pl().stops()[0].name;
    } else {
      ret["stat_name"] = "";
    }
  }
  std::string edges;
  for (auto e : circOrdering) {
    edges += util::toString(e) + ",";
  }

  ret["edge_order"] = edges;
  return ret;
}

// _____________________________________________________________________________
std::pair<OptEdge*, OptEdge*> OptGraph::isFullX(OptNode* n) const {
  if (n->getDeg() < 3) return std::pair<OptEdge*, OptEdge*>(0, 0);
  std::pair<OptEdge*, OptEdge*> ret(0, 0);

  for (auto ea : n->getAdjList()) {
    for (auto eb : n->getAdjList()) {
      if (ea == eb) continue;
      if (dirLineContains(eb, ea) && dirLineContains(ea, eb)) {
        ret.first = ea;
        ret.second = eb;
      }
    }

    if (ret.first) {
      bool nope = false;
      for (auto e : n->getAdjList()) {
        if (e == ret.second || e == ret.first) continue;
        if (dirPartialContinuedOver(ret.first, e) ||
            dirPartialContinuedOver(ret.second, e)) {
          nope = true;
          break;
        }
      }
      if (!nope) return ret;
    }
  }

  return std::pair<OptEdge*, OptEdge*>(0, 0);
}

// _____________________________________________________________________________
bool OptGraph::untangleFullX() {
  for (OptNode* n : getNds()) {
    std::pair<OptEdge*, OptEdge*> cross;
    if ((cross = isFullX(n)).first) {
      LOGTO(DEBUG, std::cerr)
          << "Found full cross at node " << n << " between " << cross.first
          << "(" << cross.first->pl().toStr() << ") and " << cross.second
          << " (" << cross.second->pl().toStr() << ")";

      auto newN = addNd(util::geo::DPoint(n->pl().getGeom()->getX() + DO,
                                          n->pl().getGeom()->getY() + DO));
      newN->pl().node = n->pl().node;

      if (cross.first->getFrom() == n) {
        addEdg(newN, cross.first->getTo(), cross.first->pl());
      } else {
        addEdg(cross.first->getFrom(), newN, cross.first->pl());
      }

      if (cross.second->getFrom() == n) {
        addEdg(newN, cross.second->getTo(), cross.second->pl());
      } else {
        addEdg(cross.second->getFrom(), newN, cross.second->pl());
      }

      auto fa = cross.first->getFrom();
      auto fb = cross.first->getTo();
      auto sa = cross.second->getFrom();
      auto sb = cross.second->getTo();

      delEdg(cross.first->getFrom(), cross.first->getTo());
      delEdg(cross.second->getFrom(), cross.second->getTo());

      updateEdgeOrder(n);
      updateEdgeOrder(newN);
      updateEdgeOrder(fa);
      updateEdgeOrder(fb);
      updateEdgeOrder(sa);
      updateEdgeOrder(sb);

      return true;
    }
  }
  return false;
}

// _____________________________________________________________________________
std::vector<OptNode*> OptGraph::explodeNodeAlong(
    OptNode* nd, const util::geo::PolyLine<double>& pl, size_t n) {
  std::vector<OptNode*> ret(n);
  for (size_t i = 0; i < n; i++) {
    double p = (n - 1 - i) / (double)n;
    ret[i] = addNd(pl.getPointAt(p).p);
    ret[i]->pl().node = nd->pl().node;
  }

  return ret;
}

// _____________________________________________________________________________
void OptGraph::untanglePartialY() {
  std::vector<OptEdge*> toUntangle;

  for (OptNode* na : getNds()) {
    if (na->getDeg() != 1) continue;  // only look at terminus nodes

    // the only outgoing edge
    OptEdge* ea = na->getAdjList().front();
    OptNode* nb = ea->getOtherNd(na);

    if (isPartialYAt(ea, nb)) {
      assert(nb->pl().node);
      assert(na->pl().node);

      LOGTO(DEBUG, std::cerr)
          << "Found partial Y at node " << nb << " with main leg " << ea << " ("
          << ea->pl().toStr() << ")";
      toUntangle.push_back(ea);
    }
  }

  for (auto ea : toUntangle) {
    // the only outgoing edge
    OptNode* na = ea->getFrom();
    if (na->getDeg() != 1) na = ea->getTo();

    assert(na->getDeg() == 1);
    OptNode* nb = ea->getOtherNd(na);

    assert(nb->pl().node);
    assert(na->pl().node);

    // the geometry of the main leg
    util::geo::PolyLine<double> pl(*nb->pl().getGeom(), *na->pl().getGeom());
    double bandW = (nb->getDeg() - 1) * (DO / (ea->pl().depth + 1));
    auto ortho = pl.getOrthoLineAtDist(pl.getLength(), bandW);

    // each leg, in clockwise fashion
    auto minLgs = partialClockwEdges(ea, nb);

    assert(minLgs.size() <= ea->pl().getCardinality());

    // for each minor leg of the Y, create a new node at the origin
    std::vector<OptNode*> origNds = explodeNodeAlong(na, ortho, minLgs.size());

    size_t offset = 0;
    for (size_t i = 0; i < minLgs.size(); i++) {
      size_t j = i;
      OptEdgePL pl;
      if (ea->getFrom() == nb) {
        if (ea->pl().lnEdgParts[0].dir) j = minLgs.size() - 1 - i;
        pl = getPartialView(ea, minLgs[j], offset);
        addEdg(nb, origNds[j], pl);
      } else {
        if (!ea->pl().lnEdgParts[0].dir) j = minLgs.size() - 1 - i;
        pl = getPartialView(ea, minLgs[j], offset);
        addEdg(origNds[j], nb, pl);
      }
      offset += pl.getLines().size();
    }

    // delete remaining stuff
    delNd(na);

    // update orderings
    for (auto n : origNds) updateEdgeOrder(n);
    updateEdgeOrder(nb);
  }
}

// _____________________________________________________________________________
const OptLO* OptGraph::isDoubleStump(OptEdge* e) const {
  // nothing to gain...
  if (e->getFrom()->getDeg() < 2 && e->getTo()->getDeg() < 2) return 0;

  for (const auto& lo : e->pl().getLines()) {
    if (terminatesAt(lo, e, e->getFrom()) && terminatesAt(lo, e, e->getTo()))
      return &lo;
  }
  return 0;
}

// _____________________________________________________________________________
std::pair<OptEdge*, bool> OptGraph::isOuterStump(OptEdge* e) const {
  std::pair<OptEdge*, bool> ret = {0, 0};
  if ((ret = isOuterStumpAt(e, e->getFrom())).first) return ret;
  if ((ret = isOuterStumpAt(e, e->getTo())).first) return ret;

  return {0, 0};
}

// _____________________________________________________________________________
std::pair<OptEdge*, bool> OptGraph::isOuterStumpAt(OptEdge* mainLeg,
                                                   OptNode* n) const {
  if (n->getDeg() < 3 || mainLeg->getOtherNd(n)->getDeg() < 2) return {0, 0};

  auto clockwEdgs = clockwEdges(mainLeg, n);

  size_t extends = 0;

  // check that at least two edges have lines that extend over n into the mainl
  // (one edge is the stump edge)
  for (auto branch : n->getAdjList()) {
    if (branch == mainLeg) continue;
    if (getCtdLinesIn(branch, mainLeg).size() > 0) extends++;
    if (extends > 1) break;
  }

  if (extends < 2) return {0, 0};
  extends = 0;

  // check that at least one edge has lines that extend over other(n) into mainl
  for (auto branch : mainLeg->getOtherNd(n)->getAdjList()) {
    if (branch == mainLeg) continue;
    if (getCtdLinesIn(branch, mainLeg).size() > 0) {
      extends++;
      break;
    }
  }

  if (extends == 0) return {0, 0};

  if (dirLineContains(mainLeg, clockwEdgs.front()) &&
      terminatesAt(clockwEdgs.front(), mainLeg, mainLeg->getOtherNd(n))) {
    bool disj = true;
    for (size_t i = 1; i < clockwEdgs.size(); i++) {
      if (!lineDisjunct({clockwEdgs.front(), clockwEdgs[i]})) disj = false;
    }
    if (disj) return {clockwEdgs.front(), 1};
  }

  if (dirLineContains(mainLeg, clockwEdgs.back()) &&
      terminatesAt(clockwEdgs.back(), mainLeg, mainLeg->getOtherNd(n))) {
    bool disj = true;
    for (size_t i = 0; i < clockwEdgs.size() - 1; i++) {
      if (!lineDisjunct({clockwEdgs.back(), clockwEdgs[i]})) disj = false;
    }
    if (disj) return {clockwEdgs.back(), 0};
  }

  return {0, 0};
}

// _____________________________________________________________________________
void OptGraph::untangleDoubleStump() {
  std::vector<OptEdge*> toUntangle;

  for (OptNode* n : getNds()) {
    for (OptEdge* mainLeg : n->getAdjList()) {
      if (mainLeg->getFrom() != n) continue;

      const OptLO* stump;
      if ((stump = isDoubleStump(mainLeg))) {
        LOGTO(DEBUG, std::cerr)
            << "Found double stump with main leg " << mainLeg << " ("
            << mainLeg->pl().toStr() << ") with stump " << stump->line->id();
        toUntangle.push_back(mainLeg);
      }
    }
  }

  for (auto mainLeg : toUntangle) {
    const OptLO* stump = isDoubleStump(mainLeg);
    OptEdgePL plMain = getPartialViewExcl(mainLeg, stump, 0);
    OptEdgePL plStump =
        getPartialView(mainLeg, stump, plMain.getLines().size());

    mainLeg->pl() = plMain;

    auto stNdA = addNd(mainLeg->getFrom()->pl());
    auto stNdB = addNd(mainLeg->getTo()->pl());
    addEdg(stNdA, stNdB, plStump);
  }
}

// _____________________________________________________________________________
void OptGraph::untangleOuterStump() {
  std::set<OptEdge*> toUntangle;

  for (OptNode* n : getNds()) {
    for (OptEdge* mainLeg : n->getAdjList()) {
      if (mainLeg->getFrom() != n) continue;

      std::pair<OptEdge*, bool> stumpEdgPair;
      if ((stumpEdgPair = isOuterStump(mainLeg)).first) {
        LOGTO(DEBUG, std::cerr)
            << "Found outer stump with main leg " << mainLeg << " ("
            << mainLeg->pl().toStr() << ") at node " << stumpEdgPair.second
            << " with stump " << stumpEdgPair.first << " ("
            << stumpEdgPair.first->pl().toStr() << ")";
        assert(mainLeg->getFrom());
        assert(mainLeg->getTo());
        toUntangle.insert(mainLeg);
      }
    }
  }

  for (auto mainLeg : toUntangle) {
    auto stumpEdgPair = isOuterStump(mainLeg);
    // check here again because the main leg may break up if there were
    // only 2 lines on it in a previous outer stump untangle, this should be
    // explicitely checked above
    if (!stumpEdgPair.first) continue;
    OptEdge* stumpEdg = stumpEdgPair.first;
    bool clockw = stumpEdgPair.second;
    OptNode* stumpN = sharedNode(mainLeg, stumpEdg);
    OptNode* notStumpN = mainLeg->getOtherNd(stumpN);

    // the geometry of the main leg
    util::geo::PolyLine<double> poly(*notStumpN->pl().getGeom(),
                                     *stumpN->pl().getGeom());
    double bandW = (DO / (mainLeg->pl().depth + 1));
    auto ortho = poly.getOrthoLineAtDist(poly.getLength(), bandW);

    // create two new nodes at the stump node
    std::vector<OptNode*> stumpNds = explodeNodeAlong(stumpN, ortho, 2);

    auto ortho2 = poly.getOrthoLineAtDist(0, bandW);

    // create two new nodes at the non-stump node
    std::vector<OptNode*> notStumpNds = explodeNodeAlong(notStumpN, ortho2, 2);

    size_t mainLegNode, stumpNode;

    // we assume that the stump edge is on the right

    if (mainLeg->getFrom() == stumpN) {
      size_t j;

      if (mainLeg->pl().lnEdgParts[0].dir)
        j = 1;
      else
        j = 0;

      mainLegNode = j;
      stumpNode = 1 - j;

      if (clockw ^ !mainLeg->pl().lnEdgParts[0].dir) {
        OptEdgePL pl = getPartialViewExcl(mainLeg, stumpEdg, 0);
        addEdg(stumpNds[mainLegNode], notStumpNds[mainLegNode], pl);

        pl = getPartialView(mainLeg, stumpEdg, pl.getLines().size());
        addEdg(stumpNds[stumpNode], notStumpNds[stumpNode], pl);
      } else {
        OptEdgePL pl = getPartialView(mainLeg, stumpEdg, 0);
        addEdg(stumpNds[stumpNode], notStumpNds[stumpNode], pl);

        pl = getPartialViewExcl(mainLeg, stumpEdg, pl.getLines().size());
        addEdg(stumpNds[mainLegNode], notStumpNds[mainLegNode], pl);
      }
    } else {
      size_t j;

      if (mainLeg->pl().lnEdgParts[0].dir)
        j = 1;
      else
        j = 0;

      mainLegNode = j;
      stumpNode = 1 - j;

      if (clockw ^ !mainLeg->pl().lnEdgParts[0].dir) {
        OptEdgePL pl = getPartialView(mainLeg, stumpEdg, 0);
        addEdg(notStumpNds[stumpNode], stumpNds[stumpNode], pl);

        pl = getPartialViewExcl(mainLeg, stumpEdg, pl.getLines().size());
        addEdg(notStumpNds[mainLegNode], stumpNds[mainLegNode], pl);
      } else {
        OptEdgePL pl = getPartialViewExcl(mainLeg, stumpEdg, 0);
        addEdg(notStumpNds[mainLegNode], stumpNds[mainLegNode], pl);

        pl = getPartialView(mainLeg, stumpEdg, pl.getLines().size());
        addEdg(notStumpNds[stumpNode], stumpNds[stumpNode], pl);
      }
    }

    if (stumpEdg->getFrom() == stumpN)
      addEdg(stumpNds[stumpNode], stumpEdg->getTo(), stumpEdg->pl());
    else
      addEdg(stumpEdg->getFrom(), stumpNds[stumpNode], stumpEdg->pl());

    for (auto e : stumpN->getAdjList()) {
      if (e == stumpEdg || e == mainLeg) continue;

      OptEdge* newE = 0;

      if (e->getFrom() == stumpN)
        newE = addEdg(stumpNds[mainLegNode], e->getTo(), e->pl());
      else
        newE = addEdg(e->getFrom(), stumpNds[mainLegNode], e->pl());

      // Important: the untangling rule assumes that if the node of the main
      // leg are split, the original node is kept. This is not the case in this
      // implementation - the original node is deleted, and a new split node
      // is inserted. But the original node may have *another* dangling outer
      // stump split in the toUntangle list, which we must replace here
      if (toUntangle.count(e)) {
        toUntangle.erase(e);
        toUntangle.insert(newE);
      }
    }

    for (auto e : notStumpN->getAdjList()) {
      if (e == mainLeg) continue;

      OptEdge* newE = 0;

      if (e->getFrom() == notStumpN)
        newE = addEdg(notStumpNds[mainLegNode], e->getTo(), e->pl());
      else
        newE = addEdg(e->getFrom(), notStumpNds[mainLegNode], e->pl());

      // see above
      if (toUntangle.count(e)) {
        toUntangle.erase(e);
        toUntangle.insert(newE);
      }
    }

    delNd(stumpN);
    delNd(notStumpN);

    // update orderings
    for (auto n : stumpNds) {
      updateEdgeOrder(n);
      for (auto e : n->getAdjList()) updateEdgeOrder(e->getOtherNd(n));
    }
    for (auto n : notStumpNds) {
      updateEdgeOrder(n);
      for (auto e : n->getAdjList()) updateEdgeOrder(e->getOtherNd(n));
    }
  }
}

// _____________________________________________________________________________
void OptGraph::untangleY() {
  std::vector<OptEdge*> toUntangle;

  for (OptNode* na : getNds()) {
    if (na->getDeg() != 1) continue;  // only look at terminus nodes

    // the only outgoing edge
    OptEdge* ea = na->getAdjList().front();
    OptNode* nb = ea->getOtherNd(na);

    if (isYAt(ea, nb)) {
      assert(nb->pl().node);
      assert(na->pl().node);
      LOGTO(DEBUG, std::cerr)
          << "Found full Y at node " << nb << " with main leg " << ea << " ("
          << ea->pl().toStr() << ")";
      toUntangle.push_back(ea);
    }
  }

  for (auto ea : toUntangle) {
    // the only outgoing edge
    OptNode* na = ea->getFrom();
    if (na->getDeg() != 1) na = ea->getTo();

    assert(na->getDeg() == 1);
    OptNode* nb = ea->getOtherNd(na);

    // the geometry of the main leg
    util::geo::PolyLine<double> pl(*nb->pl().getGeom(), *na->pl().getGeom());
    double bandW = (nb->getDeg() - 1) * (DO / (ea->pl().depth + 1));
    auto orthoPl = pl.getOrthoLineAtDist(0, bandW);
    auto orthoPlOrig = pl.getOrthoLineAtDist(pl.getLength(), bandW);

    // for each minor leg of the stump, create a new node at the origin
    std::vector<OptNode*> centerNds =
        explodeNodeAlong(nb, orthoPl, nb->getDeg() - 1);
    std::vector<OptNode*> origNds =
        explodeNodeAlong(na, orthoPlOrig, nb->getDeg() - 1);

    // each leg, in clockwise fashion
    auto minLgs = clockwEdges(ea, nb);
    std::vector<OptEdge*> minLgsN(minLgs.size());
    assert(minLgs.size() == nb->pl().circOrdering.size() - 1);

    for (size_t i = 0; i < minLgs.size(); i++) {
      if (minLgs[i]->getFrom() == nb)
        minLgsN[i] = addEdg(centerNds[i], minLgs[i]->getTo(), minLgs[i]->pl());
      else
        minLgsN[i] =
            addEdg(minLgs[i]->getFrom(), centerNds[i], minLgs[i]->pl());
    }

    size_t offset = 0;
    for (size_t i = 0; i < minLgs.size(); i++) {
      size_t j = i;
      OptEdgePL pl;
      if (ea->getFrom() == nb) {
        if (ea->pl().lnEdgParts[0].dir) j = minLgs.size() - 1 - i;
        pl = getView(ea, minLgs[j], offset);
        addEdg(centerNds[j], origNds[j], pl);
      } else {
        if (!ea->pl().lnEdgParts[0].dir) j = minLgs.size() - 1 - i;
        pl = getView(ea, minLgs[j], offset);
        addEdg(origNds[j], centerNds[j], pl);
      }
      offset += pl.getLines().size();
    }

    // delete remaining stuff
    delNd(nb);
    delNd(na);

    // update orderings
    for (auto n : origNds) updateEdgeOrder(n);  // TODO: is this redundant?
    for (auto n : centerNds) {
      updateEdgeOrder(n);
      for (auto e : n->getAdjList()) updateEdgeOrder(e->getOtherNd(n));
    }
  }
}

// _____________________________________________________________________________
OptEdgePL OptGraph::getView(OptEdge* parent, OptEdge* leg, size_t offset) {
  OptEdgePL ret(parent->pl());
  ret.depth++;

  for (auto& lnEdgPart : ret.lnEdgParts) {
    if (parent->pl().lnEdgParts[0].dir ^ lnEdgPart.dir)
      lnEdgPart.order +=
          parent->pl().getLines().size() - leg->pl().getLines().size() - offset;
    else
      lnEdgPart.order += offset;
  }

  ret.lines.clear();

  for (auto ro : leg->pl().getLines()) {
    auto lo = getCtdLineIn(ro.line, ro.dir, leg, parent);
    assert(sharedNode(leg, parent));
    if (lo) ret.lines.push_back(*lo);
  }

  std::sort(ret.lines.begin(), ret.lines.end());

  return ret;
}

// _____________________________________________________________________________
OptEdgePL OptGraph::getPartialViewExcl(OptEdge* parent, const OptLO* leg,
                                       size_t offset) {
  OptEdgePL ret(parent->pl());
  ret.depth++;

  for (auto& lnEdgPart : ret.lnEdgParts) {
    if (parent->pl().lnEdgParts[0].dir ^ lnEdgPart.dir) {
      lnEdgPart.order += 1 - offset;
    } else {
      lnEdgPart.order += offset;
    }
  }

  ret.lines.clear();

  for (auto ro : parent->pl().getLines()) {
    if (ro.line != leg->line) ret.lines.push_back(ro);
  }

  std::sort(ret.lines.begin(), ret.lines.end());

  return ret;
}

// _____________________________________________________________________________
OptEdgePL OptGraph::getPartialViewExcl(OptEdge* parent, OptEdge* leg,
                                       size_t offset) {
  OptEdgePL ret(parent->pl());
  ret.depth++;

  size_t shared = getCtdLinesIn(parent, leg).size();

  for (auto& lnEdgPart : ret.lnEdgParts) {
    if (parent->pl().lnEdgParts[0].dir ^ lnEdgPart.dir) {
      lnEdgPart.order += parent->pl().getLines().size() -
                         (parent->pl().getLines().size() - shared) - offset;
    } else {
      lnEdgPart.order += offset;
    }
  }

  ret.lines.clear();

  for (auto ro : parent->pl().getLines()) {
    auto lo = getCtdLineIn(ro.line, ro.dir, parent, leg);
    if (!lo) ret.lines.push_back(ro);
  }

  std::sort(ret.lines.begin(), ret.lines.end());

  return ret;
}

// _____________________________________________________________________________
OptEdgePL OptGraph::getPartialView(OptEdge* parent, const OptLO* leg,
                                   size_t offset) {
  OptEdgePL ret(parent->pl());
  ret.depth++;

  for (auto& lnEdgPart : ret.lnEdgParts) {
    if (parent->pl().lnEdgParts[0].dir ^ lnEdgPart.dir)
      lnEdgPart.order += parent->pl().getLines().size() - 1 - offset;
    else
      lnEdgPart.order += offset;
  }

  ret.lines.clear();

  ret.lines.push_back(*leg);

  std::sort(ret.lines.begin(), ret.lines.end());

  return ret;
}

// _____________________________________________________________________________
OptEdgePL OptGraph::getPartialView(OptEdge* parent, OptEdge* leg,
                                   size_t offset) {
  OptEdgePL ret(parent->pl());
  ret.depth++;

  size_t shared = getCtdLinesIn(parent, leg).size();

  for (auto& lnEdgPart : ret.lnEdgParts) {
    if (parent->pl().lnEdgParts[0].dir ^ lnEdgPart.dir)
      lnEdgPart.order += parent->pl().getLines().size() - shared - offset;
    else
      lnEdgPart.order += offset;
  }

  ret.lines.clear();

  for (auto ro : leg->pl().getLines()) {
    auto lo = getCtdLineIn(ro.line, ro.dir, leg, parent);
    if (lo) ret.lines.push_back(*lo);
  }

  std::sort(ret.lines.begin(), ret.lines.end());

  return ret;
}

// _____________________________________________________________________________
void OptGraph::untanglePartialDogBone() {
  std::vector<OptEdge*> toUntangle;

  for (OptNode* na : getNds()) {
    if (na->getDeg() < 3) continue;  // only look at nodes with deg > 2

    for (OptEdge* mainLeg : na->getAdjList()) {
      if (mainLeg->getFrom() != na) continue;
      OptNode* notPartN = 0;
      if ((notPartN = isPartialDogBone(mainLeg))) {
        LOGTO(DEBUG, std::cerr)
            << "Found partial dog bone with main leg " << mainLeg << " ("
            << mainLeg->pl().toStr() << ") at node " << notPartN;

        toUntangle.push_back(mainLeg);
      }
    }
  }

  for (auto mainLeg : toUntangle) {
    OptNode* notPartN = isPartialDogBone(mainLeg);
    OptNode* partN = mainLeg->getOtherNd(notPartN);

    // the geometry of the main leg
    util::geo::PolyLine<double> pl(*notPartN->pl().getGeom(),
                                   *partN->pl().getGeom());
    double bandW = (partN->getDeg() - 1) * (DO / (mainLeg->pl().depth + 1));
    auto orthoPlOrig = pl.getOrthoLineAtDist(pl.getLength(), bandW);
    orthoPlOrig.reverse();

    // each leg, in clockwise fashion
    auto minLgs = clockwEdges(mainLeg, partN);
    auto minLgsNotPart = partialClockwEdges(mainLeg, notPartN);
    std::reverse(minLgsNotPart.begin(), minLgsNotPart.end());
    std::vector<OptEdge*> minLgsNew(minLgs.size());

    assert(minLgs.size() <= mainLeg->pl().getCardinality());
    // map positions in b to positions in a
    std::vector<size_t> aToB = mapPositions(minLgsNotPart, mainLeg, minLgs);
    std::vector<size_t> iden(aToB.size());
    for (size_t i = 0; i < iden.size(); i++) iden[i] = i;

    // for each minor leg of the Y, create a new node at the origin
    std::vector<OptNode*> partNds =
        explodeNodeAlong(partN, orthoPlOrig, minLgs.size());

    for (size_t i = 0; i < minLgs.size(); i++) {
      if (minLgs[i]->getFrom() == partN)
        minLgsNew[i] = addEdg(partNds[i], minLgs[i]->getTo(), minLgs[i]->pl());
      else
        minLgsNew[i] =
            addEdg(minLgs[i]->getFrom(), partNds[i], minLgs[i]->pl());
    }

    std::vector<OptEdge*>& refLegs = minLgs;
    std::vector<size_t>& toB = iden;
    if (_scorer->getCrossingPenDiffSeg(partN) <
        _scorer->getCrossingPenDiffSeg(notPartN)) {
      refLegs = minLgsNotPart;
      toB = aToB;
    }

    size_t offset = 0;
    for (size_t i = 0; i < minLgs.size(); i++) {
      size_t j = i;
      OptEdgePL pl;
      if (mainLeg->getFrom() == partN) {
        if (mainLeg->pl().lnEdgParts[0].dir) j = minLgs.size() - 1 - i;
        pl = getPartialView(mainLeg, refLegs[j], offset);
        addEdg(partNds[toB[j]], notPartN, pl);
      } else {
        if (!mainLeg->pl().lnEdgParts[0].dir) j = minLgs.size() - 1 - i;
        pl = getPartialView(mainLeg, refLegs[j], offset);
        addEdg(notPartN, partNds[toB[j]], pl);
      }
      offset += pl.getLines().size();
    }

    // delete remaining stuff
    delNd(partN);

    // update orderings
    for (auto n : partNds) {
      updateEdgeOrder(n);
      for (auto e : n->getAdjList()) updateEdgeOrder(e->getOtherNd(n));
    }
    updateEdgeOrder(notPartN);
  }
}

// _____________________________________________________________________________
void OptGraph::untangleInnerStump() {
  std::vector<OptEdge*> toUntangle;

  for (OptNode* na : getNds()) {
    for (OptEdge* mainLeg : na->getAdjList()) {
      if (mainLeg->getFrom() != na) continue;
      if (isInnerStump(mainLeg)) {
        LOGTO(DEBUG, std::cerr) << "Found inner stump with main leg " << mainLeg
                                << " (" << mainLeg->pl().toStr() << ")";

        toUntangle.push_back(mainLeg);
      }
    }
  }

  for (auto mainLeg : toUntangle) {
    auto na = mainLeg->getFrom();
    OptNode* nb = mainLeg->getOtherNd(na);

    std::vector<OptNode*> dummies;

    // make relative to main leg
    na->pl().circOrdering = clockwEdges(mainLeg, na);
    na->pl().circOrdering.insert(na->pl().circOrdering.begin(), mainLeg);
    nb->pl().circOrdering = clockwEdges(mainLeg, nb);
    nb->pl().circOrdering.insert(nb->pl().circOrdering.begin(), mainLeg);

    // add dummy edges
    size_t i = 0;
    for (auto a : clockwEdges(mainLeg, na)) {
      if (terminatesAt(a, mainLeg, nb)) {
        auto dummyNode = addNd(nb->pl());
        dummies.push_back(dummyNode);
        auto dummyEdge = addEdg(nb, dummyNode);
        for (auto roOld : a->pl().getLines()) {
          dummyEdge->pl().lines.push_back(OptLO(roOld.line, 0));
        }
        dummyEdge->pl().lnEdgParts.push_back({0, 0, 0, 0});

        nb->pl().circOrdering.insert(nb->pl().circOrdering.end() - i,
                                     dummyEdge);
      }

      i++;
      // only necessary on this side, not the other side
      while (i < nb->pl().circOrdering.size() &&
             terminatesAt(nb->pl().circOrdering[i], mainLeg, na)) {
        i++;
      }
    }

    i = 0;
    for (auto a : clockwEdges(mainLeg, nb)) {
      if (terminatesAt(a, mainLeg, na)) {
        auto dummyNode = addNd(na->pl());
        dummies.push_back(dummyNode);
        auto dummyEdge = addEdg(na, dummyNode);
        for (auto roOld : a->pl().getLines()) {
          dummyEdge->pl().lines.push_back(OptLO(roOld.line, 0));
        }
        dummyEdge->pl().lnEdgParts.push_back({0, 0, 0, 0});
        na->pl().circOrdering.insert(na->pl().circOrdering.end() - i,
                                     dummyEdge);
      }
      i++;
    }

    // the geometry of the main leg
    util::geo::PolyLine<double> pl(*na->pl().getGeom(), *nb->pl().getGeom());
    double bandW = (nb->getDeg() - 1) * (DO / (mainLeg->pl().depth + 1));
    auto orthoPlA = pl.getOrthoLineAtDist(0, bandW);
    auto orthoPlB = pl.getOrthoLineAtDist(pl.getLength(), bandW);

    // for each minor leg at node a, create a new node
    std::vector<OptNode*> aNds =
        explodeNodeAlong(na, orthoPlA, na->getDeg() - 1);

    // each leg in a, in clockwise fashion
    auto minLgsA = clockwEdges(mainLeg, na);
    std::vector<OptEdge*> minLgsANew(minLgsA.size());
    assert(minLgsA.size() == na->pl().circOrdering.size() - 1);

    std::vector<OptNode*> bNds =
        explodeNodeAlong(nb, orthoPlB, nb->getDeg() - 1);

    // each leg in b, in counter-clockwise fashion
    auto minLgsB = clockwEdges(mainLeg, nb);
    std::vector<OptEdge*> minLgsBNew(minLgsB.size());
    assert(minLgsB.size() == nb->pl().circOrdering.size() - 1);

    std::reverse(minLgsB.begin(), minLgsB.end());

    // map positions in b to positions in a
    std::vector<size_t> bToA = mapPositions(minLgsA, mainLeg, minLgsB);
    std::vector<size_t> aToB = mapPositions(minLgsB, mainLeg, minLgsA);
    std::vector<size_t> iden(bToA.size());
    for (size_t i = 0; i < iden.size(); i++) iden[i] = i;

    for (size_t i = 0; i < minLgsA.size(); i++) {
      if (minLgsA[i]->getFrom() == na) {
        minLgsANew[i] = addEdg(aNds[i], minLgsA[i]->getTo(), minLgsA[i]->pl());
      } else {
        minLgsANew[i] =
            addEdg(minLgsA[i]->getFrom(), aNds[i], minLgsA[i]->pl());
      }
    }

    for (size_t i = 0; i < minLgsB.size(); i++) {
      if (minLgsB[i]->getFrom() == nb)
        minLgsBNew[i] = addEdg(bNds[i], minLgsB[i]->getTo(), minLgsB[i]->pl());
      else
        minLgsBNew[i] =
            addEdg(minLgsB[i]->getFrom(), bNds[i], minLgsB[i]->pl());
    }

    std::vector<OptEdge*>& refLegs = minLgsA;
    std::vector<size_t>& toA = bToA;
    std::vector<size_t>& toB = iden;

    if (dogBoneCheaper(na, nb, mainLeg->pl().getLines())) {
      refLegs = minLgsB;
      toA = iden;
      toB = aToB;
    }

    size_t offset = 0;
    for (size_t i = 0; i < minLgsA.size(); i++) {
      size_t j = i;
      OptEdgePL pl;
      if (mainLeg->getFrom() == na) {
        if (mainLeg->pl().lnEdgParts[0].dir) j = minLgsA.size() - 1 - i;
        pl = getView(mainLeg, refLegs[j], offset);
        addEdg(aNds[toB[j]], bNds[toA[j]], pl);
      } else {
        if (!mainLeg->pl().lnEdgParts[0].dir) j = minLgsA.size() - 1 - i;
        pl = getView(mainLeg, refLegs[j], offset);
        addEdg(bNds[toA[j]], aNds[toB[j]], pl);
      }
      offset += pl.getLines().size();
    }

    // delete remaining stuff
    delNd(nb);
    delNd(na);

    // delete dummy nodes
    for (auto nd : dummies) delNd(nd);

    // update orderings
    for (auto n : aNds) {
      updateEdgeOrder(n);
      for (auto e : n->getAdjList()) updateEdgeOrder(e->getOtherNd(n));
    }
    for (auto n : bNds) {
      updateEdgeOrder(n);
      for (auto e : n->getAdjList()) updateEdgeOrder(e->getOtherNd(n));
    }
  }
}

// _____________________________________________________________________________
void OptGraph::untangleDogBone() {
  std::vector<OptEdge*> toUntangle;

  for (OptNode* na : getNds()) {
    for (OptEdge* mainLeg : na->getAdjList()) {
      if (mainLeg->getFrom() != na) continue;
      if (isDogBone(mainLeg)) {
        LOGTO(DEBUG, std::cerr)
            << "Found full dog bone with main leg " << mainLeg << " ("
            << mainLeg->pl().toStr() << ")";

        toUntangle.push_back(mainLeg);
      }
    }
  }

  for (auto mainLeg : toUntangle) {
    auto na = mainLeg->getFrom();
    OptNode* nb = mainLeg->getOtherNd(na);
    // the geometry of the main leg
    util::geo::PolyLine<double> pl(*na->pl().getGeom(), *nb->pl().getGeom());
    double bandW = (nb->getDeg() - 1) * (DO / (mainLeg->pl().depth + 1));
    auto orthoPlA = pl.getOrthoLineAtDist(0, bandW);
    auto orthoPlB = pl.getOrthoLineAtDist(pl.getLength(), bandW);

    // for each minor leg at node a, create a new node
    std::vector<OptNode*> aNds =
        explodeNodeAlong(na, orthoPlA, na->getDeg() - 1);

    // each leg in a, in clockwise fashion
    auto minLgsA = clockwEdges(mainLeg, na);
    std::vector<OptEdge*> minLgsANew(minLgsA.size());
    assert(minLgsA.size() == na->pl().circOrdering.size() - 1);

    std::vector<OptNode*> bNds =
        explodeNodeAlong(nb, orthoPlB, nb->getDeg() - 1);

    // each leg in b, in counter-clockwise fashion
    auto minLgsB = clockwEdges(mainLeg, nb);
    std::vector<OptEdge*> minLgsBNew(minLgsB.size());
    assert(minLgsB.size() == nb->pl().circOrdering.size() - 1);

    std::reverse(minLgsB.begin(), minLgsB.end());

    // map positions in b to positions in a
    std::vector<size_t> bToA = mapPositions(minLgsA, mainLeg, minLgsB);
    std::vector<size_t> aToB = mapPositions(minLgsB, mainLeg, minLgsA);
    std::vector<size_t> iden(bToA.size());
    for (size_t i = 0; i < iden.size(); i++) iden[i] = i;

    for (size_t i = 0; i < minLgsA.size(); i++) {
      if (minLgsA[i]->getFrom() == na) {
        minLgsANew[i] = addEdg(aNds[i], minLgsA[i]->getTo(), minLgsA[i]->pl());
      } else {
        minLgsANew[i] =
            addEdg(minLgsA[i]->getFrom(), aNds[i], minLgsA[i]->pl());
      }
    }

    for (size_t i = 0; i < minLgsB.size(); i++) {
      if (minLgsB[i]->getFrom() == nb)
        minLgsBNew[i] = addEdg(bNds[i], minLgsB[i]->getTo(), minLgsB[i]->pl());
      else
        minLgsBNew[i] =
            addEdg(minLgsB[i]->getFrom(), bNds[i], minLgsB[i]->pl());
    }

    std::vector<OptEdge*>& refLegs = minLgsA;
    std::vector<size_t>& toA = bToA;
    std::vector<size_t>& toB = iden;

    if (dogBoneCheaper(na, nb, mainLeg->pl().getLines())) {
      refLegs = minLgsB;
      toA = iden;
      toB = aToB;
    }

    size_t offset = 0;
    for (size_t i = 0; i < minLgsA.size(); i++) {
      size_t j = i;
      OptEdgePL pl;
      if (mainLeg->getFrom() == na) {
        if (mainLeg->pl().lnEdgParts[0].dir) j = minLgsA.size() - 1 - i;
        pl = getView(mainLeg, refLegs[j], offset);
        addEdg(aNds[toB[j]], bNds[toA[j]], pl);
      } else {
        if (!mainLeg->pl().lnEdgParts[0].dir) j = minLgsA.size() - 1 - i;
        pl = getView(mainLeg, refLegs[j], offset);
        addEdg(bNds[toA[j]], aNds[toB[j]], pl);
      }
      offset += pl.getLines().size();
    }

    // delete remaining stuff
    delNd(nb);
    delNd(na);

    // update orderings
    for (auto n : aNds) {
      updateEdgeOrder(n);
      for (auto e : n->getAdjList()) updateEdgeOrder(e->getOtherNd(n));
    }
    for (auto n : bNds) {
      updateEdgeOrder(n);
      for (auto e : n->getAdjList()) updateEdgeOrder(e->getOtherNd(n));
    }
  }
}

// _____________________________________________________________________________
void OptGraph::writeEdgeOrder() {
  for (auto nd : getNds()) {
    updateEdgeOrder(nd);
  }
}

// _____________________________________________________________________________
size_t OptNodePL::circOrder(OptEdge* e) const {
  if (circOrdering.size() == 1) return 0;
  if (circOrdering.size() < 4) {
    return std::find(circOrdering.begin(), circOrdering.end(), e) -
           circOrdering.begin();
  } else {
    return circOrderMap.find(e)->second;
  }

  // TODO: throw exception
  assert(false);
}

// _____________________________________________________________________________
void OptGraph::updateEdgeOrder(OptNode* n) {
  n->pl().circOrdering.clear();

  if (n->getDeg() == 1) {
    n->pl().circOrdering.push_back(n->getAdjList().front());
    return;
  }

  for (auto e : n->getAdjList()) {
    n->pl().circOrdering.push_back(e);
  }
  std::sort(n->pl().circOrdering.begin(), n->pl().circOrdering.end(), cmpEdge);

  for (size_t i = 0; i < n->pl().circOrdering.size(); i++) {
    n->pl().circOrderMap[n->pl().circOrdering[i]] = i;
  }
}

// _____________________________________________________________________________
bool OptGraph::dirLineContains(const OptEdge* a, const OptEdge* b) {
  for (auto& to : b->pl().getLines()) {
    if (!getCtdLineIn(to.line, to.dir, b, a)) return false;
  }
  return true;
}

// _____________________________________________________________________________
bool OptGraph::dirLineEqualIn(const OptEdge* a, const OptEdge* b) {
  if (a->pl().getCardinality() != b->pl().getCardinality()) {
    return false;
  }

  for (auto& to : a->pl().getLines()) {
    if (!getSameDirLineIn(to.line, to.dir, a, b)) return false;
  }

  for (auto& to : b->pl().getLines()) {
    if (!getSameDirLineIn(to.line, to.dir, b, a)) return false;
  }

  return true;
}

// _____________________________________________________________________________
bool OptGraph::isYAt(OptEdge* eLeg, OptNode* n) const {
  if (eLeg->getFrom() == eLeg->getTo()) return 0;
  if (eLeg->getOtherNd(n)->getDeg() != 1) return false;
  return branchesAt(eLeg, n).size();
}
// _____________________________________________________________________________
bool OptGraph::isPartialYAt(OptEdge* eLeg, OptNode* n) const {
  if (eLeg->getFrom() == eLeg->getTo()) return 0;
  if (eLeg->getOtherNd(n)->getDeg() != 1) return false;

  return partiallyBranchesAt(eLeg, n).size();
}

// _____________________________________________________________________________
bool OptGraph::isInnerStump(OptEdge* leg) const {
  if (leg->getFrom() == leg->getTo()) return 0;
  auto branches = branchesAt(leg, leg->getFrom(), true);
  if (!branches.size()) return false;

  // check whether all leg lines may cross cheaper or equal at one of the nodes
  if (!(dogBoneCheaper(leg->getFrom(), leg->getTo(), leg->pl().getLines()) ||
        dogBoneCheaper(leg->getTo(), leg->getFrom(), leg->pl().getLines())))
    return false;

  if (!branchesAtIntoStump(leg, leg->getTo(), branches)) return false;

  auto minLgsA = clockwEdges(leg, leg->getFrom());
  auto minLgsB = clockwEdges(leg, leg->getTo());

  bool hasTerm = false;

  size_t b = minLgsB.size() - 1;
  for (size_t a = 0; a < minLgsA.size(); a++) {
    if (minLgsA[a] == leg) continue;
    if (terminatesAt(minLgsA[a], leg, leg->getTo())) continue;
    for (size_t bb = b + 1; bb > 0; bb--) {
      if (minLgsB[bb - 1] == leg) {
        b--;
        continue;
      }
      if (terminatesAt(minLgsB[bb - 1], leg, leg->getFrom())) {
        hasTerm = true;
        b--;
        continue;
      } else {
        break;
      }
    }
    if (dirContinuedOver(minLgsA[a], leg, minLgsB[b])) {
      b--;
    } else {
      return false;
    }
  }

  size_t a = minLgsA.size() - 1;
  for (size_t b = 0; b < minLgsB.size(); b++) {
    if (minLgsB[b] == leg) continue;
    if (terminatesAt(minLgsB[b], leg, leg->getFrom())) continue;
    for (size_t aa = a + 1; aa > 0; aa--) {
      if (minLgsA[aa - 1] == leg) {
        a--;
        continue;
      }
      if (terminatesAt(minLgsA[aa - 1], leg, leg->getTo())) {
        hasTerm = true;
        a--;
        continue;
      } else {
        break;
      }
    }
    if (dirContinuedOver(minLgsB[b], leg, minLgsA[a])) {
      a--;
    } else {
      return false;
    }
  }

  return hasTerm;
}

// _____________________________________________________________________________
bool OptGraph::isDogBone(OptEdge* leg) const {
  if (leg->getFrom() == leg->getTo()) return false;
  if (leg->getFrom()->getDeg() != leg->getTo()->getDeg()) return false;

  auto branches = branchesAt(leg, leg->getFrom());
  if (!branches.size()) return false;

  // check whether all leg lines may cross cheaper or equal at one of the nodes
  if (!(dogBoneCheaper(leg->getFrom(), leg->getTo(), leg->pl().getLines()) ||
        dogBoneCheaper(leg->getTo(), leg->getFrom(), leg->pl().getLines())))
    return false;

  return branchesAtInto(leg, leg->getTo(), branches);
}

// _____________________________________________________________________________
OptNode* OptGraph::isPartialDogBone(OptEdge* leg) const {
  if (leg->getFrom() == leg->getTo()) return 0;
  auto branchesA = branchesAt(leg, leg->getFrom());

  // check whether all leg lines may cross cheaper or equal at one of the nodes
  if (!(dogBoneCheaper(leg->getFrom(), leg->getTo(), leg->pl().getLines()) ||
        dogBoneCheaper(leg->getTo(), leg->getFrom(), leg->pl().getLines())))
    return 0;

  if (branchesA.size() && partiallyBranchesAtInto(leg, leg->getTo(), branchesA))
    return leg->getTo();

  branchesA = branchesAt(leg, leg->getTo());
  if (branchesA.size() &&
      partiallyBranchesAtInto(leg, leg->getFrom(), branchesA))
    return leg->getFrom();

  return 0;
}

// _____________________________________________________________________________
bool OptGraph::partiallyBranchesAtInto(OptEdge* leg, OptNode* n,
                                       std::vector<OptEdge*> branchesA) const {
  auto branchesB = partiallyBranchesAt(leg, n);
  if (branchesA.size() != branchesB.size()) return false;

  std::set<OptEdge*> used;

  for (auto branch : branchesA) {
    bool found = false;
    for (auto oBranch : branchesB) {
      if (!used.count(oBranch) && dirContinuedOver(branch, leg, oBranch)) {
        found = true;
        used.insert(oBranch);
        break;
      }
    }
    if (!found) return false;
  }

  return true;
}

// _____________________________________________________________________________
bool OptGraph::branchesAtIntoStump(OptEdge* leg, OptNode* n,
                                   std::vector<OptEdge*> branchesA) const {
  auto branchesB = branchesAt(leg, n, true);
  if (branchesB.size() == 0) return false;

  assert(branchesB.size() > 0);
  assert(branchesA.size() > 0);

  size_t c = 0;

  for (auto branch : branchesA) {
    bool found = false;
    if (terminatesAt(branch, leg, sharedNode(leg, branchesB[0]))) {
      c += branch->pl().getLines().size();
      found = true;
      continue;
    }
    for (auto oBranch : branchesB) {
      if (dirContinuedOver(branch, leg, oBranch)) {
        c += branch->pl().getLines().size();
        found = true;
        break;
      }
    }
    if (!found) return false;
  }

  for (auto branch : branchesB) {
    bool found = false;
    if (terminatesAt(branch, leg, sharedNode(leg, branchesA[0]))) {
      c += branch->pl().getLines().size();
      found = true;
      continue;
    }
    for (auto oBranch : branchesA) {
      if (dirContinuedOver(branch, leg, oBranch)) {
        found = true;
        break;
      }
    }
    if (!found) return false;
  }

  return c == leg->pl().getLines().size();
}

// _____________________________________________________________________________
bool OptGraph::branchesAtInto(OptEdge* leg, OptNode* n,
                              std::vector<OptEdge*> branchesA) const {
  auto branchesB = branchesAt(leg, n);
  if (branchesA.size() != branchesB.size()) return false;

  for (auto branch : branchesA) {
    bool found = false;
    for (auto oBranch : branchesB) {
      if (dirContinuedOver(branch, leg, oBranch)) {
        found = true;
        break;
      }
    }
    if (!found) return false;
  }

  for (auto branch : branchesB) {
    bool found = false;
    for (auto oBranch : branchesA) {
      if (dirContinuedOver(branch, leg, oBranch)) {
        found = true;
        break;
      }
    }
    if (!found) return false;
  }

  return true;
}

// _____________________________________________________________________________
std::vector<OptEdge*> OptGraph::branchesAt(OptEdge* e, OptNode* n) const {
  return branchesAt(e, n, false);
}

// _____________________________________________________________________________
std::vector<OptEdge*> OptGraph::branchesAt(OptEdge* e, OptNode* n,
                                           bool allowStumps) const {
  std::vector<OptEdge*> ret;
  if (e->getFrom() != n && e->getTo() != n) return ret;
  if (n->getDeg() < 3) return ret;
  if (e->pl().getCardinality() < 2) return ret;

  size_t c = 0;
  std::vector<const OptEdge*> branchEdgs;

  for (OptEdge* ea : n->getAdjList()) {
    if (ea == e) continue;
    branchEdgs.push_back(ea);
    if (!dirLineContains(e, ea)) return {};

    c += ea->pl().getCardinality();
    ret.push_back(ea);
  }

  if (!allowStumps && c != e->pl().getCardinality()) return {};

  if (!lineDisjunct(branchEdgs)) return {};

  return ret;
}

// _____________________________________________________________________________
std::vector<OptEdge*> OptGraph::partiallyBranchesAt(OptEdge* eMain,
                                                    OptNode* n) const {
  std::set<OptEdge*> ret;
  if (n->getDeg() < 3) return std::vector<OptEdge*>();
  if (eMain->pl().getCardinality() < 2) return std::vector<OptEdge*>();
  if (eMain->getFrom() != n && eMain->getTo() != n)
    return std::vector<OptEdge*>();

  // use this to filter out full Ys
  bool onlyPartial = false;
  bool branches = false;
  OptEdge* first = 0;

  for (auto ro : eMain->pl().getLines()) {
    bool found = false;

    for (OptEdge* e : n->getAdjList()) {
      if (e == eMain) continue;

      if (dirContinuedOver(ro, eMain, e)) {
        if (found) return std::vector<OptEdge*>();

        found = true;
        if (first && first != e) branches = true;
        if (!first) first = e;
        if (!onlyPartial && !dirLineContains(eMain, e)) onlyPartial = true;
        ret.insert(e);
      } else {
        onlyPartial = true;
      }
    }
    if (!found) return std::vector<OptEdge*>();
  }

  if (onlyPartial && branches)
    return std::vector<OptEdge*>(ret.begin(), ret.end());
  return std::vector<OptEdge*>();
}

// _____________________________________________________________________________
std::vector<OptEdge*> OptGraph::clockwEdges(const OptEdge* noon,
                                            const OptNode* n) {
  std::vector<OptEdge*> clockwise;

  size_t passed = n->getDeg();
  for (size_t i = 0; i < n->pl().circOrdering.size(); i++) {
    auto e = n->pl().circOrdering[i];
    if (e == noon) {
      passed = i;
      continue;
    }
    if (passed > n->getDeg() - 1)
      clockwise.push_back(e);
    else {
      clockwise.insert(clockwise.begin() + (i - 1 - passed), e);
    }
  }

  return clockwise;
}

// _____________________________________________________________________________
std::vector<OptEdge*> OptGraph::partialClockwEdges(const OptEdge* noon,
                                                   const OptNode* n) {
  std::vector<OptEdge*> clockwise = clockwEdges(noon, n);

  for (auto i = clockwise.begin(); i != clockwise.end(); i++) {
    auto e = *i;
    if (!dirPartialContinuedOver(noon, e)) {
      i = clockwise.erase(i) - 1;
      continue;
    }
  }

  return clockwise;
}

// _____________________________________________________________________________
bool OptGraph::dirContinuedOver(const OptEdge* a, const OptEdge* b,
                                const OptEdge* c) {
  OptNode* ab = sharedNode(a, b);
  OptNode* bc = sharedNode(b, c);
  if (!ab || !bc) return false;

  for (auto& to : a->pl().getLines()) {
    auto too = getCtdLineIn(to.line, to.dir, a, b);
    if (too && !getCtdLineIn(too->line, too->dir, b, c)) return false;
  }
  return true;
}

// _____________________________________________________________________________
bool OptGraph::dirPartialContinuedOver(const OptEdge* a, const OptEdge* b) {
  OptNode* ab = sharedNode(a, b);
  if (!ab) return false;

  for (auto& to : a->pl().getLines()) {
    if (getCtdLineIn(to.line, to.dir, a, b)) return true;
  }

  return false;
}

// _____________________________________________________________________________
bool OptGraph::dirContinuedOver(const OptLO& ro, const OptEdge* a,
                                const OptEdge* b) {
  OptNode* ab = sharedNode(a, b);
  if (!ab) return false;

  return getCtdLineIn(ro.line, ro.dir, a, b);
}

// _____________________________________________________________________________
bool OptGraph::dirContinuedOver(const OptLO& ro, const OptEdge* eA,
                                const OptNode* n) {
  for (auto eB : n->getAdjList()) {
    if (eA == eB) continue;
    if (dirContinuedOver(ro, eA, eB)) return true;
  }
  return false;
}

// _____________________________________________________________________________
bool OptGraph::linesBranchAt(const OptLO& roA, const OptLO& roB,
                             const OptEdge* a, const OptNode* n) {
  const OptEdge* roATo = 0;
  const OptEdge* roBTo = 0;
  for (auto e : n->getAdjList()) {
    if (e == a) continue;
    if (dirContinuedOver(roA, a, e)) {
      roATo = e;
      continue;
    }
    if (dirContinuedOver(roB, a, e)) {
      roBTo = e;
      continue;
    }
  }
  return roATo != 0 && roBTo != 0 && roATo != roBTo;
}

// _____________________________________________________________________________
bool OptGraph::linesCtnOver(const OptLO& roA, const OptLO& roB,
                            const OptEdge* a, const OptNode* n) {
  for (auto e : n->getAdjList()) {
    if (e == a) continue;
    if (dirContinuedOver(roA, a, e) && dirContinuedOver(roB, a, e)) {
      return true;
    }
  }
  return false;
}

// _____________________________________________________________________________
OptNode* OptGraph::sharedNode(const OptEdge* a, const OptEdge* b) {
  if (a->getFrom() == b->getFrom() || a->getFrom() == b->getTo())
    return a->getFrom();
  if (a->getTo() == b->getFrom() || a->getTo() == b->getTo()) return a->getTo();
  return 0;
}

// _____________________________________________________________________________
const OptLO* OptGraph::getSameDirLineIn(const Line* r, const LineNode* dir,
                                        const OptEdge* fromEdge,
                                        const OptEdge* toEdge) {
  const OptNode* n = sharedNode(fromEdge, toEdge);
  if (!n || !n->pl().node || n->getDeg() == 1) return 0;

  auto to = toEdge->pl().getLineOcc(r);
  if (!to) return 0;
  if ((to->dir == 0 && dir == 0) ||
      (to->dir == n->pl().node && dir != 0 && dir != n->pl().node) ||
      (to->dir != n->pl().node && to->dir != 0 && dir == n->pl().node)) {
    if (n->pl().node->pl().connOccurs(r, getAdjEdg(fromEdge, n),
                                      getAdjEdg(toEdge, n))) {
      return to;
    }
  }

  return 0;
}

// _____________________________________________________________________________
bool OptGraph::hasCtdLineIn(const Line* r, const LineNode* dir,
                            const OptEdge* fromEdge, const OptEdge* toEdge) {
  return getCtdLineIn(r, dir, fromEdge, toEdge);
}

// _____________________________________________________________________________
const OptLO* OptGraph::getCtdLineIn(const Line* r, const LineNode* dir,
                                    const OptEdge* fromEdge,
                                    const OptEdge* toEdge) {
  const OptNode* n = sharedNode(fromEdge, toEdge);
  if (!n || !n->pl().node || n->getDeg() == 1) return 0;

  auto to = toEdge->pl().getLineOcc(r);
  if (!to) return 0;
  if (to->dir == 0 || dir == 0 ||
      (to->dir == n->pl().node && dir != n->pl().node) ||
      (to->dir != n->pl().node && dir == n->pl().node)) {
    if (n->pl().node->pl().connOccurs(r, getAdjEdg(fromEdge, n),
                                      getAdjEdg(toEdge, n))) {
      return to;
    }
  }

  return 0;
}

// _____________________________________________________________________________
bool OptGraph::lineDisjunct(const std::vector<const OptEdge*>& edges) {
  std::set<const Line*> lines;

  for (const auto edg : edges) {
    for (const OptLO& to : edg->pl().getLines()) {
      if (lines.count(to.line)) return false;
      lines.insert(to.line);
    }
  }

  return true;
}

// _____________________________________________________________________________
std::vector<OptLO> OptGraph::getCtdLinesIn(const OptEdge* fromEdge,
                                           const OptEdge* toEdge) {
  std::vector<OptLO> ret;
  const OptNode* n = sharedNode(fromEdge, toEdge);
  if (!n || !n->pl().node) return ret;

  for (const OptLO& to : fromEdge->pl().getLines()) {
    auto r = getCtdLineIn(to.line, to.dir, fromEdge, toEdge);
    if (r) ret.push_back(*r);
  }

  return ret;
}

// _____________________________________________________________________________
std::vector<size_t> OptGraph::mapPositions(std::vector<OptEdge*> a,
                                           OptEdge* leg,
                                           std::vector<OptEdge*> b) const {
  std::vector<size_t> ret(a.size());
  assert(a.size() == b.size());

  for (size_t i = 0; i < a.size(); i++) {
    bool found = false;
    for (size_t j = 0; j < b.size(); j++) {
      if (dirContinuedOver(a[i], leg, b[j])) {
        ret[i] = j;
        found = true;
        break;
      }
    }
    assert(found);
  }

  return ret;
}

// _____________________________________________________________________________
bool OptGraph::terminatesAt(const OptEdge* from, const OptEdge* over,
                            const OptNode* nd) {
  for (const auto& lo : from->pl().getLines()) {
    auto ro = getLO(over, lo.line);
    if (ro.isNull()) return false;
    if (!terminatesAt(ro, over, nd)) return false;
  }
  return true;
}

// _____________________________________________________________________________
bool OptGraph::terminatesAt(const OptLO& lo, const OptEdge* eA,
                            const OptNode* nd) {
  for (auto eB : nd->getAdjList()) {
    if (eA == eB) continue;
    if (dirContinuedOver(lo, eA, eB)) return false;
  }
  return true;
}

// _____________________________________________________________________________
bool OptGraph::terminatesAt(const OptEdge* eA, const OptNode* nd) {
  for (const auto& lo : eA->pl().getLines()) {
    for (auto eB : nd->getAdjList()) {
      if (eA == eB) continue;
      if (dirContinuedOver(lo, eA, eB)) return false;
    }
  }
  return true;
}
