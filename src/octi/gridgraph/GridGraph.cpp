// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include "octi/gridgraph/GridGraph.h"
#include "octi/gridgraph/NodeCost.h"
#include "util/graph/Node.h"

using namespace octi::gridgraph;
using octi::gridgraph::NodeCost;
using octi::gridgraph::GridGraph;
using util::geo::DBox;
using util::geo::DPoint;
using util::geo::dist;

double INF = std::numeric_limits<double>::infinity();

// _____________________________________________________________________________
GridGraph::GridGraph(const DBox& bbox, double cellSize, double spacer,
                     const Penalties& pens)
    : _bbox(bbox),
      _c(pens),
      _grid(cellSize, cellSize, bbox),
      _cellSize(cellSize),
      _spacer(spacer) {
  assert(_c.p_0 < _c.p_135);
  assert(_c.p_135 < _c.p_90);
  assert(_c.p_90 < _c.p_45);

  // cut off illegal spacer values
  if (spacer > cellSize / 2) spacer = cellSize / 2;

  // write nodes
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      writeNd(x, y, 0, 0, 0);
    }
  }

  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      GridNode* center = getNode(x, y);
      if (!center) continue;

      for (size_t p = 0; p < 8; p++) {
        GridNode* from = center->pl().getPort(p);
        GridNode* toN = getNeighbor(x, y, p);
        if (from != 0 && toN != 0) {
          GridNode* to = toN->pl().getPort((p + 4) % 8);
          // if the edge already exists, this step is a no-op
          addEdg(from, to, GridEdgePL(0, false));
        }
      }
    }
  }

  writeInitialCosts();
}

// _____________________________________________________________________________
GridNode* GridGraph::getNode(size_t x, size_t y) const {
  if (x >= _grid.getXWidth() || y >= _grid.getYHeight()) return 0;
  std::set<GridNode*> r;
  _grid.get(x, y, &r);

  if (r.size()) return *r.begin();
  return 0;
}

// _____________________________________________________________________________
GridNode* GridGraph::getNeighbor(size_t cx, size_t cy, size_t i) const {
  int8_t x = 1;
  if (i % 4 == 0) x = 0;
  if (i > 4) x = -1;

  int8_t y = 1;
  if (i == 2 || i == 6) y = 0;
  if (i == 3 || i == 4 || i == 5) y = -1;

  return getNode(cx + x, cy + y);
}

// _____________________________________________________________________________
void GridGraph::balanceEdge(GridNode* a, GridNode* b) {
  if (a == b) return;

  // search the edge, TODO
  size_t dir = 0;
  for (; dir < 8; dir++) {
    if (getEdg(a->pl().getPort(dir), b->pl().getPort((dir + 4) % 8))) {
      break;
    }
  }

  size_t xa = a->pl().getX();
  size_t ya = a->pl().getY();

  size_t xb = b->pl().getX();
  size_t yb = b->pl().getY();

  for (size_t p = 0; p < 8; p++) {
    auto neighA = getNeighbor(xa, ya, p);
    for (auto e : getNEdges(a, neighA)) e->pl().setUsed();

    auto neighB = getNeighbor(xb, yb, p);
    for (auto e : getNEdges(b, neighB)) e->pl().setUsed();
  }

  if (dir == 1 || dir == 3 || dir == 5 || dir == 7) {
    auto na = getNeighbor(xa, ya, (dir + 7) % 8);
    auto nb = getNeighbor(xa, ya, (dir + 1) % 8);

    if (na && nb) {
      auto es = getNEdges(na, nb);

      // this closes edges b that cross this edge a - we _never_ want to
      // allow that
      for (auto e : es) e->pl().setCost(INF);
    }
  }
}

// _____________________________________________________________________________
std::vector<GridEdge*> GridGraph::getNEdges(GridNode* a, GridNode* b) {
  if (!a) return {};
  if (!b) return {};
  a = a->pl().getParent()->pl().getStepMother();
  b = b->pl().getParent()->pl().getStepMother();

  std::vector<GridEdge*> ret;

  // TODO: better implementation of this

  for (size_t i = 0; i < 8; i++) {
    for (auto pa : a->pl().getMetaPorts(i)) {
      for (auto pb : b->pl().getMetaPorts((i + 4) % 8)) {
        auto e = getEdg(pa, pb);
        if (e) ret.push_back(e);
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
void GridGraph::getSettledOutgoingEdges(GridNode* n,
                                        std::set<CombEdge*> outgoing[8]) {
  size_t x = n->pl().getX();
  size_t y = n->pl().getY();

  for (size_t i = 0; i < 8; i++) {
    auto neigh = getNeighbor(x, y, i);
    auto es = getNEdges(n, neigh);

    for (auto e : es) {
      if (!e->pl().getResEdges().size()) continue;
      outgoing[i].insert(*e->pl().getResEdges().begin());
    }
  }
}

// _____________________________________________________________________________
NodeCost GridGraph::spacingPenalty(GridNode* n, CombNode* origNode,
                                   CombEdge* e) {
  NodeCost addC;

  int origEdgeNumber = origNode->getAdjList().size();
  size_t optimDistance = (8 / origEdgeNumber) - 1;

  if (!origNode->pl().getEdgeOrdering().has(e)) {
    std::cerr << "Warning: tried to balance edge " << e << " in node "
              << origNode << ", but the edge does not appear there."
              << std::endl;
    return addC;
  }

  std::set<CombEdge*> outgoing[8];
  getSettledOutgoingEdges(n, outgoing);

  for (size_t i = 0; i < 8; i++) {
    if (!outgoing[i].size()) continue;

    // this is the number of edges that will occur between the currently checked
    // edge and the inserted edge, in clockwise and counter-clockwise dir
    int32_t dCw =
        origNode->pl().getEdgeOrdering().dist(*outgoing[i].begin(), e) - 1;
    int32_t dCCw =
        origNode->pl().getEdgeOrdering().dist(e, *outgoing[i].begin()) - 1;

    // dd and ddd are the optimal distances between outgoing[i] and e, based on
    // the total number
    // of edges in this node
    int dd = ((((dCw + 1) + dCw) % 8) * optimDistance) % 8;
    int ddd = (6 - dd) % 8;

    double pen = _c.p_45 * 2 - 1;

    for (int j = 1; dd != 0 && j <= dd + 1; j++) {
      if (addC[(i + j) % 8] < -1) continue;
      addC[(i + j) % 8] += pen * (1.0 - (j - 1.0) / (dd));
    }

    for (int j = 1; ddd != 0 && j <= ddd + 1; j++) {
      if (addC[(i + (8 - j)) % 8] < -1) continue;
      addC[(i + (8 - j)) % 8] += pen * (1.0 - (j - 1.0) / (ddd));
    }

    // addC[i] = pen;

    // for (int j = 1; j <= dCw; j++) {
      // addC[(i + j) % 8] = -1.0 * std::numeric_limits<double>::max();
    // }

    // for (int j = 1; j <= dCCw; j++) {
      // addC[(i + (8 - j)) % 8] = -1.0 * std::numeric_limits<double>::max();
    // }
  }

  return addC;
}

// _____________________________________________________________________________
NodeCost GridGraph::outDegDeviationPenalty(CombNode* origNode, CombEdge* e) {
  NodeCost ret;
  double degA = util::geo::angBetween(
      *origNode->pl().getParent()->pl().getGeom(),
      *e->getOtherNd(origNode)->pl().getParent()->pl().getGeom());

  int deg = -degA * (180.0 / M_PI);
  if (deg < 0) deg += 360;

  deg = (deg + 90) % 360;

  for (int i = 0; i < 8; i++) {
    double diff = std::min<int>(abs(deg - (45 * i)), 360 - abs(deg - (45 * i)));
    double multiplier = .1;
    ret[i] += multiplier * diff;
  }
  return ret;
}

// _____________________________________________________________________________
NodeCost GridGraph::addCostVector(GridNode* n, const NodeCost& addC) {
  NodeCost invCost;
  n = n->pl().getStepMother();
  size_t x = n->pl().getX();
  size_t y = n->pl().getY();

  for (size_t p = 0; p < 8; p++) {
    auto neigh = getNeighbor(x, y, p);
    if (!neigh) continue;

    auto es = getNEdges(n, neigh);

    for (auto e : es) {
      if (addC[p] < -1) {
        if (e->pl().closed()) {
          // already closed, so dont remove closedness in inv costs
          invCost[p] = 0;
        } else {
          // close edges
          e->pl().close();

          // close the other node to avoid "stealing" the edge
          // IMPORTANT: because we check if this edge is already closed
          // above, it is impossible for this node to be already closed -
          // then the edge would also be closed, too. So we can close this node
          // here without danger of re-opening an already closed node later on.
          closeNode(neigh);

          invCost[p] = addC[p];
        }
      } else {
        e->pl().setCost(e->pl().rawCost() + addC[p]);
        invCost[p] = addC[p];
      }
    }
  }
  return invCost;
}

// _____________________________________________________________________________
void GridGraph::removeCostVector(GridNode* n, const NodeCost& addC) {
  size_t x = n->pl().getX();
  size_t y = n->pl().getY();

  for (size_t p = 0; p < 8; p++) {
    auto neigh = getNeighbor(x, y, p);
    if (!neigh) continue;

    auto es = getNEdges(n, neigh);

    for (auto e : es) {
      if (addC[p] < -1) {
        e->pl().open();
        openNode(neigh);
      } else {
        e->pl().setCost(e->pl().rawCost() - addC[p]);
      }
    }
  }
}

// _____________________________________________________________________________
CombEdgeSet GridGraph::getResEdges(GridNode* n) const {
  CombEdgeSet ret;

  for (size_t i = 0; i < 8; i++) {
    auto port = n->pl().getPort(i);
    for (auto e : port->getAdjList()) {
      ret.insert(e->pl().getResEdges().begin(), e->pl().getResEdges().end());
    }
  }

  return ret;
}

// _____________________________________________________________________________
void GridGraph::writeInitialCosts() {
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      auto n = getNode(x, y);
      for (size_t i = 0; i < 8; i++) {
        auto port = n->pl().getPort(i);
        auto neigh = getNeighbor(x, y, i);
        if (!neigh || !port) continue;
        auto e = getEdg(port, neigh->pl().getPort((i + 4) % 8));

        if (i % 4 == 0) {
          e->pl().setCost(_c.verticalPen);
        }
        if ((i + 2) % 4 == 0) {
          e->pl().setCost(_c.horizontalPen);
        }
        if (i % 2) {
          e->pl().setCost(_c.diagonalPen);
        }
      }
    }
  }
}

// _____________________________________________________________________________
std::priority_queue<Candidate> GridGraph::getNearestCandidatesFor(
    const DPoint& p, double maxD) const {
  std::priority_queue<Candidate> ret;
  std::set<GridNode*> neigh;
  DBox b(DPoint(p.getX() - maxD, p.getY() - maxD),
         DPoint(p.getX() + maxD, p.getY() + maxD));
  _grid.get(b, &neigh);

  for (auto n : neigh) {
    if (n->pl().isClosed()) continue;
    double d = dist(*n->pl().getGeom(), p);
    if (d < maxD) ret.push(Candidate(n, d));
  }

  return ret;
}

// _____________________________________________________________________________
const Grid<GridNode*, Point, double>& GridGraph::getGrid() const {
  return _grid;
}

// _____________________________________________________________________________
double GridGraph::heurCost(int64_t xa, int64_t ya, int64_t xb,
                           int64_t yb) const {
  if (xa == xb && ya == yb) return 0;
  size_t minHops = std::max(labs(xb - xa), labs(yb - ya));

  size_t edgeCost =
      minHops *
      (std::min(_c.verticalPen, std::min(_c.horizontalPen, _c.diagonalPen)));
  size_t hopCost = (minHops - 1) * (_c.p_45 - _c.p_135);

  return edgeCost + hopCost;
}

// _____________________________________________________________________________
void GridGraph::openNode(GridNode* n) {
  if (!n->pl().isClosed()) return;
  size_t x = n->pl().getX();
  size_t y = n->pl().getY();

  for (size_t p = 0; p < 8; p++) {
    auto neigh = getNeighbor(x, y, p);
    if (!neigh || neigh->pl().isClosed()) continue;
    auto es = getNEdges(n, neigh);
    for (auto e : es) {
      e->pl().open();
    }
  }

  n->pl().setClosed(false);
}

// _____________________________________________________________________________
void GridGraph::closeNode(GridNode* n) {
  if (n->pl().isClosed()) return;
  size_t x = n->pl().getX();
  size_t y = n->pl().getY();

  for (size_t p = 0; p < 8; p++) {
    auto neigh = getNeighbor(x, y, p);
    if (!neigh) continue;
    auto es = getNEdges(n, neigh);
    for (auto e : es) {
      e->pl().close();
    }
  }

  n->pl().setClosed(true);
}

// _____________________________________________________________________________
void GridGraph::openNodeSink(GridNode* n, double cost) {
  for (size_t i = 0; i < 8; i++) {
    auto e = getEdg(n->pl().getPort(i), n);
    if (e) e->pl().setCost(cost);
  }
}

// _____________________________________________________________________________
void GridGraph::closeNodeSink(GridNode* n) {
  for (size_t i = 0; i < 8; i++) {
    auto e = getEdg(n->pl().getPort(i), n);
    if (e) e->pl().setCost(INF);
  }
}

// _____________________________________________________________________________
GridNode* GridGraph::getGridNodeFrom(CombNode* n, double maxDis) {
  if (!isSettled(n)) {
    auto cands = getNearestCandidatesFor(*n->pl().getGeom(), maxDis);

    while (!cands.empty()) {
      if (!cands.top().n->pl().isClosed() && cands.top().n->pl().getStepChilds().size() == 0) {
        return cands.top().n;
      }
      cands.pop();
    }
    return 0;
  }

  return _settled[n];
}

// _____________________________________________________________________________
std::set<GridNode*> GridGraph::getGridNodesTo(CombNode* n, GridNode* ex, double maxDis) {
  std::set<GridNode*> tos;
  if (!isSettled(n)) {
    auto cands = getNearestCandidatesFor(*n->pl().getGeom(), maxDis);

    while (!cands.empty()) {
      if (cands.top().n != ex && !cands.top().n->pl().isClosed() && cands.top().n->pl().getStepChilds().size() == 0) {
        tos.insert(cands.top().n);
      }
      cands.pop();
    }
  } else {
    tos.insert(_settled.find(n)->second);
  }

  return tos;
}

// _____________________________________________________________________________
void GridGraph::settleGridNode(GridNode* n, CombNode* cn) {
  if (_settled.find(cn) != _settled.end() && _settled.find(cn)->second == n) {
    assert(_settled[cn] == n);
    return;
  }
  for (auto kv :_settled) {
    // check for double settle
    assert(kv.second != n);
  }
  std::set<CombEdge*> outgoing[8];
  getSettledOutgoingEdges(n, outgoing);
  for (auto o : outgoing) assert(o.size() == 0);

  _settled[cn] = n;
}

// _____________________________________________________________________________
bool GridGraph::isSettled(CombNode* cn) {
  return _settled.find(cn) != _settled.end();
}

// _____________________________________________________________________________
bool GridGraph::isSettled(GridNode* gn) {
  for (auto kv :_settled) {
    // check for double settle
    if (kv.second == gn) return true;
  }
  return false;
}

// _____________________________________________________________________________
std::pair<GridNode*, std::pair<size_t, size_t>> GridGraph::splitNode(
    GridNode* nd, GridEdge* a, GridEdge* b) {
  // copy the node
  auto nn = writeNd(nd->pl().getX(), nd->pl().getY(), _spacer * 2.5,
                    _spacer * 2.5, nd);
  // node nn will have erroneous meta ports in it - but we don't care, as the
  // node will not turn up in the index

  bool in = false;

  std::pair<GridNode*, std::pair<size_t, size_t>> ret;
  ret.first = nn;

  for (size_t p = 0; p < 8; p++) {
    GridNode* fromN = nd->pl().getPort(p);

    // get neighbor only returns stepmother / parent nodes
    GridNode* toN = getNeighbor(nd->pl().getX(), nd->pl().getY(), p);
    if (!toN) continue;

    GridEdge* edg = getEdg(fromN, toN->pl().getPort((p + 4) % 8));

    if (!edg) continue;

    if ((edg == a || edg == b)) {
      if (edg == a) ret.second.first = p;
      if (edg == b) ret.second.second = p;
      if (in) break;
      if (!in) in = true;
    } else if (in) {
      addEdg(nn->pl().getPort(p), toN->pl().getPort((p + 4) % 8), edg->pl());
      delEdg(edg->getFrom(), edg->getTo());
    }
  }

  return ret;
}

// _____________________________________________________________________________
void GridGraph::splitAlong(const GrEdgList& path) {
  GridEdge* last = 0;
  std::vector<std::pair<GridNode*, std::pair<size_t, size_t>>> nodePortPairs;
  std::vector<GridEdge*> edges;

  for (auto e : path) {
    if (e->pl().isSecondary()) continue;
    edges.push_back(e);
    if (last) {
      auto nd = GridGraph::sharedParentNode(last, e);
      assert(nd);
      assert(!isSettled(nd->pl().getParent()->pl().getStepMother()));
      auto newNdWithPorts = splitNode(nd, last, e);
      nodePortPairs.push_back(newNdWithPorts);
    }
    last = e;
  }

  if (edges.size() < 2) return;

  auto firstPrt =
      nodePortPairs[0].first->pl().getPort(nodePortPairs[0].second.first);
  auto lastPrt = nodePortPairs.back().first->pl().getPort(
      nodePortPairs.back().second.second);

  addEdg(firstPrt,
         edges[0]->getOtherNd(GridGraph::sharedParentNode(edges[0], edges[1])
                                  ->pl()
                                  .getPort(nodePortPairs[0].second.first)),
         edges[0]->pl());

  // add meta port to stepmother
  nodePortPairs[0].first->pl().getStepMother()->pl().addMetaPort(
      nodePortPairs[0].second.first,
      nodePortPairs[0].first->pl().getPort(nodePortPairs[0].second.first));

  for (size_t i = 1; i < edges.size() - 1; i++) {
    assert(i < nodePortPairs.size());
    addEdg(nodePortPairs[i - 1].first->pl().getPort(
               nodePortPairs[i - 1].second.second),
           nodePortPairs[i].first->pl().getPort(nodePortPairs[i].second.first),
           edges[i]->pl());

    nodePortPairs[i - 1].first->pl().getStepMother()->pl().addMetaPort(
        nodePortPairs[i - 1].second.second,
        nodePortPairs[i - 1].first->pl().getPort(
            nodePortPairs[i - 1].second.second));
    nodePortPairs[i].first->pl().getStepMother()->pl().addMetaPort(
        nodePortPairs[i].second.first,
        nodePortPairs[i].first->pl().getPort(nodePortPairs[i].second.first));
  }

  addEdg(lastPrt,
         edges.back()->getOtherNd(
             GridGraph::sharedParentNode(edges[edges.size() - 2], edges.back())
                 ->pl()
                 .getPort(nodePortPairs.back().second.second)),
         edges.back()->pl());

  nodePortPairs.back().first->pl().getStepMother()->pl().addMetaPort(
      nodePortPairs.back().second.second,
      nodePortPairs.back().first->pl().getPort(
          nodePortPairs.back().second.second));
}

// _____________________________________________________________________________
GridNode* GridGraph::writeNd(size_t x, size_t y, double xOff, double yOff,
                             GridNode* stepMother) {
  double xPos = _bbox.getLowerLeft().getX() + x * _cellSize;
  double yPos = _bbox.getLowerLeft().getY() + y * _cellSize;

  double c_0 = _c.p_45 - _c.p_135;
  double c_135 = _c.p_45;
  double c_90 = _c.p_45 - _c.p_135 + _c.p_90;

  GridNode* n = addNd(DPoint(xPos + xOff, yPos + yOff));

  if (!stepMother) _grid.add(x, y, n);

  n->pl().setXY(x, y);
  n->pl().setParent(n);

  if (stepMother) {
    stepMother->pl().addStepChild(n);
    n->pl().setStepMother(stepMother);
  } else {
    n->pl().setStepMother(n);
  }

  for (int i = 0; i < 8; i++) {
    int xi = ((4 - (i % 8)) % 4);
    xi /= abs(abs(xi) - 1) + 1;
    int yi = ((4 - ((i + 2) % 8)) % 4);
    yi /= abs(abs(yi) - 1) + 1;
    GridNode* nn =
        addNd(DPoint(xOff + xPos + xi * _spacer, yOff + yPos + yi * _spacer));
    nn->pl().setParent(n);
    n->pl().setPort(i, nn);
    if (!stepMother) n->pl().addMetaPort(i, nn);
    addEdg(n, nn, GridEdgePL(INF, true, false));
  }

  // in-node connections
  for (size_t i = 0; i < 8; i++) {
    for (size_t j = i + 1; j < 8; j++) {
      int d = (int)(i) - (int)(j);
      size_t deg = abs((((d + 4) % 8) + 8) % 8 - 4);
      double pen = c_0;

      if (deg == 1) continue;
      if (deg == 2) pen = c_90;
      if (deg == 3) pen = c_135;
      addEdg(n->pl().getPort(i), n->pl().getPort(j), GridEdgePL(pen, true));
    }
  }

  return n;
}

// _____________________________________________________________________________
GridNode* GridGraph::sharedParentNode(const GridEdge* a, const GridEdge* b) {
  GridNode* r = 0;
  if (a->getFrom()->pl().getParent() == b->getFrom()->pl().getParent() ||
      a->getFrom()->pl().getParent() == b->getTo()->pl().getParent())
    r = a->getFrom()->pl().getParent();
  if (a->getTo()->pl().getParent() == b->getFrom()->pl().getParent() ||
      a->getTo()->pl().getParent() == b->getTo()->pl().getParent())
    r = a->getTo()->pl().getParent();
  return r;
}

// _____________________________________________________________________________
void GridGraph::addResidentEdges(CombEdge* ce, GridNode* frGrNd, const std::vector<GridEdge*>& res, GridNode* toGrNd) {
  GridEdge* last = 0;
  _settledEdges[ce].push_back(toGrNd->pl().getParent()->pl().getStepMother());
  for (auto ge : res) {
    if (ge->pl().isSecondary()) continue;
    ge->pl().addResidentEdge(ce);
    if (last) {
      _settledEdges[ce].push_back(sharedParentNode(last, ge)->pl().getParent()->pl().getStepMother());
    }
    last = ge;
  }
  _settledEdges[ce].push_back(frGrNd->pl().getParent()->pl().getStepMother());
}

// _____________________________________________________________________________
const std::vector<GridNode*>& GridGraph::getGridNodes(const CombEdge* ce) const {
  auto it = _settledEdges.find(ce);
  if (it == _settledEdges.end()) return _emptyNdVec;
  return it->second;
}
