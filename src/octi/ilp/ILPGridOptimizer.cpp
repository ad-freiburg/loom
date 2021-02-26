// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <glpk.h>
#include <fstream>
#include "octi/basegraph/BaseGraph.h"
#include "octi/ilp/ILPGridOptimizer.h"
#include "shared/optim/ILPSolvProv.h"
#include "util/log/Log.h"

using octi::basegraph::BaseGraph;
using octi::basegraph::GeoPensMap;
using octi::basegraph::GridEdge;
using octi::basegraph::GridNode;
using octi::ilp::ILPGridOptimizer;
using shared::optim::ILPSolver;
using shared::optim::StarterSol;

// _____________________________________________________________________________
double ILPGridOptimizer::optimize(BaseGraph* gg, const CombGraph& cg,
                                  combgraph::Drawing* d, double maxGrDist,
                                  bool noSolve, const GeoPensMap* geoPensMap,
                                  int timeLim, const std::string& solverStr,
                                  const std::string& path) const {
  // extract first feasible solution from gridgraph
  StarterSol sol = extractFeasibleSol(gg, cg, maxGrDist);
  gg->reset();

  for (auto nd : *gg->getNds()) {
    // if we presolve, some edges may be blocked
    for (auto e : nd->getAdjList()) {
      e->pl().open();
      e->pl().unblock();
    }
    if (!nd->pl().isSink()) continue;
    gg->openTurns(nd);
    gg->closeSinkFr(nd);
    gg->closeSinkTo(nd);
  }

  // clear drawing
  d->crumble();

  auto lp = createProblem(gg, cg, geoPensMap, maxGrDist, solverStr);
  lp->setStarter(sol);

  if (path.size()) {
    std::string basename = path;
    size_t pos = basename.find_last_of(".");
    if (pos != std::string::npos) basename = basename.substr(0, pos);

    std::string outf = basename + ".sol";
    std::string solutionF = basename + ".mst";
    lp->writeMst(solutionF, sol);
    lp->writeMps(path);
  }

  if (noSolve) {
    return std::numeric_limits<double>::infinity();
  } else {
    if (timeLim >= 0) lp->setTimeLim(timeLim);
    auto status = lp->solve();

    if (status == shared::optim::SolveType::INF) {
      throw std::runtime_error(
          "No solution found for ILP problem (most likely because of a time "
          "limit)!");
    }

    extractSolution(lp, gg, cg, d);

    double obj = lp->getObjVal();

    delete lp;

    return  obj;
  }

  return std::numeric_limits<double>::infinity();
}

// _____________________________________________________________________________
ILPSolver* ILPGridOptimizer::createProblem(BaseGraph* gg, const CombGraph& cg,
                                           const GeoPensMap* geoPensMap,
                                           double maxGrDist,
                                           const std::string& solverStr) const {
  ILPSolver* lp = shared::optim::getSolver(solverStr, shared::optim::MIN);

  // grid nodes that may potentially be a position for an
  // input station
  std::map<const CombNode*, std::set<const GridNode*>> cands;

  for (auto nd : cg.getNds()) {
    if (nd->getDeg() == 0) continue;
    std::stringstream oneAssignment;
    // must sum up to 1
    oneAssignment << "oneass(" << nd << ")";
    int rowStat = lp->addRow(oneAssignment.str(), 1, shared::optim::FIX);

    size_t i = 0;

    for (const GridNode* n : *gg->getNds()) {
      if (!n->pl().isSink()) continue;

      double gridD = dist(*n->pl().getGeom(), *nd->pl().getGeom());

      // threshold for speedup
      double maxDis = gg->getCellSize() * maxGrDist;
      if (gridD >= maxDis) continue;

      cands[nd].insert(n);

      gg->openSinkFr(const_cast<GridNode*>(n), 0);
      gg->openSinkTo(const_cast<GridNode*>(n), 0);

      auto varName = getStatPosVar(n, nd);

      int col = lp->addCol(varName, shared::optim::BIN, gg->ndMovePen(nd, n));

      lp->addColToRow(rowStat, col, 1);

      i++;
    }
  }

  // for every edge, we define a binary variable telling us whether this edge
  // is used in a path for the original edge
  for (auto nd : cg.getNds()) {
    for (auto edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      for (const GridNode* n : *gg->getNds()) {
        for (const GridEdge* e : n->getAdjList()) {
          if (e->getFrom() != n) continue;
          if (e->pl().cost() == std::numeric_limits<float>::infinity()) {
            // skip infinite edges, we cannot use them.
            // this also skips sink edges of nodes not used as
            // candidates
            continue;
          }

          if (e->getFrom()->pl().isSink() &&
              !cands[edg->getFrom()].count(e->getFrom())) {
            continue;
          }

          if (e->getTo()->pl().isSink() &&
              !cands[edg->getTo()].count(e->getTo())) {
            continue;
          }

          auto edgeVarName = getEdgeUseVar(e, edg);

          double coef;
          if (geoPensMap && !e->pl().isSecondary()) {
            // add geo pen
            coef = e->pl().cost() +
                   (*geoPensMap).find(edg)->second[e->pl().getId()];
          } else {
            coef = e->pl().cost();
          }
          lp->addCol(edgeVarName, shared::optim::BIN, coef);
        }
      }
    }
  }

  lp->update();

  // an edge can only be used a single time
  std::set<const GridEdge*> proced;
  for (const GridNode* n : *gg->getNds()) {
    for (const GridEdge* e : n->getAdjList()) {
      if (e->pl().isSecondary()) continue;
      if (proced.count(e)) continue;
      auto f = gg->getEdg(e->getTo(), e->getFrom());
      proced.insert(e);
      proced.insert(f);

      std::stringstream constName;
      constName << "ue(" << e->getFrom()->pl().getId() << ","
                << e->getTo()->pl().getId() << ")";
      int row = lp->addRow(constName.str(), 1, shared::optim::UP);

      for (auto nd : cg.getNds()) {
        for (auto edg : nd->getAdjList()) {
          if (edg->getFrom() != nd) continue;
          if (e->pl().cost() == std::numeric_limits<float>::infinity())
            continue;

          auto eVarName = getEdgeUseVar(e, edg);
          auto fVarName = getEdgeUseVar(f, edg);

          int eCol = lp->getVarByName(eVarName);
          if (eCol > -1)
            lp->addColToRow(row, eCol, 1);  // vm.addVar(row, eCol, 1);
          int fCol = lp->getVarByName(fVarName);
          if (fCol > -1)
            lp->addColToRow(row, fCol, 1);  // vm.addVar(row, fCol, 1);
        }
      }
    }
  }

  // for every node, the number of outgoing and incoming used edges must be
  // the same, except for the start and end node
  for (const GridNode* n : *gg->getNds()) {
    if (nonInfDeg(n) == 0) continue;

    for (auto nd : cg.getNds()) {
      for (auto edg : nd->getAdjList()) {
        if (edg->getFrom() != nd) continue;
        std::stringstream constName;
        constName << "as(" << n->pl().getId() << "," << edg << ")";

        // an upper bound is enough here
        int row = lp->addRow(constName.str(), 0, shared::optim::UP);

        // normally, we count an incoming edge as 1 and an outgoing edge as -1
        // later on, we make sure that each node has a some of all out and in
        // edges of 0
        int inCost = -1;
        int outCost = 1;

        // for sink nodes, we apply a trick: an outgoing edge counts as 2 here.
        // this means that a sink node cannot make up for an outgoing edge
        // with an incoming edge - it would need 2 incoming edges to achieve
        // that.
        // however, this would mean (as sink nodes are never adjacent) that 2
        // ports
        // have outgoing edges - which would mean the path "split" somewhere
        // before
        // the ports, which is impossible and forbidden by our other
        // constraints.
        // the only way a sink node can make up for in outgoin edge
        // is thus if we add -2 if the sink is marked as the start station of
        // this edge
        if (n->pl().isSink()) {
          // subtract the variable for this start node and edge, if used
          // as a candidate
          int ndColFrom = lp->getVarByName(getStatPosVar(n, edg->getFrom()));
          if (ndColFrom > -1) lp->addColToRow(row, ndColFrom, -2);

          // add the variable for this end node and edge, if used
          // as a candidate
          int ndColTo = lp->getVarByName(getStatPosVar(n, edg->getTo()));
          if (ndColTo > -1) lp->addColToRow(row, ndColTo, 1);

          outCost = 2;
        }

        for (auto e : n->getAdjListIn()) {
          int edgCol = lp->getVarByName(getEdgeUseVar(e, edg));
          if (edgCol < 0) continue;
          lp->addColToRow(row, edgCol, inCost);
        }

        for (auto e : n->getAdjListOut()) {
          int edgCol = lp->getVarByName(getEdgeUseVar(e, edg));
          if (edgCol < 0) continue;
          lp->addColToRow(row, edgCol, outCost);
        }
      }
    }
  }

  lp->update();

  // only a single sink edge can be activated per input edge and settled grid
  // node
  // THIS RULE IS REDUNDANT AND IMPLICITELY ENFORCED BY OTHER RULES,
  // BUT SEEMS TO LEAD TO FASTER SOLUTION TIMES
  for (GridNode* n : *gg->getNds()) {
    if (!n->pl().isSink()) continue;

    for (auto nd : cg.getNds()) {
      for (auto e : nd->getAdjList()) {
        if (e->getFrom() != nd) continue;

        std::stringstream constName;
        constName << "ss(" << n->pl().getId() << "," << e << ")";

        int row = lp->addRow(constName.str(), 0, shared::optim::FIX);

        if (!cands[e->getFrom()].count(n) && !cands[e->getTo()].count(n)) {
          // node does not appear as start or end cand, so the number of
          // sink edges for this node is 0

        } else {
          if (cands[e->getTo()].count(n)) {
            int ndColTo = lp->getVarByName(getStatPosVar(n, e->getTo()));
            if (ndColTo > -1) lp->addColToRow(row, ndColTo, -1);
          }

          if (cands[e->getFrom()].count(n)) {
            int ndColFr = lp->getVarByName(getStatPosVar(n, e->getFrom()));
            if (ndColFr > -1) lp->addColToRow(row, ndColFr, -1);
          }
        };

        for (size_t p = 0; p < gg->maxDeg(); p++) {
          auto varSinkTo = getEdgeUseVar(gg->getEdg(n->pl().getPort(p), n), e);
          auto varSinkFr = getEdgeUseVar(gg->getEdg(n, n->pl().getPort(p)), e);

          int ndColTo = lp->getVarByName(varSinkTo);
          if (ndColTo > -1) lp->addColToRow(row, ndColTo, 1);

          int ndColFr = lp->getVarByName(varSinkFr);
          if (ndColFr > -1) lp->addColToRow(row, ndColFr, 1);
        }
      }
    }
  }

  // a meta node can either be an activated sink, or a single pass through
  // edge is used
  for (GridNode* n : *gg->getNds()) {
    if (!n->pl().isSink()) continue;

    std::stringstream constName;
    constName << "iu(" << n->pl().getId() << ")";

    int row = lp->addRow(constName.str(), 1, shared::optim::UP);

    // a meta grid node can either be a sink for a single input node, or
    // a pass-through

    for (auto nd : cg.getNds()) {
      int ndcolto = lp->getVarByName(getStatPosVar(n, nd).c_str());
      if (ndcolto > -1) lp->addColToRow(row, ndcolto, 1);
    }

    // go over all ports
    for (size_t pf = 0; pf < gg->maxDeg(); pf++) {
      auto from = n->pl().getPort(pf);
      for (size_t pt = 0; pt < gg->maxDeg(); pt++) {
        auto to = n->pl().getPort(pt);
        if (from == to) continue;

        auto innerE = gg->getEdg(from, to);
        for (auto nd : cg.getNds()) {
          for (auto edg : nd->getAdjList()) {
            if (edg->getFrom() != nd) continue;

            int edgCol = lp->getVarByName(getEdgeUseVar(innerE, edg));
            if (edgCol < 0) continue;
            lp->addColToRow(row, edgCol, 1);
          }
        }
      }
    }
  }

  lp->update();

  // dont allow crossing edges
  size_t rowId = 0;
  for (auto edgPair : gg->getCrossEdgPairs()) {
    std::stringstream constName;
    constName << "nc(" << rowId << ")";
    rowId++;

    int row = lp->addRow(constName.str(), 1, shared::optim::UP);

    for (auto nd : cg.getNds()) {
      for (auto edg : nd->getAdjList()) {
        if (edg->getFrom() != nd) continue;

        int col = lp->getVarByName(getEdgeUseVar(edgPair.first.first, edg));
        if (col > -1) lp->addColToRow(row, col, 1);

        col = lp->getVarByName(getEdgeUseVar(edgPair.first.second, edg));
        if (col > -1) lp->addColToRow(row, col, 1);

        col = lp->getVarByName(getEdgeUseVar(edgPair.second.first, edg));
        if (col > -1) lp->addColToRow(row, col, 1);

        col = lp->getVarByName(getEdgeUseVar(edgPair.second.second, edg));
        if (col > -1) lp->addColToRow(row, col, 1);
      }
    }
  }

  lp->update();

  // for each input node N, define a var x_dirNE which tells the direction of
  // E at N
  for (auto nd : cg.getNds()) {
    if (nd->getDeg() < 2) continue;  // we don't need this for deg 1 nodes
    for (auto edg : nd->getAdjList()) {
      std::stringstream dirName;
      dirName << "d(" << nd << "," << edg << ")";
      int col = lp->addCol(dirName.str(), shared::optim::INT, 0, 0, 7);

      std::stringstream constName;
      constName << "dc(" << nd << "," << edg << ")";

      int row = lp->addRow(constName.str(), 0, shared::optim::FIX);

      lp->addColToRow(row, col, -1);

      for (GridNode* n : *gg->getNds()) {
        if (!n->pl().isSink()) continue;

        if (edg->getFrom() == nd) {
          // the 0 can be skipped here
          for (size_t i = 1; i < gg->maxDeg(); i++) {
            auto e = gg->getEdg(n, n->pl().getPort(i));
            int col = lp->getVarByName(getEdgeUseVar(e, edg));
            if (col > -1) lp->addColToRow(row, col, i);
          }
        } else {
          // the 0 can be skipped here
          for (size_t i = 1; i < gg->maxDeg(); i++) {
            auto e = gg->getEdg(n->pl().getPort(i), n);
            int col = lp->getVarByName(getEdgeUseVar(e, edg));
            if (col > -1) lp->addColToRow(row, col, i);
          }
        }
      }
    }
  }

  lp->update();

  // for each input node N, make sure that the circular ordering of the final
  // drawing matches the input ordering
  int M = gg->maxDeg();
  for (auto nd : cg.getNds()) {
    // for degree < 3, the circular ordering cannot be violated
    if (nd->getDeg() < 3) continue;

    std::stringstream vulnConstName;
    vulnConstName << "vc(" << nd << ")";
    // an upper bound would also work here, at most one
    // of the vuln vars may be 1

    int vulnRow = lp->addRow(vulnConstName.str(), 1, shared::optim::FIX);

    for (size_t i = 0; i < nd->getDeg(); i++) {
      std::stringstream n;
      n << "vuln(" << nd << "," << i << ")";
      int col = lp->addCol(n.str(), shared::optim::BIN, 0);
      lp->addColToRow(vulnRow, col, 1);
    }

    lp->update();

    auto order = nd->pl().getEdgeOrdering().getOrderedSet();
    assert(order.size() > 2);
    for (size_t i = 0; i < order.size(); i++) {
      CombEdge* edgA;
      if (i == 0) {
        edgA = nd->pl().getEdgeOrdering().getOrderedSet().back().first;
      } else {
        edgA = nd->pl().getEdgeOrdering().getOrderedSet()[i - 1].first;
      }
      auto edgB = nd->pl().getEdgeOrdering().getOrderedSet()[i].first;

      assert(edgA != edgB);

      std::stringstream colNameA;
      colNameA << "d(" << nd << "," << edgA << ")";
      int colA = lp->getVarByName(colNameA.str());
      assert(colA > -1);

      std::stringstream colNameB;
      colNameB << "d(" << nd << "," << edgB << ")";
      int colB = lp->getVarByName(colNameB.str());
      assert(colB > -1);

      std::stringstream constName;
      constName << "oc(" << nd << "," << i << ")";
      int row = lp->addRow(constName.str(), 1, shared::optim::LO);

      std::stringstream vulnColName;
      vulnColName << "vuln(" << nd << "," << i << ")";
      int vulnCol = lp->getVarByName(vulnColName.str());
      assert(vulnCol > -1);

      lp->addColToRow(row, colB, 1);
      lp->addColToRow(row, colA, -1);
      lp->addColToRow(row, vulnCol, M);
    }
  }

  lp->update();

  // for each adjacent edge pair, add variables telling the accuteness of the
  // angle between them
  for (auto nd : cg.getNds()) {
    for (size_t i = 0; i < nd->getAdjList().size(); i++) {
      auto edgA = nd->getAdjList()[i];
      for (size_t j = i + 1; j < nd->getAdjList().size(); j++) {
        auto edgB = nd->getAdjList()[j];
        assert(edgA != edgB);

        // note: we can identify pairs of edges by the edges only as we dont
        // have a multigraph - we dont need the need for uniqueness

        size_t sharedLines = 0;
        // TODO: not all lines in getChilds are equal, take the "right" end of
        // the childs here!
        for (auto lo : edgA->pl().getChilds().front()->pl().getLines()) {
          if (edgB->pl().getChilds().front()->pl().hasLine(lo.line)) {
            sharedLines++;
          }
        }

        if (!sharedLines) continue;

        std::stringstream negVar;
        negVar << "negdist(" << edgA << "," << edgB << ")";

        int colNeg = lp->addCol(negVar.str(), shared::optim::BIN, 0);

        std::stringstream constName;
        constName << "nc(" << edgA << "," << edgB << ")";

        int row1 = lp->addRow(constName.str() + "lo", 0, shared::optim::LO);
        int row2 = lp->addRow(constName.str() + "up", 7, shared::optim::UP);

        std::stringstream dirNameA;
        dirNameA << "d(" << nd << "," << edgA << ")";
        std::stringstream dirNameB;
        dirNameB << "d(" << nd << "," << edgB << ")";

        int colA = lp->getVarByName(dirNameA.str());
        assert(colA > -1);
        lp->addColToRow(row1, colA, 1);
        lp->addColToRow(row2, colA, 1);

        int colB = lp->getVarByName(dirNameB.str());
        assert(colB > -1);
        lp->addColToRow(row1, colB, -1);
        lp->addColToRow(row2, colB, -1);

        lp->addColToRow(row1, colNeg, 8);
        lp->addColToRow(row2, colNeg, 8);

        std::stringstream angConst;
        angConst << "ac(" << edgA << "," << edgB << ")";
        int rowAng = lp->addRow(angConst.str(), 0, shared::optim::FIX);

        lp->addColToRow(rowAng, colA, 1);
        lp->addColToRow(rowAng, colB, -1);
        lp->addColToRow(rowAng, colNeg, 8);

        std::stringstream sumConst;
        sumConst << "asc(" << edgA << "," << edgB << ")";

        int rowSum = lp->addRow(sumConst.str(), 1, shared::optim::UP);

        std::vector<std::string> names = {"d45",   "d90",  "d135", "d180",
                                          "d135'", "d90'", "d45'"};

        // TODO: derive from configuration!!!
        std::vector<double> pens = {3, 2.5, 2, 1, 2, 2.5, 3};

        for (int k = 0; k < 7; k++) {
          std::stringstream var;
          var << names[k] << "(" << edgA << "," << edgB << ")";

          // TODO: maybe multiply per shared lines - but this actually
          // makes the drawings look worse.
          int col = lp->addCol(var.str(), shared::optim::BIN, pens[k]);

          lp->addColToRow(rowAng, col, -(k + 1));
          lp->addColToRow(rowSum, col, 1);
        }
      }
    }
  }

  lp->update();

  return lp;
}

// _____________________________________________________________________________
std::string ILPGridOptimizer::getEdgeUseVar(const GridEdge* e,
                                            const CombEdge* cg) const {
  std::stringstream varName;
  varName << "edg(" << e->getFrom()->pl().getId() << ","
          << e->getTo()->pl().getId() << "," << cg << ")";

  return varName.str();
}

// _____________________________________________________________________________
std::string ILPGridOptimizer::getStatPosVar(const GridNode* n,
                                            const CombNode* nd) const {
  std::stringstream varName;
  varName << "sp(" << n->pl().getId() << "," << nd << ")";

  return varName.str();
}

// _____________________________________________________________________________
void ILPGridOptimizer::extractSolution(ILPSolver* lp, BaseGraph* gg,
                                       const CombGraph& cg,
                                       combgraph::Drawing* d) const {
  std::map<const CombNode*, const GridNode*> gridNds;
  std::map<const CombEdge*, std::set<const GridEdge*>> gridEdgs;

  // write solution to grid graph
  for (GridNode* n : *gg->getNds()) {
    for (GridEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;

      for (auto nd : cg.getNds()) {
        for (auto edg : nd->getAdjList()) {
          if (edg->getFrom() != nd) continue;
          auto varName = getEdgeUseVar(e, edg);

          int i = lp->getVarByName(varName);
          if (i > -1) {
            double val = lp->getVarVal(i);
            if (val > 0.5) {
              gg->addResEdg(e, edg);
              gridEdgs[edg].insert(e);
            }
          }
        }
      }
    }
  }

  for (GridNode* n : *gg->getNds()) {
    if (!n->pl().isSink()) continue;
    for (auto nd : cg.getNds()) {
      auto varName = getStatPosVar(n, nd);

      int i = lp->getVarByName(varName);
      if (i > -1) {
        double val = lp->getVarVal(i);
        if (val > 0.5) {
          n->pl().setStation();
          gridNds[nd] = n;
        }
      }
    }
  }

  // draw solution
  for (auto nd : cg.getNds()) {
    for (auto edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;

      std::vector<GridEdge*> edges(gridEdgs[edg].size());

      // get the start and end grid nodes
      auto grStart = gridNds[edg->getFrom()];
      auto grEnd = gridNds[edg->getTo()];

      auto curNode = grStart;
      GridEdge* last = 0;

      size_t i = 0;

      while (curNode != grEnd) {
        for (auto adj : curNode->getAdjList()) {
          if (adj != last && gridEdgs[edg].count(adj)) {
            last = adj;
            i++;
            edges[edges.size() - i] = adj;
            curNode = adj->getOtherNd(curNode);
            break;
          }
        }
      }

      assert(i == edges.size());

      d->draw(edg, edges, false);
    }
  }
}

// _____________________________________________________________________________
size_t ILPGridOptimizer::nonInfDeg(const GridNode* g) const {
  size_t ret = 0;
  for (auto e : g->getAdjList()) {
    if (e->pl().cost() != std::numeric_limits<float>::infinity()) ret++;
  }

  return ret;
}

// _____________________________________________________________________________
StarterSol ILPGridOptimizer::extractFeasibleSol(BaseGraph* gg,
                                                const CombGraph& cg,
                                                double maxGrDist) const {
  StarterSol sol;

  for (auto nd : cg.getNds()) {
    if (nd->getDeg() == 0) continue;
    auto settled = gg->getSettled(nd);

    for (auto gnd : *gg->getNds()) {
      if (!gnd->pl().isSink()) continue;
      double gridD = dist(*nd->pl().getGeom(), *gnd->pl().getGeom());

      // threshold for speedup
      double maxDis = gg->getCellSize() * maxGrDist;
      if (gridD >= maxDis) continue;

      auto varName = getStatPosVar(gnd, nd);
      if (gnd == settled) {
        sol[varName] = 1;

        // if settled, all bend edges are unused
        for (size_t p = 0; p < gg->maxDeg(); p++) {
          auto portNd = gnd->pl().getPort(p);
          for (auto bendEdg : portNd->getAdjList()) {
            if (!bendEdg->pl().isSecondary()) continue;
            for (auto cEdg : nd->getAdjList()) {
              if (cEdg->getFrom() != nd) continue;
              auto varName = getEdgeUseVar(bendEdg, cEdg);
              sol[varName] = 0;
            }
          }
        }
      } else {
        sol[varName] = 0;

        // if not settled, all sink edges are unused
        // for all input edges
        for (auto sinkEdg : gnd->getAdjList()) {
          assert(sinkEdg->pl().isSecondary());
          for (auto cEdg : nd->getAdjList()) {
            if (cEdg->getFrom() != nd) continue;
            auto varName = getEdgeUseVar(sinkEdg, cEdg);
            sol[varName] = 0;
          }
        }
      }
    }
  }

  for (auto grNd : *gg->getNds()) {
    for (auto grEdg : grNd->getAdjListOut()) {
      if (grEdg->pl().isSecondary()) continue;
      // we assume here that the heuristic solution is feasible!
      auto resEdg = *gg->getResEdgs(grEdg).begin();

      if (resEdg) {
        auto varName = getEdgeUseVar(grEdg, resEdg);
        sol[varName] = 1;
      } else {
        for (auto cNd : cg.getNds()) {
          for (auto cEdg : cNd->getAdjList()) {
            if (cEdg->getFrom() != cNd) continue;
            auto varName = getEdgeUseVar(grEdg, cEdg);
            sol[varName] = 0;
          }
        }
      }
    }
  }

  return sol;
}
