// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <glpk.h>
#include <fstream>
#include "octi/gridgraph/GridGraph.h"
#include "octi/ilp/ILPGridOptimizer.h"
#include "util/log/Log.h"

using octi::gridgraph::GridGraph;
using octi::gridgraph::GridNode;
using octi::gridgraph::GridEdge;
using octi::ilp::ILPGridOptimizer;
using octi::ilp::VariableMatrix;

// _____________________________________________________________________________
double ILPGridOptimizer::optimize(GridGraph* gg, const CombGraph& cg,
                                  combgraph::Drawing* d,
                                  double maxGrDist) const {
  for (auto nd : *gg->getNds()) {
    if (!nd->pl().isSink()) continue;
    gg->openNodeTurns(nd);
    gg->openNodeSink(nd, 0);
  }

  glp_prob* lp = createProblem(*gg, cg, maxGrDist);

  // TODO: make configurable
  glp_write_mps(lp, GLP_MPS_FILE, 0, "ilp_prob.mps");

  preSolve(lp);
  solveProblem(lp);

  extractSolution(lp, gg, cg, d);

  double score = glp_get_obj_val(lp);

  glp_delete_prob(lp);
  glp_free_env();

  return score;
}

// _____________________________________________________________________________
glp_prob* ILPGridOptimizer::createProblem(const GridGraph& gg,
                                          const CombGraph& cg,
                                          double maxGrDist) const {
  glp_prob* lp = glp_create_prob();

  glp_set_prob_name(lp, "griddrawing");
  glp_set_obj_dir(lp, GLP_MIN);

  VariableMatrix vm;

  glp_create_index(lp);

  for (auto nd : cg.getNds()) {
    if (nd->getDeg() == 0) continue;
    std::stringstream oneAssignment;
    // must sum up to 1
    size_t rowStat = glp_add_rows(lp, 1);
    oneAssignment << "oneass(st=" << nd << ")";
    glp_set_row_name(lp, rowStat, oneAssignment.str().c_str());
    glp_set_row_bnds(lp, rowStat, GLP_FX, 1, 1);

    size_t i = 0;

    for (const GridNode* n : gg.getNds()) {
      if (!n->pl().isSink()) continue;

      double gridD = dist(*n->pl().getGeom(), *nd->pl().getGeom());

      // threshold for speedup
      double maxDis = gg.getCellSize() * maxGrDist;
      if (gridD >= maxDis) continue;

      auto varName = getStatPosVar(n, nd);

      size_t col = glp_add_cols(lp, 1);
      glp_set_col_name(lp, col, varName.c_str());
      // binary variable € {0,1}, node is either this station, or not
      glp_set_col_kind(lp, col, GLP_BV);

      glp_set_obj_coef(lp, col, gg.ndMovePen(nd, n));

      vm.addVar(rowStat, col, 1);
      i++;
    }
  }

  // for every edge, we define a binary variable telling us whether this edge
  // is used in a path for the original edge
  for (auto nd : cg.getNds()) {
    for (auto edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      for (const GridNode* n : gg.getNds()) {
        for (const GridEdge* e : n->getAdjList()) {
          if (e->getFrom() != n) continue;
          if (e->pl().cost() == std::numeric_limits<float>::infinity())
            continue;

          auto edgeVarName = getEdgeUseVar(e, edg);

          size_t col = glp_add_cols(lp, 1);
          glp_set_col_name(lp, col, edgeVarName.c_str());
          // binary variable € {0,1}, edge is either used, or not
          glp_set_col_kind(lp, col, GLP_BV);

          glp_set_obj_coef(lp, col, e->pl().cost());
        }
      }
    }
  }

  glp_create_index(lp);

  // for every node, the number of outgoing and incoming used edges must be
  // the same, except for the start and end node
  for (const GridNode* n : gg.getNds()) {
    for (auto nd : cg.getNds()) {
      for (auto edg : nd->getAdjList()) {
        if (edg->getFrom() != nd) continue;
        std::stringstream constName;
        size_t row = glp_add_rows(lp, 1);
        constName << "adjsum(" << n << "," << edg << ")";
        glp_set_row_name(lp, row, constName.str().c_str());

        // an upper bound is enough here
        glp_set_row_bnds(lp, row, GLP_UP, 0, 0);

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
        // this
        // edge
        if (n->pl().isSink()) {
          // subtract the variable for this start node and edge, if used
          // as a candidate
          std::stringstream ndPosFromVarName;
          ndPosFromVarName << "statpos(" << n << "," << edg->getFrom() << ")";
          size_t ndColFrom = glp_find_col(lp, ndPosFromVarName.str().c_str());
          if (ndColFrom > 0) vm.addVar(row, ndColFrom, -2);

          // add the variable for this end node and edge, if used
          // as a candidate
          std::stringstream ndPosToVarName;
          ndPosToVarName << "statpos(" << n << "," << edg->getTo() << ")";
          size_t ndColTo = glp_find_col(lp, ndPosToVarName.str().c_str());
          if (ndColTo > 0) vm.addVar(row, ndColTo, 1);

          outCost = 2;
        }

        for (auto e : n->getAdjListIn()) {
          size_t edgCol = glp_find_col(lp, getEdgeUseVar(e, edg).c_str());
          if (!edgCol) continue;
          vm.addVar(row, edgCol, inCost);
        }

        for (auto e : n->getAdjListOut()) {
          size_t edgCol = glp_find_col(lp, getEdgeUseVar(e, edg).c_str());
          if (!edgCol) continue;
          vm.addVar(row, edgCol, outCost);
        }
      }
    }
  }

  glp_create_index(lp);

  // a meta node can either be an activated sink, or a single pass through
  // edge is used
  for (const GridNode* n : gg.getNds()) {
    if (!n->pl().isSink()) continue;

    std::stringstream constName;
    size_t row = glp_add_rows(lp, 1);
    constName << "inneruse(" << n << ")";
    glp_set_row_name(lp, row, constName.str().c_str());
    glp_set_row_bnds(lp, row, GLP_UP, 0, 1);

    // a meta grid node can either be a sink for a single input node, or
    // a pass-through

    for (auto nd : cg.getNds()) {
      std::stringstream ndPosToVarName;
      ndPosToVarName << "statpos(" << n << "," << nd << ")";
      size_t ndColTo = glp_find_col(lp, ndPosToVarName.str().c_str());
      if (ndColTo > 0) vm.addVar(row, ndColTo, 1);
    }

    // go over all ports
    for (size_t pf = 0; pf < 8; pf++) {
      auto from = n->pl().getPort(pf);
      for (size_t pt = 0; pt < 8; pt++) {
        auto to = n->pl().getPort(pt);
        if (from == to) continue;

        auto innerE = gg.getEdg(from, to);
        for (auto nd : cg.getNds()) {
          for (auto edg : nd->getAdjList()) {
            if (edg->getFrom() != nd) continue;
            size_t edgCol =
                glp_find_col(lp, getEdgeUseVar(innerE, edg).c_str());
            if (!edgCol) continue;
            vm.addVar(row, edgCol, 1);
          }
        }
      }
    }
  }

  glp_create_index(lp);

  // dont allow crossing edges
  for (const GridNode* n : gg.getNds()) {
    if (!n->pl().isSink()) continue;
    size_t x = n->pl().getX();
    size_t y = n->pl().getY();

    std::stringstream constName;
    size_t row = glp_add_rows(lp, 1);
    constName << "nocross(" << n << ")";
    glp_set_row_name(lp, row, constName.str().c_str());
    glp_set_row_bnds(lp, row, GLP_UP, 0, 1);

    auto eOr = gg.getNEdg(n, gg.getNeighbor(x, y, 3));
    auto fOr = gg.getNEdg(gg.getNeighbor(x, y, 3), n);

    if (!eOr || !fOr) continue;

    auto na = gg.getNeighbor(x, y, (3 + 7) % 8);
    auto nb = gg.getNeighbor(x, y, (3 + 1) % 8);

    if (!na || !nb) continue;

    auto e = gg.getNEdg(na, nb);
    auto f = gg.getNEdg(nb, na);

    for (auto nd : cg.getNds()) {
      for (auto edg : nd->getAdjList()) {
        if (edg->getFrom() != nd) continue;
        size_t col = glp_find_col(lp, getEdgeUseVar(eOr, edg).c_str());
        if (col) vm.addVar(row, col, 1);

        col = glp_find_col(lp, getEdgeUseVar(fOr, edg).c_str());
        if (col) vm.addVar(row, col, 1);

        col = glp_find_col(lp, getEdgeUseVar(e, edg).c_str());
        if (col) vm.addVar(row, col, 1);

        col = glp_find_col(lp, getEdgeUseVar(f, edg).c_str());
        if (col) vm.addVar(row, col, 1);
      }
    }
  }

  glp_create_index(lp);

  // for each input node N, define a var x_dirNE which tells the direction of
  // E at N
  for (auto nd : cg.getNds()) {
    if (nd->getDeg() < 2) continue;  // we don't need this for deg 1 nodes
    for (auto edg : nd->getAdjList()) {
      std::stringstream dirName;
      dirName << "dir(" << nd << "," << edg << ")";
      size_t col = glp_add_cols(lp, 1);
      assert(col);
      glp_set_col_name(lp, col, dirName.str().c_str());
      glp_set_col_kind(lp, col, GLP_IV);
      glp_set_col_bnds(lp, col, GLP_UP, 0, 7);

      std::stringstream constName;
      constName << "dirconst(" << nd << "," << edg << ")";
      size_t row = glp_add_rows(lp, 1);
      assert(row);
      glp_set_row_name(lp, row, constName.str().c_str());
      glp_set_row_bnds(lp, row, GLP_FX, 0, 0);

      vm.addVar(row, col, -1);

      for (GridNode* n : gg.getNds()) {
        if (!n->pl().isSink()) continue;

        if (edg->getFrom() == nd) {
          // the 0 can be skipped here
          for (size_t i = 1; i < 8; i++) {
            auto e = gg.getEdg(n, n->pl().getPort(i));
            size_t col = glp_find_col(lp, getEdgeUseVar(e, edg).c_str());
            if (col) vm.addVar(row, col, i);
          }
        } else {
          // the 0 can be skipped here
          for (size_t i = 1; i < 8; i++) {
            auto e = gg.getEdg(n->pl().getPort(i), n);
            size_t col = glp_find_col(lp, getEdgeUseVar(e, edg).c_str());
            if (col) vm.addVar(row, col, i);
          }
        }
      }
    }
  }

  glp_create_index(lp);

  // for each input node N, make sure that the circular ordering of the final
  // drawing matches the input ordering
  int M = 8;
  for (auto nd : cg.getNds()) {
    // for degree < 3, the circular ordering cannot be violated
    if (nd->getDeg() < 3) continue;

    std::stringstream vulnConstName;
    vulnConstName << "vulnconst(" << nd << ")";
    size_t vulnRow = glp_add_rows(lp, 1);
    glp_set_row_name(lp, vulnRow, vulnConstName.str().c_str());
    // an upper bound would also work here, at most one
    // of the vuln vars may be 1
    glp_set_row_bnds(lp, vulnRow, GLP_FX, 1, 1);

    for (size_t i = 0; i < nd->getDeg(); i++) {
      std::stringstream n;
      n << "vuln(" << nd << "," << i << ")";
      size_t col = glp_add_cols(lp, 1);
      glp_set_col_name(lp, col, n.str().c_str());
      glp_set_col_kind(lp, col, GLP_BV);

      vm.addVar(vulnRow, col, 1);
    }

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
      colNameA << "dir(" << nd << "," << edgA << ")";
      size_t colA = glp_find_col(lp, colNameA.str().c_str());
      assert(colA);

      std::stringstream colNameB;
      colNameB << "dir(" << nd << "," << edgB << ")";
      size_t colB = glp_find_col(lp, colNameB.str().c_str());
      assert(colB);

      std::stringstream constName;
      constName << "orderconst(" << nd << "," << i << ")";
      size_t row = glp_add_rows(lp, 1);
      glp_set_row_name(lp, row, constName.str().c_str());
      glp_set_row_bnds(lp, row, GLP_LO, 1, 1);

      std::stringstream vulnColName;
      vulnColName << "vuln(" << nd << "," << i << ")";
      size_t vulnCol = glp_find_col(lp, vulnColName.str().c_str());
      assert(vulnCol);

      vm.addVar(row, colB, 1);
      vm.addVar(row, colA, -1);
      vm.addVar(row, vulnCol, M);
    }
  }

  glp_create_index(lp);

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
        for (auto ro : edgA->pl().getChilds().front()->pl().getRoutes()) {
          if (edgB->pl().getChilds().front()->pl().hasRoute(ro.route)) {
            sharedLines++;
          }
        }

        if (!sharedLines) continue;

        std::stringstream negVar;
        negVar << "negdist(" << edgA << "," << edgB << ")";
        size_t colNeg = glp_add_cols(lp, 1);
        assert(colNeg);
        glp_set_col_name(lp, colNeg, negVar.str().c_str());
        glp_set_col_kind(lp, colNeg, GLP_BV);

        std::stringstream constName;
        constName << "negconst(" << edgA << "," << edgB << ")";
        size_t row = glp_add_rows(lp, 1);
        assert(row);
        glp_set_row_name(lp, row, constName.str().c_str());
        glp_set_row_bnds(lp, row, GLP_DB, 0, 7);

        std::stringstream dirNameA;
        dirNameA << "dir(" << nd << "," << edgA << ")";
        std::stringstream dirNameB;
        dirNameB << "dir(" << nd << "," << edgB << ")";

        size_t colA = glp_find_col(lp, dirNameA.str().c_str());
        assert(colA);
        vm.addVar(row, colA, 1);

        size_t colB = glp_find_col(lp, dirNameB.str().c_str());
        assert(colB);
        vm.addVar(row, colB, -1);

        vm.addVar(row, colNeg, 8);

        std::stringstream angConst;
        angConst << "angconst(" << edgA << "," << edgB << ")";
        size_t rowAng = glp_add_rows(lp, 1);
        assert(rowAng);
        glp_set_row_name(lp, rowAng, angConst.str().c_str());
        glp_set_row_bnds(lp, rowAng, GLP_FX, 0, 0);

        vm.addVar(rowAng, colA, 1);
        vm.addVar(rowAng, colB, -1);
        vm.addVar(rowAng, colNeg, 8);

        std::stringstream sumConst;
        sumConst << "angsumconst(" << edgA << "," << edgB << ")";
        size_t rowSum = glp_add_rows(lp, 1);
        assert(rowSum);
        glp_set_row_name(lp, rowSum, sumConst.str().c_str());
        glp_set_row_bnds(lp, rowSum, GLP_UP, 0, 1);

        std::vector<std::string> names = {"d45",   "d90",  "d135", "d180",
                                          "d135'", "d90'", "d45'"};

        // TODO: derive from configuration!!!
        std::vector<double> pens = {3, 2.5, 2, 1, 2, 2.5, 3};

        for (int k = 0; k < 7; k++) {
          std::stringstream var;
          var << names[k] << "(" << edgA << "," << edgB << ")";
          size_t col = glp_add_cols(lp, 1);
          glp_set_col_name(lp, col, var.str().c_str());
          glp_set_col_kind(lp, col, GLP_BV);

          vm.addVar(rowAng, col, -(k + 1));
          vm.addVar(rowSum, col, 1);

          // TODO: maybe multiply per shared lines - but this actually
          // makes the drawings look worse.
          glp_set_obj_coef(lp, col, pens[k]);
        }
      }
    }
  }

  glp_create_index(lp);

  int* ia = 0;
  int* ja = 0;
  double* res = 0;
  vm.getGLPKArrs(&ia, &ja, &res);

  glp_load_matrix(lp, vm.getNumVars(), ia, ja, res);

  delete[](ia);
  delete[](ja);
  delete[](res);

  return lp;
}

// _____________________________________________________________________________
void ILPGridOptimizer::preSolve(glp_prob* lp) const {
  // write temporary file
  // std::string f = std::string(std::tmpnam(0)) + ".mps";
  std::string f = "prob.mps";
  std::string outf = "sol.sol";  // std::string(std::tmpnam(0)) + ".sol";
  // glp_term_out(GLP_OFF);
  glp_write_mps(lp, GLP_MPS_FILE, 0, f.c_str());

  // std::string cmd =
  // "/home/patrick/repos/Cbc-2.9/bin/cbc {INPUT} -randomCbcSeed 0 -threads "
  // "{THREADS} -printingOptions rows -solve -solution {OUTPUT}";

  std::string cmd = "gurobi_cl ResultFile={OUTPUT} {INPUT} > ./gurobi.log";
  util::replaceAll(cmd, "{INPUT}", f);
  util::replaceAll(cmd, "{OUTPUT}", outf);
  util::replaceAll(cmd, "{THREADS}", "4");
  system(cmd.c_str());

  std::ifstream fin;
  fin.open(outf.c_str());
  std::string line;

  // skip first line
  std::getline(fin, line);

  // glp_term_out(GLP_OFF);

  while (std::getline(fin, line)) {
    std::istringstream iss(line);
    int number;
    std::string name;
    double value;

    iss >> number;

    // could not read line number, is missing, which is fine for us
    if (iss.fail()) iss.clear();

    iss >> name;
    iss >> value;

    // gurobi may sometimes write non-integer values to the solution file
    // (within a small tolerance). Values like 0.999999999 will then be
    // truncated to 0 without proper rounding
    value = round(value);

    int intVal = value;

    size_t col = glp_find_col(lp, name.c_str());
    if (col != 0) {
      glp_set_col_bnds(lp, col, GLP_FX, intVal, intVal);
    }
  }
}

// _____________________________________________________________________________
void ILPGridOptimizer::solveProblem(glp_prob* lp) const {
  glp_iocp params;
  glp_smcp sparams;
  // default initialization
  glp_init_iocp(&params);
  glp_init_smcp(&sparams);

  // params.presolve = GLP_ON;
  // params.binarize = GLP_OFF;
  // params.ps_tm_lim = 10000;
  // params.tm_lim = 30000;
  // params.fp_heur = GLP_ON;
  // params.ps_heur = GLP_ON;

  // glp_term_out(GLP_OFF);
  glp_simplex(lp, &sparams);
  glp_intopt(lp, &params);
}

// _____________________________________________________________________________
void VariableMatrix::addVar(int row, int col, double val) {
  rowNum.push_back(row);
  colNum.push_back(col);
  vals.push_back(val);
}

// _____________________________________________________________________________
void VariableMatrix::getGLPKArrs(int** ia, int** ja, double** r) const {
  assert(rowNum.size() == colNum.size());
  assert(colNum.size() == vals.size());

  *ia = new int[rowNum.size() + 1];
  *ja = new int[rowNum.size() + 1];
  *r = new double[rowNum.size() + 1];

  // glpk arrays always start at 1 for some reason
  for (size_t i = 1; i <= rowNum.size(); ++i) {
    (*ia)[i] = rowNum[i - 1];
    (*ja)[i] = colNum[i - 1];
    (*r)[i] = vals[i - 1];
  }
}

// _____________________________________________________________________________
std::string ILPGridOptimizer::getEdgeUseVar(const GridEdge* e,
                                            const CombEdge* cg) const {
  std::stringstream varName;
  varName << "edg(" << e->getFrom() << "," << e->getTo() << "," << cg << ")";

  return varName.str();
}

// _____________________________________________________________________________
std::string ILPGridOptimizer::getStatPosVar(const GridNode* n,
                                            const CombNode* nd) const {
  std::stringstream varName;
  varName << "statpos(" << n << "," << nd << ")";

  return varName.str();
}

// _____________________________________________________________________________
void ILPGridOptimizer::extractSolution(glp_prob* lp, GridGraph* gg,
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

          size_t i = glp_find_col(lp, varName.c_str());
          if (i) {
            double val = glp_mip_col_val(lp, i);
            if (val > 0) {
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

      size_t i = glp_find_col(lp, varName.c_str());
      if (i) {
        double val = glp_mip_col_val(lp, i);
        if (val > 0) {
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

      // d->draw(edg, edges, false);
    }
  }
}
