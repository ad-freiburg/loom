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
int ILPGridOptimizer::optimize(GridGraph* gg, const CombGraph& cg) const {
  std::vector<CombNode*> nodes;
  std::vector<CombEdge*> edges;

  for (auto nd : cg.getNds()) {
    nodes.push_back(nd);
    for (auto e : nd->getAdjList()) {
      if (e->getFrom() == nd) {
        edges.push_back(e);
      }
    }
  }

  // edges.push_back((*cg.getNds().begin())->getAdjList().front());
  // nodes.push_back((*cg.getNds().begin())->getAdjList().front()->getFrom());
  // nodes.push_back((*cg.getNds().begin())->getAdjList().front()->getTo());

  LOG(INFO) << nodes.size() << " nodes, " << edges.size() << " edges.";

  for (auto nd : *gg->getNds()) {
    if (!nd->pl().isSink()) continue;
    gg->openNode(nd);
    gg->openNodeSink(nd, 0);
  }

  LOG(INFO) << "Creating ILP problem... ";
  glp_prob* lp = createProblem(*gg, nodes, edges);
  LOG(INFO) << " .. done";

  LOG(INFO) << "(stats) ILP has " << glp_get_num_cols(lp) << " cols and "
            << glp_get_num_rows(lp) << " rows.";

  // TODO: make configurable
  glp_write_mps(lp, GLP_MPS_FILE, 0, "ilp_prob.mps");

  LOG(INFO) << "Solving problem...";
  preSolve(lp);
  solveProblem(lp);

  LOG(INFO) << "(stats) ILP obj = " << glp_mip_obj_val(lp);

  // write solution to grid graph
  for (GridNode* n : *gg->getNds()) {
    for (GridEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;

      for (auto edg : edges) {
        auto varName = getEdgeUseVar(e, edg);

        size_t i = glp_find_col(lp, varName.c_str());
        if (i) {
          double val = glp_mip_col_val(lp, i);
          if (val > 0) e->pl().addResidentEdge(0);
        }
      }
    }
  }

  for (GridNode* n : *gg->getNds()) {
    if (!n->pl().isSink()) continue;
    for (auto nd : nodes) {
      auto varName = getStatPosVar(n, nd);

      size_t i = glp_find_col(lp, varName.c_str());
      if (i) {
        double val = glp_mip_col_val(lp, i);
        if (val > 0) {
          n->pl().setStation();
        }
      }
    }
  }

  glp_delete_prob(lp);
  glp_free_env();

  return 0;
}

// _____________________________________________________________________________
glp_prob* ILPGridOptimizer::createProblem(
    const GridGraph& gg, const std::vector<CombNode*>& nds,
    const std::vector<CombEdge*>& edgs) const {
  glp_prob* lp = glp_create_prob();

  glp_set_prob_name(lp, "griddrawing");
  glp_set_obj_dir(lp, GLP_MIN);

  VariableMatrix vm;

  glp_create_index(lp);

  size_t i = 0;

  for (auto nd : nds) {
    double nx = nd->pl().getGeom()->getX();
    double ny = nd->pl().getGeom()->getY();

    std::stringstream oneAssignment;
    // must sum up to 1
    size_t rowStat = glp_add_rows(lp, 1);
    oneAssignment << "oneass(st=" << nd << ")";
    glp_set_row_name(lp, rowStat, oneAssignment.str().c_str());
    glp_set_row_bnds(lp, rowStat, GLP_FX, 1, 1);

    for (const GridNode* n : gg.getNds()) {
      if (!n->pl().isSink()) continue;

      double d = sqrt(
          (nx - n->pl().getGeom()->getX()) * (nx - n->pl().getGeom()->getX()) +
          (ny - n->pl().getGeom()->getY()) * (ny - n->pl().getGeom()->getY()));

      if (d > 2000) continue;

      auto varName = getStatPosVar(n, nd);

      size_t col = glp_add_cols(lp, 1);
      glp_set_col_name(lp, col, varName.c_str());
      // binary variable € {0,1}, node is either this station, or not
      glp_set_col_kind(lp, col, GLP_BV);
      i++;

      // TODO: correct distance, at the moment eucludian distance based on grid
      // coord
      // TODO: the movement penalty must be based on the penalty costs,
      // keep in mind that the cost given via the command line are CHANGED
      // later on!
      glp_set_obj_coef(lp, col, 20 * (d / 1000));

      vm.addVar(rowStat, col, 1);
    }
  }

  // for every edge, we define a binary variable telling us whether this edge
  // is used in a path for the original edge
  for (auto edg : edgs) {
    for (const GridNode* n : gg.getNds()) {
      for (const GridEdge* e : n->getAdjList()) {
        if (e->getFrom() != n) continue;
        if (e->pl().cost() == std::numeric_limits<float>::infinity()) continue;

        auto edgeVarName = getEdgeUseVar(e, edg);

        size_t col = glp_add_cols(lp, 1);
        glp_set_col_name(lp, col, edgeVarName.c_str());
        // binary variable € {0,1}, edge is either used, or not
        glp_set_col_kind(lp, col, GLP_BV);
        i++;

        glp_set_obj_coef(lp, col, e->pl().cost());
      }
    }
  }

  glp_create_index(lp);

  // for every node, the number of outgoing and incoming used edges must be
  // the same, except for the start and end node
  for (const GridNode* n : gg.getNds()) {
    for (auto edg : edgs) {
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
      // with an incoming edge - it would need 2 incoming edges to achieve that.
      // however, this would mean (as sink nodes are never adjacent) that 2
      // ports
      // have outgoing edges - which would mean the path "split" somewhere
      // before
      // the ports, which is impossible and forbidden by our other constraints.
      // the only way a sink node can make up for in outgoin edge
      // is thus if we add -2 if the sink is marked as the start station of this
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

    for (auto nd : nds) {
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
        for (auto edg : edgs) {
          size_t edgCol = glp_find_col(lp, getEdgeUseVar(innerE, edg).c_str());
          if (!edgCol) continue;
          vm.addVar(row, edgCol, 1);
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

    auto eOr = gg.getNEdge(n, gg.getNeighbor(x, y, 3));
    auto fOr = gg.getNEdge(gg.getNeighbor(x, y, 3), n);

    if (!eOr || !fOr) continue;

    auto na = gg.getNeighbor(x, y, (3 + 7) % 8);
    auto nb = gg.getNeighbor(x, y, (3 + 1) % 8);

    if (!na || !nb) continue;

    auto e = gg.getNEdge(na, nb);
    auto f = gg.getNEdge(nb, na);

    for (auto edg : edgs) {
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

  glp_create_index(lp);

  // for each grid edge e and edge input edge E, we define a variable
  // x_eEadjN which should be 1 if e is the adjacent grid edge for E at input node N
  for (auto edg : edgs) {
    for (GridNode* n : gg.getNds()) {
      if (!n->pl().isSink()) continue;
      for (size_t i = 0; i < 8; i++) {
        ///// for the original <from> node
        auto e = gg.getEdg(n, n->pl().getPort(i));

        std::stringstream ndPosFromVarName;
        ndPosFromVarName << "statpos(" << e->getFrom() << "," << edg->getFrom() << ")";
        size_t ndColFrom = glp_find_col(lp, ndPosFromVarName.str().c_str());

        // it may be that the node was not marked as a candidate at all
        if (!ndColFrom) continue;

        std::stringstream x_eEadjNNameFrom;
        x_eEadjNNameFrom << "adjedge(" << edg->getFrom() << "," << e << "," << edg << ")";

        size_t colFr = glp_add_cols(lp, 1);
        glp_set_col_name(lp, colFr, x_eEadjNNameFrom.str().c_str());
        glp_set_col_kind(lp, colFr, GLP_BV);

        std::stringstream constNameFr;
        size_t rowFr = glp_add_rows(lp, 1);
        constNameFr << "adjedgeconst(" << edg->getFrom() << "," << e << "," << edg << ")";
        glp_set_row_name(lp, rowFr, constNameFr.str().c_str());
        glp_set_row_bnds(lp, rowFr, GLP_DB, 0, 1);
        vm.addVar(rowFr, colFr, -2);
        size_t colEdgeUse = glp_find_col(lp, getEdgeUseVar(e, edg).c_str());
        assert(colEdgeUse);

        vm.addVar(rowFr, colEdgeUse, 1);
        vm.addVar(rowFr, ndColFrom, 1);
      }

      for (size_t i = 0; i < 8; i++) {
        ///// for the original <to> node, edges arriving at n
        auto e = gg.getEdg(n->pl().getPort(i), n);

        std::stringstream ndPosToVarName;
        ndPosToVarName << "statpos(" << e->getTo() << "," << edg->getTo() << ")";
        size_t ndColTo = glp_find_col(lp, ndPosToVarName.str().c_str());

        // it may be that the node was not marked as a candidate at all
        if (!ndColTo) continue;

        std::stringstream x_eEadjNNameTo;
        x_eEadjNNameTo << "adjedge(" << edg->getTo() << "," << e << "," << edg << ")";

        size_t colTo = glp_add_cols(lp, 1);
        glp_set_col_name(lp, colTo, x_eEadjNNameTo.str().c_str());
        glp_set_col_kind(lp, colTo, GLP_BV);

        std::stringstream constNameTo;
        size_t rowTo = glp_add_rows(lp, 1);
        constNameTo << "adjedgeconst(" << edg->getTo() << "," << e << "," << edg << ")";
        glp_set_row_name(lp, rowTo, constNameTo.str().c_str());
        glp_set_row_bnds(lp, rowTo, GLP_DB, 0, 1);
        vm.addVar(rowTo, colTo, -2);
        size_t colEdgeUse = glp_find_col(lp, getEdgeUseVar(e, edg).c_str());
        assert(colEdgeUse);

        vm.addVar(rowTo, colEdgeUse, 1);
        vm.addVar(rowTo, ndColTo, 1);
      }
    }
  }

  glp_create_index(lp);

  // for each input node N, define a var x_dirNE which tells the direction of
  // E at N
  for (auto nd : nds) {
    for (auto edg : edgs) {
      if (edg->getTo() != nd && edg->getFrom() != nd) continue;
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
            std::stringstream adjVar;
            adjVar << "adjedge(" << edg->getFrom() << "," << e << "," << edg << ")";
            size_t col = glp_find_col(lp, adjVar.str().c_str());
            if (col) vm.addVar(row, col, i);
          }
        } else {
          // the 0 can be skipped here
          for (size_t i = 1; i < 8; i++) {
            auto e = gg.getEdg(n->pl().getPort(i), n);
            std::stringstream adjVar;
            adjVar << "adjedge(" << edg->getTo() << "," << e << "," << edg << ")";
            size_t col = glp_find_col(lp, adjVar.str().c_str());
            if (col) vm.addVar(row, col, i);
          }
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
  std::string f = std::string(std::tmpnam(0)) + ".mps";
  std::string outf = "sol";//std::string(std::tmpnam(0)) + ".sol";
  glp_write_mps(lp, GLP_MPS_FILE, 0, f.c_str());
  LOG(INFO) << "Calling external solver...";

  std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();

  // std::string cmd =
      // "/home/patrick/gurobi/gurobi751/linux64/bin/gurobi_cl "
      // "ResultFile={OUTPUT} {INPUT}";
  std::string cmd =
  "/home/patrick/repos/Cbc-2.9/bin/cbc {INPUT} -randomCbcSeed 0 -threads "
  "{THREADS} -printingOptions rows -solve -solution {OUTPUT}";
  util::replaceAll(cmd, "{INPUT}", f);
  util::replaceAll(cmd, "{OUTPUT}", outf);
  util::replaceAll(cmd, "{THREADS}", "4");
  int r = system(cmd.c_str());

  std::chrono::high_resolution_clock::time_point t2 =
      std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  LOG(INFO) << " === External solve done (ret=" << r << ") in " << duration
            << " ms ===";
  LOG(INFO) << "Parsing solution...";

  std::ifstream fin;
  fin.open(outf.c_str());
  std::string line;

  // skip first line
  std::getline(fin, line);

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
  glp_init_iocp(&params);
  // params.presolve = GLP_ON;
  // params.binarize = GLP_OFF;
  // params.ps_tm_lim = 10000;
  // params.tm_lim = 30000;
  // params.fp_heur = GLP_ON;
  // params.ps_heur = GLP_ON;

  glp_simplex(lp, 0);
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
