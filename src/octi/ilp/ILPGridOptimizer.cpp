// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <glpk.h>
#include "octi/gridgraph/GridGraph.h"
#include "octi/ilp/ILPGridOptimizer.h"
#include "util/log/Log.h"

using octi::gridgraph::GridGraph;
using octi::gridgraph::GridNode;
using octi::gridgraph::GridEdge;
using octi::ilp::ILPGridOptimizer;
using octi::ilp::VariableMatrix;

const static size_t NUM_NDS = 5;
const static size_t NUM_EDG = 3;


// _____________________________________________________________________________
int ILPGridOptimizer::optimize(GridGraph* gg) const {

  for (auto nd : *gg->getNds()) {
    if (!nd->pl().isSink()) continue;
    gg->openNode(nd);
    gg->openNodeSink(nd, 0);
  }

  LOG(INFO) << "Creating ILP problem... ";
  glp_prob* lp = createProblem(*gg);
  LOG(INFO) << " .. done";

  LOG(INFO) << "(stats) ILP has " << glp_get_num_cols(lp) << " cols and "
            << glp_get_num_rows(lp) << " rows.";

  // TODO: make configurable
  glp_write_mps(lp, GLP_MPS_FILE, 0, "ilp_prob.mps");

  LOG(INFO) << "Solving problem...";
  solveProblem(lp);

  LOG(INFO) << "(stats) ILP obj = " << glp_mip_obj_val(lp);

  // write solution to grid graph
  for (GridNode* n : *gg->getNds()) {
    for (GridEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;

      for (size_t j = 0; j < NUM_EDG; j++) {
        std::string varName = getEdgeVar(e, j);

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
    for (size_t j = 0; j < NUM_NDS; j++) {
      std::stringstream ndPosVarName;
      ndPosVarName << "statpos(" << n << "," << j << ")";

      size_t i = glp_find_col(lp, ndPosVarName.str().c_str());
      if (i) {
        double val = glp_mip_col_val(lp, i);
        if (val > 0) {
          LOG(INFO) << "STAT " << j << " AT " << n->pl().getX() << ", " << n->pl().getY();
          LOG(INFO) << " (with cost "<< glp_get_obj_coef(lp, i) << ")";
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
glp_prob* ILPGridOptimizer::createProblem(const GridGraph& gg) const {
  glp_prob* lp = glp_create_prob();

  glp_set_prob_name(lp, "edgeorder");
  glp_set_obj_dir(lp, GLP_MIN);

  VariableMatrix vm;

  // three nodes, two edges, for testing
  size_t x[NUM_NDS] = {5, 9, 4, 0, 10};
  size_t y[NUM_NDS] = {5, 0, 10, 0, 8};

  // edges 0-1, 0-2, 3-4
  size_t edgs[NUM_EDG][2] = {{0, 1}, {0, 2}, {3, 4}};

  glp_create_index(lp);

  for (size_t i = 0; i < NUM_NDS; i++) {
    size_t nx = x[i];
    size_t ny = y[i];

    std::stringstream oneAssignment;
    // must sum up to 1
    size_t rowStat = glp_add_rows(lp, 1);
    oneAssignment << "oneass(st=" << i << ")";
    glp_set_row_name(lp, rowStat, oneAssignment.str().c_str());
    glp_set_row_bnds(lp, rowStat, GLP_FX, 1, 1);

    for (const GridNode* n : gg.getNds()) {
      if (!n->pl().isSink()) continue;
      std::stringstream ndPosVarName;
      ndPosVarName << "statpos(" << n << "," << i << ")";

      size_t col = glp_add_cols(lp, 1);
      glp_set_col_name(lp, col, ndPosVarName.str().c_str());
      // binary variable € {0,1}, node is either this station, or not
      glp_set_col_kind(lp, col, GLP_BV);

      // TODO: correct distance, at the moment eucludian distance based on grid
      // coord
      // TODO: the movement penalty must be based on the penalty costs,
      // keep in mind that the cost given via the command line are CHANGED
      // later on!
      glp_set_obj_coef(
          lp, col, 20 * sqrt((nx - n->pl().getX()) * (nx - n->pl().getX()) +
                             (ny - n->pl().getY()) * (ny - n->pl().getY())));

      vm.addVar(rowStat, col, 1);
    }
  }


  size_t i = 0;

  // for every edge, we define a binary variable telling us whether this edge
  // is used in a path for the original edge
  for (size_t i = 0; i < NUM_EDG; i++) {
    for (const GridNode* n : gg.getNds()) {
      for (const GridEdge* e : n->getAdjList()) {
        if (e->getFrom() != n) continue;
        if (e->pl().cost() == std::numeric_limits<float>::infinity()) continue;

        auto edgeVarName = getEdgeVar(e, i);

        size_t col = glp_add_cols(lp, 1);
        glp_set_col_name(lp, col, edgeVarName.c_str());
        // binary variable € {0,1}, edge is either used, or not
        glp_set_col_kind(lp, col, GLP_BV);

        // catch infinity cost, as 0 * inf = NaN
        glp_set_obj_coef(lp, col, e->pl().cost());
      }
    }
  }

  LOG(INFO) << i << " edges...";

  glp_create_index(lp);

  // for every node, the number of outgoing and incoming used edges must be
  // the same, except for the start and end node
  for (const GridNode* n : gg.getNds()) {
    for (size_t i = 0; i < NUM_EDG; i++) {
      std::stringstream constName;
      size_t row = glp_add_rows(lp, 1);
      constName << "adjsum(" << n << "," << i << ")";
      glp_set_row_name(lp, row, constName.str().c_str());
      glp_set_row_bnds(lp, row, GLP_FX, 0, 0);

      // normally, we count an incoming edge as 1 and an outgoing edge as -1
      // later on, we make sure that each node has a some of all out and in
      // edges of 0
      int inCost = -1;
      int outCost = 1;

      // for sink nodes, we apply a trick: an outgoing edge counts as 2 here.
      // this means that a sink node cannot make up for an outgoing edge
      // with an incoming edge - it would need 2 incoming edges to achieve that.
      // however, this would mean (as sink nodes are never adjacent) that 2 ports
      // have outgoing edges - which would mean the path "split" somewhere before
      // the ports, which is impossible and forbidden by our other constraints.
      // the only way a sink node can make up for in outgoin edge
      // is thus if we add -2 if the sink is marked as the start station of this
      // edge
      if (n->pl().isSink()) {
        // subtract the variable for this start node and edge
        std::stringstream ndPosFromVarName;
        ndPosFromVarName << "statpos(" << n << "," << edgs[i][0] << ")";
        size_t ndColFrom = glp_find_col(lp, ndPosFromVarName.str().c_str());
        assert(ndColFrom > 0);
        vm.addVar(row, ndColFrom, -2);

        // add the variable for this end node and edge
        std::stringstream ndPosToVarName;
        ndPosToVarName << "statpos(" << n << "," << edgs[i][1] << ")";
        size_t ndColTo = glp_find_col(lp, ndPosToVarName.str().c_str());
        assert(ndColTo > 0);
        vm.addVar(row, ndColTo, 1);

        outCost = 2;
      }

      for (auto e : n->getAdjListIn()) {
        size_t edgCol = glp_find_col(lp, getEdgeVar(e, i).c_str());
        if (!edgCol) continue;
        vm.addVar(row, edgCol, inCost);
      }

      for (auto e : n->getAdjListOut()) {
        size_t edgCol = glp_find_col(lp, getEdgeVar(e, i).c_str());
        if (!edgCol) continue;
        vm.addVar(row, edgCol, outCost);
      }
    }
  }

  glp_create_index(lp);

  // for every meta node, the number of outgoing edges must be 1
  // for every meta node, the number of incoming edges must be 1
  // this prevents both a crossing at used nodes and a double usage of edges
  for (const GridNode* n : gg.getNds()) {
    if (!n->pl().isSink()) continue;

    std::stringstream constName;
    size_t row = glp_add_rows(lp, 1);
    constName << "inneruse(" << n << ")";
    glp_set_row_name(lp, row, constName.str().c_str());
    glp_set_row_bnds(lp, row, GLP_UP, 0, 1);

    // a meta grid node can either be a sink for a single input node, or
    // a pass-through

    for (size_t i = 0; i < NUM_NDS; i++) {
      std::stringstream ndPosToVarName;
      ndPosToVarName << "statpos(" << n << "," << i << ")";
      size_t ndColTo = glp_find_col(lp, ndPosToVarName.str().c_str());
      assert(ndColTo > 0);
      vm.addVar(row, ndColTo, 1);
    }

    // go over all ports
    for (size_t pf = 0; pf < 8; pf++) {
      auto from = n->pl().getPort(pf);
      for (size_t pt = 0; pt < 8; pt++) {
        auto to = n->pl().getPort(pt);
        if (from == to) continue;

        auto innerE = gg.getEdg(from, to);
        for (size_t i = 0; i < NUM_EDG; i++) {
          size_t edgCol = glp_find_col(lp, getEdgeVar(innerE, i).c_str());
          if (!edgCol) continue;
          vm.addVar(row, edgCol, 1);
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
void ILPGridOptimizer::solveProblem(glp_prob* lp) const {
  glp_iocp params;
  glp_init_iocp(&params);
  params.presolve = GLP_ON;
  params.binarize = GLP_ON;
  // params.ps_tm_lim = _cfg->glpkPSTimeLimit;
  // params.tm_lim = _cfg->glpkTimeLimit;
  params.fp_heur = GLP_ON;
  params.ps_heur = GLP_ON;

  glp_intopt(lp, &params);
}

// _____________________________________________________________________________
// void ILPOptimizer::preSolveCoinCbc(glp_prob* lp) const {
// // write temporary file
// std::string f = std::string(std::tmpnam(0)) + ".mps";
// std::string outf = std::string(std::tmpnam(0)) + ".sol";
// glp_write_mps(lp, GLP_MPS_FILE, 0, f.c_str());
// LOG(INFO) << "Calling external solver..." << std::endl;

// std::chrono::high_resolution_clock::time_point t1 =
// std::chrono::high_resolution_clock::now();

// std::string cmd = _cfg->externalSolver;
// util::replaceAll(cmd, "{INPUT}", f);
// util::replaceAll(cmd, "{OUTPUT}", outf);
// util::replaceAll(cmd, "{THREADS}",
// util::toString(std::thread::hardware_concurrency()));
// LOG(INFO) << "Cmd: '" << cmd << "'" << std::endl;
// int r = system(cmd.c_str());

// std::chrono::high_resolution_clock::time_point t2 =
// std::chrono::high_resolution_clock::now();
// auto duration =
// std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
// LOG(INFO) << " === External solve done (ret=" << r << ") in " << duration
// << " ms ===" << std::endl;
// LOG(INFO) << "Parsing solution..." << std::endl;

// std::ifstream fin;
// fin.open(outf.c_str());
// std::string line;

// // skip first line
// std::getline(fin, line);

// while (std::getline(fin, line)) {
// std::istringstream iss(line);
// int number;
// string name;
// double value;

// iss >> number;

// // could not read line number, is missing, which is fine for us
// if (iss.fail()) iss.clear();

// iss >> name;
// iss >> value;

// int intVal = value;

// size_t col = glp_find_col(lp, name.c_str());
// if (col != 0) {
// glp_set_col_bnds(lp, col, GLP_FX, intVal, intVal);
// }
// }
// }

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
std::string ILPGridOptimizer::getEdgeVar(const GridEdge* e, size_t i) const {
  std::stringstream varName;
  varName << "edg(" << e->getFrom() << "," << e->getTo() << ","
          << std::to_string(i) << ")";

  return varName.str();
}
