// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <glpk.h>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <thread>
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/optim/ILPOptimizer.h"
#include "transitmap/optim/OptGraph.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using namespace transitmapper;
using namespace optim;
using namespace transitmapper::graph;

// _____________________________________________________________________________
int ILPOptimizer::optimize(TransitGraph* tg) const {
  // create optim graph
  OptGraph g(tg);

  if (_cfg->createCoreOptimGraph) {
    g.simplify();
    g.untangle();
    g.split();
  }

  if (_cfg->outOptGraph) {
    util::geo::output::GeoGraphJsonOutput out;
    std::ofstream fstr(_cfg->dbgPath + "/optgraph.json");
    out.print(g, fstr);
  }

  HierarchOrderingConfig hc;
  OrderingConfig c;

  optimize(*g.getNds(), &hc);

  hc.writeFlatCfg(&c);

  Optimizer::expandRelatives(tg, &c);

  tg->setConfig(c);

  return 0;
}

// _____________________________________________________________________________
int ILPOptimizer::optimize(const std::set<OptNode*>& g,
                           HierarchOrderingConfig* hc) const {
  LOG(DEBUG) << "Creating ILP problem... " << std::endl;
  glp_prob* lp = createProblem(g);
  LOG(DEBUG) << " .. done" << std::endl;

  LOG(INFO) << "(stats) ILP has " << glp_get_num_cols(lp) << " cols and "
            << glp_get_num_rows(lp) << " rows." << std::endl;

  if (!_cfg->glpkHOutputPath.empty()) {
    LOG(DEBUG) << "Writing human readable ILP to '" << _cfg->glpkHOutputPath
               << "'" << std::endl;
    printHumanReadable(lp, _cfg->glpkHOutputPath.c_str());
  }

  if (!_cfg->glpkMPSOutputPath.empty()) {
    LOG(DEBUG) << "Writing ILP as .mps to '" << _cfg->glpkMPSOutputPath << "'"
               << std::endl;
    glp_write_mps(lp, GLP_MPS_FILE, 0, _cfg->glpkMPSOutputPath.c_str());
  }

  LOG(DEBUG) << "Solving problem..." << std::endl;

  if (!_cfg->externalSolver.empty()) {
    preSolveCoinCbc(lp);
  }

  std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();
  solveProblem(lp);
  std::chrono::high_resolution_clock::time_point t2 =
      std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (_cfg->externalSolver.empty()) {
    LOG(INFO) << " === Solve done in " << duration << " ms ===" << std::endl;
  }

  LOG(INFO) << "(stats) ILP obj = " << glp_mip_obj_val(lp) << std::endl;

  if (!_cfg->glpkSolutionOutputPath.empty()) {
    LOG(DEBUG) << "Writing ILP full solution to '"
               << _cfg->glpkSolutionOutputPath << "'" << std::endl;
    glp_print_mip(lp, _cfg->glpkSolutionOutputPath.c_str());
  }

  getConfigurationFromSolution(lp, hc, g);

  glp_delete_prob(lp);
  glp_free_env();

  return 0;
}

// _____________________________________________________________________________
int ILPOptimizer::getCrossingPenaltySameSeg(const OptNode* n) const {
  return _scorer->getCrossingPenaltySameSeg(n->pl().node);
}

// _____________________________________________________________________________
int ILPOptimizer::getCrossingPenaltyDiffSeg(const OptNode* n) const {
  return _scorer->getCrossingPenaltyDiffSeg(n->pl().node);
}

// _____________________________________________________________________________
int ILPOptimizer::getSplittingPenalty(const OptNode* n) const {
  // double the value because we only count a splitting once for each pair!
  return _scorer->getSplittingPenalty(n->pl().node) * 1;
}

// _____________________________________________________________________________
double ILPOptimizer::getConstraintCoeff(glp_prob* lp, int constraint,
                                        int col) const {
  int indices[glp_get_num_cols(lp)];
  double values[glp_get_num_cols(lp)];
  glp_get_mat_col(lp, col, indices, values);

  for (int i = 1; i < glp_get_num_cols(lp); i++) {
    if (indices[i] == constraint) return values[i];
  }

  return 0;
}

// _____________________________________________________________________________
void ILPOptimizer::getConfigurationFromSolution(
    glp_prob* lp, HierarchOrderingConfig* hc, const std::set<OptNode*>& g) const {
  // build name index for faster lookup
  glp_create_index(lp);

  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      for (auto etgp : e->pl().etgs) {
        if (etgp.wasCut) continue;
        for (size_t tp = 0; tp < etgp.etg->getCardinality(true); tp++) {
          bool found = false;
          for (size_t p = 0; p < etgp.etg->getCardinality(); p++) {
            auto r = (*etgp.etg->getRoutes())[p];
            if (r.route->relativeTo()) continue;
            std::string varName = getILPVarName(e, r.route, tp);

            size_t i = glp_find_col(lp, varName.c_str());
            assert(i > 0);
            double val = glp_mip_col_val(lp, i);

            if (val > 0.5) {
              if (!(etgp.dir ^ e->pl().etgs.front().dir)) {
                (*hc)[etgp.etg][etgp.order].insert((*hc)[etgp.etg][etgp.order].begin(), p);
              } else {
                (*hc)[etgp.etg][etgp.order].push_back(p);
              }
              assert(!found);  // should be assured by ILP constraints
              found = true;
            }
          }
          assert(found);
        }
      }
    }
  }
}


// _____________________________________________________________________________
glp_prob* ILPOptimizer::createProblem(const std::set<OptNode*>& g) const {
  glp_prob* lp = glp_create_prob();

  glp_set_prob_name(lp, "edgeorder");
  glp_set_obj_dir(lp, GLP_MIN);

  VariableMatrix vm;

  // for every segment s, we define |L(s)|^2 decision variables x_slp
  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      // get string repr of etg

      size_t newCols = e->pl().getCardinality() * e->pl().getCardinality();
      size_t cols = glp_add_cols(lp, newCols);
      size_t i = 0;
      size_t rowA = glp_add_rows(lp, e->pl().getCardinality());

      for (size_t p = 0; p < e->pl().getCardinality(); p++) {
        std::stringstream varName;

        varName << "sum(" << e->pl().getStrRepr() << ",p=" << p << ")";
        glp_set_row_name(lp, rowA + p, varName.str().c_str());
        glp_set_row_bnds(lp, rowA + p, GLP_FX, 1, 1);
      }

      for (auto r : e->pl().getRoutes()) {
        if (r.route->relativeTo()) continue;
        // constraint: the sum of all x_slp over p must be 1 for equal sl
        size_t row = glp_add_rows(lp, 1);
        std::stringstream varName;
        varName << "sum(" << e->pl().getStrRepr() << ",l=" << r.route << ")";
        glp_set_row_name(lp, row, varName.str().c_str());
        glp_set_row_bnds(lp, row, GLP_FX, 1, 1);

        for (size_t p = 0; p < e->pl().getCardinality(); p++) {
          std::string varName = getILPVarName(e, r.route, p);
          size_t curCol = cols + i;
          glp_set_col_name(lp, curCol, varName.c_str());

          // binary variable € {0,1}
          glp_set_col_kind(lp, curCol, GLP_BV);

          vm.addVar(row, curCol, 1);
          vm.addVar(rowA + p, curCol, 1);

          i++;
        }
      }
    }
  }

  glp_create_index(lp);

  writeSameSegConstraints(g, &vm, lp);
  writeDiffSegConstraints(g, &vm, lp);

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
void ILPOptimizer::writeSameSegConstraints(const std::set<OptNode*>& g,
                                           VariableMatrix* vm,
                                           glp_prob* lp) const {
  // go into nodes and build crossing constraints for adjacent
  for (OptNode* node : g) {
    std::set<OptEdge*> processed;
    for (OptEdge* segmentA : node->getAdjList()) {
      processed.insert(segmentA);
      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segmentA)) {
        // iterate over all edges this
        // pair traverses to _TOGETHER_
        // (its possible that there are multiple edges if a line continues
        //  in more then 1 segment)
        for (OptEdge* segmentB : getEdgePartners(node, segmentA, linepair)) {
          if (processed.find(segmentB) != processed.end()) continue;
          // try all position combinations
          size_t decisionVar = glp_add_cols(lp, 1);

          // introduce dec var
          std::stringstream ss;
          ss << "x_dec(" << segmentA->pl().getStrRepr() << ","
             << segmentB->pl().getStrRepr() << "," << linepair.first << "("
             << linepair.first->getId() << ")," << linepair.second << "("
             << linepair.second->getId() << ")," << node << ")";
          glp_set_col_name(lp, decisionVar, ss.str().c_str());
          glp_set_col_kind(lp, decisionVar, GLP_BV);
          glp_set_obj_coef(
              lp, decisionVar,
              getCrossingPenaltySameSeg(node)
                  // multiply the penalty with the number of collapsed lines!
                  * (linepair.first->getNumCollapsedPartners() + 1) *
                  (linepair.second->getNumCollapsedPartners() + 1));

          for (PosComPair poscomb :
               getPositionCombinations(segmentA, segmentB)) {
            if (crosses(node, segmentA, segmentB, poscomb)) {
              size_t lineAinAatP = glp_find_col(
                  lp,
                  getILPVarName(segmentA, linepair.first, poscomb.first.first)
                      .c_str());
              size_t lineBinAatP = glp_find_col(
                  lp,
                  getILPVarName(segmentA, linepair.second, poscomb.second.first)
                      .c_str());
              size_t lineAinBatP = glp_find_col(
                  lp,
                  getILPVarName(segmentB, linepair.first, poscomb.first.second)
                      .c_str());
              size_t lineBinBatP =
                  glp_find_col(lp, getILPVarName(segmentB, linepair.second,
                                                 poscomb.second.second)
                                       .c_str());

              assert(lineAinAatP > 0);
              assert(lineAinBatP > 0);
              assert(lineBinAatP > 0);
              assert(lineBinBatP > 0);

              size_t row = glp_add_rows(lp, 1);
              std::stringstream ss;
              ss << "dec_sum(" << segmentA->pl().getStrRepr() << ","
                 << segmentB->pl().getStrRepr() << "," << linepair.first << ","
                 << linepair.second << "pa=" << poscomb.first.first
                 << ",pb=" << poscomb.second.first
                 << ",pa'=" << poscomb.first.second
                 << ",pb'=" << poscomb.second.second << ",n=" << node << ")";
              glp_set_row_name(lp, row, ss.str().c_str());
              glp_set_row_bnds(lp, row, GLP_UP, 0, 3);

              vm->addVar(row, lineAinAatP, 1);
              vm->addVar(row, lineBinAatP, 1);
              vm->addVar(row, lineAinBatP, 1);
              vm->addVar(row, lineBinBatP, 1);
              vm->addVar(row, decisionVar, -1);
            }
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
void ILPOptimizer::writeDiffSegConstraints(const std::set<OptNode*>& g,
                                           VariableMatrix* vm,
                                           glp_prob* lp) const {
  // go into nodes and build crossing constraints for adjacent
  for (OptNode* node : g) {
    std::set<OptEdge*> processed;
    for (OptEdge* segmentA : node->getAdjList()) {
      processed.insert(segmentA);
      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segmentA)) {
        for (EdgePair segments :
             getEdgePartnerPairs(node, segmentA, linepair)) {
          // try all position combinations
          size_t decisionVar = glp_add_cols(lp, 1);

          // introduce dec var
          std::stringstream ss;
          ss << "x_dec(" << segmentA->pl().getStrRepr() << ","
             << segments.first->pl().getStrRepr()
             << segments.second->pl().getStrRepr() << "," << linepair.first
             << "(" << linepair.first->getId() << ")," << linepair.second << "("
             << linepair.second->getId() << ")," << node << ")";
          glp_set_col_name(lp, decisionVar, ss.str().c_str());
          glp_set_col_kind(lp, decisionVar, GLP_BV);
          glp_set_obj_coef(
              lp, decisionVar,
              getCrossingPenaltyDiffSeg(node)
                  // multiply the penalty with the number of collapsed lines!
                  * (linepair.first->getNumCollapsedPartners() + 1) *
                  (linepair.second->getNumCollapsedPartners() + 1));

          for (PosCom poscomb : getPositionCombinations(segmentA)) {
            if (crosses(node, segmentA, segments, poscomb)) {
              size_t lineAinAatP = glp_find_col(
                  lp, getILPVarName(segmentA, linepair.first, poscomb.first)
                          .c_str());
              size_t lineBinAatP = glp_find_col(
                  lp, getILPVarName(segmentA, linepair.second, poscomb.second)
                          .c_str());

              assert(lineAinAatP > 0);
              assert(lineBinAatP > 0);

              size_t row = glp_add_rows(lp, 1);
              std::stringstream ss;
              ss << "dec_sum(" << segmentA->pl().getStrRepr() << ","
                 << segments.first->pl().getStrRepr()
                 << segments.second->pl().getStrRepr() << "," << linepair.first
                 << "," << linepair.second << "pa=" << poscomb.first
                 << ",pb=" << poscomb.second << ",n=" << node << ")";
              glp_set_row_name(lp, row, ss.str().c_str());
              glp_set_row_bnds(lp, row, GLP_UP, 0, 1);

              vm->addVar(row, lineAinAatP, 1);
              vm->addVar(row, lineBinAatP, 1);
              vm->addVar(row, decisionVar, -1);
            }
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
std::vector<PosComPair> ILPOptimizer::getPositionCombinations(
    OptEdge* a, OptEdge* b) const {
  std::vector<PosComPair> ret;
  for (size_t posLineAinA = 0; posLineAinA < a->pl().getCardinality();
       posLineAinA++) {
    for (size_t posLineBinA = 0; posLineBinA < a->pl().getCardinality();
         posLineBinA++) {
      if (posLineAinA == posLineBinA) continue;

      for (size_t posLineAinB = 0; posLineAinB < b->pl().getCardinality();
           posLineAinB++) {
        for (size_t posLineBinB = 0; posLineBinB < b->pl().getCardinality();
             posLineBinB++) {
          if (posLineAinB == posLineBinB) continue;

          ret.push_back(PosComPair(PosCom(posLineAinA, posLineAinB),
                                   PosCom(posLineBinA, posLineBinB)));
        }
      }
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::vector<PosCom> ILPOptimizer::getPositionCombinations(OptEdge* a) const {
  std::vector<PosCom> ret;
  for (size_t posLineAinA = 0; posLineAinA < a->pl().getCardinality();
       posLineAinA++) {
    for (size_t posLineBinA = 0; posLineBinA < a->pl().getCardinality();
         posLineBinA++) {
      if (posLineAinA == posLineBinA) continue;
      ret.push_back(PosCom(posLineAinA, posLineBinA));
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::string ILPOptimizer::getILPVarName(OptEdge* seg, const Route* r,
                                        size_t p) const {
  std::stringstream varName;
  varName << "x_(" << seg->pl().getStrRepr() << ",l=" << r << ",p=" << p << ")";
  return varName.str();
}

// _____________________________________________________________________________
std::vector<OptEdge*> ILPOptimizer::getEdgePartners(
    OptNode* node, OptEdge* segmentA, const LinePair& linepair) const {
  std::vector<OptEdge*> ret;

  graph::Edge* fromEtg = OptGraph::getAdjEdg(segmentA, node);
  const Node* dirA = fromEtg->getRoute(linepair.first)->direction;
  const Node* dirB = fromEtg->getRoute(linepair.second)->direction;

  for (OptEdge* segmentB : node->getAdjList()) {
    if (segmentB == segmentA) continue;
    graph::Edge* e = OptGraph::getAdjEdg(segmentB, node);

    if (e->getCtdRoutesIn(node->pl().node, linepair.first, dirA, fromEtg)
            .size() &&
        e->getCtdRoutesIn(node->pl().node, linepair.second, dirB, fromEtg)
            .size()) {
      ret.push_back(segmentB);
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::vector<EdgePair> ILPOptimizer::getEdgePartnerPairs(
    OptNode* node, OptEdge* segmentA, const LinePair& linepair) const {
  std::vector<EdgePair> ret;

  graph::Edge* fromEtg = OptGraph::getAdjEdg(segmentA, node);
  const Node* dirA = fromEtg->getRoute(linepair.first)->direction;
  const Node* dirB = fromEtg->getRoute(linepair.second)->direction;

  for (OptEdge* segmentB : node->getAdjList()) {
    if (segmentB == segmentA) continue;
    graph::Edge* etg = OptGraph::getAdjEdg(segmentB, node);

    if (etg->getCtdRoutesIn(node->pl().node, linepair.first, dirA, fromEtg)
            .size()) {
      EdgePair curPair;
      curPair.first = segmentB;
      for (OptEdge* segmentC : node->getAdjList()) {
        if (segmentC == segmentA || segmentC == segmentB) continue;
        graph::Edge* e = OptGraph::getAdjEdg(segmentC, node);

        if (e->getCtdRoutesIn(node->pl().node, linepair.second, dirB, fromEtg)
                .size()) {
          curPair.second = segmentC;
          ret.push_back(curPair);
        }
      }
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::vector<LinePair> ILPOptimizer::getLinePairs(OptEdge* segment) const {
  return getLinePairs(segment, false);
}

// _____________________________________________________________________________
std::vector<LinePair> ILPOptimizer::getLinePairs(OptEdge* segment,
                                                 bool unique) const {
  std::set<const Route*> processed;
  std::vector<LinePair> ret;
  for (auto& toA : segment->pl().getRoutes()) {
    if (toA.route->relativeTo()) continue;
    processed.insert(toA.route);
    for (auto& toB : segment->pl().getRoutes()) {
      if (unique && processed.find(toB.route) != processed.end()) continue;
      if (toB.route->relativeTo()) continue;
      if (toA.route == toB.route) continue;

      // this is to make sure that we always get the same line pairs in unique
      // mode -> take the smaller pointer first
      if (!unique || toA.route < toB.route)
        ret.push_back(LinePair(toA.route, toB.route));
      else
        ret.push_back(LinePair(toB.route, toA.route));
    }
  }
  return ret;
}

// _____________________________________________________________________________
void ILPOptimizer::solveProblem(glp_prob* lp) const {
  glp_iocp params;
  glp_init_iocp(&params);
  params.presolve = GLP_ON;
  params.binarize = GLP_ON;
  params.ps_tm_lim = _cfg->glpkPSTimeLimit;
  params.tm_lim = _cfg->glpkTimeLimit;
  if (_cfg->externalSolver.empty()) {
    params.fp_heur = _cfg->useGlpkFeasibilityPump ? GLP_ON : GLP_OFF;
    params.ps_heur = _cfg->useGlpkProximSearch ? GLP_ON : GLP_OFF;
  }

  glp_intopt(lp, &params);
}

// _____________________________________________________________________________
void ILPOptimizer::preSolveCoinCbc(glp_prob* lp) const {
  // write temporary file
  std::string f = std::string(std::tmpnam(0)) + ".mps";
  std::string outf = std::string(std::tmpnam(0)) + ".sol";
  glp_write_mps(lp, GLP_MPS_FILE, 0, f.c_str());
  LOG(INFO) << "Calling external solver..." << std::endl;

  std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();

  std::string cmd = _cfg->externalSolver;
  util::replaceAll(cmd, "{INPUT}", f);
  util::replaceAll(cmd, "{OUTPUT}", outf);
  util::replaceAll(cmd, "{THREADS}",
                   util::toString(std::thread::hardware_concurrency()));
  LOG(INFO) << "Cmd: '" << cmd << "'" << std::endl;
  int r = system(cmd.c_str());

  std::chrono::high_resolution_clock::time_point t2 =
      std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  LOG(INFO) << " === External solve done (ret=" << r << ") in " << duration
            << " ms ===" << std::endl;
  LOG(INFO) << "Parsing solution..." << std::endl;

  std::ifstream fin;
  fin.open(outf.c_str());
  std::string line;

  // skip first line
  std::getline(fin, line);

  while (std::getline(fin, line)) {
    std::istringstream iss(line);
    int number;
    string name;
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
bool ILPOptimizer::crosses(OptNode* node, OptEdge* segmentA, OptEdge* segmentB,
                           PosComPair poscomb) const {
  bool otherWayA =
      (segmentA->getFrom() != node) ^ segmentA->pl().etgs.front().dir;
  bool otherWayB =
      (segmentB->getFrom() != node) ^ segmentB->pl().etgs.front().dir;

  size_t cardA = segmentA->pl().etgs.front().etg->getCardinality(true);
  size_t cardB = segmentB->pl().etgs.front().etg->getCardinality(true);

  size_t posAinA =
      otherWayA ? cardA - 1 - poscomb.first.first : poscomb.first.first;
  size_t posAinB =
      otherWayB ? cardB - 1 - poscomb.first.second : poscomb.first.second;
  size_t posBinA =
      otherWayA ? cardA - 1 - poscomb.second.first : poscomb.second.first;
  size_t posBinB =
      otherWayB ? cardB - 1 - poscomb.second.second : poscomb.second.second;

  bool pCrossing = false;

  DPoint aInA = getPos(node, segmentA, posAinA);
  DPoint bInA = getPos(node, segmentA, posBinA);
  DPoint aInB = getPos(node, segmentB, posAinB);
  DPoint bInB = getPos(node, segmentB, posBinB);

  DLine a;
  a.push_back(aInA);
  a.push_back(aInB);

  DLine b;
  b.push_back(bInA);
  b.push_back(bInB);

  if (util::geo::intersects(aInA, aInB, bInA, bInB) ||
      util::geo::dist(a, b) < 1)
    pCrossing = true;

  // pCrossing = (posAinA < posBinA && posAinB < posBinB);

  return pCrossing;
}

// _____________________________________________________________________________
bool ILPOptimizer::crosses(OptNode* node, OptEdge* segmentA, EdgePair segments,
                           PosCom postcomb) const {
  bool otherWayA =
      (segmentA->getFrom() != node) ^ segmentA->pl().etgs.front().dir;
  bool otherWayB = (segments.first->getFrom() != node) ^
                   segments.first->pl().etgs.front().dir;
  bool otherWayC = (segments.second->getFrom() != node) ^
                   segments.second->pl().etgs.front().dir;

  size_t cardA = segmentA->pl().etgs.front().etg->getCardinality(true);
  size_t cardB = segments.first->pl().etgs.front().etg->getCardinality(true);
  size_t cardC = segments.second->pl().etgs.front().etg->getCardinality(true);

  size_t posAinA = otherWayA ? cardA - 1 - postcomb.first : postcomb.first;
  size_t posBinA = otherWayA ? cardA - 1 - postcomb.second : postcomb.second;

  DPoint aInA = getPos(node, segmentA, posAinA);
  DPoint bInA = getPos(node, segmentA, posBinA);

  for (size_t i = 0;
       i < segments.first->pl().etgs.front().etg->getCardinality(true); ++i) {
    for (size_t j = 0;
         j < segments.second->pl().etgs.front().etg->getCardinality(true);
         ++j) {
      size_t posAinB = otherWayB ? cardB - 1 - i : i;
      size_t posBinC = otherWayC ? cardC - 1 - j : j;

      DPoint aInB = getPos(node, segments.first, posAinB);
      DPoint bInC = getPos(node, segments.second, posBinC);

      DLine a;
      a.push_back(aInA);
      a.push_back(aInB);

      DLine b;
      b.push_back(bInA);
      b.push_back(bInC);

      if (util::geo::intersects(aInA, aInB, bInA, bInC) ||
          util::geo::dist(a, b) < 1)
        return true;
    }
  }
  return false;
}

// _____________________________________________________________________________
DPoint ILPOptimizer::getPos(OptNode* n, OptEdge* segment, size_t p) const {
  // look for correct nodefront
  const NodeFront* nf = 0;
  for (auto etg : segment->pl().etgs) {
    const NodeFront* test = n->pl().node->getNodeFrontFor(etg.etg);
    if (test) {
      nf = test;
      break;
    }
  }

  if (!nf)
    std::cerr << n->pl().node->getId()
              << " node fronts: " << n->pl().node->getMainDirs().size()
              << std::endl;
  assert(nf);

  return nf->getTripPos(segment->pl().etgs.front().etg, p, false);
}

// _____________________________________________________________________________
bool ILPOptimizer::printHumanReadable(glp_prob* lp,
                                      const std::string& path) const {
  double maxDblDst = 0.000001;
  std::ofstream file;
  file.open(path);
  if (!file) return false;

  // get objective function
  std::stringstream obj;

  for (int i = 0; i < glp_get_num_cols(lp) - 1; ++i) {
    std::string colName = i > 0 ? glp_get_col_name(lp, i) : "";
    double coef = glp_get_obj_coef(lp, i);
    if (fabs(coef) > maxDblDst) {
      if (coef > maxDblDst && obj.str().size() > 0)
        obj << " + ";
      else if (coef < 0 && obj.str().size() > 0)
        obj << " - ";
      if (fabs(fabs(coef) - 1) > maxDblDst) {
        obj << (obj.str().size() > 0 ? fabs(coef) : coef) << " ";
      }
      obj << colName;
    }
  }

  if (glp_get_obj_dir(lp) == GLP_MIN) {
    file << "min ";
  } else {
    file << "max ";
  }
  file << obj.str() << std::endl;

  for (int j = 1; j < glp_get_num_rows(lp) - 1; ++j) {
    std::stringstream row;
    for (int i = 1; i < glp_get_num_cols(lp) - 1; ++i) {
      std::string colName = glp_get_col_name(lp, i);
      double coef = getConstraintCoeff(lp, j, i);
      if (fabs(coef) > maxDblDst) {
        if (coef > maxDblDst && row.str().size() > 0)
          row << " + ";
        else if (coef < 0 && row.str().size() > 0)
          row << " - ";
        if (fabs(fabs(coef) - 1) > maxDblDst) {
          row << (row.str().size() > 0 ? fabs(coef) : coef) << " ";
        }
        row << colName;
      }
    }

    switch (glp_get_row_type(lp, j)) {
      case GLP_FR:
        break;
      case GLP_LO:
        file << std::endl << row.str() << " >= " << glp_get_row_lb(lp, j);
        break;
      case GLP_UP:
        file << std::endl << row.str() << " <= " << glp_get_row_ub(lp, j);
        break;
      case GLP_DB:
        file << std::endl << row.str() << " >= " << glp_get_row_lb(lp, j);
        file << std::endl << row.str() << " <= " << glp_get_row_ub(lp, j);
        break;
      case GLP_FX:
        file << std::endl << row.str() << " = " << glp_get_row_lb(lp, j);
        break;
    }
  }

  file.close();
  return true;
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
