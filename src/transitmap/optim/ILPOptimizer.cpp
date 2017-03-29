// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <glpk.h>
#include <fstream>
#include "./../graph/OrderingConfiguration.h"
#include "./../output/OgrOutput.h"
#include "./ILPOptimizer.h"
#include "./OptGraph.h"
#include "pbutil/geo/Geo.h"
#include "pbutil/log/Log.h"

using namespace transitmapper;
using namespace optim;
using namespace graph;

// _____________________________________________________________________________
int ILPOptimizer::optimize() const {
  // create optim graph
  OptGraph g(_g);
  g.simplify();

  // output::OgrOutput ogrOut("/home/patrick/optimgraph", _cfg);
  // ogrOut.print(g);

  LOG(DEBUG) << "Creating ILP problem... " << std::endl;
  glp_prob* lp = createProblem(g);
  LOG(DEBUG) << " .. done" << std::endl;

  // printHumanReadable(lp, "/home/patrick/optim.ilp");

  // write problem for debugging...
  // glp_write_mps(lp, GLP_MPS_FILE, 0, "/home/patrick/ilp");

  solveProblem(lp);

  LOG(DEBUG) << "ILP obj = " << glp_mip_obj_val(lp) << std::endl;

  glp_print_mip(lp, "/home/patrick/ilp.sol");

  Configuration c;
  getConfigurationFromSolution(lp, &c, g);
  _g->setConfig(c);

  glp_delete_prob(lp);
  glp_free_env();

  return 0;
}

// _____________________________________________________________________________
double ILPOptimizer::getConstraintCoeff(glp_prob* lp, size_t constraint,
                                        size_t col) const {
  int indices[glp_get_num_cols(lp)];
  double values[glp_get_num_cols(lp)];
  glp_get_mat_col(lp, col, indices, values);

  for (size_t i = 1; i < glp_get_num_cols(lp); i++) {
    if (indices[i] == constraint) return values[i];
  }

  return 0;
}

// _____________________________________________________________________________
bool ILPOptimizer::printHumanReadable(glp_prob* lp,
                                      const std::string& path) const {
  double maxDblDst = 0.000001;
  std::ofstream myfile;
  myfile.open(path);
  if (!myfile) return false;

  // get objective function
  std::stringstream obj;

  for (size_t i = 0; i < glp_get_num_cols(lp) - 1; ++i) {
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
    myfile << "min ";
  } else {
    myfile << "max ";
  }
  myfile << obj.str() << std::endl;

  for (size_t j = 1; j < glp_get_num_rows(lp) - 1; ++j) {
    std::stringstream row;
    for (size_t i = 1; i < glp_get_num_cols(lp) - 1; ++i) {
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
        myfile << std::endl << row.str() << " >= " << glp_get_row_lb(lp, j);
        break;
      case GLP_UP:
        myfile << std::endl << row.str() << " <= " << glp_get_row_ub(lp, j);
        break;
      case GLP_DB:
        myfile << std::endl << row.str() << " >= " << glp_get_row_lb(lp, j);
        myfile << std::endl << row.str() << " <= " << glp_get_row_ub(lp, j);
        break;
      case GLP_FX:
        myfile << std::endl << row.str() << " = " << glp_get_row_lb(lp, j);
        break;
    }
  }

  myfile.close();
}

// _____________________________________________________________________________
void ILPOptimizer::getConfigurationFromSolution(glp_prob* lp, Configuration* c,
                                                const OptGraph& g) const {
  // build name index for faster lookup
  glp_create_index(lp);

  for (OptNode* n : g.getNodes()) {
    for (OptEdge* e : n->adjListOut) {
      for (auto etgp : e->etgs) {
        for (size_t tp = 0; tp < etgp.etg->getCardinality(); tp++) {
          bool found = false;
          for (size_t p = 0; p < etgp.etg->getCardinality(); p++) {
            auto r = (*etgp.etg->getTripsUnordered())[p];
            std::stringstream varName;
            varName << "x_(" << e->getStrRepr() << ",l=" << r.route
                    << ",p=" << tp << ")";

            size_t i = glp_find_col(lp, varName.str().c_str());
            assert(i > 0);
            double val = glp_mip_col_val(lp, i);

            if (val > 0.5) {
              if (!(etgp.dir ^ e->etgs.front().dir)) {
                (*c)[etgp.etg].insert((*c)[etgp.etg].begin(), p);
              } else {
                (*c)[etgp.etg].push_back(p);
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
glp_prob* ILPOptimizer::createProblem(const OptGraph& g) const {
  glp_prob* lp = glp_create_prob();

  glp_set_prob_name(lp, "edgeorder");
  glp_set_obj_dir(lp, GLP_MIN);

  size_t c = 0;

  // TODO: array sizes
  int* ia = new int[1000000];
  int* ja = new int[1000000];
  double* res = new double[1000000];

  // for every segment s, we define |L(s)|^2 decision variables x_slp
  for (OptNode* n : g.getNodes()) {
    for (OptEdge* e : n->adjListOut) {
      // the first stored etg is always the ref
      graph::Edge* etg = e->etgs[0].etg;

      // get string repr of etg

      size_t newCols = etg->getCardinality() * etg->getCardinality();
      size_t cols = glp_add_cols(lp, newCols);
      size_t i = 0;
      size_t rowA = glp_add_rows(lp, etg->getCardinality());

      for (size_t p = 0; p < etg->getCardinality(); p++) {
        std::stringstream varName;

        varName << "sum(" << e->getStrRepr() << ",p=" << p << ")";
        glp_set_row_name(lp, rowA + p, varName.str().c_str());
        glp_set_row_bnds(lp, rowA + p, GLP_FX, 1, 1);
      }

      for (auto r : *etg->getTripsUnordered()) {
        // constraint: the sum of all x_slp over p must be 1 for equal sl
        size_t row = glp_add_rows(lp, 1);
        std::stringstream varName;
        varName << "sum(" << e->getStrRepr() << ",l=" << r.route << ")";
        glp_set_row_name(lp, row, varName.str().c_str());
        glp_set_row_bnds(lp, row, GLP_FX, 1, 1);

        for (size_t p = 0; p < etg->getCardinality(); p++) {
          std::stringstream varName;
          varName << "x_(" << e->getStrRepr() << ",l=" << r.route << ",p=" << p
                  << ")";
          size_t curCol = cols + i;
          glp_set_col_name(lp, curCol, varName.str().c_str());

          // binary variable € {0,1}
          glp_set_col_kind(lp, curCol, GLP_BV);
          c++;

          ia[c] = row;
          ja[c] = curCol;
          res[c] = 1;

          c++;

          ia[c] = rowA + p;
          ja[c] = curCol;
          res[c] = 1;

          i++;
        }
      }
    }
  }

  glp_create_index(lp);

  writeSameSegConstraints(g, ia, ja, res, &c, lp);
  writeDiffSegConstraints(g, ia, ja, res, &c, lp);

  glp_load_matrix(lp, c, ia, ja, res);

  delete[](ia);
  delete[](ja);
  delete[](res);

  return lp;
}

// _____________________________________________________________________________
void ILPOptimizer::writeSameSegConstraints(const OptGraph& g, int* ia, int* ja,
                                           double* res, size_t* c,
                                           glp_prob* lp) const {
  // go into nodes and build crossing constraints for adjacent
  for (OptNode* node : g.getNodes()) {
    std::set<OptEdge*> processed;
    for (OptEdge* segmentA : node->adjList) {
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
          ss << "x_dec(" << segmentA->getStrRepr() << ","
             << segmentB->getStrRepr() << "," << linepair.first << "("
             << linepair.first->id << ")," << linepair.second << "("
             << linepair.second->id << ")," << node << ")";
          glp_set_col_name(lp, decisionVar, ss.str().c_str());
          glp_set_col_kind(lp, decisionVar, GLP_BV);
          glp_set_obj_coef(lp, decisionVar, 1);

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
              ss << "dec_sum(" << segmentA->getStrRepr() << ","
                 << segmentB->getStrRepr() << "," << linepair.first << ","
                 << linepair.second << "pa=" << poscomb.first.first
                 << ",pb=" << poscomb.second.first
                 << ",pa'=" << poscomb.first.second
                 << ",pb'=" << poscomb.second.second << ",n=" << node << ")";
              glp_set_row_name(lp, row, ss.str().c_str());
              glp_set_row_bnds(lp, row, GLP_UP, 0, 3);

              (*c)++;
              ia[*c] = row;
              ja[*c] = lineAinAatP;
              res[*c] = 1;

              (*c)++;
              ia[*c] = row;
              ja[*c] = lineBinAatP;
              res[*c] = 1;

              (*c)++;
              ia[*c] = row;
              ja[*c] = lineAinBatP;
              res[*c] = 1;

              (*c)++;
              ia[*c] = row;
              ja[*c] = lineBinBatP;
              res[*c] = 1;

              (*c)++;
              ia[*c] = row;
              ja[*c] = decisionVar;
              res[*c] = -1;
            }
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
void ILPOptimizer::writeDiffSegConstraints(const OptGraph& g, int* ia, int* ja,
                                           double* res, size_t* c,
                                           glp_prob* lp) const {
  // go into nodes and build crossing constraints for adjacent
  for (OptNode* node : g.getNodes()) {
    std::set<OptEdge*> processed;
    for (OptEdge* segmentA : node->adjList) {
      processed.insert(segmentA);
      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segmentA)) {
        for (EdgePair segments :
             getEdgePartnerPairs(node, segmentA, linepair)) {
          // if (processed.find(segmentB) != processed.end()) continue;

          // try all position combinations
          size_t decisionVar = glp_add_cols(lp, 1);

          // introduce dec var
          std::stringstream ss;
          ss << "x_dec(" << segmentA->getStrRepr() << ","
             << segments.first->getStrRepr() << segments.second->getStrRepr()
             << "," << linepair.first << "(" << linepair.first->id << "),"
             << linepair.second << "(" << linepair.second->id << ")," << node
             << ")";
          glp_set_col_name(lp, decisionVar, ss.str().c_str());
          glp_set_col_kind(lp, decisionVar, GLP_BV);
          glp_set_obj_coef(lp, decisionVar, 1);

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
              ss << "dec_sum(" << segmentA->getStrRepr() << ","
                 << segments.first->getStrRepr()
                 << segments.second->getStrRepr() << "," << linepair.first
                 << "," << linepair.second << "pa=" << poscomb.first
                 << ",pb=" << poscomb.second << ",n=" << node << ")";
              glp_set_row_name(lp, row, ss.str().c_str());
              glp_set_row_bnds(lp, row, GLP_UP, 0, 1);

              (*c)++;
              ia[*c] = row;
              ja[*c] = lineAinAatP;
              res[*c] = 1;

              (*c)++;
              ia[*c] = row;
              ja[*c] = lineBinAatP;
              res[*c] = 1;

              (*c)++;
              ia[*c] = row;
              ja[*c] = decisionVar;
              res[*c] = -1;
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
  graph::Edge* etgA = a->etgs[0].etg;
  graph::Edge* etgB = b->etgs[0].etg;
  for (size_t posLineAinA = 0; posLineAinA < etgA->getCardinality();
       posLineAinA++) {
    for (size_t posLineBinA = 0; posLineBinA < etgA->getCardinality();
         posLineBinA++) {
      if (posLineAinA == posLineBinA) continue;

      for (size_t posLineAinB = 0; posLineAinB < etgB->getCardinality();
           posLineAinB++) {
        for (size_t posLineBinB = 0; posLineBinB < etgB->getCardinality();
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
  graph::Edge* etgA = a->etgs[0].etg;
  for (size_t posLineAinA = 0; posLineAinA < etgA->getCardinality();
       posLineAinA++) {
    for (size_t posLineBinA = 0; posLineBinA < etgA->getCardinality();
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
  varName << "x_(" << seg->getStrRepr() << ",l=" << r << ",p=" << p << ")";
  return varName.str();
}

// _____________________________________________________________________________
std::vector<OptEdge*> ILPOptimizer::getEdgePartners(
    OptNode* node, OptEdge* segmentA, const LinePair& linepair) const {
  std::vector<OptEdge*> ret;

  const Node* routeDirA = 0;
  const Node* routeDirB = 0;

  graph::Edge* fromEtg = 0;

  // outfactor (see below)
  if (segmentA->from == node) {
    fromEtg = segmentA->getFirstEdge().etg;
  } else {
    fromEtg = segmentA->getLastEdge().etg;
  }

  routeDirA = fromEtg->getTripsForRoute(linepair.first)->direction;
  routeDirB = fromEtg->getTripsForRoute(linepair.second)->direction;
  //

  for (OptEdge* segmentB : node->adjList) {
    if (segmentB == segmentA) continue;
    graph::Edge* etg = 0;
    graph::Edge* fromEtg = 0;

    // outfactor (see below)
    if (segmentB->from == node) {
      etg = segmentB->getFirstEdge().etg;
    } else {
      etg = segmentB->getLastEdge().etg;
    }
    //

    if (etg->getContinuedRoutesIn(node->node, linepair.first, routeDirA,
                                  fromEtg)
            .size() &&
        etg->getContinuedRoutesIn(node->node, linepair.second, routeDirB,
                                  fromEtg)
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

  const Node* routeDirA = 0;
  const Node* routeDirB = 0;

  graph::Edge* fromEtg = 0;

  // outfactor (see below)
  if (segmentA->from == node) {
    fromEtg = segmentA->getFirstEdge().etg;
  } else {
    fromEtg = segmentA->getLastEdge().etg;
  }

  routeDirA = fromEtg->getTripsForRoute(linepair.first)->direction;
  routeDirB = fromEtg->getTripsForRoute(linepair.second)->direction;
  //

  for (OptEdge* segmentB : node->adjList) {
    if (segmentB == segmentA) continue;

    graph::Edge* etg = 0;

    // outfactor (see below and above)
    if (segmentB->from == node) {
      etg = segmentB->getFirstEdge().etg;
    } else {
      etg = segmentB->getLastEdge().etg;
    }
    //

    if (etg->getContinuedRoutesIn(node->node, linepair.first, routeDirA,
                                  fromEtg)
            .size()) {
      EdgePair curPair;
      curPair.first = segmentB;
      for (OptEdge* segmentC : node->adjList) {
        if (segmentC == segmentA || segmentC == segmentB) continue;

        graph::Edge* etg = 0;

        // outfactor (see above)
        if (segmentC->from == node) {
          etg = segmentC->getFirstEdge().etg;
        } else {
          etg = segmentC->getLastEdge().etg;
        }
        //

        if (etg->getContinuedRoutesIn(node->node, linepair.second, routeDirB,
                                      fromEtg)
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
  std::set<const Route*> processed;
  std::vector<LinePair> ret;
  for (auto& toA : *segment->etgs[0].etg->getTripsUnordered()) {
    processed.insert(toA.route);
    for (auto& toB : *segment->etgs[0].etg->getTripsUnordered()) {
      if (toA.route == toB.route) continue;
      ret.push_back(LinePair(toA.route, toB.route));
    }
  }
  return ret;
}

// _____________________________________________________________________________
void ILPOptimizer::solveProblem(glp_prob* lp) const {
  glp_iocp params;
  glp_init_iocp(&params);
  params.presolve = GLP_ON;
  // params.tm_lim = 30000;
  // params.mip_gap = 0.5;
  params.binarize = 1;
  // params.bt_tech = GLP_BT_BPH;
  // params.fp_heur = GLP_ON;
  // params.ps_heur = GLP_ON;
  params.ps_tm_lim = 120000;

  glp_intopt(lp, &params);
}

// _____________________________________________________________________________
bool ILPOptimizer::crosses(OptNode* node, OptEdge* segmentA, OptEdge* segmentB,
                           PosComPair postcomb) const {
  bool otherWayA = (segmentA->from != node) ^ segmentA->etgs.front().dir;
  bool otherWayB = (segmentB->from != node) ^ segmentB->etgs.front().dir;

  size_t cardA = segmentA->etgs.front().etg->getCardinality();
  size_t cardB = segmentB->etgs.front().etg->getCardinality();

  size_t posAinA =
      otherWayA ? cardA - 1 - postcomb.first.first : postcomb.first.first;
  size_t posAinB =
      otherWayB ? cardB - 1 - postcomb.first.second : postcomb.first.second;
  size_t posBinA =
      otherWayA ? cardA - 1 - postcomb.second.first : postcomb.second.first;
  size_t posBinB =
      otherWayB ? cardB - 1 - postcomb.second.second : postcomb.second.second;

  bool pCrossing = (posAinA < posBinA && posAinB < posBinB);

  return pCrossing;
}

// _____________________________________________________________________________
bool ILPOptimizer::crosses(OptNode* node, OptEdge* segmentA, EdgePair segments,
                           PosCom postcomb) const {
  bool otherWayA = (segmentA->from != node) ^ segmentA->etgs.front().dir;
  bool otherWayB =
      (segments.first->from != node) ^ segments.first->etgs.front().dir;
  bool otherWayC =
      (segments.second->from != node) ^ segments.second->etgs.front().dir;

  size_t cardA = segmentA->etgs.front().etg->getCardinality();
  size_t cardB = segments.first->etgs.front().etg->getCardinality();
  size_t cardC = segments.second->etgs.front().etg->getCardinality();

  size_t posAinA = otherWayA ? cardA - 1 - postcomb.first : postcomb.first;
  size_t posBinA = otherWayA ? cardA - 1 - postcomb.second : postcomb.second;

  Point aInA = getPos(node, segmentA, posAinA);
  Point bInA = getPos(node, segmentA, posBinA);

  for (size_t i = 0; i < segments.first->etgs.front().etg->getCardinality();
       ++i) {
    for (size_t j = 0; j < segments.second->etgs.front().etg->getCardinality();
         ++j) {
      size_t posAinB = otherWayB ? cardB - 1 - i : i;
      size_t posBinC = otherWayC ? cardC - 1 - j : j;

      Point aInB = getPos(node, segments.first, posAinB);
      Point bInB = getPos(node, segments.second, posBinC);

      Line a;
      a.push_back(aInA);
      a.push_back(aInB);

      Line b;
      b.push_back(bInA);
      b.push_back(bInB);

      if (pbutil::geo::intersects(aInA, aInB, bInA, bInA) ||
          bgeo::distance(a, b) < 1)
        return true;
    }
  }
  return false;
}

// _____________________________________________________________________________
Point ILPOptimizer::getPos(OptNode* n, OptEdge* segment, size_t p) const {
  // look for correct nodefront
  const NodeFront* nf = 0;
  for (auto etg : segment->etgs) {
    const NodeFront* test = n->node->getNodeFrontFor(etg.etg);
    if (test) {
      nf = test;
      break;
    }
  }

  assert(nf);

  return nf->getTripPos(segment->etgs.front().etg, p, false);
}
