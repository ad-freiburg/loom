// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "./../../log/Log.h"
#include "./ILPEdgeOrderOptimizer.h"
#include "./OptGraph.h"
#include "./../output/OgrOutput.h"
#include "./../graph/OrderingConfiguration.h"
#include "./../util/Geo.h"
#include <glpk.h>

using namespace transitmapper;
using namespace optim;
using namespace graph;

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::optimize() {
  // create optim graph
  OptGraph g(_g);
  g.simplify();

  output::OgrOutput ogrOut("/home/patrick/optimgraph", _cfg);
  //ogrOut.print(g);

  glp_prob* lp = createProblem(g);

  // write problem for debugging...
  glp_write_mps(lp, GLP_MPS_FILE, 0, "/home/patrick/ilp");

  solveProblem(lp);

  std::cout << "ILP obj = " << glp_mip_obj_val(lp) << std::endl;

  glp_print_mip(lp, "/home/patrick/ilp.sol");

  Configuration c;
  getConfigurationFromSoluation(lp, &c, g);
  _g->setConfig(c);

  glp_delete_prob(lp);
  glp_free_env();
}

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::getConfigurationFromSoluation(glp_prob* lp,
    Configuration* c, const OptGraph& g) const {
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
            varName << "x_(" << e->getStrRepr() << ",l="
              << r.route << ",p=" << tp << ")";

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
glp_prob* ILPEdgeOrderOptimizer::createProblem(const OptGraph& g) const {
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
      graph::EdgeTripGeom* etg = e->etgs[0].etg;

      // get string repr of etg

      size_t newCols = etg->getCardinality() * etg->getCardinality();
      size_t cols = glp_add_cols(lp, newCols);
      size_t i = 0;
      size_t rowA = glp_add_rows(lp, etg->getCardinality());

      for (size_t p = 0; p < etg->getCardinality(); p++) {
        std::stringstream varName;

        varName << "sum(" << e->getStrRepr() << ",p="
          << p << ")";
        glp_set_row_name(lp, rowA + p, varName.str().c_str());
        glp_set_row_bnds(lp, rowA + p, GLP_FX, 1, 1);
      }

      for (auto r : *etg->getTripsUnordered()) {

        // constraint: the sum of all x_slp over p must be 1 for equal sl
        size_t row = glp_add_rows(lp, 1);
        std::stringstream varName;
        varName << "sum(" << e->getStrRepr() << ",l="
          << r.route << ")";
        glp_set_row_name(lp, row, varName.str().c_str());
        glp_set_row_bnds(lp, row, GLP_FX, 1, 1);

        for (size_t p = 0; p < etg->getCardinality(); p++) {
          std::stringstream varName;
          varName << "x_(" << e->getStrRepr() << ",l="
            << r.route << ",p=" << p << ")";
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
void ILPEdgeOrderOptimizer::writeSameSegConstraints(const OptGraph& g,
    int* ia, int* ja, double* res, size_t* c, glp_prob* lp)
const {
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
          ss << "x_dec(" << segmentA->getStrRepr() << "," << segmentB->getStrRepr()
            << "," << linepair.first << "(" << linepair.first->getShortName() << "),"
            << linepair.second << "(" << linepair.second->getShortName() << ")," << node << ")";
          glp_set_col_name(lp, decisionVar, ss.str().c_str());
          glp_set_col_kind(lp, decisionVar, GLP_BV);
          glp_set_obj_coef(lp, decisionVar, 1);

          for (PosComPair poscomb : getPositionCombinations(segmentA, segmentB)) {
            if (crosses(node, segmentA, segmentB, poscomb)) {
							size_t lineAinAatP = glp_find_col(lp, getILPVarName(segmentA, linepair.first, poscomb.first.first).c_str());
							size_t lineBinAatP = glp_find_col(lp, getILPVarName(segmentA, linepair.second, poscomb.second.first).c_str());
							size_t lineAinBatP = glp_find_col(lp, getILPVarName(segmentB, linepair.first, poscomb.first.second).c_str());
							size_t lineBinBatP = glp_find_col(lp, getILPVarName(segmentB, linepair.second, poscomb.second.second).c_str());

							assert(lineAinAatP > 0);
							assert(lineAinBatP > 0);
							assert(lineBinAatP > 0);
							assert(lineBinBatP > 0);

							size_t row = glp_add_rows(lp, 1);
							std::stringstream ss;
							ss << "dec_sum(" << segmentA->getStrRepr() << "," << segmentB->getStrRepr()
								<< "," << linepair.first << "," << linepair.second
								<< "pa=" << poscomb.first.first << ",pb=" << poscomb.second.first
							 << ",pa'=" << poscomb.first.second << ",pb'=" << poscomb.second.second
								<< ",n=" << node << ")";
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
void ILPEdgeOrderOptimizer::writeDiffSegConstraints(const OptGraph& g,
    int* ia, int* ja, double* res, size_t* c, glp_prob* lp)
const {
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
        for (EdgePair segments : getEdgePartnerPairs(node, segmentA, linepair)) {
          // if (processed.find(segmentB) != processed.end()) continue;

          // try all position combinations
          size_t decisionVar = glp_add_cols(lp, 1);

          // introduce dec var
          std::stringstream ss;
          ss << "x_dec(" << segmentA->getStrRepr() << ","
            << segments.first->getStrRepr()
            << segments.second->getStrRepr()
            << "," << linepair.first << "(" << linepair.first->getShortName() << "),"
            << linepair.second << "(" << linepair.second->getShortName() << ")," << node << ")";
          glp_set_col_name(lp, decisionVar, ss.str().c_str());
          glp_set_col_kind(lp, decisionVar, GLP_BV);
          glp_set_obj_coef(lp, decisionVar, 1);

          for (PosCom poscomb : getPositionCombinations(segmentA)) {
            std::cout << "checking: "
              << linepair.first->getShortName() << " and " 
              << linepair.second->getShortName() << " @ "
              << poscomb.first << " and " << poscomb.second
              << " in " << (node->node->getStops().size() ? (*node->node->getStops().begin())->getName() : "?")
              << std::endl;
            if (crosses(node, segmentA, segments, poscomb)) {
              std::cout << "crossing: "
                << linepair.first->getShortName() << " and " 
                << linepair.second->getShortName() << " @ "
                << poscomb.first << " and " << poscomb.second
                << " in " << (node->node->getStops().size() ? (*node->node->getStops().begin())->getName() : "?")
                << std::endl;
							size_t lineAinAatP = glp_find_col(lp, getILPVarName(segmentA, linepair.first, poscomb.first).c_str());
							size_t lineBinAatP = glp_find_col(lp, getILPVarName(segmentA, linepair.second, poscomb.second).c_str());

							assert(lineAinAatP > 0);
							assert(lineBinAatP > 0);

							size_t row = glp_add_rows(lp, 1);
							std::stringstream ss;
							ss << "dec_sum(" << segmentA->getStrRepr() << ","
                << segments.first->getStrRepr()
                << segments.second->getStrRepr()
								<< "," << linepair.first << "," << linepair.second
								<< "pa=" << poscomb.first << ",pb=" << poscomb.second
								<< ",n=" << node << ")";
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
std::vector<PosComPair> ILPEdgeOrderOptimizer::getPositionCombinations(OptEdge* a,
    OptEdge* b) const {
  std::vector<PosComPair> ret;
  graph::EdgeTripGeom* etgA = a->etgs[0].etg;
  graph::EdgeTripGeom* etgB = b->etgs[0].etg;
  for (size_t posLineAinA = 0; posLineAinA < etgA->getCardinality(); posLineAinA++) {
    for (size_t posLineBinA = 0; posLineBinA < etgA->getCardinality(); posLineBinA++) {
      if (posLineAinA == posLineBinA) continue;

      for (size_t posLineAinB = 0; posLineAinB < etgB->getCardinality(); posLineAinB++) {
        for (size_t posLineBinB = 0; posLineBinB < etgB->getCardinality(); posLineBinB++) {
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
std::vector<PosCom> ILPEdgeOrderOptimizer::getPositionCombinations(OptEdge* a)
const {
  std::vector<PosCom> ret;
  graph::EdgeTripGeom* etgA = a->etgs[0].etg;
  for (size_t posLineAinA = 0; posLineAinA < etgA->getCardinality(); posLineAinA++) {
    for (size_t posLineBinA = 0; posLineBinA < etgA->getCardinality(); posLineBinA++) {
      if (posLineAinA == posLineBinA) continue;
      ret.push_back(PosCom(posLineAinA, posLineBinA));
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::string ILPEdgeOrderOptimizer::getILPVarName(OptEdge* seg,
    const gtfs::Route* r, size_t p) const {
  std::stringstream varName;
  varName << "x_(" << seg->getStrRepr() << ",l="
    << r << ",p=" << p << ")";
  return varName.str();
}

// _____________________________________________________________________________
std::vector<OptEdge*> ILPEdgeOrderOptimizer::getEdgePartners(OptNode* node,
  OptEdge* segmentA, const LinePair& linepair) const {
  std::vector<OptEdge*> ret;
  for (OptEdge* segmentB : node->adjList) {
    if (segmentB == segmentA) continue;
    graph::EdgeTripGeom* etg = segmentB->etgs[0].etg;
    if (etg->getTripsForRoute(linepair.first) &&
        etg->getTripsForRoute(linepair.second)) {
      ret.push_back(segmentB);
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::vector<EdgePair> ILPEdgeOrderOptimizer::getEdgePartnerPairs(OptNode* node,
  OptEdge* segmentA, const LinePair& linepair) const {
  std::vector<EdgePair> ret;
  for (OptEdge* segmentB : node->adjList) {
    if (segmentB == segmentA) continue;
    graph::EdgeTripGeom* etg = segmentB->etgs[0].etg;
    if (etg->getTripsForRoute(linepair.first)) {
      EdgePair curPair;
      curPair.first = segmentB;
      for (OptEdge* segmentC : node->adjList) {
        if (segmentC == segmentA || segmentC == segmentB) continue;
        graph::EdgeTripGeom* etg = segmentC->etgs[0].etg;
        if (etg->getTripsForRoute(linepair.second)) {
          curPair.second = segmentC;
          ret.push_back(curPair);
        }
      }
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::vector<LinePair> ILPEdgeOrderOptimizer::getLinePairs(OptEdge* segment)
const {
  std::set<Route*> processed;
  std::vector<LinePair> ret;
  for (auto& toA : *segment->etgs[0].etg->getTripsUnordered()) {
    processed.insert(toA.route);
    for (auto& toB : *segment->etgs[0].etg->getTripsUnordered()) {
      if (processed.find(toB.route) != processed.end()) continue;
      ret.push_back(LinePair(toA.route, toB.route));
    }
  }
  return ret;
}

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::solveProblem(glp_prob* lp) const {
  glp_iocp params;
  glp_init_iocp(&params);
  params.presolve = GLP_ON;
  params.tm_lim = 10000;
  glp_intopt(lp, &params);
}

// _____________________________________________________________________________
bool ILPEdgeOrderOptimizer::crosses(OptNode* node, OptEdge* segmentA,
      OptEdge* segmentB, PosComPair postcomb) const {
  Point aInA = getPos(node, segmentA, postcomb.first.first);
  Point aInB = getPos(node, segmentB, postcomb.first.second);

  Point bInA = getPos(node, segmentA, postcomb.second.first);
  Point bInB = getPos(node, segmentB, postcomb.second.second);

  Line a;
  a.push_back(aInA);
  a.push_back(aInB);

  Line b;
  b.push_back(bInA);
  b.push_back(bInB);

  if (bgeo::distance(a, b) < 1)  return true;
  return util::geo::intersects(aInA, aInB, bInA, bInA);
}

// _____________________________________________________________________________
bool ILPEdgeOrderOptimizer::crosses(OptNode* node, OptEdge* segmentA,
    EdgePair segments, PosCom postcomb) const {
  Point aInA = getPos(node, segmentA, postcomb.first);
  Point bInA = getPos(node, segmentA, postcomb.second);

  for (size_t i = 0; i < segments.first->etgs.front().etg->getCardinality(); ++i) {
    for (size_t j = 0; j < segments.second->etgs.front().etg->getCardinality(); ++j) {
      Point aInB = getPos(node, segments.first, i);
      Point bInB = getPos(node, segments.second, j);

      Line a;
      a.push_back(aInA);
      a.push_back(aInB);

      Line b;
      b.push_back(bInA);
      b.push_back(bInB);

      if (util::geo::intersects(aInA, aInB, bInA, bInA) ||
          bgeo::distance(a, b) < 1) return true;
    }
  }
  return false;
}

// _____________________________________________________________________________
Point ILPEdgeOrderOptimizer::getPos(OptNode* n, OptEdge* segment, size_t p) const {
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

  bool otherWay = (segment->from != n) ^ segment->etgs.front().dir;

  return nf->getTripPos(*segment->etgs.front().etg, p, otherWay);
}
