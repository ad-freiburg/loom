// Copyright 2022
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "shared/linegraph/LineGraph.h"
#include "topoeval/DirLineGraph.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/EDijkstra.h"
#include "util/log/Log.h"

std::vector<shared::linegraph::LineNode*> gtNds;
std::vector<std::string> gtStations;

util::geo::Grid<topoeval::DirLineNode*, Point, double> testGrid;
util::geo::Grid<topoeval::DirLineNode*, Point, double> gtGrid;

std::unordered_map<std::string, std::set<topoeval::DirLineNode*>> stationsGt;
std::unordered_map<std::string, std::set<topoeval::DirLineNode*>> stationsTest;

struct CostFunc
    : public util::graph::EDijkstra::CostFunc<topoeval::DirLineNodePL,
                                              topoeval::DirLineEdgePL, double> {
  CostFunc(const shared::linegraph::Line* r) : _line(r) {}
  double inf() const { return std::numeric_limits<double>::infinity(); };
  double operator()(const topoeval::DirLineEdge* from,
                    const topoeval::DirLineNode* n,
                    const topoeval::DirLineEdge* to) const {
    // don't count start edge
    if (!from) return 0;

    // if an edge does not contain the line we are routing for, set
    // cost to inf
    if (!from->pl().lines.count(_line)) return inf();
    if (!to->pl().lines.count(_line)) return inf();

    if (n) {
      auto lRestrs = n->pl().restrs.find(_line);
      if (lRestrs != n->pl().restrs.end() && lRestrs->second.count(from)) {
        if (lRestrs->second.find(from)->second.count(to)) return inf();
      }

      // don't allow going back the same edge
      if (from->getOtherNd(n) == to->getOtherNd(n)) return inf();
    }

    return util::geo::len(util::geo::Line<double>{
        *to->getFrom()->pl().getGeom(), *to->getTo()->pl().getGeom()});
  };

  const shared::linegraph::Line* _line;
};

// _____________________________________________________________________________
void getDirLineGraph(const shared::linegraph::LineGraph* g,
                     topoeval::DirLineGraph* ret) {
  std::unordered_map<const shared::linegraph::LineNode*, topoeval::DirLineNode*>
      nMap;
  std::unordered_map<const shared::linegraph::LineEdge*,
                     std::vector<topoeval::DirLineEdge*>>
      eMap;
  for (auto nd : g->getNds()) {
    if (nd->pl().stops().size())
      nMap[nd] =
          ret->addNd({nd->pl().stops().front().name, *nd->pl().getGeom()});
    else
      nMap[nd] = ret->addNd({"", *nd->pl().getGeom()});
  }

  for (auto nd : g->getNds()) {
    for (auto edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      const auto& pl = edg->pl().getPolyline();
      eMap[edg] = {
          ret->addEdg(nMap[edg->getFrom()], nMap[edg->getTo()], pl),
          ret->addEdg(nMap[edg->getTo()], nMap[edg->getFrom()], pl.reversed())};

      for (auto r : edg->pl().getLines()) {
        if (r.direction == 0 || r.direction == edg->getTo()) {
          eMap[edg][0]->pl().lines.insert(r.line);
        }

        if (r.direction == 0 || r.direction == edg->getFrom()) {
          eMap[edg][1]->pl().lines.insert(r.line);
        }
      }
    }
  }

  // copy turn restrictions from original graph
  for (auto nd : g->getNds()) {
    for (auto ex : nd->pl().getConnExc()) {
      auto line = ex.first;
      for (auto exPair : ex.second) {
        const auto* edgeFr = exPair.first;

        for (auto* rEdgeFr : eMap[edgeFr]) {
          for (auto edgeTo : exPair.second) {
            for (auto* rEdgeTo : eMap[edgeTo]) {
              nMap[nd]->pl().restrs[line][rEdgeFr].insert(rEdgeTo);
              nMap[nd]->pl().restrs[line][rEdgeTo].insert(rEdgeFr);
            }
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
void denseSample(shared::linegraph::LineGraph* g, double d) {
  std::vector<shared::linegraph::LineEdge*> edgs;
  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      edgs.push_back(e);
    }
  }

  for (auto e : edgs) {
    auto denseL = util::geo::densify(*e->pl().getGeom(), d);

    auto fr = e->getFrom();
    auto to = e->getTo();
    auto pl = e->pl();
    g->delEdg(e->getFrom(), e->getTo());

    auto prev = fr;

    for (size_t i = 1; i < denseL.size() - 1; i++) {
      auto supNd = g->addNd(denseL[i]);

      auto eA = g->addEdg(prev, supNd, pl);

      shared::linegraph::LineGraph::nodeRpl(eA, fr, prev);
      shared::linegraph::LineGraph::nodeRpl(eA, to, supNd);

      eA->pl().setGeom(
          {*eA->getFrom()->pl().getGeom(), *eA->getTo()->pl().getGeom()});

      if (i == 1) shared::linegraph::LineGraph::edgeRpl(fr, e, eA);

      prev = supNd;
    }

    auto eA = g->addEdg(prev, to, pl);

    shared::linegraph::LineGraph::nodeRpl(eA, fr, prev);

    eA->pl().setGeom(
        {*eA->getFrom()->pl().getGeom(), *eA->getTo()->pl().getGeom()});
    shared::linegraph::LineGraph::edgeRpl(to, e, eA);
  }
}

// _____________________________________________________________________________
std::pair<std::pair<std::set<topoeval::DirLineNode*>,
                    std::set<topoeval::DirLineNode*>>,
          std::pair<std::set<topoeval::DirLineNode*>,
                    std::set<topoeval::DirLineNode*>>>
getFromToCandsStat() {
  std::pair<std::pair<std::set<topoeval::DirLineNode*>,
                      std::set<topoeval::DirLineNode*>>,
            std::pair<std::set<topoeval::DirLineNode*>,
                      std::set<topoeval::DirLineNode*>>>
      ret;
  auto fromStat = gtStations[rand() % gtStations.size()];
  auto toStat = gtStations[rand() % gtStations.size()];

  ret.first.first = stationsGt[fromStat];
  ret.first.second = stationsGt[toStat];
  ret.second.first = stationsTest[fromStat];
  ret.second.second = stationsTest[toStat];

  return ret;
}

// _____________________________________________________________________________
std::pair<std::pair<std::set<topoeval::DirLineNode*>,
                    std::set<topoeval::DirLineNode*>>,
          std::pair<std::set<topoeval::DirLineNode*>,
                    std::set<topoeval::DirLineNode*>>>
getFromToCands(double d) {
  std::pair<std::pair<std::set<topoeval::DirLineNode*>,
                      std::set<topoeval::DirLineNode*>>,
            std::pair<std::set<topoeval::DirLineNode*>,
                      std::set<topoeval::DirLineNode*>>>
      ret;
  auto from = gtNds[rand() % gtNds.size()];
  auto to = gtNds[rand() % gtNds.size()];

  std::set<topoeval::DirLineNode*> testNeighsFr;
  std::set<topoeval::DirLineNode*> testNeighsTo;
  std::set<topoeval::DirLineNode*> testFr;
  std::set<topoeval::DirLineNode*> testTo;

  std::set<topoeval::DirLineNode*> gtNeighsFr;
  std::set<topoeval::DirLineNode*> gtNeighsTo;
  std::set<topoeval::DirLineNode*> gtFr;
  std::set<topoeval::DirLineNode*> gtTo;

  auto fromGeom = *from->pl().getGeom();
  auto toGeom = *to->pl().getGeom();

  // auto fromGeom = util::geo::Point<double>{875686, 6.11379e+06};
  // auto toGeom = util::geo::Point<double>{876380, 6.11301e+06};

  testGrid.get(fromGeom, d, &testNeighsFr);
  testGrid.get(toGeom, d, &testNeighsTo);

  for (auto neigh : testNeighsFr) {
    if (util::geo::dist(*neigh->pl().getGeom(), fromGeom) <= d)
      testFr.insert(neigh);
  }

  for (auto neigh : testNeighsTo) {
    if (util::geo::dist(*neigh->pl().getGeom(), toGeom) <= d)
      testTo.insert(neigh);
  }

  gtGrid.get(fromGeom, d, &gtNeighsFr);
  gtGrid.get(toGeom, d, &gtNeighsTo);

  for (auto neigh : gtNeighsFr) {
    if (util::geo::dist(*neigh->pl().getGeom(), fromGeom) <= d)
      gtFr.insert(neigh);
  }

  for (auto neigh : gtNeighsTo) {
    if (util::geo::dist(*neigh->pl().getGeom(), toGeom) <= d)
      gtTo.insert(neigh);
  }

  ret.first.first = gtFr;
  ret.first.second = gtTo;
  ret.second.first = testFr;
  ret.second.second = testTo;

  return ret;
}

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  double d = 150;
  size_t SAMPLES = 10000;

  bool fromStations = false;

  std::string gtPath, testPath;

  for (int i = 1; i < argc; i++) {
    std::string cur = argv[i];
    if (cur == "-h" || cur == "--help") {
      std::cerr << "Usage: " << argv[0]
                << "[-d <maxdist=150>] [-s <numsamples=10000>] "
                   "[--sample-stations] <ground truth "
                   "graph> <test graph>"
                << std::endl;
      exit(0);
    } else if (cur == "-d") {
      if (++i >= argc) {
        LOG(ERROR) << "Missing argument for maxdist (-d).";
        exit(1);
      }
      d = atof(argv[i]);
    } else if (cur == "--sample-stations") {
      fromStations = true;
    } else if (cur == "-s") {
      if (++i >= argc) {
        LOG(ERROR) << "Missing argument for samples (-s).";
        exit(1);
      }
      SAMPLES = atoi(argv[i]);
    } else {
      if (gtPath.empty())
        gtPath = cur;
      else if (testPath.empty())
        testPath = cur;
      else {
        std::cerr << "Usage: " << argv[0]
                  << "-d <maxdist=150> <ground truth graph> <test graph>"
                  << std::endl;
        exit(1);
      }
    }
  }

  if (gtPath.empty()) {
    std::cerr << "Missing ground truth graph path." << std::endl;
    exit(1);
  }

  if (testPath.empty()) {
    std::cerr << "Missing test graph path." << std::endl;
    exit(1);
  }

  shared::linegraph::LineGraph gtGraph;
  shared::linegraph::LineGraph testGraph;

  std::ifstream ifs;

  ifs.open(gtPath);
  gtGraph.readFromJson(&ifs, false);
  ifs.close();

  ifs.open(testPath);
  testGraph.readFromJson(&ifs, false);
  ifs.close();

  LOG(DEBUG) << "Ground truth graph: " << gtGraph.getNds()->size() << " nodes";
  LOG(DEBUG) << "Test graph: " << testGraph.getNds()->size() << " nodes";

  denseSample(&gtGraph, 10);
  denseSample(&testGraph, 10);

  topoeval::DirLineGraph gtDirGraph;
  topoeval::DirLineGraph testDirGraph;

  getDirLineGraph(&gtGraph, &gtDirGraph);
  getDirLineGraph(&testGraph, &testDirGraph);

  // collect stations
  for (auto nd : *gtDirGraph.getNds()) {
    if (nd->pl().getStatLabel().size()) {
      if (!stationsGt.count(nd->pl().getStatLabel()))
        gtStations.push_back(nd->pl().getStatLabel());
      stationsGt[nd->pl().getStatLabel()].insert(nd);
    }
  }
  for (auto nd : *testDirGraph.getNds()) {
    if (nd->pl().getStatLabel().size()) {
      stationsTest[nd->pl().getStatLabel()].insert(nd);
    }
  }

  testGrid = util::geo::Grid<topoeval::DirLineNode*, Point, double>(
      120, 120, testGraph.getBBox());
  for (auto nd : *testDirGraph.getNds()) {
    testGrid.add(*nd->pl().getGeom(), nd);
  }

  gtGrid = util::geo::Grid<topoeval::DirLineNode*, Point, double>(
      120, 120, gtGraph.getBBox());
  for (auto nd : *gtDirGraph.getNds()) {
    gtGrid.add(*nd->pl().getGeom(), nd);
  }

  gtNds = std::vector<shared::linegraph::LineNode*>(gtGraph.getNds()->begin(),
                                                    gtGraph.getNds()->end());

  std::map<const shared::linegraph::Line*, const shared::linegraph::Line*>
      lineMap;

  for (auto nd : *gtGraph.getNds()) {
    for (auto e : nd->getAdjList()) {
      for (auto l : e->pl().getLines()) {
        lineMap[l.line] = 0;
      }
    }
  }

  for (auto i : lineMap) {
    for (auto nd : *testGraph.getNds()) {
      for (auto e : nd->getAdjList()) {
        for (auto l : e->pl().getLines()) {
          if (l.line->id() == i.first->id()) {
            lineMap[i.first] = l.line;
          }
        }
      }
    }
  }

  double match = 0, unmatch = 0;

  double minFrech = std::numeric_limits<double>::max();
  double maxFrech = 0;
  double frechAvg = 0;

  for (size_t i = 0; i < SAMPLES; i++) {
    std::vector<const shared::linegraph::Line*> lines;
    std::set<const shared::linegraph::Line*> linesSet;

    std::pair<std::pair<std::set<topoeval::DirLineNode*>,
                        std::set<topoeval::DirLineNode*>>,
              std::pair<std::set<topoeval::DirLineNode*>,
                        std::set<topoeval::DirLineNode*>>>
        cands;

    if (fromStations)
      cands = getFromToCandsStat();
    else
      cands = getFromToCands(d);

    auto gtFr = cands.first.first;
    auto gtTo = cands.first.second;

    auto testFr = cands.second.first;
    auto testTo = cands.second.second;

    for (auto from : gtFr) {
      for (auto e : from->getAdjList()) {
        for (auto l : e->pl().lines) {
          if (linesSet.count(l)) continue;
          lines.push_back(l);
          linesSet.insert(l);
        }
      }
    }

    for (auto to : gtTo) {
      for (auto e : to->getAdjList()) {
        for (auto l : e->pl().lines) {
          if (linesSet.count(l)) continue;
          lines.push_back(l);
          linesSet.insert(l);
        }
      }
    }

    auto gtLine = lines[rand() % lines.size()];
    CostFunc cFuncGt(gtLine);

    util::graph::EList<topoeval::DirLineNodePL, topoeval::DirLineEdgePL>*
        resEdges = 0;
    util::graph::NList<topoeval::DirLineNodePL, topoeval::DirLineEdgePL>
        resNodesGt;
    util::graph::NList<topoeval::DirLineNodePL, topoeval::DirLineEdgePL>
        resNodesTest;

    util::graph::EList<topoeval::DirLineNodePL, topoeval::DirLineEdgePL>
        resEdgesTest;

    auto cGt = util::graph::EDijkstra::shortestPath(gtFr, gtTo, cFuncGt,
                                                    resEdges, &resNodesGt);

    if (!lineMap[gtLine]) {
      LOG(ERROR) << "Input line " << gtLine->id() << " (" << gtLine->label()
                 << ") not found in test data";
      exit(1);
    }
    CostFunc cFuncTest(lineMap[gtLine]);
    auto cTest = util::graph::EDijkstra::shortestPath(
        testFr, testTo, cFuncTest, &resEdgesTest, &resNodesTest);

    util::geo::Line<double> lineTest, lineGt;

    for (auto nd : resNodesGt) lineGt.push_back(*nd->pl().getGeom());
    for (auto nd : resNodesTest) lineTest.push_back(*nd->pl().getGeom());

    if (cGt > std::numeric_limits<double>::max() &&
        cTest > std::numeric_limits<double>::max()) {
      continue;
    }

    if ((cGt > std::numeric_limits<double>::max()) ^
        (cTest > std::numeric_limits<double>::max())) {
      unmatch += 1;
      continue;
    }

    double frechetDist = util::geo::frechetDist(lineTest, lineGt, 15);
    LOG(DEBUG) << " for line " << gtLine->label() << " : " << cGt << " vs "
               << cTest << ": fr " << frechetDist;

    frechAvg += frechetDist;

    if (frechetDist > maxFrech) maxFrech = frechetDist;
    if (frechetDist < minFrech) minFrech = frechetDist;

    if (frechetDist < d)
      match += 1;
    else {
      unmatch += 1;
    }
  }

  std::cout << match / (match + unmatch) << "\t" << minFrech << "\t" << maxFrech
            << "\t" << frechAvg / (match + unmatch) << std::endl;

  return (0);
}
