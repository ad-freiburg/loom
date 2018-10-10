// Copyright 2016
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "transitmap/config/ConfigReader.cpp"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/GraphBuilder.h"
#include "transitmap/graph/Node.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/optim/ILPEdgeOrderOptimizer.h"
#include "transitmap/optim/Scorer.h"
#include "transitmap/graph/Penalties.h"
#include "transitmap/output/SvgOutput.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using namespace transitmapper;
using namespace util::geo;
using std::string;

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  config::Config cfg;

  config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  if (true) {
    LOG(INFO) << "reading graph ";
    transitmapper::graph::TransitGraph g(cfg.name, cfg.projectionString);
    transitmapper::graph::GraphBuilder b(&cfg);
    if (!b.build(&(std::cin), &g)) {
      exit(1);
    }

    if (cfg.collapseLinePartners) {
      b.combinePartnerRoutes(&g);
    }


    LOG(INFO) << "Creating node fronts...";
    b.writeMainDirs(&g);

    if (cfg.dontExpandStations) {
      b.writeStationGeoms(&g);
    }

    if (cfg.expandFronts) {
      b.expandOverlappinFronts(&g);
    }

    if (!cfg.dontExpandStations) {
      b.writeStationGeoms(&g);
    }

    b.createMetaNodes(&g);

    LOG(INFO) << "Writing initial ordering configuration...";
    b.writeInitialConfig(&g);

    LOG(INFO) << "Optimizing...";

    double maxCrossPen = g.getMaxDegree() * (cfg.crossPenMultiSameSeg > cfg.crossPenMultiDiffSeg ? cfg.crossPenMultiSameSeg : cfg.crossPenMultiDiffSeg);
    double maxSplitPen = g.getMaxDegree() * cfg.splitPenWeight;

    // TODO move this into configuration, at least partially
    transitmapper::graph::Penalties pens{
      maxCrossPen,
      maxSplitPen,
      cfg.crossPenMultiSameSeg,
      cfg.crossPenMultiDiffSeg,
      cfg.splitPenWeight,
      cfg.stationCrossWeightSameSeg,
      cfg.stationCrossWeightDiffSeg,
      cfg.stationSplitWeight,
      true,
      true
    };

    optim::Scorer scorer(&g, pens);

    if (cfg.outputStats) {
      LOG(INFO) << "(stats) Stats for graph '" << g.getName();
      LOG(INFO) << "(stats)   Total node count: " << g.getNumNodes() << " ("
                << g.getNumNodes(true) << " topo, " << g.getNumNodes(false)
                << " non-topo)";
      LOG(INFO) << "(stats)   Total edge count: " << g.getNumEdges()
                ;
      LOG(INFO) << "(stats)   Total unique route count: " << g.getNumRoutes()
                ;
      LOG(INFO) << "(stats)   Max edge route cardinality: "
                << g.getMaxCardinality();
      LOG(INFO) << "(stats)   Number of poss. solutions: "
                << scorer.getNumPossSolutions();
      LOG(INFO) << "(stats)   Highest node degree: " << g.getMaxDegree();
    }

    LOG(INFO) << "(stats) Max crossing pen: " << maxCrossPen;
    LOG(INFO) << "(stats) Max splitting pen: " << maxSplitPen;

    LOG(INFO) << "(stats) Total graph score BEFORE optim is -- "
              << scorer.getScore()
              << " -- (incl. unavoidable crossings!)";
    LOG(INFO) << "(stats)   Per node graph score: "
              << scorer.getScore() /
                     g.getNodes()->size();
    LOG(INFO) << "(stats)   Crossings: " << scorer.getNumCrossings() << " (score: "
              << scorer.getCrossScore()
              << ")";
    LOG(INFO) << "(stats)   Separations: " << scorer.getNumSeparations()
              << " (score: "
              << scorer.getSeparationScore()
              << ")";

    if (!cfg.noOptim) {
      if (cfg.optimMethod == "ilp_impr") {
        optim::ILPEdgeOrderOptimizer ilpEoOptim(&g, &cfg, &scorer);
        ilpEoOptim.optimize();
      } else if (cfg.optimMethod == "ilp") {
        optim::ILPOptimizer ilpEoOptim(&g, &cfg, &scorer);
        ilpEoOptim.optimize();
      }
    }

    LOG(INFO) << "(stats) Total graph score AFTER optim is -- "
              << scorer.getScore()
              << " -- (incl. unavoidable crossings!)";
    LOG(INFO) << "(stats)   Per node graph score: "
              << scorer.getScore() /
                     g.getNodes()->size()
              ;
    LOG(INFO) << "(stats)   Crossings: " << scorer.getNumCrossings() << " (score: "
              << scorer.getCrossScore()
              << ")";
    LOG(INFO) << "(stats)   Separations: " << scorer.getNumSeparations()
              << " (score: "
              << scorer.getSeparationScore()
              << ")";

    if (cfg.renderMethod == "svg") {
      std::string path = cfg.outputPath;
      LOG(INFO) << "Outputting to SVG " << path << " ..." << std::endl;
      std::ofstream o;
      o.open(path);
      output::SvgOutput svgOut(&o, &cfg, &scorer);
      svgOut.print(g);
    }

    if (cfg.renderMethod == "svg_sep") {
      {
        std::string path = cfg.outputPath + "/edges.svg";
        cfg.renderEdges = 1;
        cfg.renderStations = 0;
        cfg.renderNodeConnections = 0;
        LOG(INFO) << "Outputting edge SVG to " << path << " ..." << std::endl;
        std::ofstream o;
        o.open(path);
        output::SvgOutput svgOut(&o, &cfg, &scorer);
        svgOut.print(g);
      }

      {
        std::string path = cfg.outputPath + "/nodes.svg";
        cfg.renderEdges = 0;
        cfg.renderStations = 0;
        cfg.renderNodeConnections = 1;
        LOG(INFO) << "Outputting node SVG to " << path << " ..." << std::endl;
        std::ofstream o;
        o.open(path);
        output::SvgOutput svgOut(&o, &cfg, &scorer);
        svgOut.print(g);
      }

      {
        std::string path = cfg.outputPath + "/stations.svg";
        cfg.renderEdges = 0;
        cfg.renderStations = 1;
        cfg.renderNodeConnections = 0;
        LOG(INFO) << "Outputting station SVG to " << path << " ..." << std::endl;
        std::ofstream o;
        o.open(path);
        output::SvgOutput svgOut(&o, &cfg, &scorer);
        svgOut.print(g);
      }
    }

    if (!cfg.worldFilePath.empty()) {
      LOG(INFO) << "Writing world file for this map to " << cfg.worldFilePath;

      std::ofstream file;
      file.open(cfg.worldFilePath);
      if (file) {
        file << 1 / cfg.outputResolution << std::endl
                << 0 << std::endl
                << 0 << std::endl
                << -1 / cfg.outputResolution << std::endl
                << std::fixed << g.getBoundingBox(cfg.outputPadding).getLowerLeft().getX()
                << std::endl
                << g.getBoundingBox(cfg.outputPadding).getUpperRight().getY() << std::endl;
        file.close();
      }
    }
  } else if (false) {
    PolyLine<double> a;
    PolyLine<double> b;

    a << DPoint(1, 1);
    a << DPoint(1, 2);
    a << DPoint(1, 3);
    a << DPoint(1, 4);
    a << DPoint(1, 5);

    b << DPoint(2, 3);
    b << DPoint(2, 4);
    b << DPoint(2, 5);
    b << DPoint(2, 6);
    b << DPoint(2, 7);

    SharedSegments<double> ss = a.getSharedSegments(b, 3);
    LOG(INFO) << "RESULT: " << ss.segments.size();

  } else {
    // just testing parallel drawing...
    PolyLine<double> p;
    PolyLine<double> a;
    PolyLine<double> b;
    PolyLine<double> c;

    for (double i = 0; i < 2000; i += 20) {
      // sin
      a << DPoint(i, 1000 + sin(0.006 * i) * 500);

      // cos
      b << DPoint(i, 1000 + cos(0.003 * i) * 500);

      // linear
      // a << Point(i, 1000);

      // triangle
      if (i < 1000) {
        p << DPoint(i, 1000 - (i / 2));
      } else {
        p << DPoint(i, (i / 2));
      }
    }

    std::ofstream o;
    o.open("/home/patrick/test.svg");
    cfg.outputResolution = 1;

    output::SvgOutput svgOut(&o, &cfg, 0);

    o << "<svg width=\"2000\" height=\"2000\">";

    DPoint pp(500, 470);
    DPoint ppp(550, 500);
    c << DPoint(400, 400) << DPoint(450, 450) << DPoint(550, 450)
      << DPoint(600, 700);

    c = c.getSegment(pp, ppp);

    DPoint projected = c.projectOn(pp).p;
    DPoint projectedB = c.projectOn(ppp).p;

    transitmapper::output::RenderParams rmp;
    rmp.height = 2000;
    rmp.width = 2000;
    rmp.xOff = 0;
    rmp.yOff = 0;

    svgOut.printLine(c, "fill:none;stroke:red;stroke-width:5", rmp);
    svgOut.printPoint(pp, "fill:none;stroke:red;stroke-width:5", rmp);
    svgOut.printPoint(ppp, "fill:none;stroke:red;stroke-width:5", rmp);
    svgOut.printPoint(projected, "fill:none;stroke:red;stroke-width:5", rmp);
    svgOut.printPoint(projectedB, "fill:none;stroke:red;stroke-width:5", rmp);

    // svgOut.printLine(p, "fill:none;stroke:black;stroke-width:10");

    std::vector<const PolyLine<double>*> ls;
    ls.push_back(&a);
    ls.push_back(&p);
    ls.push_back(&b);
    PolyLine<double> avg = PolyLine<double>::average(ls);

    svgOut.printLine(avg, "fill:none;stroke:#8844AA;stroke-width:15", rmp);

    svgOut.printLine(a, "fill:none;stroke:black;stroke-width:5", rmp);
    svgOut.printLine(p, "fill:none;stroke:red;stroke-width:5", rmp);
    svgOut.printLine(b, "fill:none;stroke:red;stroke-width:5", rmp);

    PolyLine<double> ptest;

    ptest << DPoint(90, 100);
    ptest << DPoint(100, 100);
    ptest << DPoint(95, 102);
    ptest << DPoint(101, 103);
    ptest << DPoint(121, 103);

    ptest.smoothenOutliers(10);
    svgOut.printLine(ptest, "fill:none;stroke:red;stroke-width:1", rmp);

    PolyLine<double> ptest2;

    ptest2 << DPoint(130, 100);
    ptest2 << DPoint(190, 100);
    ptest2 << DPoint(192, 180);
    ptest2 << DPoint(210, 103);
    ptest2 << DPoint(280, 193);
    ptest2 << DPoint(290, 183);

    svgOut.printLine(ptest2, "fill:none;stroke:red;stroke-width:1", rmp);

    ptest2.offsetPerp(18);
    svgOut.printLine(ptest2, "fill:none;stroke:red;stroke-width:1", rmp);

    for (size_t i = 1; i < 5; i++) {
      PolyLine<double> pl = avg;

      pl.offsetPerp(18 * i);
      if (i % 2)
        svgOut.printLine(pl, "fill:none;stroke:green;stroke-width:15", rmp);
      else
        svgOut.printLine(pl, "fill:none;stroke:red;stroke-width:15", rmp);
    }

    for (int i = 1; i < 5; i++) {
      PolyLine<double> pl = avg;

      pl.offsetPerp(-18 * i);
      if (i % 2)
        svgOut.printLine(pl, "fill:none;stroke:blue;stroke-width:15", rmp);
      else
        svgOut.printLine(pl, "fill:none;stroke:orange;stroke-width:15", rmp);
    }

    PolyLine<double> ptest3;
    PolyLine<double> ptest4;

    ptest3 << DPoint(130, 300);
    ptest3 << DPoint(150, 300);
    ptest3 << DPoint(160, 305);
    ptest3 << DPoint(165, 310);
    ptest3 << DPoint(170, 320);
    ptest3 << DPoint(170, 353);

    ptest4 = ptest3;

    svgOut.printLine(ptest3, "fill:none;stroke:red;stroke-width:1", rmp);

    ptest3.applyChaikinSmooth(1);
    ptest3.offsetPerp(39);
    svgOut.printLine(ptest3, "fill:none;stroke:red;stroke-width:1", rmp);

    PolyLine<double> ptest5;

    ptest5 << DPoint(210, 300);
    ptest5 << DPoint(230, 300);
    ptest5 << DPoint(250, 300);
    ptest5 << DPoint(260, 305);
    ptest5 << DPoint(265, 310);
    ptest5 << DPoint(270, 320);
    ptest5 << DPoint(270, 353);
    ptest5 << DPoint(240, 353);
    ptest5 << DPoint(240, 253);
    std::cout << ptest5.getLine().size() << std::endl;
    ptest5.fixTopology(100);
    std::cout << ptest5.getLine().size() << std::endl;
    svgOut.printLine(ptest5, "fill:none;stroke:red;stroke-width:1", rmp);

    std::cout << "Intersects: " << intersects(DPoint(4.914, 8.505), DPoint(7.316, 9.094),
                                              DPoint(12.198, 10.008), DPoint(14.676, 10.332))
              << std::endl;

    PolyLine<double> ptest6;
    PolyLine<double> ptest7;

    ptest6 << DPoint(500, 500) << DPoint(600, 600);
    ptest7 << DPoint(500, 480) << DPoint(600, 620);

    svgOut.printLine(ptest6, "fill:none;stroke:red;stroke-width:1", rmp);
    svgOut.printLine(ptest7, "fill:none;stroke:blue;stroke-width:1", rmp);

    SharedSegments<double> ss = ptest6.getSharedSegments(ptest7, 2);
    for (auto seg : ss.segments) {
      svgOut.printPoint(seg.first.first.p, "fill:green;r:4", rmp);
      svgOut.printPoint(seg.second.first.p, "fill:green;r:4", rmp);
    }

    PolyLine<double> ptest8;
    PolyLine<double> ptest9;

    ptest8 << DPoint(700, 300) << DPoint(800, 400);
    ptest9 << DPoint(700, 280) << DPoint(750, 350) << DPoint(800, 380);

    svgOut.printLine(ptest8, "fill:none;stroke:red;stroke-width:1", rmp);
    svgOut.printLine(ptest9, "fill:none;stroke:blue;stroke-width:1", rmp);

    ss = ptest8.getSharedSegments(ptest9, 2);
    for (auto seg : ss.segments) {
      svgOut.printPoint(seg.first.first.p, "fill:green;r:4", rmp);
      svgOut.printPoint(seg.second.first.p, "fill:green;r:4", rmp);
    }

    o << "</svg>";
  }

  return (0);
}
