// Copyright 2016
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "./config/ConfigReader.cpp"
#include "./config/TransitMapConfig.h"
#include "./graph/GraphBuilder.h"
#include "./graph/Node.h"
#include "./graph/TransitGraph.h"
#include "./optim/ILPEdgeOrderOptimizer.h"
#include "./output/OgrOutput.h"
#include "./output/SvgOutput.h"
#include "pbutil/geo/PolyLine.h"
#include "pbutil/log/Log.h"

using namespace transitmapper;
using namespace pbutil::geo;
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
    LOG(INFO) << "reading graph " << std::endl;
    graph::TransitGraph g(cfg.name, cfg.projectionString);
    graph::GraphBuilder b(&cfg);
    if (!b.build(&(std::cin), &g)) {
      exit(1);
    }

    if (cfg.collapseLinePartners) {
      b.combinePartnerRoutes(&g);
    }

    if (cfg.outputStats) {
      LOG(INFO) << "(stats) Stats for graph '" << g.getName() << std::endl;
      LOG(INFO) << "(stats)   Total node count: " << g.getNumNodes() << " ("
                << g.getNumNodes(true) << " topo, " << g.getNumNodes(false)
                << " non-topo)" << std::endl;
      LOG(INFO) << "(stats)   Total edge count: " << g.getNumEdges()
                << std::endl;
      LOG(INFO) << "(stats)   Total unique route count: " << g.getNumRoutes()
                << std::endl;
      LOG(INFO) << "(stats)   Max edge route cardinality: "
                << g.getMaxCardinality() << std::endl;
      LOG(INFO) << "(stats)   Number of poss. solutions: "
                << g.getNumPossSolutions() << std::endl;
    }

    LOG(INFO) << "Creating node fronts..." << std::endl;
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

    LOG(INFO) << "Writing initial ordering configuration..." << std::endl;
    b.writeInitialConfig(&g);

    LOG(INFO) << "Optimizing..." << std::endl;

    LOG(INFO) << "(stats) Total graph score BEFORE optim is -- "
              << g.getScore(cfg.inStationCrossPenalty, cfg.crossPenMultiSameSeg,
                            cfg.crossPenMultiDiffSeg, cfg.splitPenWeight)
              << " --" << std::endl;
    LOG(INFO) << "(stats)   Per node graph score: "
              << g.getScore(cfg.inStationCrossPenalty, cfg.crossPenMultiSameSeg,
                            cfg.crossPenMultiDiffSeg, cfg.splitPenWeight) /
                     g.getNodes()->size()
              << std::endl;
    LOG(INFO) << "(stats)   Crossings: " << g.getNumCrossings() << " (score: "
              << g.getCrossScore(cfg.inStationCrossPenalty,
                                 cfg.crossPenMultiSameSeg,
                                 cfg.crossPenMultiDiffSeg)
              << ")" << std::endl;
    LOG(INFO) << "(stats)   Separations: " << g.getNumSeparations()
              << " (score: "
              << g.getSeparationScore(cfg.inStationCrossPenalty,
                                      cfg.splitPenWeight)
              << ")" << std::endl;

    if (cfg.renderMethod != "ogr" && !cfg.noOptim) {
      if (cfg.optimMethod == "ilp_impr") {
        optim::ILPEdgeOrderOptimizer ilpEoOptim(&g, &cfg);
        ilpEoOptim.optimize();
      } else if (cfg.optimMethod == "ilp") {
        optim::ILPOptimizer ilpEoOptim(&g, &cfg);
        ilpEoOptim.optimize();
      }
    }

    LOG(INFO) << "(stats) Total graph score AFTER optim is -- "
              << g.getScore(cfg.inStationCrossPenalty, cfg.crossPenMultiSameSeg,
                            cfg.crossPenMultiDiffSeg, cfg.splitPenWeight)
              << " --" << std::endl;
    LOG(INFO) << "(stats)   Per node graph score: "
              << g.getScore(cfg.inStationCrossPenalty, cfg.crossPenMultiSameSeg,
                            cfg.crossPenMultiDiffSeg, cfg.splitPenWeight) /
                     g.getNodes()->size()
              << std::endl;
    LOG(INFO) << "(stats)   Crossings: " << g.getNumCrossings() << " (score: "
              << g.getCrossScore(cfg.inStationCrossPenalty,
                                 cfg.crossPenMultiSameSeg,
                                 cfg.crossPenMultiDiffSeg)
              << ")" << std::endl;
    LOG(INFO) << "(stats)   Separations: " << g.getNumSeparations()
              << " (score: "
              << g.getSeparationScore(cfg.inStationCrossPenalty,
                                      cfg.splitPenWeight)
              << ")" << std::endl;

    if (cfg.renderMethod == "ogr") {
      std::string path = cfg.outputPath;
      LOG(INFO) << "Outputting to OGR " << path << " ..." << std::endl;
      output::OgrOutput ogrOut(path, &cfg);
      ogrOut.print(g);
    }

    if (cfg.renderMethod == "svg") {
      std::string path = cfg.outputPath;
      LOG(INFO) << "Outputting to SVG " << path << " ..." << std::endl;
      std::ofstream o;
      o.open(path);
      output::SvgOutput svgOut(&o, &cfg);
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
        output::SvgOutput svgOut(&o, &cfg);
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
        output::SvgOutput svgOut(&o, &cfg);
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
        output::SvgOutput svgOut(&o, &cfg);
        svgOut.print(g);
      }
    }

    if (!cfg.worldFilePath.empty()) {
      LOG(INFO) << "Writing world file for this map to " << cfg.worldFilePath << std::endl;

      std::ofstream file;
      file.open(cfg.worldFilePath);
      if (file) {
        file << 1 / cfg.outputResolution << std::endl
                << 0 << std::endl
                << 0 << std::endl
                << -1 / cfg.outputResolution << std::endl
                << std::fixed << g.getBoundingBox(cfg.outputPadding).min_corner().get<0>()
                << std::endl
                << g.getBoundingBox(cfg.outputPadding).max_corner().get<1>() << std::endl;
        file.close();
      }
    }
  } else if (false) {
    PolyLine a;
    PolyLine b;

    a << Point(1, 1);
    a << Point(1, 2);
    a << Point(1, 3);
    a << Point(1, 4);
    a << Point(1, 5);

    b << Point(2, 3);
    b << Point(2, 4);
    b << Point(2, 5);
    b << Point(2, 6);
    b << Point(2, 7);

    SharedSegments ss = a.getSharedSegments(b, 3);
    LOG(INFO) << "RESULT: " << ss.segments.size() << std::endl;

  } else {
    // just testing parallel drawing...
    PolyLine p;
    PolyLine a;
    PolyLine b;
    PolyLine c;

    for (double i = 0; i < 2000; i += 20) {
      // sin
      a << Point(i, 1000 + sin(0.006 * i) * 500);

      // cos
      b << Point(i, 1000 + cos(0.003 * i) * 500);

      // linear
      // a << Point(i, 1000);

      // triangle
      if (i < 1000) {
        p << Point(i, 1000 - (i / 2));
      } else {
        p << Point(i, (i / 2));
      }
    }

    std::ofstream o;
    o.open("/home/patrick/test.svg");
    cfg.outputResolution = 1;
    output::SvgOutput svgOut(&o, &cfg);

    o << "<svg width=\"2000\" height=\"2000\">";

    Point pp(500, 470);
    Point ppp(550, 500);
    c << Point(400, 400) << Point(450, 450) << Point(550, 450)
      << Point(600, 700);

    c = c.getSegment(pp, ppp);

    Point projected = c.projectOn(pp).p;
    Point projectedB = c.projectOn(ppp).p;

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

    std::vector<const PolyLine*> ls;
    ls.push_back(&a);
    ls.push_back(&p);
    ls.push_back(&b);
    PolyLine avg = PolyLine::average(ls);

    svgOut.printLine(avg, "fill:none;stroke:#8844AA;stroke-width:15", rmp);

    svgOut.printLine(a, "fill:none;stroke:black;stroke-width:5", rmp);
    svgOut.printLine(p, "fill:none;stroke:red;stroke-width:5", rmp);
    svgOut.printLine(b, "fill:none;stroke:red;stroke-width:5", rmp);

    PolyLine ptest;

    ptest << Point(90, 100);
    ptest << Point(100, 100);
    ptest << Point(95, 102);
    ptest << Point(101, 103);
    ptest << Point(121, 103);

    ptest.smoothenOutliers(10);
    svgOut.printLine(ptest, "fill:none;stroke:red;stroke-width:1", rmp);

    PolyLine ptest2;

    ptest2 << Point(130, 100);
    ptest2 << Point(190, 100);
    ptest2 << Point(192, 180);
    ptest2 << Point(210, 103);
    ptest2 << Point(280, 193);
    ptest2 << Point(290, 183);

    svgOut.printLine(ptest2, "fill:none;stroke:red;stroke-width:1", rmp);

    ptest2.offsetPerp(18);
    svgOut.printLine(ptest2, "fill:none;stroke:red;stroke-width:1", rmp);

    for (size_t i = 1; i < 5; i++) {
      PolyLine pl = avg;

      pl.offsetPerp(18 * i);
      if (i % 2)
        svgOut.printLine(pl, "fill:none;stroke:green;stroke-width:15", rmp);
      else
        svgOut.printLine(pl, "fill:none;stroke:red;stroke-width:15", rmp);
    }

    for (int i = 1; i < 5; i++) {
      PolyLine pl = avg;

      pl.offsetPerp(-18 * i);
      if (i % 2)
        svgOut.printLine(pl, "fill:none;stroke:blue;stroke-width:15", rmp);
      else
        svgOut.printLine(pl, "fill:none;stroke:orange;stroke-width:15", rmp);
    }

    PolyLine ptest3;
    PolyLine ptest4;

    ptest3 << Point(130, 300);
    ptest3 << Point(150, 300);
    ptest3 << Point(160, 305);
    ptest3 << Point(165, 310);
    ptest3 << Point(170, 320);
    ptest3 << Point(170, 353);

    ptest4 = ptest3;

    svgOut.printLine(ptest3, "fill:none;stroke:red;stroke-width:1", rmp);

    ptest3.applyChaikinSmooth(1);
    ptest3.offsetPerp(39);
    svgOut.printLine(ptest3, "fill:none;stroke:red;stroke-width:1", rmp);

    PolyLine ptest5;

    ptest5 << Point(210, 300);
    ptest5 << Point(230, 300);
    ptest5 << Point(250, 300);
    ptest5 << Point(260, 305);
    ptest5 << Point(265, 310);
    ptest5 << Point(270, 320);
    ptest5 << Point(270, 353);
    ptest5 << Point(240, 353);
    ptest5 << Point(240, 253);
    std::cout << ptest5.getLine().size() << std::endl;
    ptest5.fixTopology(100);
    std::cout << ptest5.getLine().size() << std::endl;
    svgOut.printLine(ptest5, "fill:none;stroke:red;stroke-width:1", rmp);

    std::cout << "Intersects: " << intersects(4.914, 8.505, 7.316, 9.094,
                                              12.198, 10.008, 14.676, 10.332)
              << std::endl;

    PolyLine ptest6;
    PolyLine ptest7;

    ptest6 << Point(500, 500) << Point(600, 600);
    ptest7 << Point(500, 480) << Point(600, 620);

    svgOut.printLine(ptest6, "fill:none;stroke:red;stroke-width:1", rmp);
    svgOut.printLine(ptest7, "fill:none;stroke:blue;stroke-width:1", rmp);

    SharedSegments ss = ptest6.getSharedSegments(ptest7, 2);
    for (auto seg : ss.segments) {
      svgOut.printPoint(seg.first.first.p, "fill:green;r:4", rmp);
      svgOut.printPoint(seg.second.first.p, "fill:green;r:4", rmp);
    }

    PolyLine ptest8;
    PolyLine ptest9;

    ptest8 << Point(700, 300) << Point(800, 400);
    ptest9 << Point(700, 280) << Point(750, 350) << Point(800, 380);

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
