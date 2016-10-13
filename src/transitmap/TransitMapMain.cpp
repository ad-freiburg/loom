// Copyright 2016
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <stdio.h>
#include "log/Log.h"
#include "gtfsparser/Parser.h"
#include "./graph/TransitGraph.h"
#include "./graph/GraphBuilder.h"
#include "./graph/Node.h"
#include "./util/XmlWriter.cpp"
#include "./output/SvgOutput.h"
#include "./output/OgrOutput.h"
#include "./geo/PolyLine.h"
#include "./gtfsparser/gtfs/Service.h"
#include "./optim/EdgeOrderOptimizer.h"
#include "./config/ConfigReader.cpp"
#include "./config/TransitMapConfig.h"

using namespace transitmapper;
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

  // parse an example feed
  gtfsparser::Parser parser;
  gtfsparser::gtfs::Feed feed;

  LOG(INFO) << "Resolution is " << cfg.outputResolution << std::endl;
  LOG(INFO) << cfg.lineWidth << std::endl;
  LOG(INFO) << cfg.inputFeedPath << std::endl;
  LOG(INFO) << cfg.lineSpacing << std::endl;
  LOG(INFO) << cfg.renderMethod << std::endl;
  LOG(INFO) << cfg.optimIterations << std::endl;
  LOG(INFO) << cfg.renderStations << std::endl;
  LOG(INFO) << cfg.renderStationNames << std::endl;

  if (!cfg.inputFeedPath.empty()) {
    LOG(INFO) << "reading feed at " << cfg.inputFeedPath << std::endl;
    parser.parse(&feed, cfg.inputFeedPath);

    graph::TransitGraph g("shinygraph", cfg.projectionString);
    graph::GraphBuilder b(&g, &cfg);

    LOG(INFO) << "Building graph..." << std::endl;
    b.consume(feed);

    LOG(INFO) << "Simplyfing..." << std::endl;
    b.simplify();

    LOG(INFO) << "Building topological nodes..." << std::endl;
    while (b.createTopologicalNodes()) {}

    LOG(INFO) << "Averaging node positions" << std::endl;
    b.averageNodePositions();

    b.removeArtifacts();

    LOG(INFO) << "Creating node fronts..." << std::endl;
    b.writeMainDirs();
    b.expandOverlappinFronts();

    b.removeArtifacts();

    LOG(INFO) << "Writing initial ordering configuration..." << std::endl;
    b.writeInitialConfig();

    LOG(INFO) << "Optimizing..." << std::endl;
    optim::EdgeOrderOptimizer eoOptim(&g);
    eoOptim.optimize(cfg.optimIterations);
    LOG(INFO) << "Total graph score is -- " << g.getScore() << " --" << std::endl;
    LOG(INFO) << "Per node graph score is -- "
      << g.getScore() / g.getNodes()->size() << " --" << std::endl;

    if (cfg.renderMethod == "ogr") {
      std::string path = cfg.outputPath;
      LOG(INFO) << "Outputting to OGR " << path << " ..."
        << std::endl;
      output::OgrOutput ogrOut(path, &cfg);
      ogrOut.print(g);
    }

    if (cfg.renderMethod == "svg") {
      std::string path = cfg.outputPath;
      LOG(INFO) << "Outputting to SVG " << path << " ..."
        << std::endl;
      std::ofstream o;
      o.open(path);
      output::SvgOutput svgOut(&o, &cfg);
      svgOut.print(g);
    }

    LOG(INFO) << "World file for this map: " << std::endl << std::endl;

    std::cout << 1/cfg.outputResolution << std::endl
      << 0 << std::endl << 0 << std::endl
      << -1/cfg.outputResolution << std::endl
      << std::fixed << g.getBoundingBox().min_corner().get<0>() << std::endl
      << g.getBoundingBox().max_corner().get<1>() << std::endl;
  } else if (false) {
    transitmapper::geo::PolyLine a;
    transitmapper::geo::PolyLine b;

    a << util::geo::Point(1, 1);
    a << util::geo::Point(1, 2);
    a << util::geo::Point(1, 3);
    a << util::geo::Point(1, 4);
    a << util::geo::Point(1, 5);

    b << util::geo::Point(2, 3);
    b << util::geo::Point(2, 4);
    b << util::geo::Point(2, 5);
    b << util::geo::Point(2, 6);
    b << util::geo::Point(2, 7);

    transitmapper::geo::SharedSegments ss = a.getSharedSegments(b, 3);
    LOG(INFO) << "RESULT: " << ss.segments.size() << std::endl;

  } else {

    // just testing parallel drawing...
    transitmapper::geo::PolyLine p;
    transitmapper::geo::PolyLine a;
    transitmapper::geo::PolyLine b;
    transitmapper::geo::PolyLine c;

    for (double i = 0; i < 2000; i += 20) {
      // sin
      a << util::geo::Point(i, 1000 + sin(0.006 * i) * 500);

      // cos
      b << util::geo::Point(i, 1000 + cos(0.003 * i) * 500);

      // linear
      //a << util::geo::Point(i, 1000);

      // triangle
      if (i < 1000) {
        p << util::geo::Point(i, 1000 - (i/2));
      } else {
        p << util::geo::Point(i, (i/2));
      }
    }

    std::ofstream o;
    o.open("/home/patrick/test.svg");
    cfg.outputResolution = 1;
    output::SvgOutput svgOut(&o, &cfg);

    o << "<svg width=\"2000\" height=\"2000\">";


    util::geo::Point pp(500,470);
    util::geo::Point ppp(550,500);
    c << util::geo::Point(400, 400) << util::geo::Point(450, 450) << util::geo::Point(550, 450) << util::geo::Point( 600, 700);

    c = c.getSegment(pp, ppp);

    util::geo::Point projected = c.projectOn(pp).p;
    util::geo::Point projectedB = c.projectOn(ppp).p;

    svgOut.printLine(c, "fill:none;stroke:red;stroke-width:5", 2000, 2000, 0, 0);
    svgOut.printPoint(pp, "fill:none;stroke:red;stroke-width:5", 2000, 2000, 0, 0);
    svgOut.printPoint(ppp, "fill:none;stroke:red;stroke-width:5", 2000, 2000, 0, 0);
    svgOut.printPoint(projected, "fill:none;stroke:red;stroke-width:5", 2000, 2000, 0, 0);
    svgOut.printPoint(projectedB, "fill:none;stroke:red;stroke-width:5", 2000, 2000, 0, 0);

    //svgOut.printLine(p, "fill:none;stroke:black;stroke-width:10");

    std::vector<const transitmapper::geo::PolyLine*> ls;
    ls.push_back(&a);
    ls.push_back(&p);
    ls.push_back(&b);
    transitmapper::geo::PolyLine avg = transitmapper::geo::PolyLine::average(ls);

    svgOut.printLine(avg, "fill:none;stroke:#8844AA;stroke-width:15", 2000, 2000, 0, 0);

    svgOut.printLine(a, "fill:none;stroke:black;stroke-width:5", 2000, 2000, 0, 0);
    svgOut.printLine(p, "fill:none;stroke:red;stroke-width:5", 2000, 2000, 0, 0);
    svgOut.printLine(b, "fill:none;stroke:red;stroke-width:5", 2000, 2000, 0, 0);

    transitmapper::geo::PolyLine ptest;

    ptest << util::geo::Point(90, 100);
    ptest << util::geo::Point(100, 100);
    ptest << util::geo::Point(95, 102);
    ptest << util::geo::Point(101, 103);
    ptest << util::geo::Point(121, 103);

    ptest.smoothenOutliers(10);
    svgOut.printLine(ptest, "fill:none;stroke:red;stroke-width:1", 2000, 2000, 0, 0);

    transitmapper::geo::PolyLine ptest2;

    ptest2 << util::geo::Point(130, 100);
    ptest2 << util::geo::Point(190, 100);
    ptest2 << util::geo::Point(192, 180);
    ptest2 << util::geo::Point(210, 103);
    ptest2 << util::geo::Point(280, 193);
    ptest2 << util::geo::Point(290, 183);

    svgOut.printLine(ptest2, "fill:none;stroke:red;stroke-width:1", 2000, 2000, 0, 0);

    ptest2.offsetPerp(18);
    svgOut.printLine(ptest2, "fill:none;stroke:red;stroke-width:1", 2000, 2000, 0, 0);

    for (size_t i = 1; i < 5; i++) {
      transitmapper::geo::PolyLine pl = avg;

      pl.offsetPerp(18 * i);
      if (i % 2) svgOut.printLine(pl, "fill:none;stroke:green;stroke-width:15", 2000, 2000, 0, 0);
      else svgOut.printLine(pl, "fill:none;stroke:red;stroke-width:15", 2000, 2000, 0, 0);
    }

    for (int i = 1; i < 5; i++) {
      transitmapper::geo::PolyLine pl = avg;

      pl.offsetPerp(-18 * i);
      if (i % 2) svgOut.printLine(pl, "fill:none;stroke:blue;stroke-width:15", 2000, 2000, 0, 0);
      else svgOut.printLine(pl, "fill:none;stroke:orange;stroke-width:15", 2000, 2000, 0, 0);
    }

    transitmapper::geo::PolyLine ptest3;
    transitmapper::geo::PolyLine ptest4;

    ptest3 << util::geo::Point(130, 300);
    ptest3 << util::geo::Point(150, 300);
    ptest3 << util::geo::Point(160, 305);
    ptest3 << util::geo::Point(165, 310);
    ptest3 << util::geo::Point(170, 320);
    ptest3 << util::geo::Point(170, 353);


    ptest4 = ptest3;

    svgOut.printLine(ptest3, "fill:none;stroke:red;stroke-width:1", 2000, 2000, 0, 0);

    ptest3.applyChaikinSmooth(1);
    ptest3.offsetPerp(39);
    svgOut.printLine(ptest3, "fill:none;stroke:red;stroke-width:1", 2000, 2000, 0, 0);

    transitmapper::geo::PolyLine ptest5;

    ptest5 << util::geo::Point(210, 300);
    ptest5 << util::geo::Point(230, 300);
    ptest5 << util::geo::Point(250, 300);
    ptest5 << util::geo::Point(260, 305);
    ptest5 << util::geo::Point(265, 310);
    ptest5 << util::geo::Point(270, 320);
    ptest5 << util::geo::Point(270, 353);
    ptest5 << util::geo::Point(240, 353);
    ptest5 << util::geo::Point(240, 253);
    std::cout << ptest5.getLine().size() << std::endl;
    ptest5.fixTopology(100);
    std::cout << ptest5.getLine().size() << std::endl;
    svgOut.printLine(ptest5, "fill:none;stroke:red;stroke-width:1", 2000, 2000, 0, 0);

    std::cout << "Intersects: " << util::geo::intersects(4.914, 8.505, 7.316, 9.094, 12.198, 10.008, 14.676, 10.332) << std::endl;

    o << "</svg>";

  }


  return(0);
}
