// Copyright 2016
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#define _ELPP_NO_DEFAULT_LOG_FILE
#define _ELPP_DISABLE_LOG_FILE_FROM_ARG

#include <unistd.h>
#include <iostream>
#include <string>
#include <set>
#include <stdio.h>
#include "easylogging/easylogging.h"
#include "gtfsparser/Parser.h"
#include "graph/TransitGraph.h"
#include "graph/GraphBuilder.h"
#include "graph/Node.h"
#include "util/XmlWriter.cpp"
#include "output/SvgOutput.h"
#include "geo/PolyLine.h"
#include "gtfsparser/gtfs/Service.h"
#include "optim/EdgeOrderOptimizer.h"

INITIALIZE_EASYLOGGINGPP;

using namespace transitmapper;
using std::string;

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  el::Loggers::reconfigureAllLoggers(
    el::ConfigurationType::Format, "[%datetime] %level: %msg");

  el::Loggers::reconfigureAllLoggers(
    el::ConfigurationType::ToFile, "false");

  // initialize easylogging lib
  START_EASYLOGGINGPP(argc, argv);

  // parse an example feed
  gtfsparser::Parser parser;

  gtfsparser::gtfs::Feed feed;

  if (argc > 1) {
    LOG(INFO) << "reading feed at " << argv[1];
    parser.parse(&feed, argv[1]);

    graph::TransitGraph g("shinygraph", "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs");

    graph::GraphBuilder b(&g);

    LOG(INFO) << "Building graph...";
    b.consume(feed);

    LOG(INFO) << "Simplyfing...";
    b.simplify();

    LOG(INFO) << "Building topological nodes...";
    b.createTopologicalNodes();

    LOG(INFO) << "Averaging node positions";
    b.averageNodePositions();

    b.setLineWidth(5, 20);

    LOG(INFO) << "Creating node fronts...";
    b.writeMainDirs();
    b.freeNodes();


    LOG(INFO) << "Optimizing...";
    optim::EdgeOrderOptimizer eoOptim(&g);
    //eoOptim.optimize();

    LOG(INFO) << "Graph score is -- " << g.getScore() << " --";

    LOG(INFO) << "Outputting to SVG...";
    std::ofstream o;
    o.open("/home/patrick/test.svg");

    LOG(INFO) << "outputting to SVG...";
    output::SvgOutput svgOut(&o, .3);
    svgOut.print(g);
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
    output::SvgOutput svgOut(&o, 1);

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
    o << "</svg>";

  }


  return(0);
}
