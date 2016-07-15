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
#include "gtfsparser/parser.h"
#include "graph/transitgraph.h"
#include "graph/node.h"
#include "util/XmlWriter.cpp"
#include "output/svgoutput.h"
#include "geo/PolyLine.h"

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
    try {
      LOG(INFO) << "reading feed at " << argv[1];
      parser.parse(&feed, argv[1]);

      graph::TransitGraph g("shinygraph");

      g.addEdge(g.addNode(new graph::Node(50, 30)), g.addNode(new graph::Node(100, 80)));

      graph::Node* a = g.addNode(new graph::Node(244, 600));
      g.addEdge(a, a);

      std::ofstream o;
      o.open("/home/patrick/test.svg");
      output::SvgOutput svgOut(&o);
      svgOut.print(g);
    } catch (gtfsparser::ParserException &e) {
      LOG(ERROR) << "Feed could not be parsed. " << e.what();
    }
  } else {
    transitmapper::geo::PolyLine p;

    for (double i = 0; i < 2000; i += 10) {
      // sin
      p << util::geo::Point(i, 1000 + sin(0.005 * i) * 500);

      // linear
      // p << util::geo::Point(i, 1000);
    }

    std::ofstream o;
    o.open("/home/patrick/test.svg");
    output::SvgOutput svgOut(&o);

    o << "<svg width=\"2000\" height=\"2000\">";

    svgOut.printLine(p, "fill:none;stroke:black;stroke-width:10");

    for (size_t i = 24; i < 25; i++) {
      transitmapper::geo::PolyLine pl = p;

      pl.offsetPerp(12 * i);
      if (i % 2) svgOut.printLine(pl, "fill:none;stroke:green;stroke-width:10");
      else svgOut.printLine(pl, "fill:none;stroke:red;stroke-width:10");
    }

    /**
    for (int i = 1; i < 15; i++) {
      transitmapper::geo::PolyLine pl = p;

      pl.offsetPerp(-12 * i);
      if (i % 2) svgOut.printLine(pl, "fill:none;stroke:blue;stroke-width:10");
      else svgOut.printLine(pl, "fill:none;stroke:orange;stroke-width:10");
    }
    **/
    o << "</svg>";

  }


  return(0);
}
