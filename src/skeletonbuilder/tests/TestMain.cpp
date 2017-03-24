// Copyright 2016
// Author: Patrick Brosi
//

#include <string>
#include "./lest.h"
#include "./../graph/Graph.h"
#include "./../graph/Node.h"
#include "./../graph/Edge.h"

using lest::approx;
using skeletonbuilder::graph::Graph;
using skeletonbuilder::graph::Edge;
using skeletonbuilder::graph::Node;


const std::string proj = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0"
  " +lon_0=0.0 +x_0=0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs";


// define LEST cases
const lest::test specification[] = {

// ___________________________________________________________________________
CASE("node combine") {
  EXPECT(true);

  skeletonbuilder::graph::Graph g("testgraph", proj);

  Node* a = g.addNode(new Node(0, 0));
  Node* b = g.addNode(new Node(1, 1));

  Node* f = g.addNode(new Node(0, 0));

  Edge* connEdge = g.addEdge(a, b);
  Edge* af = g.addEdge(a, f);

}

};

// _____________________________________________________________________________
int main(int argc, char** argv) {
  return(lest::run(specification, argc, argv));
}
