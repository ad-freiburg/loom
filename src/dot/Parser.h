// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef DOT_PARSER_H_
#define DOT_PARSER_H_

#include <fstream>
#include <map>
#include <vector>
#include <string>

namespace dot {
namespace parser {

enum State {
  NONE,
  KW_STRICT,
  KW_GRAPH,
  AW_KW_GRAPH,
  KW_DIRGRAPH,
  AW_GRAPH_ID,
  IN_GRAPH_ID,
  AW_LL_ID,
  IN_ID,
  IN_QUOTED_ID,
  AW_OPEN,
  AW_STMT_DEC,
  AW_ID,
  AW_EDGE_OP,
  AW_ATTR_KEY,
  IN_QUOTED_GRAPH_ID,
  IN_ATTR_KEY,
  AW_ATTR_ASSIGN,
  IN_QUOTED_ATTR_KEY,
  AW_ATTR_VAL,
  IN_ATTR_VAL,
  IN_QUOTED_ATTR_VAL,
  AW_CONT_ATTR_KEY
};

enum EntityType {
  EMPTY,
  EDGE,
  NODE,
  ATTR,
  ATTR_NODE,
  ATTR_EDGE,
  ATTR_GRAPH
};

enum GraphType {
  STRICT_GRAPH,
  STRICT_DIGRAPH,
  GRAPH,
  DIGRAPH
};

struct Entity {
  EntityType type;
  std::vector<std::string> ids;
  std::map<std::string, std::string> attrs;

  GraphType graphType;
  std::string graphName;
  size_t level;
};

class Parser {
 public:
  Parser(std::istream* s);
  const Entity& get();

  bool has();

 private:
  std::istream* _is;
  State _s;
  Entity _ret;
  char _c;
  size_t _level;
  bool _has;
  std::string tmp;
  std::string tmp2;

  bool isIDChar(const char& c) const;
};

}  // graph
}  // octi

#endif  // DOT_PARSER_H_
