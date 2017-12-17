// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include "dot/Parser.h"

using namespace dot;
using namespace parser;

// _____________________________________________________________________________
Parser::Parser(std::istream* is)
    : _is(is), _s(NONE), _level(std::numeric_limits<size_t>::max()) {}

// _____________________________________________________________________________
bool Parser::has() { return _level > 0; }

// _____________________________________________________________________________
const Entity& Parser::get() {
  _ret.ids.clear();
  _ret.attrs.clear();
  _ret.type = EMPTY;
  _level = _level == std::numeric_limits<size_t>::max() ? 0 : _level;

  while (_is->get(_c)) {
    // std::cout << "state: " << _s << " level: " << _level << " char: " << _c << std::endl;
    _ret.level = _level;
    switch (_s) {
      case NONE:
        if (std::isspace(_c)) continue;
        tmp.clear();
        if (_c == 'd' || _c == 'D') {
          _s = KW_DIRGRAPH;
          continue;
        }
        if (_c == 'g' || _c == 'G') {
          _s = KW_GRAPH;
          continue;
        }
        if (_c == 's' || _c == 'S') {
          _s = KW_STRICT;
          continue;
        }

        std::cerr << "Expected keywords 'strict', 'graph' or 'digraph'"
                  << std::endl;
        exit(1);
      case KW_STRICT:
        tmp += _c;
        while (_is->get(_c)) {
          if (std::isalpha(_c)) {
            tmp += std::tolower(_c);
          } else
            break;
        }

        if (tmp == "trict" && std::isspace(_c)) {
          _ret.graphType = STRICT_GRAPH;
          tmp.clear();
          _s = AW_KW_GRAPH;
          continue;
        }

        std::cerr << "Expected keyword 'strict'" << std::endl;
        exit(1);
      case KW_GRAPH:
        tmp += _c;
        while (_is->get(_c)) {
          if (std::isalpha(_c)) {
            tmp += std::tolower(_c);
            if (tmp == "raph") break;
          } else
            break;
        }

        if (tmp == "raph") {
          if (_ret.graphType == STRICT_GRAPH) {
            _ret.graphType = STRICT_GRAPH;
          } else {
            _ret.graphType = GRAPH;
          }
          tmp.clear();
          _s = AW_GRAPH_ID;
          continue;
        }

        std::cerr << "Expected keyword 'graph'" << std::endl;
        exit(1);
      case KW_DIRGRAPH:
        tmp += _c;
        while (_is->get(_c)) {
          if (std::isalpha(_c)) {
            tmp += std::tolower(_c);
            if (tmp == "igraph") break;
          } else
            break;
        }

        if (tmp == "igraph") {
          if (_ret.graphType == STRICT_GRAPH) {
            _ret.graphType = STRICT_DIGRAPH;
          } else {
            _ret.graphType = DIGRAPH;
          }
          tmp.clear();
          _s = AW_GRAPH_ID;
          continue;
        }

        std::cerr << tmp << std::endl;

        std::cerr << "Expected keyword 'digraph'" << std::endl;
        exit(1);
      case AW_KW_GRAPH:
        if (std::isspace(_c)) continue;
        if (_c == 'g' || _c == 'G') {
          _s = KW_GRAPH;
          continue;
        }
        if (_c == 'd' || _c == 'D') {
          _s = KW_DIRGRAPH;
          continue;
        }
      case AW_GRAPH_ID:
        if (std::isspace(_c)) continue;
        if (isIDChar(_c)) {
          tmp += _c;
          _s = IN_GRAPH_ID;
          continue;
        } else if (_c == '"') {
          _s = IN_QUOTED_GRAPH_ID;
          continue;
        } else if (_c == '{') {
          _s = AW_LL_ID;
          continue;
        }

        std::cerr << "Expected graph id or opening {" << std::endl;
        exit(1);
      case IN_QUOTED_GRAPH_ID:
        if (_c != '"' || tmp.back() == '\\') {
          tmp += _c;
          continue;
        } else {
          _ret.graphName = tmp;
          tmp.clear();
          _s = AW_OPEN;
          continue;
        }
      case IN_GRAPH_ID:
        if (isIDChar(_c)) {
          tmp += _c;
          continue;
        }

        if (_c == '{') {
          _ret.graphName = tmp;
          tmp.clear();
          _s = AW_LL_ID;
          continue;
        }

        if (std::isspace(_c)) {
          _ret.graphName = tmp;
          tmp.clear();
          _s = AW_OPEN;
          continue;
        }

        std::cerr << "Expected ID characters" << std::endl;
        exit(1);
      case AW_OPEN:
        if (std::isspace(_c)) continue;
        if (_c == '{') {
          _level += 1;
          tmp.clear();
          _s = AW_LL_ID;
          continue;
        }
        std::cerr << "Expected opening {" << std::endl;
        exit(1);
      case AW_LL_ID:
        if (std::isspace(_c)) continue;
        if (isIDChar(_c)) {
          tmp += _c;
          _s = IN_ID;
          continue;
        }
        if (_c == '"') {
          _s = IN_QUOTED_ID;
          continue;
        }
        if (_c == '}') {
          _level -= 1;
          return _ret;
        }
        if (_c == '{') {
          _level += 1;
          continue;
        }
        if (_c == ';') {
          continue;
        }

        std::cerr << "Expected statement" << std::endl;
        exit(1);

      case AW_ID:
        if (std::isspace(_c)) continue;
        if (isIDChar(_c)) {
          tmp += _c;
          _s = IN_ID;
          continue;
        }
        if (_c == '"') {
          _s = IN_QUOTED_ID;
          continue;
        }

        if (_c == '{') {
          std::cerr << "Subgraphs in edge statements not yet supported."
                    << std::endl;
          exit(1);
        }

        std::cerr << "Expected ID" << std::endl;
        exit(1);

      case IN_ID:
        if (isIDChar(_c)) {
          tmp += _c;
          continue;
        }

        if (tmp == "subgraph") {
          if (_ret.type != EMPTY) {
            std::cerr << "Subgraphs in edge statements not yet supported."
                      << std::endl;
            exit(1);
          }
          _s = AW_GRAPH_ID;
          tmp.clear();
          continue;
        }

        if (tmp == "graph") {
          if (_ret.type != EMPTY) {
            std::cerr << "'graph' is a reserved keyword"
                      << std::endl;
            exit(1);
          }
          _ret.type = ATTR_GRAPH;
        }

        if (tmp == "node") {
          if (_ret.type != EMPTY) {
            std::cerr << "'node' is a reserved keyword"
                      << std::endl;
            exit(1);
          }
          _ret.type = ATTR_NODE;
        }

        if (tmp == "edge") {
          if (_ret.type != EMPTY) {
            std::cerr << "'edge' is a reserved keyword"
                      << std::endl;
            exit(1);
          }
          _ret.type = ATTR_EDGE;
        }

        _ret.ids.push_back(tmp);
        tmp.clear();

        _s = AW_STMT_DEC;
        // slip down


      case AW_STMT_DEC:
        if (std::isspace(_c)) continue;
        if (_c == '=') {
          if (_ret.type != EMPTY) {
            std::cerr << "Syntax error." << std::endl;
            exit(1);
          }
          tmp.clear();
          _ret.type = ATTR;
          _s = AW_ID;
          continue;
        }

        if (_c == '-') {
          if (_ret.type != EMPTY && _ret.type != EDGE) {
            std::cerr << "Syntax error." << std::endl;
            exit(1);
          }
          _s = AW_EDGE_OP;
          continue;
        }

        if (_c == '[') {
          if (_ret.type == EMPTY) _ret.type = NODE;
          _s = AW_ATTR_KEY;
          continue;
        }

        if (_c == '{') {
          _level += 1;
          _s = AW_LL_ID;
          if (_ret.type == EMPTY && _ret.ids.size()) _ret.type = NODE;
          return _ret;
        }

        if (_c == '"') {
          _s = IN_QUOTED_ID;
          if (_ret.type == EMPTY && _ret.ids.size()) _ret.type = NODE;
          return _ret;
        }

        if (isIDChar(_c)) {
          tmp.clear();
          tmp += _c;
          _s = IN_ID;
          return _ret;
        }

        if (_c == ';') {
          _s = AW_LL_ID;
          return _ret;
        }

        if (_c == '}') {
          _level -= 1;
          _s = AW_LL_ID;
          return _ret;
        }

        std::cerr << "Expected edge operator, id, = or [" << std::endl;
        exit(1);

      case AW_ATTR_KEY:
        if (std::isspace(_c)) continue;
        if (isIDChar(_c)) {
          tmp += _c;
          _s = IN_ATTR_KEY;
          continue;
        }
        if (_c == '"') {
          _s = IN_QUOTED_ATTR_KEY;
          continue;
        }
        if (_c == ']') {
          _s = AW_LL_ID;
          return _ret;
        }

        std::cerr << "Expected attribute key" << std::endl;
        exit(1);

      case IN_ATTR_KEY:
        if (isIDChar(_c)) {
          tmp += _c;
          continue;
        }

        _s = AW_ATTR_ASSIGN;

        if (_c == '=') {
          _s = AW_ATTR_VAL;
        }

        _ret.attrs[tmp] = "";
        continue;

      case IN_QUOTED_ATTR_KEY:
        if (_c != '"' || tmp.back() == '\\') {
          tmp += _c;
          continue;
        }

        _ret.attrs[tmp] = "";
        _s = AW_ATTR_ASSIGN;
        continue;

      case IN_ATTR_VAL:
        if (isIDChar(_c)) {
          tmp2 += _c;
          continue;
        }

        _ret.attrs[tmp] = tmp2;
        tmp.clear();
        tmp2.clear();
        _s = AW_CONT_ATTR_KEY;
        // slip through

      case AW_CONT_ATTR_KEY:
        if (std::isspace(_c)) continue;
        if (_c == ';' || _c == ',') {
          _s = AW_ATTR_KEY;
          continue;
        }
        if (isIDChar(_c)) {
          tmp += _c;
          _s = IN_ATTR_KEY;
          continue;
        }
        if (_c == '"') {
          _s = IN_QUOTED_ATTR_KEY;
          continue;
        }
        if (_c == ']') {
          _s = AW_LL_ID;
          return _ret;
        }

        std::cerr << "Expected attribute key" << std::endl;
        exit(1);

      case IN_QUOTED_ATTR_VAL:
        if (_c != '"' || tmp.back() == '\\') {
          tmp2 += _c;
          continue;
        }

        _ret.attrs[tmp] = tmp2;
        tmp.clear();
        tmp2.clear();
        _s = AW_CONT_ATTR_KEY;
        continue;

      case AW_ATTR_ASSIGN:
        if (std::isspace(_c)) continue;
        if (_c == '=') {
          _s = AW_ATTR_VAL;
          continue;
        }

        std::cerr << "Expected '='" << std::endl;
        exit(1);

      case AW_ATTR_VAL:
        if (std::isspace(_c)) continue;
        if (isIDChar(_c)) {
          tmp2 += _c;
          _s = IN_ATTR_VAL;
          continue;
        }
        if (_c == '"') {
          _s = IN_QUOTED_ATTR_VAL;
          continue;
        }

        std::cerr << "Expected attribute value" << std::endl;
        exit(1);

      case IN_QUOTED_ID:
        if (_c != '"' || tmp.back() == '\\') {
          tmp += _c;
          continue;
        } else {
          _ret.ids.push_back(tmp);
          tmp.clear();
          _s = AW_STMT_DEC;
          continue;
        }

      case AW_EDGE_OP:
        if (_c == '>') {
          if (_ret.graphType == STRICT_GRAPH || _ret.graphType == GRAPH) {
            std::cerr << "No directed edges allowed in undirected graph."
                      << std::endl;
            exit(1);
          }
          _ret.type = EDGE;
          _s = AW_ID;
          continue;
        }

        if (_c == '-') {
          if (_ret.graphType == STRICT_DIGRAPH || _ret.graphType == DIGRAPH) {
            std::cerr << "No undirected edges allowed in directed graph."
                      << std::endl;
            exit(1);
          }
          _ret.type = EDGE;
          _s = AW_ID;
          continue;
        }

        std::cerr << "Expected edge operator" << std::endl;
        exit(1);
    }
  }

  std::cerr << "Syntax error." << std::endl;
  exit(1);
}

// _____________________________________________________________________________
bool Parser::isIDChar(const char& c) const {
  return std::isalnum(c) || c == '_' || c == '.';
}
