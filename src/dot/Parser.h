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
  IN_TAG_NAME,
  IN_TAG_NAME_META,
  IN_TAG,
  IN_TAG_CLOSE,
  IN_TAG_NAME_CLOSE,
  IN_TAG_TENTATIVE,
  IN_ATTRKEY,
  AFTER_ATTRKEY,
  AW_IN_ATTRVAL,
  IN_ATTRVAL_SQ,
  IN_ATTRVAL_DQ,
  IN_TEXT,
  AW_CLOSING,
  WS_SKIP
};

enum EntityType {
  EDGE,
  NODE
};

typedef std::pair<std::string, std::string> Attribute;

struct Entity {
  EntityType type;
  std::vector<std::string> ids;
  std::vector<Attribute> attrs;
};

class Parser {
 public:
  Parser(std::istream* s);

  bool has();

 private:
  State _s;
};

}  // graph
}  // octi

#endif  // DOT_PARSER_H_
