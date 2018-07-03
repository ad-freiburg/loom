// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_JSON_JSONWRITER_H_
#define UTIL_JSON_JSONWRITER_H_

#include <map>
#include <ostream>
#include <stack>
#include <string>

namespace util {
namespace json {

class JsonWriterException : public std::exception {
 public:
  JsonWriterException(std::string msg) : _msg(msg) {}
  ~JsonWriterException() throw() {}

  virtual const char* what() const throw() { return _msg.c_str(); };

 private:
  std::string _msg;
};

// simple JSON writer class without much overhead
class JsonWriter {
 public:
  explicit JsonWriter(std::ostream* out);
  JsonWriter(std::ostream* out, size_t prec);
  JsonWriter(std::ostream* out, size_t prec, bool pretty);
  JsonWriter(std::ostream* out, size_t prec, bool pretty, size_t indent);
  ~JsonWriter(){};

  void obj();
  void obj(const std::map<std::string, std::string> kvs);
  void arr();
  void key(const std::string& k);
  void val(const std::string& v);
  void val(int v);
  void val(double v);
  void keyVal(const std::string& k, const std::string& v);

  void close();
  void closeAll();

 private:
  std::ostream* _out;

  enum JSON_NODE_T { OBJ, ARR, KEY };

  struct JsonNode {
    JSON_NODE_T type;
    bool empty;
  };

  std::stack<JsonNode> _stack;

  bool _pretty;
  size_t _indent;

  void valCheck();
  void prettor();
};

}  // namespace json
}  // namespace util

#endif  // UTIL_JSON_JSONWRITER_H_
