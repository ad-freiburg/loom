// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_JSON_WRITER_H_
#define UTIL_JSON_WRITER_H_

#include <map>
#include <ostream>
#include <stack>
#include <string>
#include <vector>

namespace util {
namespace json {

class WriterException : public std::exception {
 public:
  WriterException(std::string msg) : _msg(msg) {}
  ~WriterException() throw() {}

  virtual const char* what() const throw() { return _msg.c_str(); };

 private:
  std::string _msg;
};

struct Null {};

struct Val {
  enum VAL_T { UINT, INT, FLOAT, STRING, ARRAY, DICT, BOOL, JSNULL };
  VAL_T type;
  int i = 0;
  uint64_t ui = 0;
  double f = 0;
  std::string str;
  std::vector<Val> arr;
  std::map<std::string, Val> dict;

  Val() { type = DICT; }
  Val(Null) { type = JSNULL; }
  Val(const std::vector<Val>& arrC) { arr = arrC, type = ARRAY; }
  Val(const std::map<std::string, Val>& dC) { dict = dC, type = DICT; }
  Val(const std::string& strC) { str = strC, type = STRING; }
  Val(const char* strC) { str = strC, type = STRING; }
  Val(double fC) { f = fC, type = FLOAT; }
  Val(size_t iC) { ui = iC, type = UINT; }
  Val(uint32_t iC) { ui = iC, type = UINT; }
  Val(int iC) { i = iC, type = INT; }
  Val(bool fC) { i = fC, type = BOOL; }
};

typedef int Int;
typedef double Float;
typedef bool Bool;
typedef std::string String;
typedef std::vector<Val> Array;
typedef std::map<std::string, Val> Dict;

// simple JSON writer class without much overhead
class Writer {
 public:
  explicit Writer(std::ostream* out);
  Writer(std::ostream* out, size_t prec);
  Writer(std::ostream* out, size_t prec, bool pretty);
  Writer(std::ostream* out, size_t prec, bool pretty, size_t indent);
  ~Writer(){};

  void obj();
  void arr();
  void key(const std::string& k);
  void val(const std::string& v);
  void val(const char* v);
  void val(double v);
  void val(int v);
  void val(uint64_t v);
  void val(bool v);
  void val(Null);
  void val(const Val& v);
  template <typename V>
  void keyVal(const std::string& k, const V& v) {
    key(k);
    val(v);
  }

  void close();
  void closeAll();

 private:
  std::ostream* _out;

  enum NODE_T { OBJ, ARR, KEY };

  struct Node {
    NODE_T type;
    bool empty;
  };

  std::stack<Node> _stack;

  bool _pretty;
  size_t _indent;
  size_t _floatPrec;

  void valCheck();
  void prettor();
};

}  // namespace json
}  // namespace util

#endif  // UTIL_JSON_WRITER_H_
