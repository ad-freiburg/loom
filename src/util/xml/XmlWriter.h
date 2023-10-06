// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_XML_XMLWRITER_H_
#define UTIL_XML_XMLWRITER_H_

#ifdef ZLIB_FOUND
#include <zlib.h>
#endif
#ifdef BZLIB_FOUND
#include <bzlib.h>
#endif
#include <map>
#include <ostream>
#include <fstream>
#include <stack>
#include <string>

namespace util {
namespace xml {

static const size_t BUFFER_S = 32 * 1024 * 1024;

class XmlWriterException : public std::exception {
 public:
  XmlWriterException(std::string msg) : _msg(msg) {}
  ~XmlWriterException() throw() {}

  virtual const char* what() const throw() { return _msg.c_str(); };

 private:
  std::string _msg;
};

// simple XML writer class without much overhead
class XmlWriter {
 public:
  explicit XmlWriter(std::ostream* out);
  XmlWriter(std::ostream* out, bool pretty);
  XmlWriter(std::ostream* out, bool pretty, size_t indent);

  explicit XmlWriter(const std::string& file);
  XmlWriter(const std::string& file, bool pretty);
  XmlWriter(const std::string& file, bool pretty, size_t indent);
  ~XmlWriter(){
#ifdef ZLIB_FOUND
    if (_gzfile) gzclose(_gzfile);
#endif
#ifdef BZLIB_FOUND
    int err;
    if (_bzfile) {
      flushBzip();
      BZ2_bzWriteClose(&err, _bzfile, 0, 0, 0);
    }
#endif
  };

  // open tag without attributes
  void openTag(const std::string& tag);

  // open tag with single attribute (for convenience...)
  void openTag(const std::string& tag, const std::string& key,
               const std::string& val);

  // open tag with attribute list
  void openTag(const std::string& tag,
               const std::map<std::string, std::string>& attrs);

  // open comment
  void openComment();

  // write text
  void writeText(const std::string& text);

  // close tag
  void closeTag();

  // close all open tags, essentially closing the document
  void closeTags();

  // pushes XML escaped text to stream
  void putEsced(const std::string& str, char quot);
  void put(const std::string& str);
  void put(const char c);

 private:
  enum XML_NODE_T { TAG, TEXT, COMMENT };

  struct XmlNode {
    XmlNode(XML_NODE_T t, const std::string& pload, bool hanging)
        : t(t), pload(pload), hanging(hanging) {}
    XML_NODE_T t;
    std::string pload;
    bool hanging;
  };

  void flushBzip();

  std::ostream* _out;
  std::ofstream _outs;
  std::stack<XmlNode> _nstack;

  bool _pretty;
  size_t _indent;

#ifdef ZLIB_FOUND
  gzFile _gzfile;
#else
  int _gzfile;
#endif

  char* _bzbuf;
  size_t _bzbufpos = 0;
#ifdef BZLIB_FOUND
  BZFILE* _bzfile;
#else
  int _bzfile;
#endif

  // handles indentation
  void doIndent();

  // close "hanging" tags
  void closeHanging();


  // checks tag names for validiy
  void checkTagName(const std::string& str) const;
};

}  // namespace xml
}  // namespace util

#endif  // UTIL_XML_XMLWRITER_H_
