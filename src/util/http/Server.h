// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <condition_variable>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <unistd.h>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifndef UTIL_HTTP_SERVER_H_
#define UTIL_HTTP_SERVER_H_

namespace util {
namespace http {

// socket backlog size
const static size_t BLOG = 128;
// socket read buffer size
const static size_t BSIZE = 16 * 1024;
// zlib compression buffer size
const size_t BSIZE_C = 128 * 1024;

// states for HTTP header parser
enum HeaderState { NONE, I_COM, I_URL, I_VER, A_KEY, I_KEY, A_VAL, I_VAL };

/*
 * HTTP Error
 */
class HttpErr : public std::exception {
 public:
  HttpErr(std::string msg) : _msg(msg) {}
  ~HttpErr() throw() {}

  virtual const char* what() const throw() { return _msg.c_str(); }

 private:
  std::string _msg;
};

/*
 * HTTP Request
 */
struct Req {
  std::string cmd, url, ver, payload;
  std::unordered_map<std::string, std::string> params;
  bool gzip = false;
};

/*
 * HTTP Answer
 */
struct Answer {
  Answer() : status(""), pl(""), gzip(false) {}
  Answer(const std::string& status, const std::string& pl)
      : status(status), pl(pl), gzip(false) {}
  Answer(const std::string& status, const std::string& pl, bool gz)
      : status(status), pl(pl), gzip(gz) {}
  std::string status, pl;
  bool gzip = false;
  bool raw = false;
  std::unordered_map<std::string, std::string> params;
};

/*
 * Virtual handler provider class
 */
class Handler {
 public:
  virtual Answer handle(const Req& request, int connection) const = 0;
};

/*
 * Queue of connections to handle
 */
class Queue {
 public:
  void add(int c);
  int get();

 private:
  std::mutex _mut;
  std::queue<int> _jobs;
  std::condition_variable _hasNew;
};

/*
 * Socket wrapper
 */
class Socket {
 public:
  Socket(int port);
  ~Socket();
  int wait();

 private:
  int _sock;
};

/*
 * Simple HTTP server, must provide a pointer to a class instance implementing
 * virtual class Handler.
 */
class HttpServer {
 public:
  HttpServer(int port, const Handler* h) : HttpServer(port, h, 0) {}
  HttpServer(int port, const Handler* h, size_t threads)
      : _port(port), _handler(h), _threads(threads) {
    if (!_threads) _threads = 8 * std::thread::hardware_concurrency();
  }
  void run();

 private:
  int _port;
  Queue _jobs;
  const Handler* _handler;
  size_t _threads;

  void handle();

  static void send(int sock, Answer* aw);
  static Req getReq(int connection);
  static std::string compress(const std::string& str, std::string* enc);
  static bool gzipSupport(const Req& req);
};
}  // http
}  // util

#endif  // UTIL_HTTP_SERVER_H_
