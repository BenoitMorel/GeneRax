#pragma once
#include <string>

class LibpllException: public std::exception {
public:
  LibpllException(const std::string &s): msg_(s) {}
  LibpllException(const std::string &s1, 
      const std::string s2): msg_(s1 + s2) {}
  virtual const char* what() const noexcept { return msg_.c_str(); }
  void append(const std::string &str) {msg_ += str;}
  virtual ~LibpllException() {}
private:
  std::string msg_;
};

