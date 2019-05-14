#pragma once

#include <ParallelContext.hpp>
#include <fstream>
#include <string>



class ParallelOfstream {
public:
  ParallelOfstream(const std::string &fileName, bool masterRankOnly = true);
  void close();
  ~ParallelOfstream();
private:
  template<typename T> friend std::ostream& operator<<(ParallelOfstream&, T);
  std::ostream *_os;
};

template<typename T> 
std::ostream& operator<<(ParallelOfstream& log, T op) {
  *log._os << op;
  return *log._os;
}

