#pragma once

#include <parallelization/ParallelContext.hpp>
#include <fstream>
#include <string>
#include <memory>


class ParallelOfstream {
public:
  ParallelOfstream(const std::string &fileName, bool masterRankOnly = true);
  void close();
  ~ParallelOfstream();
private:
  template<typename T> friend std::ostream& operator<<(ParallelOfstream&, T);
  std::unique_ptr<std::ostream> _os;
};

template<typename T> 
std::ostream& operator<<(ParallelOfstream& log, T op) {
  *log._os << op;
  return *log._os;
}

