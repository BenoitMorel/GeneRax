#pragma once

#include <ParallelContext.hpp>
#include <fstream>
#include <string>

using namespace std;

class ParallelOfstream {
public:
  ParallelOfstream(const string &fileName, bool masterRankOnly = true);
  void close();
  ~ParallelOfstream();
private:
  template<typename T> friend ostream& operator<<(ParallelOfstream&, T);
  ostream *_os;
};

template<typename T> 
ostream& operator<<(ParallelOfstream& log, T op) {
  *log._os << op;
  return *log._os;
}

