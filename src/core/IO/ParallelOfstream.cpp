#include <IO/ParallelOfstream.hpp>
  
ParallelOfstream::ParallelOfstream(const std::string &fileName, bool masterRankOnly): _os(0)
{
  if (!ParallelContext::getRank() || !masterRankOnly) {
    _os = new std::ofstream(fileName);
  } else {
    _os = new std::ostream(0);
  }
}
  
void ParallelOfstream::close()
{
  delete _os;
  _os = 0;
}

ParallelOfstream::~ParallelOfstream()
{
  close();
}
