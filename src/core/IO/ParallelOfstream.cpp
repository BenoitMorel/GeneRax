#include <IO/ParallelOfstream.hpp>
  
ParallelOfstream::ParallelOfstream(const std::string &fileName, bool masterRankOnly): _os(nullptr)
{
  if (!ParallelContext::getRank() || !masterRankOnly) {
    _os = std::make_unique<std::ofstream>(fileName);
  } else {
    _os = std::make_unique<std::ostream>(nullptr);
  }
}
  
void ParallelOfstream::close()
{
  // force the _os destruction if it holds an ofstream
  _os = std::make_unique<std::ostream>(nullptr);
}

ParallelOfstream::~ParallelOfstream()
{
  close();
}
