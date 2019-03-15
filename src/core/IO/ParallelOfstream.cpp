#include <IO/ParallelOfstream.hpp>
  
ParallelOfstream::ParallelOfstream(const string &fileName): _os(0)
{
  if (!ParallelContext::getRank()) {
    _os = new ofstream(fileName);
  } else {
    _os = new ostream(0);
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
