#include <ccp/ConditionalClades.hpp>
#include <cassert>


int main(int argc, char** argv)
{
  assert(argc == 2);
  ConditionalClades ccp(argv[1]);
  ccp.printStats();
  return 1;
}


