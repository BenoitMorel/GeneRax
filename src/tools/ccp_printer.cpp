#include <ccp/ConditionalClades.hpp>
#include <cassert>
#include <util/enums.hpp>


int main(int argc, char** argv)
{
  assert(argc == 2);
  ConditionalClades ccp(argv[1],  "", CCPRooting::UNIFORM);
  ccp.printStats();
  return 1;
}


