#include <ccp/ConditionalClades.hpp>
#include <cassert>
#include <util/enums.hpp>

void test(const std::string &newickFile)
{
  ConditionalClades cc(newickFile, "", CCPRooting::UNIFORM);
  cc.printContent();
}

int main(int argc, char** argv)
{
  assert(argc == 2);
  test(argv[1]);
  return 1;
}

