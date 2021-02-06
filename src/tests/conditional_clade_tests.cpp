#include <ccp/ConditionalClades.hpp>
#include <cassert>


void test(const std::string &newickFile)
{
  ConditionalClades cc(newickFile);
  cc.printContent();
}

int main(int argc, char** argv)
{
  assert(argc == 2);
  test(argv[1]);
}

