#pragma once

class JointTree;
#include <string>

using namespace std;


class DTLOptimizer {
public:
  static void optimizeDTLRates(JointTree &jointTree, bool transfers);
};


